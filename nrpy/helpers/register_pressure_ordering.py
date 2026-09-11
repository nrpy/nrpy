# nrpy/helpers/register_pressure_ordering.py
"""
Reorder the statements of a generated finite-difference kernel body to reduce register pressure.

c_codegen() emits a kernel body in phases: every gridfunction read, then every finite-difference
derivative, then every common subexpression, then every output. At the end of the read phase all
stencil values are live at once (646 doubles in the BSSN right-hand-side kernel against a 127-double
CUDA register file), and the compiler does not recover from that order. This module re-emits the
same statements in a list-schedule order: at each step the ready statement whose operands die is
emitted first, so a load is issued next to the derivative that consumes it and each output is
written as soon as it is complete.

Every statement keeps its exact text, so the arithmetic is unchanged operation for operation. All
definitions in such a body are single-assignment `const` values, which makes any topological order
legal. Compile-time constants are kept first and stores last, so no store precedes a load from an
array it may alias. On the standard 64x64x128 SinhCylindrical BSSN evolution the reordering alone,
with the finite-difference helpers still marked noinline, removed about 40% of the CUDA
right-hand-side kernel's spill instructions and measured 6% (CUDA) and 3% (OpenMP AVX-512)
faster; with the helpers inlined the combined step measured 13% and 9% faster. Only values that are
zero to roundoff change, because the compiler contracts different multiply/add pairs into fused
multiply-adds.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from typing import Dict, List, Optional, Set, Tuple

_STORE_PREFIXES = ("WriteSIMD", "WriteCUDA")
_LOAD_PREFIXES = ("ReadSIMD", "ReadCUDA")
_CONST_PREFIXES = ("ConstSIMD", "ConstCUDA")


def _tokenize(text: str) -> List[str]:
    """
    Split C text into identifiers, numbers, and single-character punctuation.

    :param text: C source text of one statement.
    :return: List of tokens; whitespace is dropped.

    DocTests:
    >>> _tokenize("const REAL a_b1 = MulSIMD(x, 1.0e-3);")
    ['const', 'REAL', 'a_b1', '=', 'MulSIMD', '(', 'x', ',', '1.0e-3', ')', ';']
    >>> _tokenize("rfmstruct->f0_of_xx0[i0]")
    ['rfmstruct', '-', '>', 'f0_of_xx0', '[', 'i0', ']']
    """
    tokens: List[str] = []
    i = 0
    n = len(text)
    while i < n:
        c = text[i]
        if c.isspace():
            i += 1
        elif c.isalpha() or c == "_":
            j = i + 1
            while j < n and (text[j].isalnum() or text[j] == "_"):
                j += 1
            tokens.append(text[i:j])
            i = j
        elif c.isdigit():
            j = i + 1
            while j < n and (
                text[j].isalnum()
                or text[j] == "."
                or (text[j] in "+-" and text[j - 1] in "eE")
            ):
                j += 1
            tokens.append(text[i:j])
            i = j
        else:
            tokens.append(c)
            i += 1
    return tokens


def _split_statements(body: str) -> Optional[List[str]]:
    r"""
    Drop comments and split a body at semicolons outside parentheses and brackets.

    The phase comments c_codegen() writes describe an order that no longer holds, so they are
    removed rather than moved.

    :param body: C source text.
    :return: Whitespace-normalized statements without the trailing semicolon, or None when the body
             contains braces or preprocessor lines, which this scheduler does not handle.

    DocTests:
    >>> _split_statements("/* reads */ const REAL a = f(x,\n  y); // end\nb[IDX4(0, i0, i1, i2)] = a;")
    ['const REAL a = f(x, y)', 'b[IDX4(0, i0, i1, i2)] = a']
    >>> _split_statements("{ const REAL a = 1.0; }") is None
    True
    >>> _split_statements("const REAL a = 1.0; b = a") is None
    True
    """
    code: List[str] = []
    i = 0
    n = len(body)
    while i < n:
        if body.startswith("/*", i):
            end = body.find("*/", i + 2)
            i = n if end < 0 else end + 2
        elif body.startswith("//", i):
            end = body.find("\n", i)
            i = n if end < 0 else end
        else:
            code.append(body[i])
            i += 1
    text = "".join(code)
    if "{" in text or "}" in text or "#" in text:
        return None
    statements: List[str] = []
    depth = 0
    buf: List[str] = []
    for c in text:
        if c in "([":
            depth += 1
        elif c in ")]":
            depth -= 1
        if c == ";" and depth == 0:
            statements.append(" ".join("".join(buf).split()))
            buf = []
        else:
            buf.append(c)
    if "".join(buf).strip():
        return None
    return [s for s in statements if s]


def _classify(
    statement: str, defined: Set[str]
) -> Optional[Tuple[str, str, List[str]]]:
    """
    Identify a statement as a constant, load, definition or store, and list the defined names it reads.

    :param statement: One statement without its semicolon.
    :param defined: Names defined by earlier statements.
    :return: (kind, defined name or "", operands). kind is "constant" (a `const` definition that
             reads no other statement and no memory, or a ConstSIMD/ConstCUDA broadcast), "load"
             (a `const` definition that reads memory), "definition" (any other `const` definition)
             or "store". None when the statement is none of these, including any non-`const`
             assignment and any left-hand side that is not just specifiers and a name.

    DocTests:
    >>> _classify("MAYBE_UNUSED const REAL u = AddSIMD(a, b)", {"a", "b", "c"})
    ('definition', 'u', ['a', 'b'])
    >>> _classify("static const double dblHalf = 1.0 / 2.0", set())
    ('constant', 'dblHalf', [])
    >>> _classify("const REAL_SIMD_ARRAY Half = ConstSIMD(dblHalf)", {"dblHalf"})
    ('constant', 'Half', ['dblHalf'])
    >>> _classify("const REAL u_i0m1 = ReadSIMD(&in[IDX4(0, i0 - 1, i1, i2)])", set())
    ('load', 'u_i0m1', [])
    >>> _classify("const REAL u_i0m1 = in[IDX4(0, i0 - 1, i1, i2)]", set())
    ('load', 'u_i0m1', [])
    >>> _classify("out[IDX4(0, i0, i1, i2)] = AddSIMD(a, c)", {"a", "b", "c"})
    ('store', '', ['a', 'c'])
    >>> _classify("WriteSIMD(&out[IDX4(0, i0, i1, i2)], a)", {"a"})
    ('store', '', ['a'])
    >>> _classify("if (c) y = 1", {"c"}) is None
    True
    >>> _classify("if (c) const REAL y = 1.0", {"c"}) is None
    True
    >>> _classify("printf(a)", {"a"}) is None
    True
    """
    tokens = _tokenize(statement)
    if tokens and tokens[0] in _STORE_PREFIXES:
        return ("store", "", [t for t in tokens[1:] if t in defined])
    if "=" not in tokens:
        return None
    eq = tokens.index("=")
    lhs, rhs = tokens[:eq], tokens[eq + 1 :]
    operands = [t for t in rhs if t in defined]
    if "[" in lhs:
        return ("store", "", operands)
    if (
        not rhs
        or not ("const" in lhs or "constexpr" in lhs)
        or not all(t[0].isalpha() or t[0] == "_" for t in lhs)
    ):
        return None
    if "[" in rhs or any(t in _LOAD_PREFIXES for t in rhs):
        kind = "load"
    elif rhs[0] in _CONST_PREFIXES or not operands:
        kind = "constant"
    else:
        kind = "definition"
    return (kind, lhs[-1], operands)


def order_statements_for_register_pressure(body: str) -> str:
    """
    Re-emit a generated kernel body with definitions scheduled to shorten live ranges.

    Compile-time constants come first in their original order (this also keeps the SIMD
    `UPWIND_ALG` macro's implicit `upwind_Integer_*` operands ahead of its uses), then loads and
    definitions in greedy free-most-first list-schedule order, then stores in their original order.
    Ties break toward the original order, so the result is deterministic. The body is returned
    unchanged when it contains anything other than single-assignment `const` definitions and
    stores, or when a load follows a store (the load might read what the store wrote).

    :param body: Kernel body produced by c_codegen() with enable_fd_codegen=True.
    :return: The reordered body, or the input when it cannot be scheduled safely.

    DocTests:
    >>> body = '''
    ... /* reads */
    ... const REAL u_i0m1 = ReadSIMD(&in[IDX4(0, i0 - 1, i1, i2)]);
    ... const REAL u_i0p1 = ReadSIMD(&in[IDX4(0, i0 + 1, i1, i2)]);
    ... const REAL v_i0m1 = ReadSIMD(&in[IDX4(1, i0 - 1, i1, i2)]);
    ... const REAL v_i0p1 = ReadSIMD(&in[IDX4(1, i0 + 1, i1, i2)]);
    ... static const double dblHalf = 0.5;
    ... const REAL Half = ConstSIMD(dblHalf);
    ... const REAL u_dD0 = MulSIMD(Half, SubSIMD(u_i0p1, u_i0m1));
    ... const REAL v_dD0 = MulSIMD(Half, SubSIMD(v_i0p1, v_i0m1));
    ... WriteSIMD(&out[IDX4(0, i0, i1, i2)], u_dD0);
    ... WriteSIMD(&out[IDX4(1, i0, i1, i2)], v_dD0);
    ... '''
    >>> print(order_statements_for_register_pressure(body))
    static const double dblHalf = 0.5;
    const REAL Half = ConstSIMD(dblHalf);
    const REAL u_i0m1 = ReadSIMD(&in[IDX4(0, i0 - 1, i1, i2)]);
    const REAL u_i0p1 = ReadSIMD(&in[IDX4(0, i0 + 1, i1, i2)]);
    const REAL u_dD0 = MulSIMD(Half, SubSIMD(u_i0p1, u_i0m1));
    const REAL v_i0m1 = ReadSIMD(&in[IDX4(1, i0 - 1, i1, i2)]);
    const REAL v_i0p1 = ReadSIMD(&in[IDX4(1, i0 + 1, i1, i2)]);
    const REAL v_dD0 = MulSIMD(Half, SubSIMD(v_i0p1, v_i0m1));
    WriteSIMD(&out[IDX4(0, i0, i1, i2)], u_dD0);
    WriteSIMD(&out[IDX4(1, i0, i1, i2)], v_dD0);
    <BLANKLINE>

    A broadcast of a loaded value is not a compile-time constant and stays behind its operand:

    >>> print(order_statements_for_register_pressure(
    ...     "const REAL a = in[IDX4(0, i0, i1, i2)]; const REAL b = ConstSIMD(a); out[IDX4(1, i0, i1, i2)] = b;"))
    const REAL a = in[IDX4(0, i0, i1, i2)];
    const REAL b = ConstSIMD(a);
    out[IDX4(1, i0, i1, i2)] = b;
    <BLANKLINE>

    Bodies the scheduler does not understand come back unchanged:

    >>> unchanged = [
    ...     "x = 1; x = x + 1;",
    ...     "const REAL c = in[IDX4(0, i0, i1, i2)]; if (c) y = 1; const REAL z = y; out[IDX4(1, i0, i1, i2)] = z;",
    ...     "const REAL a = in[X]; out[X] = a; const REAL b = out[X]; out[Y] = b;",
    ...     "const REAL a = 1.0; const REAL a = 2.0;",
    ... ]
    >>> all(order_statements_for_register_pressure(s) == s for s in unchanged)
    True
    """
    statements = _split_statements(body)
    if statements is None:
        return body

    kinds: List[str] = []
    names: List[str] = []
    operands: List[List[str]] = []
    defined: Set[str] = set()
    constant_names: Set[str] = set()
    seen_store = False
    for statement in statements:
        classified = _classify(statement, defined)
        if classified is None or (classified[1] and classified[1] in defined):
            return body
        kind, name, ops = classified
        if kind == "store":
            seen_store = True
        elif kind == "load" and seen_store:
            return body
        elif kind == "constant" and not all(op in constant_names for op in ops):
            kind = "definition"
        kinds.append(kind)
        names.append(name)
        operands.append(ops)
        if name:
            defined.add(name)
        if kind == "constant":
            constant_names.add(name)

    # Remaining consumers of each defined name; a definition is ready once all its operands exist.
    remaining_uses: Dict[str, int] = {name: 0 for name in names if name}
    consumers: Dict[str, List[int]] = {name: [] for name in names if name}
    for k, ops in enumerate(operands):
        for op in set(ops):
            remaining_uses[op] += 1
            consumers[op].append(k)

    constants = [k for k, kind in enumerate(kinds) if kind == "constant"]
    stores = [k for k, kind in enumerate(kinds) if kind == "store"]
    scheduled: Set[int] = set(constants)
    waiting = {k: len(set(operands[k])) for k in range(len(statements))}
    for k in constants:
        for consumer in consumers[names[k]]:
            waiting[consumer] -= 1
    ready = [
        k
        for k in range(len(statements))
        if k not in scheduled and kinds[k] in ("load", "definition") and waiting[k] == 0
    ]
    order: List[int] = list(constants)
    while ready:
        best = min(
            ready,
            key=lambda k: (
                -(
                    sum(1 for op in set(operands[k]) if remaining_uses[op] == 1)
                    - (1 if remaining_uses[names[k]] > 0 else 0)
                ),
                k,
            ),
        )
        ready.remove(best)
        scheduled.add(best)
        order.append(best)
        for op in set(operands[best]):
            remaining_uses[op] -= 1
        for consumer in consumers[names[best]]:
            waiting[consumer] -= 1
            if waiting[consumer] == 0 and kinds[consumer] in ("load", "definition"):
                ready.append(consumer)
    order += stores
    if len(order) != len(statements):
        return body
    return "".join(f"{statements[k]};\n" for k in order)


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
