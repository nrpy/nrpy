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

from typing import Dict, List, Set

_STORE_PREFIXES = ("WriteSIMD", "WriteCUDA")
_LOAD_PREFIXES = ("ReadSIMD", "ReadCUDA")
_CONST_PREFIXES = ("ConstSIMD", "ConstCUDA")


def order_statements_for_register_pressure(body: str) -> str:
    r"""
    Re-emit a generated kernel body with definitions scheduled to shorten live ranges.

    Compile-time constants come first in their original order (this also keeps the SIMD
    `UPWIND_ALG` macro's implicit `upwind_Integer_*` operands ahead of its uses), then loads and
    definitions in greedy free-most-first list-schedule order, then stores in their original order.
    Ties break toward the original order, so the result is deterministic. The body is returned
    unchanged when it contains anything other than single-assignment `const` definitions and
    stores, or when a load follows a store (the load might read what the store wrote).

    :param body: Kernel body produced by c_codegen() with enable_fd_codegen=True.
    :return: The reordered body, or the input when it cannot be scheduled safely.

    Doctests:
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

    Comments are dropped, a statement may span lines, and a `MAYBE_UNUSED` or `static` specifier,
    a number with an exponent, a struct member read and a plain array store are all understood:

    >>> print(order_statements_for_register_pressure(
    ...     "/* reads */ const REAL f0 = rfmstruct->f0_of_xx0[i0]; // rfm\n"
    ...     "static const double dblEps = 1.0e-3; MAYBE_UNUSED const REAL a = MulSIMD(f0,\n  dblEps);\n"
    ...     "out[IDX4(0, i0, i1, i2)] = a;"))
    static const double dblEps = 1.0e-3;
    const REAL f0 = rfmstruct->f0_of_xx0[i0];
    MAYBE_UNUSED const REAL a = MulSIMD(f0, dblEps);
    out[IDX4(0, i0, i1, i2)] = a;
    <BLANKLINE>

    A broadcast of a loaded value is not a compile-time constant and stays behind its operand:

    >>> print(order_statements_for_register_pressure(
    ...     "const REAL a = in[IDX4(0, i0, i1, i2)]; const REAL b = ConstSIMD(a); out[IDX4(1, i0, i1, i2)] = b;"))
    const REAL a = in[IDX4(0, i0, i1, i2)];
    const REAL b = ConstSIMD(a);
    out[IDX4(1, i0, i1, i2)] = b;
    <BLANKLINE>

    Bodies the scheduler does not understand come back unchanged: braces or preprocessor lines, a
    statement without its semicolon, a non-`const` assignment, a left-hand side that is not just
    specifiers and a name, a call statement, a redefinition, or a load after a store:

    >>> unchanged = [
    ...     "{ const REAL a = 1.0; }",
    ...     "const REAL a = 1.0; b = a",
    ...     "x = 1; x = x + 1;",
    ...     "const REAL c = in[IDX4(0, i0, i1, i2)]; if (c) y = 1; const REAL z = y; out[IDX4(1, i0, i1, i2)] = z;",
    ...     "if (c) const REAL y = 1.0; out[X] = y;",
    ...     "printf(a);",
    ...     "const REAL a = 1.0; const REAL a = 2.0;",
    ...     "const REAL a = in[X]; out[X] = a; const REAL b = out[X]; out[Y] = b;",
    ... ]
    >>> all(order_statements_for_register_pressure(s) == s for s in unchanged)
    True
    """
    # Step 1: Drop comments and split at semicolons outside parentheses and brackets. The phase
    #         comments c_codegen() writes describe an order that no longer holds, so they are
    #         removed rather than moved. Braces and preprocessor lines are not handled.
    code: List[str] = []
    i = 0
    while i < len(body):
        if body.startswith("/*", i):
            end = body.find("*/", i + 2)
            i = len(body) if end < 0 else end + 2
        elif body.startswith("//", i):
            end = body.find("\n", i)
            i = len(body) if end < 0 else end
        else:
            code.append(body[i])
            i += 1
    text = "".join(code)
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
    if any(ch in text for ch in "{}#") or "".join(buf).strip():
        return body
    statements = [s for s in statements if s]

    # Step 2: Classify each statement as a constant (a `const` definition that reads no other
    #         statement and no memory, or a ConstSIMD/ConstCUDA broadcast), a load (a `const`
    #         definition that reads memory), a definition (any other `const` definition) or a
    #         store, and record the defined names it reads. Anything else, including a
    #         non-`const` assignment or a left-hand side that is not just specifiers and a name,
    #         leaves the body unchanged, as does a redefinition or a load after a store.
    kinds: List[str] = []
    names: List[str] = []
    operands: List[List[str]] = []
    defined: Set[str] = set()
    constant_names: Set[str] = set()
    seen_store = False
    for statement in statements:
        # Step 2.a: Split into identifiers, numbers and single-character punctuation.
        tokens: List[str] = []
        i = 0
        while i < len(statement):
            c = statement[i]
            if c.isspace():
                i += 1
            elif c.isalpha() or c == "_":
                j = i + 1
                while j < len(statement) and (
                    statement[j].isalnum() or statement[j] == "_"
                ):
                    j += 1
                tokens.append(statement[i:j])
                i = j
            elif c.isdigit():
                j = i + 1
                while j < len(statement) and (
                    statement[j].isalnum()
                    or statement[j] == "."
                    or (statement[j] in "+-" and statement[j - 1] in "eE")
                ):
                    j += 1
                tokens.append(statement[i:j])
                i = j
            else:
                tokens.append(c)
                i += 1
        # Step 2.b: Name the statement's kind, the name it defines and the names it reads.
        kind, name, ops = "", "", []
        if tokens and tokens[0] in _STORE_PREFIXES:
            kind, ops = "store", [t for t in tokens[1:] if t in defined]
        elif "=" in tokens:
            eq = tokens.index("=")
            lhs, rhs = tokens[:eq], tokens[eq + 1 :]
            ops = [t for t in rhs if t in defined]
            if "[" in lhs:
                kind = "store"
            elif (
                rhs
                and ("const" in lhs or "constexpr" in lhs)
                and all(t[0].isalpha() or t[0] == "_" for t in lhs)
            ):
                name = lhs[-1]
                if "[" in rhs or any(t in _LOAD_PREFIXES for t in rhs):
                    kind = "load"
                elif rhs[0] in _CONST_PREFIXES or not ops:
                    kind = "constant"
                else:
                    kind = "definition"
        if kind == "store":
            seen_store = True
        elif kind == "constant" and not all(op in constant_names for op in ops):
            kind = "definition"
        if not kind or (name and name in defined) or (kind == "load" and seen_store):
            return body
        kinds.append(kind)
        names.append(name)
        operands.append(ops)
        if name:
            defined.add(name)
        if kind == "constant":
            constant_names.add(name)

    # Step 3: Remaining consumers of each defined name; a definition is ready once all its
    #         operands exist.
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
