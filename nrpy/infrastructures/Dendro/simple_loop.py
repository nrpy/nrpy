# nrpy/infrastructures/Dendro/simple_loop.py
"""
Dendro interior point-loop emission built on the generic NRPy loop helper.

The interior point loop is x-fastest (``i0`` innermost) over the padded block
interior.  It emits the interior base index ``pp``, which
:meth:`gri.DendroGridFunction.read_gf_from_memory_Ccode_onept` uses for all
one-point reads, and the interior coordinates ``xx0``, ``xx1``, ``xx2``.
"""

import nrpy.helpers.loop as lp
import nrpy.params as par


def _dendro_scalar_alias() -> str:
    """
    Return the live registered Dendro scalar type alias.

    Both scalar spellings come from frozen NRPy parameters (whitepaper
    sections 4.1/6.6); the loop emitter reads the live registry at build
    time and freeze asserts equality with the snapshot.

    :return: The scalar alias (e.g., ``"DendroScalar"``).
    :raises ValueError: If the alias is missing or not a C identifier.
    """
    scalar_type = par.parval_from_str("Dendro_scalar_type")
    if not isinstance(scalar_type, str) or not scalar_type.isidentifier():
        raise ValueError(f"Invalid Dendro scalar type: {scalar_type!r}")
    return scalar_type


def require_serial_parallelization() -> None:
    """
    Require the qualified serial point-loop profile (whitepaper section 7.2).

    The generated point kernel runs inside Dendro's own block traversal, so an
    inner OpenMP pragma would nest parallelism.  Builders assert this rather
    than overwriting the registered value: silently discarding a caller's
    request would produce an unqualified configuration whose manifest
    disagrees with the invocation (section 9.3).

    :raises ValueError: If ``parallelization`` is not ``"none"``.
    """
    parallelization = par.parval_from_str("parallelization")
    if parallelization != "none":
        raise ValueError(
            f"Dendro generation requires parallelization='none', got "
            f"{parallelization!r}: the generated point loop runs inside "
            "Dendro's own block traversal, and nested parallelism is not "
            "qualified (whitepaper section 7.2)."
        )


def interior_loop(
    loop_body: str,
    nx: str = "nx",
    ny: str = "ny",
    nz: str = "nz",
    padding: str = "padding",
    pmin_padded: str = "pmin_padded",
    dx: str = "dx",
) -> str:
    """
    Emit a Dendro interior point loop (x-fastest) around a loop body.

    :param loop_body: C code evaluated at every interior point; has access to
        ``i0``, ``i1``, ``i2``, ``pp``, the stride locals ``nx``, ``ny``,
        ``nxy``, ``padding``, ``xx0``, ``xx1``, ``xx2``, and the spacing
        inverses ``invdxx0``, ``invdxx1``, ``invdxx2``.
    :param nx: C expression for the block's x extent (including padding).
    :param ny: C expression for the block's y extent (including padding).
    :param nz: C expression for the block's z extent (including padding).
    :param padding: C expression for the block's per-axis padding.
    :param pmin_padded: C expression for the coordinate of local padded index
        zero, an array of three.
    :param dx: C expression for the spacing, an array of three.
    :return: The generated nested-loop C code string.
    :raises ValueError: If ``parallelization`` is not ``"none"`` or ``"openmp"``.

    Doctests:
    >>> import nrpy.params as par
    >>> par.set_parval_from_str("parallelization", "none")
    >>> import nrpy.infrastructures.Dendro.generation_parameters  # noqa: F401
    >>> print(interior_loop("f(xx0, xx1, xx2, invdxx1, invdxx2)", nx="NX", ny="NY", nz="NZ"))  # doctest: +ELLIPSIS
    const std::ptrdiff_t nx = ...
    ...
    for (int i2 = static_cast<int>(padding); i2 < static_cast<int>(nz - padding); i2++) {
    ...
    } // END LOOP: for i2 over ...
    <BLANKLINE>
    """
    parallelization = par.parval_from_str("parallelization")
    if parallelization not in ("none", "openmp"):
        raise ValueError(
            "Unsupported Dendro CPU parallelization: "
            f"{parallelization!r} (expected 'none' or 'openmp')."
        )
    pragma = (
        "#pragma omp parallel for collapse(3)" if parallelization == "openmp" else ""
    )
    scalar_type = _dendro_scalar_alias()
    # Hoisted loop invariants (Appendix A): block extents, strides, padding,
    # and spacing inverses are computed once per block, not per point.
    # Bounds use these locals; `nxy` is maybe-unused when a kernel only
    # reaches along x.
    outer_preamble = (
        f"const std::ptrdiff_t nx = static_cast<std::ptrdiff_t>({nx});\n"
        f"const std::ptrdiff_t ny = static_cast<std::ptrdiff_t>({ny});\n"
        f"const std::ptrdiff_t nz = static_cast<std::ptrdiff_t>({nz});\n"
        "[[maybe_unused]] const std::ptrdiff_t nxy = nx * ny;\n"
        f"const std::ptrdiff_t padding = static_cast<std::ptrdiff_t>({padding});\n"
        f"[[maybe_unused]] const {scalar_type} invdxx0 = static_cast<{scalar_type}>(1) / {dx}[0];\n"
        f"[[maybe_unused]] const {scalar_type} invdxx1 = static_cast<{scalar_type}>(1) / {dx}[1];\n"
        f"[[maybe_unused]] const {scalar_type} invdxx2 = static_cast<{scalar_type}>(1) / {dx}[2];\n"
    )
    point_body = (
        "const std::ptrdiff_t pp = i0 + nx * (i1 + ny * i2);\n"
        f"[[maybe_unused]] const {scalar_type} xx0 = {pmin_padded}[0] + static_cast<{scalar_type}>(i0) * {dx}[0];\n"
        f"[[maybe_unused]] const {scalar_type} xx1 = {pmin_padded}[1] + static_cast<{scalar_type}>(i1) * {dx}[1];\n"
        f"[[maybe_unused]] const {scalar_type} xx2 = {pmin_padded}[2] + static_cast<{scalar_type}>(i2) * {dx}[2];\n"
        + loop_body
    )
    return outer_preamble + str(
        lp.loop(
            idx_var=["i2", "i1", "i0"],
            lower_bound=[
                "static_cast<int>(padding)",
                "static_cast<int>(padding)",
                "static_cast<int>(padding)",
            ],
            upper_bound=[
                "static_cast<int>(nz - padding)",
                "static_cast<int>(ny - padding)",
                "static_cast<int>(nx - padding)",
            ],
            increment=["1", "1", "1"],
            pragma=[pragma, "", ""],
            loop_body=point_body,
        )
    )


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    print(f"Doctest passed: All {results.attempted} test(s) passed")
