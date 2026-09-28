# nrpy/infrastructures/Dendro/simple_loop.py
"""
Dendro interior point-loop emission built on the generic NRPy loop helper.

The interior point loop is x-fastest (``i0`` innermost) over the padded block
interior.  It emits the interior base index ``pp``, which
:meth:`gri.DendroGridFunction.read_gf_from_memory_Ccode_onept` uses for all
one-point reads, and the interior coordinates ``xx0``, ``xx1``, ``xx2`` in
scalar mode.

Deliberate divergence from the established ``simple_loop`` signature: BHaH and
ETLegacy take ``loop_region``, ``enable_OpenMP``, ``OMP_custom_pragma`` and
``OMP_collapse``.  This form takes none of them and takes the block extents,
padding, padded origin and spacing instead, because a Dendro point loop runs
inside Dendro's own block traversal over one padded block: the interior is the
only region a generated kernel writes, and an inner OpenMP pragma would nest
parallelism.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

import nrpy.grid as gri
import nrpy.helpers.loop as lp
import nrpy.params as par


def simple_loop(
    loop_body: str,
    nx: str = "geom.nx",
    ny: str = "geom.ny",
    nz: str = "geom.nz",
    padding: str = "geom.padding",
    pmin_padded: str = "geom.pmin_padded",
    dx: str = "geom.dx",
    enable_intrinsics: bool = False,
) -> str:
    """
    Emit a Dendro interior point loop (x-fastest) around a loop body.

    With ``enable_intrinsics``, ``i0`` advances by ``SIMD_WIDTH`` and each
    vector starts at ``min(i0_vector, max(0, nx - padding - SIMD_WIDTH))``.
    When a row has at least ``SIMD_WIDTH`` interior points, its final vector
    ends at the last interior point, recomputing a few interior points
    already written by the previous vector.  A shorter row has one vector,
    which covers the interior and padding points of the same row.  Every
    load stays inside the block and every store inside its own row, so no
    remainder loop is needed.

    :param loop_body: C code evaluated at every interior point; has access to
        ``i0``, ``i1``, ``i2``, ``pp``, the stride locals ``nx``, ``ny``,
        ``nxy``, ``padding``, and the spacing inverses ``invdxx0``,
        ``invdxx1``, ``invdxx2``. Scalar mode also provides coordinates.
    :param nx: C expression for the block's x extent (including padding).
    :param ny: C expression for the block's y extent (including padding).
    :param nz: C expression for the block's z extent (including padding).
    :param padding: C expression for the block's per-axis padding.
    :param pmin_padded: C expression for the coordinate of local padded index
        zero, an array of three.
    :param dx: C expression for the spacing, an array of three.
    :param enable_intrinsics: Emit ``SIMD_WIDTH``-point vector iterations with
        vector spacing inverses and no coordinates.
    :return: The generated nested-loop C code string.
    :raises ValueError: If nested point-loop parallelism is requested.

    Doctests:
    >>> import nrpy.params as par
    >>> par.set_parval_from_str("parallelization", "none")
    >>> print(simple_loop("f(xx0, xx1, xx2, invdxx1, invdxx2)", nx="NX", ny="NY", nz="NZ", padding="PAD", pmin_padded="PMIN", dx="DX"))
    const std::ptrdiff_t nx = static_cast<std::ptrdiff_t>(NX);
    const std::ptrdiff_t ny = static_cast<std::ptrdiff_t>(NY);
    const std::ptrdiff_t nz = static_cast<std::ptrdiff_t>(NZ);
    [[maybe_unused]] const std::ptrdiff_t nxy = nx * ny;
    const std::ptrdiff_t padding = static_cast<std::ptrdiff_t>(PAD);
    [[maybe_unused]] const DendroScalar invdxx0 = static_cast<DendroScalar>(1) / DX[0];
    [[maybe_unused]] const DendroScalar invdxx1 = static_cast<DendroScalar>(1) / DX[1];
    [[maybe_unused]] const DendroScalar invdxx2 = static_cast<DendroScalar>(1) / DX[2];
    for (int i2 = static_cast<int>(padding); i2 < static_cast<int>(nz - padding); i2++) {
    for (int i1 = static_cast<int>(padding); i1 < static_cast<int>(ny - padding); i1++) {
    for (int i0 = static_cast<int>(padding); i0 < static_cast<int>(nx - padding); i0++) {
    const std::ptrdiff_t pp = i0 + nx * (i1 + ny * i2);
    [[maybe_unused]] const DendroScalar xx0 = PMIN[0] + static_cast<DendroScalar>(i0) * DX[0];
    [[maybe_unused]] const DendroScalar xx1 = PMIN[1] + static_cast<DendroScalar>(i1) * DX[1];
    [[maybe_unused]] const DendroScalar xx2 = PMIN[2] + static_cast<DendroScalar>(i2) * DX[2];
    f(xx0, xx1, xx2, invdxx1, invdxx2)
    } // END LOOP: for i0 over [static_cast<int>(padding), static_cast<int>(nx - padding))
    } // END LOOP: for i1 over [static_cast<int>(padding), static_cast<int>(ny - padding))
    } // END LOOP: for i2 over [static_cast<int>(padding), static_cast<int>(nz - padding))
    <BLANKLINE>
    >>> print(simple_loop("WriteSIMD(&rhs[pp], value);", enable_intrinsics=True))
    const std::ptrdiff_t nx = static_cast<std::ptrdiff_t>(geom.nx);
    const std::ptrdiff_t ny = static_cast<std::ptrdiff_t>(geom.ny);
    const std::ptrdiff_t nz = static_cast<std::ptrdiff_t>(geom.nz);
    [[maybe_unused]] const std::ptrdiff_t nxy = nx * ny;
    const std::ptrdiff_t padding = static_cast<std::ptrdiff_t>(geom.padding);
    [[maybe_unused]] const REAL_SIMD_ARRAY invdxx0 = ConstSIMD(static_cast<DendroScalar>(1) / geom.dx[0]);
    [[maybe_unused]] const REAL_SIMD_ARRAY invdxx1 = ConstSIMD(static_cast<DendroScalar>(1) / geom.dx[1]);
    [[maybe_unused]] const REAL_SIMD_ARRAY invdxx2 = ConstSIMD(static_cast<DendroScalar>(1) / geom.dx[2]);
    if (nx < SIMD_WIDTH) throw std::invalid_argument("Dendro block row is shorter than SIMD_WIDTH");
    const std::ptrdiff_t last_vector_start = nx - padding - SIMD_WIDTH > 0 ? nx - padding - SIMD_WIDTH : 0;
    for (int i2 = static_cast<int>(padding); i2 < static_cast<int>(nz - padding); i2++) {
    for (int i1 = static_cast<int>(padding); i1 < static_cast<int>(ny - padding); i1++) {
    for (int i0_vector = static_cast<int>(padding); i0_vector < static_cast<int>(nx - padding); i0_vector += SIMD_WIDTH) {
    const std::ptrdiff_t i0 = i0_vector < last_vector_start ? i0_vector : last_vector_start;
    const std::ptrdiff_t pp = i0 + nx * (i1 + ny * i2);
    WriteSIMD(&rhs[pp], value);
    } // END LOOP: for i0_vector over [static_cast<int>(padding), static_cast<int>(nx - padding))
    } // END LOOP: for i1 over [static_cast<int>(padding), static_cast<int>(ny - padding))
    } // END LOOP: for i2 over [static_cast<int>(padding), static_cast<int>(nz - padding))
    <BLANKLINE>
    """
    parallelization = par.parval_from_str("parallelization")
    if parallelization != "none":
        raise ValueError(
            f"Dendro generation requires parallelization='none', got "
            f"{parallelization!r}: the generated point loop runs inside "
            "Dendro's own block traversal, and nested parallelism is not "
            "qualified."
        )
    # One spelling of the scalar alias, the core constant every Dendro emitter
    # reads.
    scalar_type = gri.DENDRO_SCALAR_TYPE
    # Hoisted loop invariants: block extents, strides, padding,
    # and spacing inverses are computed once per block, not per point.
    # Bounds use these locals; `nxy` is maybe-unused when a kernel only
    # reaches along x.
    spacing_type = "REAL_SIMD_ARRAY" if enable_intrinsics else scalar_type
    outer_preamble = (
        f"const std::ptrdiff_t nx = static_cast<std::ptrdiff_t>({nx});\n"
        f"const std::ptrdiff_t ny = static_cast<std::ptrdiff_t>({ny});\n"
        f"const std::ptrdiff_t nz = static_cast<std::ptrdiff_t>({nz});\n"
        "[[maybe_unused]] const std::ptrdiff_t nxy = nx * ny;\n"
        f"const std::ptrdiff_t padding = static_cast<std::ptrdiff_t>({padding});\n"
    )
    for direction in range(3):
        inverse = f"static_cast<{scalar_type}>(1) / {dx}[{direction}]"
        if enable_intrinsics:
            inverse = f"ConstSIMD({inverse})"
        outer_preamble += (
            f"[[maybe_unused]] const {spacing_type} invdxx{direction} = {inverse};\n"
        )
    loop_bounds = {
        "lower_bound": ["static_cast<int>(padding)"] * 3,
        "upper_bound": [
            "static_cast<int>(nz - padding)",
            "static_cast<int>(ny - padding)",
            "static_cast<int>(nx - padding)",
        ],
        "pragma": ["", "", ""],
    }
    point_index = "const std::ptrdiff_t pp = i0 + nx * (i1 + ny * i2);\n"
    if enable_intrinsics:
        outer_preamble += (
            "if (nx < SIMD_WIDTH) throw std::invalid_argument("
            '"Dendro block row is shorter than SIMD_WIDTH");\n'
            "const std::ptrdiff_t last_vector_start = nx - padding - SIMD_WIDTH > 0 "
            "? nx - padding - SIMD_WIDTH : 0;\n"
        )
        return outer_preamble + str(
            lp.loop(
                idx_var=["i2", "i1", "i0_vector"],
                increment=["1", "1", "SIMD_WIDTH"],
                loop_body=(
                    "const std::ptrdiff_t i0 = i0_vector < last_vector_start "
                    "? i0_vector : last_vector_start;\n" + point_index + loop_body
                ),
                **loop_bounds,
            )
        )
    point_body = (
        point_index
        + f"[[maybe_unused]] const {scalar_type} xx0 = {pmin_padded}[0] + static_cast<{scalar_type}>(i0) * {dx}[0];\n"
        f"[[maybe_unused]] const {scalar_type} xx1 = {pmin_padded}[1] + static_cast<{scalar_type}>(i1) * {dx}[1];\n"
        f"[[maybe_unused]] const {scalar_type} xx2 = {pmin_padded}[2] + static_cast<{scalar_type}>(i2) * {dx}[2];\n"
        + loop_body
    )
    return outer_preamble + str(
        lp.loop(
            idx_var=["i2", "i1", "i0"],
            increment=["1", "1", "1"],
            loop_body=point_body,
            **loop_bounds,
        )
    )


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
