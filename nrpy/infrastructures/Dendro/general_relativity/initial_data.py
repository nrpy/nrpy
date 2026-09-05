# nrpy/infrastructures/Dendro/general_relativity/initial_data.py
"""
PR 7/8: initial-data CFunctions for the Dendro backend.

The generated Minkowski fill writes every EVOL field to its asymptotic value
through exact-name ``rhs_<name>`` bindings derived from the frozen EVOL
registry (whitepaper section 14.1: the first end-to-end proof of state
allocation, field order, loops, and CFunctions).  No field name, count, or
asymptotic value is hardcoded in the builder: all three come from the
registered gridfunction records.

PR 8 adds the smooth ADM-to-fCCZ4 conversion (whitepaper section 14.2) and
the separate connection (``lambdaU``) initialization pass (section 14.3).
The conversion reuses the established
:class:`nrpy.equations.general_relativity.ADM_to_BSSN.ADM_to_BSSN` map, so
the registered ``EvolvedConformalFactor_cf`` choice governs the conformal
convention in exactly one place; the connection pass writes
``lambdaU^i = DeltaGamma^i / ReU^i``, which is the statement that the
connection constraint ``C^i = 0`` holds.
"""

from typing import Dict, List, Tuple

import sympy as sp

import nrpy.grid as gri
import nrpy.indexedexp as ixp
import nrpy.params as par
import nrpy.reference_metric as refmetric
from nrpy.c_codegen import c_codegen
from nrpy.equations.general_relativity.ADM_to_BSSN import ADM_to_BSSN
from nrpy.equations.general_relativity.BSSN_quantities import BSSN_quantities
from nrpy.infrastructures.Dendro import access_capture as cap
from nrpy.infrastructures.Dendro import generation_parameters
from nrpy.infrastructures.Dendro import gridfunction_output as gfo
from nrpy.infrastructures.Dendro import naming
from nrpy.infrastructures.Dendro import registration as reg
from nrpy.infrastructures.Dendro.block_loop import block_loop
from nrpy.infrastructures.Dendro.simple_loop import (
    interior_loop,
    require_serial_parallelization,
)

# CFunction names (Dendro scheduling role keys).
MINKOWSKI_BLOCK_CFUNCTION = "fccz4_minkowski_initial_data_block"
MINKOWSKI_GLOBAL_CFUNCTION = "fccz4_minkowski_initial_data"

# Smooth analytic perturbation used by the lifecycle gates (section 16.8): a
# spatially varying state is what makes the generated derivative stencils
# observable at all.
PERTURBATION_BLOCK_CFUNCTION = "fccz4_smooth_perturbation_block"
PERTURBATION_GLOBAL_CFUNCTION = "fccz4_smooth_perturbation"

# PR 8 CFunction names (whitepaper section 7.1).
ADM_TO_EVOLVED_BLOCK_CFUNCTION = "fccz4_adm_to_evolved_block"
INITIALIZE_LAMBDA_BLOCK_CFUNCTION = "fccz4_initialize_lambda_block"


def _block_pointer_bindings(evol_order: Tuple[str, ...], scalar_type: str) -> str:
    """
    Emit the per-field output pointer bindings for the block layout.

    The bindings are rendered by the single shared emitter from the registry
    order, so no field name is hardcoded; every binding adds
    ``geom.component_offset`` (section 8.2).  The role is ``out_`` because
    this writer produces *state*, not a right-hand side (section 5.2).

    :param evol_order: The EVOL names, in registry order.
    :param scalar_type: The registered Dendro scalar alias.
    :return: The binding statements.
    """
    return gfo.render_component_bindings(
        evol_order,
        scalar_type,
        array="out_gfs",
        role=naming.out_pointer,
        const_pointee=False,
        index_expression=lambda _name, position: str(position),
    )


def build_minkowski_initial_data() -> Tuple[str, str, str, str]:
    """
    Build the Minkowski initial data block and all-block CFunction bodies.

    Each EVOL field is written to its registered asymptotic value
    (``f_infinity``) at every cell of the padded block.  The point loop is
    emitted through the NRPy Dendro loop helper; the all-block entry point
    wraps it in the NRPy numerical block loop (sections 7.2/7.3/14.1).

    :return: (block_body, block_params, global_body, global_params).
    :raises ValueError: If Infrastructure is not Dendro.
    """
    if par.parval_from_str("Infrastructure") != "Dendro":
        raise ValueError(
            "Infrastructure must be 'Dendro' to build the Minkowski initial "
            f"data, got {par.parval_from_str('Infrastructure')!r}."
        )
    require_serial_parallelization()
    scalar_type = str(par.parval_from_str("Dendro_scalar_type"))
    evol_order = reg.registered_evol_order()

    # The fill writes every EVOL field to its registered asymptotic value.
    # Access capture is not applicable: no neighbor reads occur, and the
    # leaf CFunction has no calls, so freeze's capture gate would reject it
    # - the writer is therefore registered with an empty capture through the
    # explicit ``non_numerical`` opt-out below.
    fill_lines: List[str] = []
    for name in evol_order:
        gf = gri.glb_gridfcs_dict[name]
        value = str(gf.f_infinity)
        # Render a double literal (registered values are float).
        fill_lines.append(f"{naming.out_pointer(name)}[pp] = {scalar_type}{{{value}}};")
    fill_body = "\n".join(fill_lines)

    # Minkowski is spatially constant: every cell of the padded block (ghost
    # cells included) equals the registered asymptotic value.  Filling the
    # whole block (padding 0) makes the state a true fixed point, so the
    # finite-difference RHS vanishes everywhere, including the interior
    # cells adjacent to a block boundary that read ghost cells.
    point_loop_body = interior_loop(
        fill_body,
        nx="geom.nx",
        ny="geom.ny",
        nz="geom.nz",
        padding="0",
        pmin_padded="geom.pmin_padded",
        dx="geom.dx",
    )
    block_body = (
        _block_pointer_bindings(evol_order, scalar_type) + "\n" + point_loop_body
    )
    block_params = f"const BlockGeometry& geom, {scalar_type}* const* out_gfs"
    global_body = block_loop(
        f"{MINKOWSKI_BLOCK_CFUNCTION}(world.geom[blk_id], out_gfs);",
        count="world.num_blocks",
    )
    global_params = f"const MockWorld& world, {scalar_type}* const* out_gfs"
    return block_body, block_params, global_body, global_params


def build_smooth_perturbation() -> Tuple[str, str, str, str]:
    """
    Build the smooth analytic perturbation CFunction bodies.

    The perturbation adds one bounded, infinitely differentiable analytic
    profile to every evolved component, so a lifecycle test can exercise the
    generated derivative stencils on a state that is not spatially constant
    (a Minkowski state has an identically vanishing RHS whatever the stencil
    coefficients).  The profile is authored in NRPy and lowered by
    ``c_codegen``: it is formulation content, so it belongs in a registered
    CFunction and never in a fixed template (whitepaper section 9.16).

    :return: (block_body, block_params, global_body, global_params).
    :raises ValueError: If Infrastructure is not Dendro.
    """
    if par.parval_from_str("Infrastructure") != "Dendro":
        raise ValueError(
            "Infrastructure must be 'Dendro' to build the Dendro perturbation, "
            f"got {par.parval_from_str('Infrastructure')!r}."
        )
    require_serial_parallelization()
    scalar_type = str(par.parval_from_str("Dendro_scalar_type"))
    evol_order = reg.registered_evol_order()

    amplitude, wavelength = sp.symbols("amplitude wavelength", real=True)
    xx0, xx1, xx2 = sp.symbols("xx0 xx1 xx2", real=True)
    wavenumber = 2 * sp.pi / wavelength
    profile = (
        amplitude
        * sp.sin(wavenumber * xx0)
        * sp.cos(wavenumber * xx1)
        * sp.sin(wavenumber * xx2)
    )
    profile_code = c_codegen(
        [profile],
        [f"const {scalar_type} smooth_profile"],
        include_braces=False,
        enable_simd=False,
        fp_type=str(par.parval_from_str("fp_type")),
        fp_type_alias=scalar_type,
        verbose=False,
    )
    fill_lines: List[str] = [profile_code.strip()]
    for name in evol_order:
        fill_lines.append(f"{naming.out_pointer(name)}[pp] += smooth_profile;")
    # The perturbation is applied over the whole padded block, ghost cells
    # included, so the interior RHS sees a consistent field on every stencil.
    point_loop_body = interior_loop(
        "\n".join(fill_lines),
        nx="geom.nx",
        ny="geom.ny",
        nz="geom.nz",
        padding="0",
        pmin_padded="geom.pmin_padded",
        dx="geom.dx",
    )
    block_body = (
        _block_pointer_bindings(evol_order, scalar_type) + "\n" + point_loop_body
    )
    block_params = (
        f"const BlockGeometry& geom, {scalar_type}* const* out_gfs, "
        f"const {scalar_type} amplitude, const {scalar_type} wavelength"
    )
    global_body = block_loop(
        f"{PERTURBATION_BLOCK_CFUNCTION}"
        "(world.geom[blk_id], out_gfs, amplitude, wavelength);",
        count="world.num_blocks",
    )
    global_params = (
        f"const MockWorld& world, {scalar_type}* const* out_gfs, "
        f"const {scalar_type} amplitude, const {scalar_type} wavelength"
    )
    return block_body, block_params, global_body, global_params


def register_perturbation_CFunctions() -> None:
    """
    Register the smooth-perturbation CFunctions (with Dendro roles).

    The writer touches only the current point, so it is registered with an
    explicit empty capture record (section 8.6).
    """
    block_body, block_params, global_body, global_params = build_smooth_perturbation()
    cap.record_empty_capture(PERTURBATION_BLOCK_CFUNCTION)
    reg.register_dendro_CFunction(
        role="perturbation_block",
        entry_point=False,
        lifecycle_hook="diagnostic",
        name=PERTURBATION_BLOCK_CFUNCTION,
        desc=(
            "Per-block smooth analytic perturbation of every evolved field "
            "(lifecycle-test state; NRPy-authored profile)."
        ),
        subdirectory="generated/src/initial_data",
        params=block_params,
        body=block_body,
    )
    reg.register_dendro_CFunction(
        role="perturbation",
        entry_point=True,
        calls=(PERTURBATION_BLOCK_CFUNCTION,),
        lifecycle_hook="diagnostic",
        name=PERTURBATION_GLOBAL_CFUNCTION,
        desc="All-block smooth analytic perturbation (NRPy block loop).",
        subdirectory="generated/src/initial_data",
        params=global_params,
        body=global_body,
    )


def register_minkowski_CFunctions() -> None:
    """
    Register the Minkowski initial data CFunctions (with Dendro roles).

    The block writer reads no neighbors, so it is registered with an explicit
    empty capture record (padding zero is correct: the fill touches only the
    current point).

    """
    block_body, block_params, global_body, global_params = (
        build_minkowski_initial_data()
    )
    # The fill reads no neighbors: record an explicit empty capture so
    # freeze's leaf-capture gate is satisfied with the correct (zero)
    # padding (section 8.6).
    cap.record_empty_capture(MINKOWSKI_BLOCK_CFUNCTION)
    reg.register_dendro_CFunction(
        role="minkowski_initial_data_block",
        entry_point=False,
        lifecycle_hook="initial_data",
        name=MINKOWSKI_BLOCK_CFUNCTION,
        desc="PR 7 per-block Minkowski initial data fill (all EVOL fields to their asymptotic values).",
        subdirectory="generated/src/initial_data",
        params=block_params,
        body=block_body,
    )
    reg.register_dendro_CFunction(
        role="minkowski_initial_data",
        entry_point=True,
        calls=(MINKOWSKI_BLOCK_CFUNCTION,),
        lifecycle_hook="initial_data",
        name=MINKOWSKI_GLOBAL_CFUNCTION,
        desc="PR 7 all-block Minkowski initial data fill (NRPy block loop).",
        subdirectory="generated/src/initial_data",
        params=global_params,
        body=global_body,
    )


def register_adm_source_gridfunctions() -> (
    Tuple[List[List[sp.Expr]], List[List[sp.Expr]], List[sp.Expr], List[sp.Expr]]
):
    """
    Register the ADM source fields the converter reads, as AUXEVOL.

    The physical ADM data (``gammaDD``, ``KDD``, ``betaU``, ``BU``) is supplied
    by the host (Dendro or TwoPunctures) and is not evolved, so it is
    registered in the AUXEVOL group under exact NRPy names (whitepaper section
    14.2).  Registration is idempotent: a second call returns the symbols of
    the already-registered fields.

    :return: (gammaDD, KDD, betaU, BU) symbol containers.
    """
    if "gammaDD00" in gri.glb_gridfcs_dict:
        return (
            ixp.declarerank2("gammaDD", symmetry="sym01"),
            ixp.declarerank2("KDD", symmetry="sym01"),
            ixp.declarerank1("betaU"),
            ixp.declarerank1("BU"),
        )
    gammaDD = gri.register_gridfunctions_for_single_rank2(
        "gammaDD", symmetry="sym01", group="AUXEVOL"
    )
    KDD = gri.register_gridfunctions_for_single_rank2(
        "KDD", symmetry="sym01", group="AUXEVOL"
    )
    betaU = gri.register_gridfunctions_for_single_rank1("betaU", group="AUXEVOL")
    BU = gri.register_gridfunctions_for_single_rank1("BU", group="AUXEVOL")
    return gammaDD, KDD, betaU, BU


def build_adm_to_evolved() -> Tuple[str, str]:
    """
    Build the smooth ADM-to-fCCZ4 conversion CFunction body (section 14.2).

    The conversion is pointwise (no neighbour reads), so it runs over the whole
    padded block and its padding is zero.  The three connection components are
    deliberately absent: they depend on metric derivatives and are written by
    the separate pass below.

    :return: (block_body, block_params).
    :raises ValueError: If Infrastructure is not Dendro, if the conversion
        leaves an evolved field undefined, or if it reads outside the
        registered ADM source set.

    Doctests:
    >>> import contextlib
    >>> import io
    >>> from nrpy.helpers.generic import validate_strings
    >>> from nrpy.equations.general_relativity.fCCZ4_system import (
    ...     build_fccz4_expression_bundle,
    ... )
    >>> par.set_parval_from_str("Infrastructure", "Dendro")
    >>> par.set_parval_from_str("parallelization", "none")
    >>> par.set_parval_from_str("fd_order", 4)
    >>> par.set_parval_from_str("EvolvedConformalFactor_cf", "chi")
    >>> with contextlib.redirect_stdout(io.StringIO()):
    ...     _bundle = build_fccz4_expression_bundle()
    >>> _body, _params = build_adm_to_evolved()

    Section 14.2: the conversion defines every evolved field except the three
    connection components, which the separate pass below owns.

    >>> _written = [
    ...     name
    ...     for name in reg.registered_evol_order()
    ...     if naming.out_pointer(name) + "[pp] =" in _body
    ... ]
    >>> len(_written)
    22
    >>> [n for n in ("lambdaU0", "lambdaU1", "lambdaU2") if n in _written]
    []

    The ADM source data is registered as AUXEVOL under exact NRPy names.

    >>> reg.registered_auxevol_order()[:3]
    ('betaU0', 'betaU1', 'betaU2')
    >>> len(reg.registered_auxevol_order())
    18

    >>> validate_strings(_body, "adm_to_evolved_block", file_ext="cpp")
    """
    if par.parval_from_str("Infrastructure") != "Dendro":
        raise ValueError(
            "Infrastructure must be 'Dendro' to build the ADM conversion, got "
            f"{par.parval_from_str('Infrastructure')!r}."
        )
    require_serial_parallelization()
    generation_parameters.validate_generation_parameters()
    scalar_type = str(par.parval_from_str("Dendro_scalar_type"))
    fp_type = str(par.parval_from_str("fp_type"))
    CoordSystem = str(par.parval_from_str("Dendro_fccz4_CoordSystem"))
    evol_order = reg.registered_evol_order()

    # Step 1: Map every evolved field the conversion defines to its expression,
    # reading each target name back from the registered BSSN gridfunction
    # symbols so no field name is spelled out.
    Bq = BSSN_quantities[CoordSystem]
    gammaDD, KDD, betaU, BU = register_adm_source_gridfunctions()
    adm = ADM_to_BSSN(
        gammaDD, KDD, betaU, BU, CoordSystem=CoordSystem, enable_rfm_precompute=False
    )
    # Insertion order is the emission order; Python dicts preserve it.
    targets: Dict[str, sp.Expr] = {}
    # The lapse is not part of the ADM source set this profile registers, so it
    # takes its registered asymptotic value (the same single authority the
    # Minkowski fill reads).
    targets[str(Bq.alpha)] = sp.sympify(gri.glb_gridfcs_dict[str(Bq.alpha)].f_infinity)
    targets[str(Bq.cf)] = adm.cf
    targets[str(Bq.trK)] = adm.trK
    for i in range(3):
        for j in range(i, 3):
            targets[str(Bq.hDD[i][j])] = adm.hDD[i][j]
    for i in range(3):
        for j in range(i, 3):
            targets[str(Bq.aDD[i][j])] = adm.aDD[i][j]
    for i in range(3):
        targets[str(Bq.vetU[i])] = adm.vetU[i]
        targets[str(Bq.betU[i])] = adm.betU[i]

    # Step 2: The one evolved field the conversion neither defines nor leaves
    # to the connection pass is the Z4 scalar, which constraint-satisfying
    # source data sets to zero (section 14.2).  It is identified by set
    # difference rather than by name.
    connection = {str(Bq.lambdaU[i]) for i in range(3)}
    residual = [
        name for name in evol_order if name not in targets and name not in connection
    ]
    if len(residual) != 1:
        raise ValueError(
            "The ADM conversion must leave exactly one evolved field (the Z4 "
            f"scalar) to be zeroed, found {residual}."
        )
    targets[residual[0]] = sp.sympify(0)
    missing = sorted(set(evol_order) - set(targets) - connection)
    if missing:
        raise ValueError(f"ADM conversion leaves evolved fields undefined: {missing}.")

    # Step 3: Lower the conversion inside an access capture.
    auxevol_order = reg.registered_auxevol_order()
    written = tuple(name for name in evol_order if name in targets)
    with cap.capture_gridfunction_accesses(ADM_TO_EVOLVED_BLOCK_CFUNCTION):
        kernel = c_codegen(
            [targets[name] for name in written],
            [f"{naming.out_pointer(name)}[pp]" for name in written],
            include_braces=False,
            enable_fd_codegen=True,
            enable_fd_functions=False,
            enable_simd=False,
            fp_type=fp_type,
            fp_type_alias=scalar_type,
            verbose=False,
        )
    accessed = set(cap.accessed_gridfunction_names(ADM_TO_EVOLVED_BLOCK_CFUNCTION))
    read_names = tuple(name for name in auxevol_order if name in accessed)
    unexpected = sorted(accessed - set(auxevol_order))
    if unexpected:
        raise ValueError(
            "The ADM conversion may read only registered ADM source fields, "
            f"but it also read {unexpected}."
        )
    bindings = "\n".join(
        [
            gfo.render_component_bindings(
                read_names,
                scalar_type,
                array="adm_gfs",
                role=naming.input_pointer,
                const_pointee=True,
                index_expression=lambda name, _p: str(auxevol_order.index(name)),
            ),
            gfo.render_component_bindings(
                written,
                scalar_type,
                array="out_gfs",
                role=naming.out_pointer,
                const_pointee=False,
                index_expression=lambda name, _p: str(evol_order.index(name)),
            ),
        ]
    )
    point_loop_body = interior_loop(
        kernel,
        nx="geom.nx",
        ny="geom.ny",
        nz="geom.nz",
        padding="0",
        pmin_padded="geom.pmin_padded",
        dx="geom.dx",
    )
    block_params = (
        f"const BlockGeometry& geom, const {scalar_type}* const* adm_gfs, "
        f"{scalar_type}* const* out_gfs"
    )
    return bindings + "\n" + point_loop_body, block_params


def build_lambda_initialization() -> Tuple[str, str]:
    """
    Build the separate connection-initialization CFunction body (section 14.3).

    ``lambdaU^i = DeltaGamma^i / ReU^i`` is exactly the statement that the
    connection constraint ``C^i = LambdatildeU^i - DeltaGamma^i`` vanishes, so
    the pass is initialization rather than a copy of a BSSN slot by position.
    It reads metric derivatives, so it runs over the interior and its input and
    output arrays are distinct.

    :return: (block_body, block_params).
    :raises ValueError: If Infrastructure is not Dendro or the pass is circular.

    Doctests:
    >>> import contextlib
    >>> import io
    >>> from nrpy.helpers.generic import validate_strings
    >>> from nrpy.equations.general_relativity.fCCZ4_system import (
    ...     build_fccz4_expression_bundle,
    ... )
    >>> par.set_parval_from_str("Infrastructure", "Dendro")
    >>> par.set_parval_from_str("parallelization", "none")
    >>> par.set_parval_from_str("fd_order", 4)
    >>> par.set_parval_from_str("EvolvedConformalFactor_cf", "chi")
    >>> with contextlib.redirect_stdout(io.StringIO()):
    ...     _bundle = build_fccz4_expression_bundle()
    >>> _body, _params = build_lambda_initialization()

    Section 14.3: the pass writes exactly the three connection components, and
    reads a distinct input array.

    >>> [
    ...     name
    ...     for name in reg.registered_evol_order()
    ...     if naming.out_pointer(name) + "[pp] =" in _body
    ... ]
    ['lambdaU0', 'lambdaU1', 'lambdaU2']
    >>> "const DendroScalar* const* in_gfs" in _params
    True

    >>> validate_strings(_body, "initialize_lambda_block", file_ext="cpp")
    """
    if par.parval_from_str("Infrastructure") != "Dendro":
        raise ValueError(
            "Infrastructure must be 'Dendro' to build the connection pass, got "
            f"{par.parval_from_str('Infrastructure')!r}."
        )
    require_serial_parallelization()
    scalar_type = str(par.parval_from_str("Dendro_scalar_type"))
    fp_type = str(par.parval_from_str("fp_type"))
    CoordSystem = str(par.parval_from_str("Dendro_fccz4_CoordSystem"))
    Bq = BSSN_quantities[CoordSystem]
    rfm = refmetric.reference_metric[CoordSystem]
    written = tuple(str(Bq.lambdaU[i]) for i in range(3))
    exprs = [Bq.DGammaU[i] / rfm.ReU[i] for i in range(3)]
    # The pass must not read what it writes: DeltaGamma^i is built from the
    # conformal metric and its first derivatives only.
    circular = sorted(
        {str(symbol) for expr in exprs for symbol in expr.free_symbols} & set(written)
    )
    if circular:
        raise ValueError(f"The connection pass reads the fields it writes: {circular}.")
    evol_order = reg.registered_evol_order()
    with cap.capture_gridfunction_accesses(INITIALIZE_LAMBDA_BLOCK_CFUNCTION):
        kernel = c_codegen(
            exprs,
            [f"{naming.out_pointer(name)}[pp]" for name in written],
            include_braces=False,
            enable_fd_codegen=True,
            enable_fd_functions=False,
            enable_simd=False,
            fp_type=fp_type,
            fp_type_alias=scalar_type,
            verbose=False,
        )
    accessed = set(cap.accessed_gridfunction_names(INITIALIZE_LAMBDA_BLOCK_CFUNCTION))
    read_names = tuple(name for name in evol_order if name in accessed)
    bindings = "\n".join(
        [
            gfo.render_component_bindings(
                read_names,
                scalar_type,
                array="in_gfs",
                role=naming.input_pointer,
                const_pointee=True,
                index_expression=lambda name, _p: str(evol_order.index(name)),
            ),
            gfo.render_component_bindings(
                written,
                scalar_type,
                array="out_gfs",
                role=naming.out_pointer,
                const_pointee=False,
                index_expression=lambda name, _p: str(evol_order.index(name)),
            ),
        ]
    )
    point_loop_body = interior_loop(
        kernel,
        nx="geom.nx",
        ny="geom.ny",
        nz="geom.nz",
        padding="geom.padding",
        pmin_padded="geom.pmin_padded",
        dx="geom.dx",
    )
    block_params = (
        f"const BlockGeometry& geom, const {scalar_type}* const* in_gfs, "
        f"{scalar_type}* const* out_gfs"
    )
    return bindings + "\n" + point_loop_body, block_params


def register_initial_data_conversion_CFunctions() -> None:
    """Register the ADM conversion and the connection-initialization CFunctions."""
    adm_body, adm_params = build_adm_to_evolved()
    reg.register_dendro_CFunction(
        role="adm_to_evolved_block",
        entry_point=True,
        lifecycle_hook="initial_data",
        name=ADM_TO_EVOLVED_BLOCK_CFUNCTION,
        desc=(
            "Per-block smooth ADM-to-fCCZ4 conversion (section 14.2); the "
            "connection components are written by the separate pass."
        ),
        subdirectory="generated/src/initial_data",
        params=adm_params,
        body=adm_body,
    )
    lam_body, lam_params = build_lambda_initialization()
    reg.register_dendro_CFunction(
        role="initialize_lambda_block",
        entry_point=True,
        lifecycle_hook="initial_data",
        name=INITIALIZE_LAMBDA_BLOCK_CFUNCTION,
        desc=(
            "Per-block connection initialization: lambdaU^i = DeltaGamma^i / "
            "ReU^i, so the connection constraint C^i vanishes (section 14.3)."
        ),
        subdirectory="generated/src/initial_data",
        params=lam_params,
        body=lam_body,
    )


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    print(f"Doctest passed: All {results.attempted} test(s) passed")
