"""
C function to save data from GRHayL structs into NRPy gridfunctions, and do basis transforms when needed.

Author: Terrence Pierre Jacques
        terrencepierrej **at** gmail **dot* com
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Union, cast

import sympy as sp

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.helpers.jacobians as jcb
import nrpy.helpers.parallel_codegen as pcg
import nrpy.indexedexp as ixp
import nrpy.reference_metric as refmetric


def register_CFunction_basis_transform_Cartesian_to_rfm_basis(
    CoordSystem: str,
    enable_GoldenKernels: bool = False,
    evolving_temperature: bool = False,
    evolving_entropy: bool = False,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register function for basis transforming conserved momentum and velocities from Cartesian coordinates to rfm basis, and fill in structs for GRHayL.

    :param CoordSystem: The coordinate system (e.g., "Cartesian", "Spherical").
    :param enable_GoldenKernels: Boolean to enable Golden Kernels.
    :param evolving_temperature: whether we're using grhayl to evolve temperature or not.
    :param evolving_entropy: whether we're using grhayl to evolve entropy or not.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    # Step 1: Register the function with NRPy's parallel codegen infrastructure
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    # Step 2: Set up C function headers, signature, and parameters
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    cfunc_type = "void"
    desc = (
        """
    Convert fluid velocity and Stilde from Cartesian coordinates to """
        + CoordSystem
        + """coordinates."""
    )
    name = "basis_transform_Cartesian_to_rfm_basis"
    params = "const commondata_struct *restrict commondata, const params_struct *restrict params, const ghl_primitive_quantities *restrict prims, "
    params += "const ghl_conservative_quantities *restrict cons, "
    params += "const int i0, const int i1, const int i2, REAL *restrict xx[3], "
    params += "REAL *restrict auxevol_gfs, REAL *restrict evol_gfs"

    # Step 3: Define reference metric and symbolic variables
    # GRHayL typically works with Cartesian components,
    # so we define symbols for Cartesian velocity (v) and momentum (S).
    rfm = refmetric.reference_metric[CoordSystem]

    vx, vy, vz = sp.symbols("vx vy vz")
    S_x, S_y, S_z = sp.symbols("S_x S_y S_z")
    vU_Cart = [vx, vy, vz]
    StildeD_Cart = [S_x, S_y, S_z]

    # Step 4: Perform Symbolic Basis Transformations
    # Convert vectors from Cartesian basis to the specific reference metric basis (e.g., Spherical)
    # using Jacobian matrices.
    vU = jcb.basis_transform_vectorU_from_Cartesian_to_rfmbasis(CoordSystem, vU_Cart)

    StildeD = jcb.basis_transform_vectorD_from_Cartesian_to_rfmbasis(
        CoordSystem, StildeD_Cart
    )

    # Step 5: Rescale vectors
    # NRPy+ typically rescales vectors by reference metric scale factors (rfm.ReU)
    # to handle coordinate singularities smoothly.
    # v^i_rescaled = v^i / ReU[i]
    # S_i_rescaled = S_i * ReU[i] (Covariant vectors scale inversely to Contravariant)
    rescaledvU = ixp.zerorank1()
    rescaledStildeD = ixp.zerorank1()
    for i in range(3):
        rescaledvU[i] = vU[i] / rfm.ReU[i]
        rescaledStildeD[i] = StildeD[i] * rfm.ReU[i]

    # Step 6: Generate C code to read coordinates and copy Scalar quantities
    # Scalars (Density, Pressure, Entropy, etc.) do not require basis transforms.
    read_rfm_xx_arrays = r"""
    const REAL xx0 = xx[0][i0];
    const REAL xx1 = xx[1][i1];
    const REAL xx2 = xx[2][i2];
    """

    pre_body = read_rfm_xx_arrays + r"""

    auxevol_gfs[IDX4(RHOBGF, i0, i1, i2)] = prims->rho;
    auxevol_gfs[IDX4(PGF, i0, i1, i2)] = prims->press;
    """

    if evolving_temperature:
        pre_body += r"""
    auxevol_gfs[IDX4(YEGF, i0, i1, i2)] = prims->Y_e;
    auxevol_gfs[IDX4(TEMPERATUREGF, i0, i1, i2)] = prims->temperature;
    """

    if evolving_entropy:
        pre_body += r"""
    auxevol_gfs[IDX4(SGF, i0, i1, i2)] = prims->entropy;
    """

    pre_body += r"""
    evol_gfs[IDX4(RHO_STARGF, i0, i1, i2)] = cons->rho;
    evol_gfs[IDX4(TAU_TILDEGF, i0, i1, i2)] = cons->tau;
    """

    if evolving_temperature:
        pre_body += r"""
    evol_gfs[IDX4(YE_STARGF, i0, i1, i2)] = cons->Y_e;
"""

    if evolving_entropy:
        pre_body += r"""
    evol_gfs[IDX4(S_STARGF, i0, i1, i2)] = cons->entropy;
"""

    # Read Cartesian components from GRHayL structs into local C variables
    pre_body += r"""

    const REAL vx = prims->vU[0];
    const REAL vy = prims->vU[1];
    const REAL vz = prims->vU[2];
    const REAL S_x = cons->SD[0];
    const REAL S_y = cons->SD[1];
    const REAL S_z = cons->SD[2];

"""

    # Step 7: Generate C code for vector quantities
    # Map the symbolically transformed and rescaled variables to the output grid functions.
    lhs = [
        "auxevol_gfs[IDX4(RESCALEDVU0GF, i0, i1, i2)]",
        "auxevol_gfs[IDX4(RESCALEDVU1GF, i0, i1, i2)]",
        "auxevol_gfs[IDX4(RESCALEDVU2GF, i0, i1, i2)]",
        "evol_gfs[IDX4(RESCALEDSTILDED0GF, i0, i1, i2)]",
        "evol_gfs[IDX4(RESCALEDSTILDED1GF, i0, i1, i2)]",
        "evol_gfs[IDX4(RESCALEDSTILDED2GF, i0, i1, i2)]",
    ]
    rhs = [
        rescaledvU[0],
        rescaledvU[1],
        rescaledvU[2],
        rescaledStildeD[0],
        rescaledStildeD[1],
        rescaledStildeD[2],
    ]

    expr_body = ccg.c_codegen(
        rhs,
        lhs,
        enable_simd=False,
        enable_fd_codegen=True,
        enable_GoldenKernels=enable_GoldenKernels,
    )

    # Step 8: Register the final C function
    cfc.register_CFunction(
        include_CodeParameters_h=True,
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        CoordSystem_for_wrapper_func=CoordSystem,
        name=name,
        params=params,
        body=pre_body + expr_body,
        enable_simd=False,
    )
    return pcg.NRPyEnv()
