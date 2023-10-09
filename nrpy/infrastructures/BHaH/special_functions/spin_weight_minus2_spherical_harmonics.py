"""
Register C functions for computing the spin-weight -2 spherical harmonics.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot** com
"""
from typing import cast, Union
from inspect import currentframe as cf
from types import FrameType as FT
import sympy as sp

import nrpy.params as par
import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.helpers.parallel_codegen as pcg
import nrpy.equations.special_functions.spin_weighted_spherical_harmonics as SWSH

par.register_param(
    int,
    __name__,
    "swm2sh_maximum_l_mode_generated",
    2,
)
par.register_CodeParameter(
    "int", __name__, "swm2sh_maximum_l_mode_to_compute", 2, commondata=True
)


def register_CFunction_spin_weight_minus2_sph_harmonics() -> (
    Union[None, pcg.NRPyEnv_type]
):
    """
    Register C function for computing arbitrary spin-weight -2 spherical harmonics.
    These spherical harmonics are computed at any (l, m) pair at a specific point (theta, phi).
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cf()).f_code.co_name}", locals())
        return None

    includes = ["BHaH_defines.h"]
    # Set up the C function for computing the spin-weight -2 spherical harmonic at theta,phi: Y_{s=-2, l,m}(theta,phi)
    desc = r"""// Compute at a single point (th,ph) the spin-weight -2 spherical harmonic Y_{s=-2, l,m}(th,ph)
// Manual "inline void" of this function results in compilation error with clang.
"""
    c_type = "void"
    name = "spin_weight_minus2_sph_harmonics"
    params = "const int l, const int m, const REAL th, const REAL ph, REAL *restrict reYlmswm2_l_m, REAL *restrict imYlmswm2_l_m"
    # Construct body:
    # real=True for th, ph is ESSENTIAL, so that sp.re and sp.im function properly below.
    th, ph = sp.symbols("th ph", real=True)
    body = "switch(l) {\n"

    # maximum_l: The maximum value of the angular momentum number 'l' up to which
    #                  the function will compute the spin-weight -2 spherical harmonics.
    swm2sh_maximum_l_mode_generated = par.parval_from_str(
        "swm2sh_maximum_l_mode_generated"
    )
    for l in range(
        swm2sh_maximum_l_mode_generated + 1
    ):  # Output values up to and including l=8.
        body += f"  case {l}:\n"
        body += "    switch(m) {\n"
        for m in range(-l, l + 1):
            body += f"    case {m}:\n"
            body += "       {\n"
            body += f"      // l = {l}, m = {m}:\n"
            Y_m2_lm = SWSH.Y(-2, l, m, th, ph)
            body += ccg.c_codegen(
                [sp.re(Y_m2_lm), sp.im(Y_m2_lm)],
                ["*reYlmswm2_l_m", "*imYlmswm2_l_m"],
                verbose=False,
                include_braces=False,
            )
            body += "       }\n"
            body += "      return;\n"
        body += f"    }}  // END switch(l == {l})\n"
    body += "  } // END switch blocks\n"
    body += rf"""
  fprintf(stderr, "ERROR: SpinWeight_minus2_SphHarmonics handles only l=[0,swm2sh_maximum_l_mode_generated={swm2sh_maximum_l_mode_generated}] and only m=[-l,+l] is defined.\n");
  fprintf(stderr, "       You chose l=%d and m=%d, which is out of these bounds.\n",l,m);
  exit(1);
"""
    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        c_type=c_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
    )
    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())
