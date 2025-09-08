"""
Set up C function library for the SEOBNR aligned spin expressions when using pre-computed waveform coefficients.
The waveform coefficients entering the factorized resummed waveform expressions are functions of the masses and spins.
The individual spins do not evolve when the spins are aligned and we can compute and store the flux coefficients
at the start of the code to accelerate the computation of fluxes and waveforms.


Authors: Siddharth Mahesh
        sm0193 **at** mix **dot** wvu **dot** edu
        Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Union, cast

import sympy as sp

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.equations.seobnr.SEOBNRv5_aligned_spin_precomputed_waveform_quantities as SEOBNRv5_precomp_wf
import nrpy.helpers.parallel_codegen as pcg
import nrpy.params as par


# Register waveform coefficients
def register_CFunction_SEOBNRv5_aligned_spin_waveform_coefficients() -> (
    Union[None, pcg.NRPyEnv_type]
):
    """
    Register CFunction for evaluating the SEOBNRv5 waveform coefficients.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    includes = ["BHaH_defines.h"]
    desc = """
Evaluates the SEOBNRv5 inspiral waveform coefficients.

@param commondata - Common data structure containing the model parameters."""
    cfunc_type = "void"
    name = "SEOBNRv5_aligned_spin_waveform_coefficients"
    params = "commondata_struct *restrict commondata"
    wf = SEOBNRv5_precomp_wf.SEOBNRv5_aligned_spin_waveform_quantities()

    # Compile a list of non-zero waveform coefficients from the
    # waveform quantities class. Since some of the variables are stored
    # as sympy arrays, we want to convert them to a list for passing
    # into nrpy.c_codegen. To this end, for every coefficient x, we define three objects:
    #   'x': The name of the variable as it appears in the commondata struct
    #   wf.x: The sympy expression for the variable
    #   'commondata->x': The 'label' of the variable when it is being assigned inside a function
    #                   with commondata as a pointer argument
    # Additionally, the waveform coefficients for odd-m terms
    # have different values depending on whether the mass ratio is equal.
    # For odd m's, the Newtonian coefficient c is proportional to delta = (m_1 - m_2),
    # while the PN spinning coefficient f_spin has some terms proportional to 1/delta
    # Thus, when the mass ratios are approximately equal:
    # h ~ lim_(delta->0) (c * f_spin)
    # => c = lim_(delta -> 0) (c / delta) ; f_spin = lim(delta -> 0) (f_spin * delta)
    # We will end up using a single if statement to store the values of the coefficients in either limit.

    waveform_coefficients = []
    waveform_coefficients_exprs = []
    waveform_coefficients_labels = []
    # first passage, don't assume
    for l in range(2, 9):
        for m in range(1, l + 1):
            if f"c{l + (l + m) % 2}" not in waveform_coefficients:
                waveform_coefficients.append(f"c{l + (l + m) % 2}")
                waveform_coefficients_labels.append(f"commondata->c{l + (l + m) % 2}")
                waveform_coefficients_exprs.append(wf.c[l + (l + m) % 2])
            waveform_coefficients.append(f"n_abs_{l}_{m}")
            waveform_coefficients_labels.append(f"commondata->n_abs_{l}_{m}")
            waveform_coefficients_exprs.append(sp.Abs(wf.n_complex[l, m]))
            if l == 2 and m == 2:
                waveform_coefficients.append(f"n_complex_{l}_{m}")
                waveform_coefficients_labels.append(f"commondata->n_complex_{l}_{m}")
                waveform_coefficients_exprs.append(wf.n_complex[l, m])
            for j in range(wf.max_vomega_order + 1):
                if (m % 2 and l < 5) or (l == 5 and m == 5):
                    if sp.sympify(wf.fspin[l, m, j]).is_nonzero is not False:
                        waveform_coefficients.append(f"fspin_{l}_{m}_{j}")
                        waveform_coefficients_labels.append(
                            f"commondata->fspin_{l}_{m}_{j}"
                        )
                        fspin_piecewise = (
                            wf.fspin[l, m, j] * wf.noneqcond
                            + wf.fspin_limit[l, m, j] * wf.eqcond
                        )
                        waveform_coefficients_exprs.append(fspin_piecewise)
                for k in range(m + 1):
                    if sp.sympify(wf.rho[l, m, j, k]).is_nonzero is not False:
                        waveform_coefficients.append(f"rho_{l}_{m}_{j}_{k}")
                        waveform_coefficients_labels.append(
                            f"commondata->rho_{l}_{m}_{j}_{k}"
                        )
                        waveform_coefficients_exprs.append(wf.rho[l, m, j, k])
                for k in range(wf.max_vh3_order + 1):
                    if (
                        l == 2
                        and m == 2
                        and sp.sympify(wf.deltalm[l, m, k, j]).is_nonzero is not False
                    ):
                        waveform_coefficients.append(f"deltalm_{l}_{m}_{k}_{j}")
                        waveform_coefficients_labels.append(
                            f"commondata->deltalm_{l}_{m}_{k}_{j}"
                        )
                        waveform_coefficients_exprs.append(wf.deltalm[l, m, k, j])

    # Register the coefficients as part of the commondata structure
    par.register_CodeParameters(
        "REAL",
        __name__,
        waveform_coefficients,
        commondata=True,
        add_to_parfile=False,
    )
    # Generate C code for the coefficients
    coeffs_ccode = ccg.c_codegen(
        waveform_coefficients_exprs,
        waveform_coefficients_labels,
        include_braces=False,
        verbose=False,
    )

    body = f"""
const REAL m1 = commondata->m1;
const REAL m2 = commondata->m2;
const REAL chi1 = commondata->chi1;
const REAL chi2 = commondata->chi2;
{coeffs_ccode}
"""
    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
    )
    return pcg.NRPyEnv()


def register_CFunction_SEOBNRv5_aligned_spin_waveform() -> (
    Union[None, pcg.NRPyEnv_type]
):
    """
    Register CFunction for calculating the (2,2) mode of the SEOBNRv5 inspiral when the waveform coefficients have been precomputed.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    wf = SEOBNRv5_precomp_wf.SEOBNRv5_aligned_spin_waveform_quantities()
    hlms = wf.strain()
    h22 = hlms["(2 , 2)"]
    # We are going to be doing this twice;
    # once for the fine dynamics and once for the coarse.
    h22_code = (
        ccg.c_codegen(
            h22,
            ["const double complex h22"],
            verbose=False,
            include_braces=False,
        )
        .replace("REAL", "double complex")
        .replace("exp", "cexp")
    )
    khat2_code = ccg.c_codegen(
        [wf.khat[2]],
        [
            "const REAL khat2",
        ],
        verbose=False,
        include_braces=False,
        cse_varprefix="khat",
    )

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = """
Calculates the SEOBNRv5 (2,2) mode for a single timestep.

@param dynamics - Array of dynamical variables.
@param commondata - Common data structure containing the model parameters.
@return - The (2,2) mode of the SEOBNRv5 inspiral.
"""
    cfunc_type = "double complex"
    prefunc = "#include<complex.h>"
    name = "SEOBNRv5_aligned_spin_waveform"
    params = "REAL *restrict dynamics, commondata_struct *restrict commondata"
    body = """
double complex gamma_22;
const REAL m1 = commondata->m1;
const REAL m2 = commondata->m2;
const REAL phi = dynamics[PHI];
const REAL Hreal = dynamics[H];
const REAL Omega = dynamics[OMEGA];
const REAL Omega_circ = dynamics[OMEGA_CIRC];
//compute
"""
    body += khat2_code
    body += """
  gamma_22 = SEOBNRv5_aligned_spin_gamma_wrapper(3.,-2.*khat2);
"""
    body += h22_code
    body += """
return h22;
"""
    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        prefunc=prefunc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
    )
    return pcg.NRPyEnv()


def register_CFunction_SEOBNRv5_aligned_spin_flux() -> Union[None, pcg.NRPyEnv_type]:
    """
    Register CFunction for evaluating the SEOBNRv5 aligned spin flux when the waveform coefficients have been precomputed.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    includes = ["BHaH_defines.h"]
    desc = """
Evaluates the SEOBNRv5 aligned spin flux.

@param y - Array of dynamical variables.
@param Hreal - The real EOB Hamiltonian.
@param Omega - The instantaneous angular frequency of the EOB perturber.
@param Omega_circ - The circular angular frequency of the EOB perturber.
@param f - Array to store the flux.
@param params - Pointer to the commondata structure containing the model parameters.
"""
    cfunc_type = "void"
    name = "SEOBNRv5_aligned_spin_flux"
    params = "const REAL *restrict y, const REAL Hreal, const REAL Omega, const REAL Omega_circ, REAL *restrict f, void *restrict params"
    wf = SEOBNRv5_precomp_wf.SEOBNRv5_aligned_spin_waveform_quantities()
    flux = wf.flux()
    body = """
commondata_struct *restrict commondata = (commondata_struct *restrict) params;
const REAL m1 = commondata->m1;
const REAL m2 = commondata->m2;
const REAL prstar = y[2];
const REAL pphi = y[3];
"""
    flux_ccode = ccg.c_codegen(
        [flux],
        [
            "const REAL flux",
        ],
        verbose=False,
        include_braces=False,
    )
    body += flux_ccode

    body += """
const REAL flux_over_omega = flux / Omega;
f[0] = prstar * flux_over_omega / pphi;
f[1] = flux_over_omega;
"""
    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
    )
    return pcg.NRPyEnv()
