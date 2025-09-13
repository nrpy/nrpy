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
        subdirectory="inspiral_waveform_precomputed",
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
    )
    return pcg.NRPyEnv()
