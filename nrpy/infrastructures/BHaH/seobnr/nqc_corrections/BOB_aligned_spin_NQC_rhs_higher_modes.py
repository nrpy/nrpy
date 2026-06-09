"""
Register CFunction for BOB-informed NQC right hand sides.

Authors: Siddharth Mahesh
        sm0193 **at** mix **dot** wvu **dot** edu
        Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Union, cast

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.equations.seobnr.BOB_aligned_spin_waveform_quantities_higher_modes as BOB_HM_wf
import nrpy.helpers.parallel_codegen as pcg


def register_CFunction_BOB_aligned_spin_NQC_rhs_HM() -> Union[None, pcg.NRPyEnv_type]:
    """
    Register CFunction for calculating the NQC amplitudes and phase from BOB.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    wf = BOB_HM_wf.BOB_aligned_spin_waveform_quantities_higher_modes()
    # We are going to be doing this twice;
    # once for the fine dynamics and once for the coarse.
    modes = wf.modes
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = """
Calculates the BOB-informed Non Quasi-Circular (NQC) right-hand side terms.

@param commondata - Common data structure containing the model parameters.
@param amps - Array to store the calculated amplitudes.
@param omegas - Array to store the calculated angular frequencies.
"""
    cfunc_type = "void"
    name = "BOB_aligned_spin_NQC_rhs_HM"
    params = (
        "commondata_struct *restrict commondata , REAL (*amps)[3] , REAL (*omegas)[2]"
    )
    body = """
//allocate memory for hNR and wNR
REAL *hNR = malloc(NUMMODES * sizeof(REAL));
REAL *wNR = malloc(NUMMODES * sizeof(REAL));


if (hNR == NULL){
    fprintf(stderr, "Error: in BOB_aligned_spin_NQC_rhs_higher_modes(), malloc() failed for hNR\\n");
    exit(EXIT_FAILURE);
}
if (wNR == NULL){
    fprintf(stderr, "Error: in BOB_aligned_spin_NQC_rhs_higher_modes(), malloc() failed for wNR\\n");
    exit(EXIT_FAILURE);
}

//compute hNR and wNR at t_attach
SEOBNRv5_aligned_spin_hNR_fits_at_t_attach(commondata, hNR);
SEOBNRv5_aligned_spin_omegaNR_fits_at_t_attach(commondata, wNR);

const REAL m1 = commondata->m1;
const REAL m2 = commondata->m2;
const REAL M = m1 + m2;
const REAL nu = m1 * m2/ (M * M);

//compute
"""
    declare_HM_input_symbols = ""
    output_statement = ""
    nqc_rhs_quantities = []
    nqc_rhs_labels = []
    for lm in modes:
        # amplitude and derivatives at t_0
        declare_HM_input_symbols += f"const REAL {wf.modewise_input_symbols[lm][0]} = nu * hNR[HNR{lm[0]}{lm[1]}];\n"
        nqc_rhs_quantities.append(wf.hdot_t_attach_lm[lm])
        nqc_rhs_labels.append(f"const REAL hdot_{lm[0]}{lm[1]}_t_attach")
        nqc_rhs_quantities.append(wf.hddot_t_attach_lm[lm])
        nqc_rhs_labels.append(f"const REAL hddot_{lm[0]}{lm[1]}_t_attach")

        # frequency at t_0
        declare_HM_input_symbols += f"const REAL {wf.modewise_input_symbols[lm][1]} = fabs(wNR[HNR{lm[0]}{lm[1]}]);\n"
        nqc_rhs_quantities.append(wf.wdot_t_attach_lm[lm])
        nqc_rhs_labels.append(f"const REAL wdot_{lm[0]}{lm[1]}_t_attach")

        # qnm frequency
        declare_HM_input_symbols += f"const REAL {wf.modewise_input_symbols[lm][2]} = commondata->omega_qnm_l{lm[0]}m{lm[1]};\n"
        # qnm damping time
        declare_HM_input_symbols += f"const REAL {wf.modewise_input_symbols[lm][3]} = commondata->tau_qnm_l{lm[0]}m{lm[1]};\n"

        # outputs for this mode
        output_statement += f"""
amps[HNR{lm[0]}{lm[1]}][0] = {wf.modewise_input_symbols[lm][0]};
amps[HNR{lm[0]}{lm[1]}][1] = hdot_{lm[0]}{lm[1]}_t_attach;
amps[HNR{lm[0]}{lm[1]}][2] = hddot_{lm[0]}{lm[1]}_t_attach;
omegas[HNR{lm[0]}{lm[1]}][0] = fabs({wf.modewise_input_symbols[lm][1]});
omegas[HNR{lm[0]}{lm[1]}][1] = wdot_{lm[0]}{lm[1]}_t_attach;
"""
    BOB_code = ccg.c_codegen(
        nqc_rhs_quantities,
        nqc_rhs_labels,
        verbose=False,
        include_braces=False,
    )
    body += declare_HM_input_symbols
    body += BOB_code
    body += output_statement
    body += """
free(hNR);
free(wNR);
"""
    cfc.register_CFunction(
        subdirectory="nqc_corrections",
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
    )
    return pcg.NRPyEnv()
