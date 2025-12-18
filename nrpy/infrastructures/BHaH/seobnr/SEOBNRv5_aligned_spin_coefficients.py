"""
Set up C function library for SEOBNR initial conditions.

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
import nrpy.equations.seobnr.SEOBNRv5_aligned_spin_constants as SEOBNRv5_const
import nrpy.helpers.parallel_codegen as pcg
import nrpy.params as par


def register_CFunction_SEOBNRv5_aligned_spin_coefficients(
    calibration_no_spin: bool = False,
    calibration_spin: bool = False,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register CFunction for evaluating the masses and SEOBNRv5 coefficients.
    The inputs needed to generate an EOB approximant are mass ratio, spins and
    an initial orbital frequency. Therefore, one needs to compute the individual
    masses from the mass ratio and the Hamiltonian coeffients which are a function
    of mass ratio and spins. The Hamiltonian calibration coefficients can be pre-computed
    or added to the parfile for calibrating the SEOBNRv5 approximant.

    :param calibration_no_spin: If True, the non-spinning calibration coefficients are added to the parfile.
                                pySEOBNR v5 calibration coefficients are used if False.
    :param calibration_spin: If True, the spin-dependent calibration coefficients are added to the parfile.
                                pySEOBNR v5 calibration coefficients are used if False.
    :raises ValueError: If both calibration_no_spin and calibration_spin are True.
    :return: None if in registration phase, else the updated NRPy environment.
    """
    # The calibration process for the SEOBNRv5 is done in two steps:
    # 1. Calibration of the non-spinning coefficients
    # 2. Calibration of the spin-dependent coefficients
    # Therefore, the C code can only be generated for one of the above calibration options.
    if calibration_no_spin and calibration_spin:
        raise ValueError(
            "calibration_no_spin and calibration_spin cannot both be True."
        )
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    # Needed during integration and derived from other quantities; do not set!
    par.register_CodeParameters(
        "REAL",
        __name__,
        [
            "r",
            "phi",
            "m1",
            "m2",
            "prstar",
            "pphi",
            "t_stepback",
        ],
        [
            20,
            0,
            0.5,
            0.5,
            0.0,
            3.3,
            250.0,
        ],  # r, phi, m1, m2, prstar, pphi
        commondata=True,
        add_to_parfile=False,
    )

    par.register_CodeParameters(
        "REAL",
        __name__,
        [
            "Delta_t",
            "dT",
            "t_ISCO",
            "t_attach",
            "omega_qnm",
            "tau_qnm",
            "a_f",
            "M_f",
            "nr_amp_1",
            "nr_amp_2",
            "nr_amp_3",
            "nr_omega_1",
            "nr_omega_2",
            "a_1_NQC",
            "a_2_NQC",
            "a_3_NQC",
            "b_1_NQC",
            "b_2_NQC",
            "r_stop",
            "r_ISCO",
        ],
        [
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ],
        commondata=True,
        add_to_parfile=False,
    )

    par.register_CodeParameters(
        "REAL *restrict",
        __name__,
        [
            "dynamics_low",
            "dynamics_fine",
            "dynamics_raw",
        ],
        commondata=True,
        add_to_parfile=False,
        add_to_set_CodeParameters_h=False,
    )

    par.register_CodeParameters(
        "double complex *restrict",
        __name__,
        [
            "waveform_low",
            "waveform_fine",
            "waveform_inspiral",
            "waveform_IMR",
        ],
        commondata=True,
        add_to_parfile=False,
        add_to_set_CodeParameters_h=False,
    )

    par.register_CodeParameters(
        "size_t",
        __name__,
        ["nsteps_raw", "nsteps_low", "nsteps_fine", "nsteps_inspiral", "nsteps_IMR"],
        commondata=True,
        add_to_parfile=False,
        add_to_set_CodeParameters_h=False,
    )

    par.register_CodeParameters(
        "REAL",
        __name__,
        [
            "dHreal_dr",
            "dHreal_dprstar",
            "dHreal_dpphi",
            "dHreal_dr_dr",
            "dHreal_dr_dpphi",
            "dHreal_dr_circ",
            "dHreal_dpphi_circ",
            "Hreal",
            "xi",
            "flux",
            "Omega_circ",
        ],
        commondata=True,
        add_to_parfile=False,
    )

    # This is sufficient for initial conditions. Order is the same as pySEOBNR.
    par.register_CodeParameters(
        "REAL",
        __name__,
        ["mass_ratio", "chi1", "chi2", "initial_omega", "total_mass", "dt"],
        [
            1,
            0.4,
            -0.3,
            0.01118,
            50,
            2.4627455127717882e-05,
        ],  # mass_ratio convention is m_greater/m_lesser, initial_omega chosen for r ~ 20M
        commondata=True,
        add_to_parfile=True,
    )

    # register the calibration coefficients
    if calibration_no_spin:
        par.register_CodeParameters(
            "REAL",
            __name__,
            [
                "Delta_t_NS",
                "a6",
            ],
            [
                0.0,
                0.0,
            ],
            commondata=True,
            add_to_parfile=True,
        )
        par.register_CodeParameters(
            "REAL",
            __name__,
            [
                "Delta_t_S",
                "dSO",
            ],
            [
                0.0,
                0.0,
            ],
            commondata=True,
            add_to_parfile=False,
        )
    elif calibration_spin:
        par.register_CodeParameters(
            "REAL",
            __name__,
            [
                "Delta_t_NS",
                "a6",
            ],
            [
                0.0,
                0.0,
            ],
            commondata=True,
            add_to_parfile=False,
        )
        par.register_CodeParameters(
            "REAL",
            __name__,
            [
                "Deltat_t_S",
                "dSO",
            ],
            [
                0.0,
                0.0,
            ],
            commondata=True,
            add_to_parfile=True,
        )
    else:
        par.register_CodeParameters(
            "REAL",
            __name__,
            [
                "a6",
                "dSO",
            ],
            [
                0.0,
                0.0,
            ],
            commondata=True,
            add_to_parfile=False,
        )
    includes = ["BHaH_defines.h"]
    desc = """
Evaluate and store the SEOBNRv5 calibration coefficients and remnant properties.

@param commondata - The Common data structure containing the model parameters.
"""
    cfunc_type = "void"
    name = "SEOBNRv5_aligned_spin_coefficients"
    params = "commondata_struct *restrict commondata"
    v5_const = SEOBNRv5_const.SEOBNR_aligned_spin_constants(
        calibration_no_spin, calibration_spin
    )
    body = """
REAL q = commondata->mass_ratio;
REAL eta = q / (1.0 + q) / (1.0 + q);
if (eta > 0.25){
  if (fabs(q - 1.) < 1e-13){
    q = 1.;
  } else{
    printf("mass ratio = %.15e causes eta = %.15e > 0.25\\n",q,eta);
    exit(EXIT_FAILURE);
  }
}
commondata->m1 = q / (1.0 + q);
commondata->m2 = 1.0 / (1.0 + q);
const REAL m1 = commondata->m1;
const REAL m2 = commondata->m2;
const REAL chi1 = commondata->chi1;
const REAL chi2 = commondata->chi2;
commondata->dT = commondata->dt / commondata->total_mass / 4.925490947641266978197229498498379006e-6;
"""
    if calibration_no_spin:
        body += """
const REAL Delta_t_NS = commondata->Delta_t_NS;
"""
        body += ccg.c_codegen(
            [
                v5_const.Delta_t,
                v5_const.M_f,
                v5_const.a_f,
                v5_const.rISCO,
                v5_const.rstop,
            ],
            [
                "commondata->Delta_t",
                "commondata->M_f",
                "commondata->a_f",
                "commondata->r_ISCO",
                "commondata->r_stop",
            ],
            verbose=False,
            include_braces=False,
        )
    elif calibration_spin:
        body += """
const REAL Delta_t_S = commondata->Delta_t_S;
"""
        body += ccg.c_codegen(
            [
                v5_const.pyseobnr_a6,
                v5_const.Delta_t,
                v5_const.M_f,
                v5_const.a_f,
                v5_const.rISCO,
                v5_const.rstop,
            ],
            [
                "commondata->a6",
                "commondata->Delta_t",
                "commondata->M_f",
                "commondata->a_f",
                "commondata->r_ISCO",
                "commondata->r_stop",
            ],
            verbose=False,
            include_braces=False,
        )
    else:
        body += ccg.c_codegen(
            [
                v5_const.pyseobnr_a6,
                v5_const.pyseobnr_dSO,
                v5_const.Delta_t,
                v5_const.M_f,
                v5_const.a_f,
                v5_const.rISCO,
                v5_const.rstop,
            ],
            [
                "commondata->a6",
                "commondata->dSO",
                "commondata->Delta_t",
                "commondata->M_f",
                "commondata->a_f",
                "commondata->r_ISCO",
                "commondata->r_stop",
            ],
            verbose=False,
            include_braces=False,
        )
    body += """
const REAL afinallist[107] = { -0.9996, -0.9995, -0.9994, -0.9992, -0.999, -0.9989, -0.9988,
  -0.9987, -0.9986, -0.9985, -0.998, -0.9975, -0.997, -0.996, -0.995, -0.994, -0.992, -0.99, -0.988,
  -0.986, -0.984, -0.982, -0.98, -0.975, -0.97, -0.96, -0.95, -0.94, -0.92, -0.9, -0.88, -0.86, -0.84,
  -0.82, -0.8, -0.78, -0.76, -0.74, -0.72, -0.7, -0.65, -0.6, -0.55, -0.5, -0.45, -0.4, -0.35, -0.3,
  -0.25, -0.2, -0.15, -0.1, -0.05, 0., 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6,
  0.65, 0.7, 0.72, 0.74, 0.76, 0.78, 0.8, 0.82, 0.84, 0.86, 0.88, 0.9, 0.92, 0.94, 0.95, 0.96, 0.97,
  0.975, 0.98, 0.982, 0.984, 0.986, 0.988, 0.99, 0.992, 0.994, 0.995, 0.996, 0.997, 0.9975, 0.998,
  0.9985, 0.9986, 0.9987, 0.9988, 0.9989, 0.999, 0.9992, 0.9994, 0.9995, 0.9996
  };

const REAL reomegaqnm22[107] = { 0.2915755,0.291581,0.2915866,0.2915976,0.2916086,0.2916142,0.2916197,0.2916252,
  0.2916307,0.2916362,0.2916638,0.2916915,0.2917191,0.2917744,0.2918297,0.291885,0.2919958,0.2921067,0.2922178,0.2923289,
  0.2924403,0.2925517,0.2926633,0.292943,0.2932235,0.2937871,0.2943542,0.2949249,0.2960772,0.2972442,0.2984264,0.299624,
  0.3008375,0.3020672,0.3033134,0.3045767,0.3058573,0.3071558,0.3084726,0.3098081,0.3132321,0.316784,0.3204726,0.3243073,
  0.3282986,0.3324579,0.336798,0.3413329,0.3460786,0.3510526,0.3562748,0.3617677,0.3675569,0.3736717,0.3801456,0.3870175,
  0.394333,0.4021453,0.4105179,0.4195267,0.4292637,0.4398419,0.4514022,0.464123,0.4782352,0.4940448,0.5119692,0.5326002,
  0.5417937,0.5516303,0.5622007,0.5736164,0.586017,0.5995803,0.6145391,0.631206,0.6500179,0.6716143,0.6969947,0.7278753,
  0.74632,0.7676741,0.7932082,0.8082349,0.8254294,0.8331,0.8413426,0.8502722,0.8600456,0.8708927,0.8830905,0.8969183,0.9045305
  ,0.912655,0.9213264,0.9258781,0.9305797,0.9354355,0.9364255,0.937422,0.9384248,0.9394341,0.9404498,0.9425009,0.9445784,
  0.9456271,0.9466825 };

const REAL imomegaqnm22[107] = { 0.0880269,0.0880272,0.0880274,0.088028,0.0880285,0.0880288,0.088029,0.0880293,
  0.0880296,0.0880298,0.0880311,0.0880325,0.0880338,0.0880364,0.0880391,0.0880417,0.088047,0.0880523,0.0880575,0.0880628,0.088068,
  0.0880733,0.0880785,0.0880915,0.0881045,0.0881304,0.088156,0.0881813,0.0882315,0.0882807,0.0883289,0.0883763,0.0884226,0.0884679,
  0.0885122,0.0885555,0.0885976,0.0886386,0.0886785,0.0887172,0.0888085,0.0888917,0.0889663,0.0890315,0.0890868,0.0891313,0.0891643,
  0.0891846,0.0891911,0.0891825,0.0891574,0.0891138,0.0890496,0.0889623,0.0888489,0.0887057,0.0885283,0.0883112,0.0880477,0.0877293,
  0.0873453,0.086882,0.0863212,0.0856388,0.0848021,0.0837652,0.0824618,0.0807929,0.0799908,0.0790927,0.0780817,0.0769364,0.0756296,
  0.0741258,0.072378,0.0703215,0.0678642,0.0648692,0.0611186,0.0562313,0.053149,0.0494336,0.0447904,0.0419586,0.0386302,0.0371155,
  0.0354677,0.033659,0.0316517,0.0293904,0.0268082,0.0238377,0.0221857,0.0204114,0.0185063,0.0175021,0.016462,0.015385,0.0151651,
  0.0149437,0.0147207,0.0144962,0.0142701,0.0138132,0.0133501,0.0131161,0.0128806
  };

gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, 107);
if (spline == NULL){
  fprintf(stderr,"Error: in SEOBNRv5_aligned_spin_coefficients(), gsl_spline_alloc failed to initialize\\n");
  exit(1);
}
gsl_interp_accel *acc = gsl_interp_accel_alloc();
if (acc == NULL){
  fprintf(stderr,"Error: in SEOBNRv5_aligned_spin_coefficients(), gsl_interp_acc_alloc failed to initialize\\n");
  exit(1);
}

gsl_spline_init(spline, afinallist, reomegaqnm22, 107);
commondata->omega_qnm = gsl_spline_eval(spline, commondata->a_f, acc) / commondata->M_f;

gsl_spline_init(spline, afinallist, imomegaqnm22, 107);
gsl_interp_accel_reset(acc);
commondata->tau_qnm = 1./(gsl_spline_eval(spline, commondata->a_f, acc) / commondata->M_f) ;

gsl_spline_free(spline);
gsl_interp_accel_free(acc);
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
