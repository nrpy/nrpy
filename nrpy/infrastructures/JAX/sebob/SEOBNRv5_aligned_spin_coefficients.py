"""
Set up JAX function library for SEOBNR initial conditions.

Authors: Siddharth Mahesh
        sm0193 **at** mix **dot** wvu **dot** edu
        Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Union, cast

import nrpy.equations.seobnr.SEOBNRv5_aligned_spin_constants as SEOBNRv5_const
import nrpy.helpers.parallel_codegen as pcg
import nrpy.py_codegen as pycg
import nrpy.py_function as pyfc


def register_PyFunction_SEOBNRv5_aligned_spin_coefficients() -> (
    Union[None, pcg.NRPyEnv_type]
):
    """
    Register JAX function for evaluating the masses and SEOBNRv5 coefficients.
    The inputs needed to generate an EOB approximant are mass ratio, spins and
    an initial orbital frequency. Therefore, one needs to compute the individual
    masses from the mass ratio and the Hamiltonian coeffients which are a function
    of mass ratio and spins.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None
    imports = [
        "import jax",
        "import jax.numpy as jnp",
    ]
    desc = """
Evaluate and store the SEOBNRv5 calibration coefficients and remnant properties.

:param mass_ratio: The mass ratio of the system.
:param chi1: The dimensionless spin of the first object.
:param chi2: The dimensionless spin of the second object.
:param initial_omega: The initial orbital frequency of the system.
:param dt: The time step of the system.
:param total_mass: The total mass of the system.
:return: A commondata object with SEOBNRv5 calibration coefficients and remnant properties.
"""
    name = "SEOBNRv5_aligned_spin_coefficients"
    params = "mass_ratio, chi1, chi2, initial_omega, dt, total_mass"
    v5_const = SEOBNRv5_const.SEOBNR_aligned_spin_constants()
    body = """
q = mass_ratio
eta = q / (1.0 + q) / (1.0 + q)
if (eta > 0.25):
    if (jnp.abs(q - 1.) < 1e-13):
        q = 1.
    else:
        raise ValueError(f"mass ratio = {q} causes eta = {eta} > 0.25")

m1 = q / (1.0 + q)
m2 = 1.0 / (1.0 + q)
dT = dt / total_mass / 4.925490947641266978197229498498379006e-6
"""
    body += pycg.py_codegen(
        [
            v5_const.pyseobnr_a6,
            v5_const.pyseobnr_dSO,
            v5_const.Delta_t,
            v5_const.M_f,
            v5_const.a_f,
            v5_const.rISCO,
        ],
        [
            "a6",
            "dSO",
            "Delta_t",
            "M_f",
            "a_f",
            "r_ISCO",
        ],
        verbose=False,
    )
    body += """
rstop = 0.98*rISCO if Delta_t > 0 else -1.
afinallist = jnp.array([
    -0.9996, -0.9995, -0.9994, -0.9992, -0.999, -0.9989, -0.9988,
    -0.9987, -0.9986, -0.9985, -0.998, -0.9975, -0.997, -0.996, -0.995, -0.994, -0.992, -0.99, -0.988,
    -0.986, -0.984, -0.982, -0.98, -0.975, -0.97, -0.96, -0.95, -0.94, -0.92, -0.9, -0.88, -0.86, -0.84,
    -0.82, -0.8, -0.78, -0.76, -0.74, -0.72, -0.7, -0.65, -0.6, -0.55, -0.5, -0.45, -0.4, -0.35, -0.3,
    -0.25, -0.2, -0.15, -0.1, -0.05, 0., 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6,
    0.65, 0.7, 0.72, 0.74, 0.76, 0.78, 0.8, 0.82, 0.84, 0.86, 0.88, 0.9, 0.92, 0.94, 0.95, 0.96, 0.97,
    0.975, 0.98, 0.982, 0.984, 0.986, 0.988, 0.99, 0.992, 0.994, 0.995, 0.996, 0.997, 0.9975, 0.998,
    0.9985, 0.9986, 0.9987, 0.9988, 0.9989, 0.999, 0.9992, 0.9994, 0.9995, 0.9996
])

reomegaqnm22 = jnp.array([
    0.2915755,0.291581,0.2915866,0.2915976,0.2916086,0.2916142,0.2916197,0.2916252,
    0.2916307,0.2916362,0.2916638,0.2916915,0.2917191,0.2917744,0.2918297,0.291885,0.2919958,0.2921067,0.2922178,0.2923289,
    0.2924403,0.2925517,0.2926633,0.292943,0.2932235,0.2937871,0.2943542,0.2949249,0.2960772,0.2972442,0.2984264,0.299624,
    0.3008375,0.3020672,0.3033134,0.3045767,0.3058573,0.3071558,0.3084726,0.3098081,0.3132321,0.316784,0.3204726,0.3243073,
    0.3282986,0.3324579,0.336798,0.3413329,0.3460786,0.3510526,0.3562748,0.3617677,0.3675569,0.3736717,0.3801456,0.3870175,
    0.394333,0.4021453,0.4105179,0.4195267,0.4292637,0.4398419,0.4514022,0.464123,0.4782352,0.4940448,0.5119692,0.5326002,
    0.5417937,0.5516303,0.5622007,0.5736164,0.586017,0.5995803,0.6145391,0.631206,0.6500179,0.6716143,0.6969947,0.7278753,
    0.74632,0.7676741,0.7932082,0.8082349,0.8254294,0.8331,0.8413426,0.8502722,0.8600456,0.8708927,0.8830905,0.8969183,0.9045305,
    0.912655,0.9213264,0.9258781,0.9305797,0.9354355,0.9364255,0.937422,0.9384248,0.9394341,0.9404498,0.9425009,0.9445784,
    0.9456271,0.9466825
])

imomegaqnm22 = jnp.array([ 
    0.0880269,0.0880272,0.0880274,0.088028,0.0880285,0.0880288,0.088029,0.0880293,
    0.0880296,0.0880298,0.0880311,0.0880325,0.0880338,0.0880364,0.0880391,0.0880417,0.088047,0.0880523,0.0880575,0.0880628,0.088068,
    0.0880733,0.0880785,0.0880915,0.0881045,0.0881304,0.088156,0.0881813,0.0882315,0.0882807,0.0883289,0.0883763,0.0884226,0.0884679,
    0.0885122,0.0885555,0.0885976,0.0886386,0.0886785,0.0887172,0.0888085,0.0888917,0.0889663,0.0890315,0.0890868,0.0891313,0.0891643,
    0.0891846,0.0891911,0.0891825,0.0891574,0.0891138,0.0890496,0.0889623,0.0888489,0.0887057,0.0885283,0.0883112,0.0880477,0.0877293,
    0.0873453,0.086882,0.0863212,0.0856388,0.0848021,0.0837652,0.0824618,0.0807929,0.0799908,0.0790927,0.0780817,0.0769364,0.0756296,
    0.0741258,0.072378,0.0703215,0.0678642,0.0648692,0.0611186,0.0562313,0.053149,0.0494336,0.0447904,0.0419586,0.0386302,0.0371155,
    0.0354677,0.033659,0.0316517,0.0293904,0.0268082,0.0238377,0.0221857,0.0204114,0.0185063,0.0175021,0.016462,0.015385,0.0151651,
    0.0149437,0.0147207,0.0144962,0.0142701,0.0138132,0.0133501,0.0131161,0.0128806
])

omega_qnm = jnp.interp(a_f, afinallist, reomegaqnm22)/M_f
tau_qnm = M_f/jnp.interp(a_f, afinallist, imomegaqnm22)
return CommonData(
    m1=m1,
    m2=m2,
    chi1=chi1,
    chi2=chi2,
    initial_omega=initial_omega,
    dT=dT,
    a6=a6,
    dSO=dSO,
    rISCO=rISCO,
    rstop=rstop,
    omega_qnm=omega_qnm,
    tau_qnm=tau_qnm,
    M_f=M_f,
    a_f=a_f,
)
"""

    pyfc.register_PyFunction(
        imports=imports,
        desc=desc,
        name=name,
        params=params,
        body=body,
    )
    return pcg.NRPyEnv()
