"""
Python module for managing initial parameters for binary black hole (BBH) simulations.

License: Lesser GNU Public License, version 2.0+

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""
from collections import namedtuple  # Standard Python: Enable namedtuple data type
from typing import Dict

bbh_params = namedtuple(
    "bbh_params",
    "q d M_chix M_chiy M_chiz m_chix m_chiy m_chiz p_t p_r bare_mass_m_plus bare_mass_m_minus NA NB Nphi",
)
bbh_params_dict: Dict[str, bbh_params] = {}


def add_to_bbh_params_dict(
    name: str,
    q: float = 1.0,
    d: float = 2.0,
    M_chix: float = 0.0,
    M_chiy: float = 0.0,
    M_chiz: float = 0.0,
    m_chix: float = 0.0,
    m_chiy: float = 0.0,
    m_chiz: float = 0.0,
    p_t: float = 0.0,
    p_r: float = 0.0,
    bare_mass_m_plus: float = -1.0,
    bare_mass_m_minus: float = -1.0,
    NA: int = -1,
    NB: int = -1,
    Nphi: int = -1,
) -> None:
    """
    Add a new entry to the bbh_params_dict with the given parameters.

    :param name: The name identifier for the bbh parameters.
    :param q: Mass ratio.
    :param d: Separation.
    :param M_chix: Spin of the larger black hole in the x-direction.
    :param M_chiy: Spin of the larger black hole in the y-direction.
    :param M_chiz: Spin of the larger black hole in the z-direction.
    :param m_chix: Spin of the smaller black hole in the x-direction.
    :param m_chiy: Spin of the smaller black hole in the y-direction.
    :param m_chiz: Spin of the smaller black hole in the z-direction.
    :param p_t: Tangential momentum.
    :param p_r: Radial momentum.
    :param bare_mass_m_plus: Bare mass of the larger black hole.
    :param bare_mass_m_minus: Bare mass of the smaller black hole.
    :param NA: Number of grid points in A region.
    :param NB: Number of grid points in B region.
    :param Nphi: Number of azimuthal grid points.
    """
    bbh_params_dict[name] = bbh_params(
        q,
        d,
        M_chix,
        M_chiy,
        M_chiz,
        m_chix,
        m_chiy,
        m_chiz,
        p_t,
        p_r,
        bare_mass_m_plus,
        bare_mass_m_minus,
        NA,
        NB,
        Nphi,
    )


add_to_bbh_params_dict(
    "GW150914ET",
    q=36.0 / 29.0,
    d=10.0,
    M_chiy=+0.31,
    m_chiy=-0.46,
    p_t=0.09530152296974252,
    p_r=0.00084541526517121,
    NA=30,
    NB=30,
    Nphi=16,
)
add_to_bbh_params_dict(
    "BHB_q1_chi_m_0.3__chi_M_0_sep_11p8",
    q=1.0,
    d=5.9 * 2.0,
    m_chiy=0.075 / (0.5 * 0.5),  # J=0.075 is given, and chi = J/m^2, with m=0.5
    p_t=0.0852865,
    p_r=0.000515022,
    bare_mass_m_plus=0.48811120218500131,
    bare_mass_m_minus=0.4698442439908046,
    NA=64,
    NB=64,
    Nphi=44,
)
add_to_bbh_params_dict(
    "q1sep10",
    q=1.0,
    d=10.0,
    p_t=0.0962578089026658,
    p_r=0.00100787185295814,  # from NRPyPN
    NA=44,
    NB=44,
    Nphi=28,
)
add_to_bbh_params_dict(
    "q1sep8",
    q=1.0,
    d=8.0,
    p_t=0.112845235097096,
    p_r=0.00228434381143799,  # from NRPyPN
    NA=44,
    NB=44,
    Nphi=28,
)
add_to_bbh_params_dict(
    "q4sep11RITBBH0088",
    q=4.0,
    d=11.0,
    p_t=0.0578913282748551,
    p_r=0.000302008048634052,  # from NRPyPN
    NA=48,
    NB=48,
    Nphi=20,
)
add_to_bbh_params_dict(
    "q6sep13_SXSBBH0166",
    q=6.0,
    d=13.0,
    p_t=0.0396619891606381,
    p_r=0.000101374427664241,  # from NRPyPN
    NA=66,
    NB=66,
    Nphi=28,
)
add_to_bbh_params_dict(
    "q4sep12p5_SXSBBH0167",
    q=4.0,
    d=12.5,
    p_t=0.0531333569625095,
    p_r=0.000195987020545860,  # from NRPyPN
    NA=48,
    NB=48,
    Nphi=20,
)
add_to_bbh_params_dict(
    "q1sep12-1810.00036",
    # FROM TABLE V OF https://arxiv.org/pdf/1810.00036.pdf
    q=1.0,
    d=12.0,
    p_t=0.850686e-1,
    p_r=0.468113e-3,
    NA=44,
    NB=44,
    Nphi=24,
)
add_to_bbh_params_dict(
    "QC0_mclachlan",
    # From einsteinanalysis/Extract/test/qc0-mclachlan.par:
    q=1.0,
    d=1.168642873 * 2,
    p_t=0.3331917498,
    p_r=0.0,
    bare_mass_m_plus=0.453,
    bare_mass_m_minus=0.453,
    NA=24,
    NB=24,
    Nphi=18,
)
add_to_bbh_params_dict(
    "QC0",
    q=1.0,
    d=1.168642873 * 2,
    p_t=0.3331917498,
    p_r=0.0,
    bare_mass_m_plus=0.453,
    bare_mass_m_minus=0.453,
)
add_to_bbh_params_dict(
    "TP_BL",
    # TwoPunctures version of Brill-Lindquist
    # (trivial solve, but beware orbital separation=0 won't work)
    q=3.0,
    d=14.5,
    p_t=0.0,
    p_r=0.0,
    NA=10,
    NB=10,
    Nphi=6,
)

#
# def output_bbh_parfile_library():
#     # Loop through the bbh_params_dict items to populate the body string
#     for key, item in bbh_params_dict.items():
#         body += f""" else if(strcmp(TP_ID_type, "{key}")==0) {{
#         commondata->mass_ratio = {item.q};
#         commondata->initial_orbital_separation = {item.d};
#     """
#         if item.M_chix != 0.0:
#             body += f"    commondata->chi_BH_M[0] = {item.M_chix};\n"
#         if item.M_chiy != 0.0:
#             body += f"    commondata->chi_BH_M[1] = {item.M_chiy};\n"
#         if item.M_chiz != 0.0:
#             body += f"    commondata->chi_BH_M[2] = {item.M_chiz};\n"
#
#         if item.m_chix != 0.0:
#             body += f"    commondata->chi_BH_m[0] = {item.m_chix};\n"
#         if item.m_chiy != 0.0:
#             body += f"    commondata->chi_BH_m[1] = {item.m_chiy};\n"
#         if item.m_chiz != 0.0:
#             body += f"    commondata->chi_BH_m[2] = {item.m_chiz};\n"
#
#         body += f"    commondata->initial_p_t = {item.p_t};\n"
#         body += f"    commondata->initial_p_r = {item.p_r};\n"
#
#         if item.bare_mass_m_plus != -1.0 and item.bare_mass_m_minus != -1.0:
#             body += f"    par->give_bare_mass = true; //User provides bare masses rather than target ADM masses\n"
#             body += f"    par->par_m_plus = {item.bare_mass_m_plus};\n"
#             body += f"    par->par_m_minus = {item.bare_mass_m_minus};\n"
#
#         if item.NA > 0:
#             body += f"    par->npoints_A   = {item.NA};\n"
#         if item.NB > 0:
#             body += f"    par->npoints_B   = {item.NB};\n"
#         if item.Nphi > 0:
#             body += f"    par->npoints_phi = {item.Nphi};\n"
#
#         body += "  }"
#
#     body += r""" else {
#         fprintf(stderr, "Error: did not recognize TwoPunctures ID type = %s\n", TP_ID_type);
#         fprintf(stderr, "       Please pick one of the following library cases, or choose your own parameters:\n");
#         fprintf(stderr,"""
#
#     keys_str = ",".join([f"{key}\\n" for key in bbh_params_dict.keys()])
#     body += f'"{keys_str}\\n"'
#     body += r""");
#         exit(1);
#       }"""
