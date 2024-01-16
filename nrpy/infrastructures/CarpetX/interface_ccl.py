"""
Module for constructing interface.ccl for Cactus thorns.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""
from typing import cast, Iterator, Tuple
from pathlib import Path
import nrpy.grid as gri
from nrpy.helpers.safewrite import SafeWrite

def carpetx_gfs()->Iterator[Tuple[str,gri.CarpetXGridFunction]]:
    for gfname, gf in gri.glb_gridfcs_dict.items():
        assert type(gf) == gri.CarpetXGridFunction
        yield (gfname, gf)

def construct_interface_ccl(
    project_dir: str,
    thorn_name: str,
    inherits: str,
    USES_INCLUDEs: str,
    is_evol_thorn: bool = False,
    enable_NewRad: bool = False,
) -> None:
    """
    Generate `interface.ccl` file required for the specified Thorn.

    :param thorn_name: The name of the thorn for which the interface.ccl file is generated. Defaults to "BaikalETK".
    :param enable_stress_energy_source_terms: Boolean flag to determine whether to include stress-energy source terms. Defaults to False.
    :return: None
    """
    outstr = rf"""
# This interface.ccl file was automatically generated by NRPy+.
#   You are advised against modifying it directly; instead
#   modify the Python code that generates it.

# With "implements", we give our thorn its unique name.
implements: {thorn_name}

# By "inheriting" other thorns, we tell the Toolkit that we
#   will rely on variables/function that exist within those
#   functions.
inherits: {inherits}

# Needed functions and #include's:
{USES_INCLUDEs}
"""
    if enable_NewRad:
        outstr += r"""
# Note: we don't have NewRad yet
# Needed for NewRad outer boundary condition driver:
#CCTK_INT FUNCTION                         \
#    NewRad_Apply                          \
#        (CCTK_POINTER_TO_CONST IN cctkGH, \
#         CCTK_REAL ARRAY IN var,          \
#         CCTK_REAL ARRAY INOUT rhs,       \
#         CCTK_REAL IN var0,               \
#         CCTK_REAL IN v0,                 \
#         CCTK_INT IN radpower)
#REQUIRES FUNCTION NewRad_Apply
"""
    outstr += """
# Tell the Toolkit that we want all gridfunctions
#    to be visible to other thorns by using
#    the keyword "public". Note that declaring these
#    gridfunctions *does not* allocate memory for them;
#    that is done by the schedule.ccl file.

public:
"""
    if is_evol_thorn:
        evol_parities = ""
        evol_gfs = []
        # First, EVOL type:
        for gfname, gf in carpetx_gfs():
            if gf.group == "EVOL":
                evol_parities += f"{gf.parity}  "
                evol_gfs += [f"{gfname}GF"]
        if evol_gfs:
            outstr += f"CCTK_REAL evol_variables type = GF Timelevels=1 TAGS=\'rhs=\"evol_variables_rhs\" parities={{{evol_parities}}}\'\n{{\n  "
            outstr += ", ".join(evol_gfs) + "\n"
            outstr += """} "Evolved gridfunctions."

"""

            # Second EVOL right-hand-sides
            outstr += 'CCTK_REAL evol_variables_rhs type = GF Timelevels=1 TAGS=\'InterpNumTimelevels=1 prolongation="none" checkpoint="no"\'\n{\n  '
            rhs_gfs = [
                f"{gfname}_rhsGF"
                for gfname, gf in carpetx_gfs() 
                if gf.group == "EVOL"
            ]
            outstr += ", ".join(rhs_gfs) + "\n"
            outstr += """} "Right-hand-side gridfunctions."

"""
            # Then AUXEVOL type:
            auxevol_parities = ""
            auxevol_gfs = []
            for gfname, gf in carpetx_gfs():
                if gf.group == "AUXEVOL":
                    auxevol_parities += f"{gf.parity}  "
                    auxevol_gfs += [f"{gfname}GF"]
            if auxevol_gfs:
                outstr += f'CCTK_REAL auxevol_variables type = GF Timelevels=1 TAGS=\'InterpNumTimelevels=1 prolongation="none" checkpoint="no" parities={{{auxevol_parities}}}\'\n{{\n  '
                outstr += ", ".join(auxevol_gfs) + "\n"
                outstr += """} "Auxiliary gridfunctions needed for evaluating the RHSs."

"""
        # Then AUX type:
        aux_parities = ""
        aux_gfs = []
        for gfname, gf in carpetx_gfs():
            if gf.group == "AUX":
                aux_parities += f"{gf.parity}  "
                aux_gfs += [f"{gfname}GF"]
        if aux_gfs:
            outstr += f"CCTK_REAL aux_variables type = GF Timelevels=1 TAGS=\'parities={{{aux_parities}}}\'\n{{\n  "
            outstr += ", ".join(aux_gfs) + "\n"
            outstr += """} "Auxiliary gridfunctions for e.g., diagnostics."

"""
    output_Path = Path(project_dir) / thorn_name
    output_Path.mkdir(parents=True, exist_ok=True)
    with SafeWrite(output_Path / "interface.ccl", encoding="utf-8") as file:
        file.write(outstr)
