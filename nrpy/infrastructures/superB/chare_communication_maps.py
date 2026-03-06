"""
C function for setting up chare communcation struct for maps between indices on the local chare grid and the global grid, converting idx3 of grid point to the idx3 of the chare evolving it, etc...

Authors: Nishita Jadoo
         njadoo **at** uidaho **dot* edu
"""

# Step P1: Import needed NRPy core modules:
from typing import Set

import nrpy.c_function as cfc
from nrpy.infrastructures import BHaH


def register_CFunction_charecommstruct_set_up(CoordSystem: str) -> None:
    """
    Register C function for setting up charecommstruct.

    :param CoordSystem: Coordinate system
    """
    includes = [
        "BHaH_defines.h",
        "BHaH_function_prototypes.h",
    ]
    desc = r"""Setup charecommstruct"""
    cfunc_type = "void"
    name = "charecommstruct_set_up"
    params = "const commondata_struct *restrict commondata, const params_struct *restrict params, const params_struct *restrict params_chare, charecomm_struct* charecommstruct, const int thischareindex[3]"
    body = r"""
  const int Nxx_plus_2NGHOSTS0chare = params_chare->Nxx_plus_2NGHOSTS0;
  const int Nxx_plus_2NGHOSTS1chare = params_chare->Nxx_plus_2NGHOSTS1;
  const int Nxx_plus_2NGHOSTS2chare = params_chare->Nxx_plus_2NGHOSTS2;
  const int Nxx_plus_2NGHOSTS0 = params->Nxx_plus_2NGHOSTS0;
  const int Nxx_plus_2NGHOSTS1 = params->Nxx_plus_2NGHOSTS1;
  const int Nxx0chare = params_chare->Nxx0;
  const int Nxx1chare = params_chare->Nxx1;
  const int Nxx2chare = params_chare->Nxx2;

  // Global->owner/local maps are now computed on demand from index arithmetic.
  charecommstruct->globalidx3pt_to_chareidx3 = NULL;
  charecommstruct->globalidx3pt_to_localidx3pt = NULL;

  // local to this chare
  const int ntotchare = Nxx_plus_2NGHOSTS0chare * Nxx_plus_2NGHOSTS1chare * Nxx_plus_2NGHOSTS2chare;
  charecommstruct->localidx3pt_to_globalidx3pt = (int *restrict)malloc(sizeof(int) * ntotchare);
  const int startlocali = 0;
  const int startlocalj = 0;
  const int startlocalk = 0;
  const int endlocali = Nxx_plus_2NGHOSTS0chare - 1;
  const int endlocalj = Nxx_plus_2NGHOSTS1chare - 1;
  const int endlocalk = Nxx_plus_2NGHOSTS2chare - 1;
  for (int localk = startlocalk; localk <= endlocalk; localk++) {
    for (int localj = startlocalj; localj <= endlocalj; localj++) {
      for (int locali = startlocali; locali <= endlocali; locali++) {
        int localidx3 = IDX3GENERAL(locali, localj, localk, Nxx_plus_2NGHOSTS0chare, Nxx_plus_2NGHOSTS1chare);
        int globali = MAP_LOCAL_TO_GLOBAL_IDX0(thischareindex[0], locali, Nxx0chare);
        int globalj = MAP_LOCAL_TO_GLOBAL_IDX1(thischareindex[1], localj, Nxx1chare);
        int globalk = MAP_LOCAL_TO_GLOBAL_IDX2(thischareindex[2], localk, Nxx2chare);
        int globalidx3 = IDX3GENERAL(globali, globalj, globalk, Nxx_plus_2NGHOSTS0, Nxx_plus_2NGHOSTS1);
        charecommstruct->localidx3pt_to_globalidx3pt[localidx3] = globalidx3;
      }
    }
  }
"""

    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        CoordSystem_for_wrapper_func=CoordSystem,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
    )


def chare_comm_register_C_functions(
    set_of_CoordSystems: Set[str],
) -> None:
    """
    Register C functions for chare communication.

    :param set_of_CoordSystems: Set of coordinate systems to use.
    :return None
    """
    for CoordSystem in set_of_CoordSystems:
        # Register C function to set up the chare communication struct.
        register_CFunction_charecommstruct_set_up(CoordSystem=CoordSystem)

    # Register charecomm_struct's contribution to griddata_struct:
    BHaH.griddata_commondata.register_griddata_commondata(
        __name__,
        "charecomm_struct charecommstruct",
        "maps that convert between index of a pt in chare's local grid to the index on the global grid, etc",
    )


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
