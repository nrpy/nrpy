"""

Authors: Nishita Jadoo
         njadoo **at** uidaho **dot* edu
"""

# Step P1: Import needed NRPy+ core modules:
from typing import List, Dict, Tuple, Union
import nrpy.c_function as cfc
from nrpy.infrastructures.BHaH import griddata_commondata
from nrpy.infrastructures.BHaH import BHaH_defines_h

def register_CFunction_charecommstruct_set_up(CoordSystem: str) -> None:
    """
    Register C function for setting up charecommstruct.
    
    """
    includes = [
        "BHaH_defines.h",
        "BHaH_function_prototypes.h",
    ]
    desc = r"""Setup charecommstruct"""
    c_type = "void"
    name = "charecommstruct_set_up"
    params = "const commondata_struct *restrict commondata, const params_struct *restrict params, const params_struct *restrict params_chare, charecomm_struct* charecommstruct, const int thischareindex[3]"
    body = r"""
  const int Nxx_plus_2NGHOSTS0 = params->Nxx_plus_2NGHOSTS0;
  const int Nxx_plus_2NGHOSTS1 = params->Nxx_plus_2NGHOSTS1;
  const int Nxx_plus_2NGHOSTS2 = params->Nxx_plus_2NGHOSTS2;  
  const int Nxx0 = params->Nxx0;
  const int Nxx1 = params->Nxx1;
  const int Nxx2 = params->Nxx2;  
  const int Nxx_plus_2NGHOSTS0chare = params_chare->Nxx_plus_2NGHOSTS0;
  const int Nxx_plus_2NGHOSTS1chare = params_chare->Nxx_plus_2NGHOSTS1;
  const int Nxx_plus_2NGHOSTS2chare = params_chare->Nxx_plus_2NGHOSTS2;
  const int Nxx0chare = params_chare->Nxx0;
  const int Nxx1chare = params_chare->Nxx1;
  const int Nxx2chare = params_chare->Nxx2;

  const int ntot = Nxx_plus_2NGHOSTS0 * Nxx_plus_2NGHOSTS1 * Nxx_plus_2NGHOSTS2;
  charecommstruct->globalidx3pt_to_chareidx3 = (int *)malloc(sizeof(int) * ntot); 
  charecommstruct->globalidx3pt_to_localidx3pt = (int *)malloc(sizeof(int) * ntot);

  int startlocali;
  int startlocalj;
  int startlocalk;
  int endlocali;
  int endlocalj;
  int endlocalk;
  for (int charek = 0; charek < Nchare2; charek++) {
    if (charek == 0) {
      startlocalk = 0;
    } else {
      startlocalk = NGHOSTS;
    }
    if (charek == Nchare2 - 1) {
      endlocalk = Nxx_plus_2NGHOSTS2chare - 1;
    } else  {
      endlocalk = Nxx2chare + NGHOSTS - 1;
    }
    for (int charej = 0; charej < Nchare1; charej++) {    
      if (charej == 0) {
        startlocalj = 0;
      } else {
        startlocalj = NGHOSTS;
      }
      if (charej == Nchare1 - 1) {
        endlocalj = Nxx_plus_2NGHOSTS1chare - 1;
      } else {
        endlocalj = Nxx1chare + NGHOSTS - 1;
      }
      for (int charei = 0; charei < Nchare0; charei++) {
        if (charei == 0) {
          startlocali = 0;
        } else {
          startlocali = NGHOSTS;
        }
        if (charei == Nchare0 - 1) {
          endlocali = Nxx_plus_2NGHOSTS0chare - 1;
        } else  {
          endlocali = Nxx0chare + NGHOSTS - 1;
        }
        for (int localk = startlocalk; localk <= endlocalk; localk++) {
          for (int localj = startlocalj; localj <= endlocalj; localj++) {
            for (int locali = startlocali; locali <= endlocali; locali++) {
              
              int localidx3 =  IDX3GENERAL(locali, localj, localk, Nxx0chare, Nxx1chare);
              
              int globali = mapLocalToGlobalIdx0(charei, locali, Nxx0chare);
              int globalj = mapLocalToGlobalIdx1(charej, localj, Nxx1chare);
              int globalk = mapLocalToGlobalIdx2(charek, localk, Nxx2chare);                            
              int globalidx3 =  IDX3GENERAL(globali, globalj, globalk, Nxx0, Nxx1);
                            
              charecommstruct->globalidx3pt_to_chareidx3[globalidx3] = IDX3_OF_CHARE(charei, charej, charek);              
              
              charecommstruct->globalidx3pt_to_localidx3pt[globalidx3] = localidx3;
            }
          }
        }
      }// end charei
    }// end charej
  }// end charek
  
  // local to this chare
  const int ntotchare = Nxx_plus_2NGHOSTS0chare * Nxx_plus_2NGHOSTS1chare * Nxx_plus_2NGHOSTS2chare;
  charecommstruct->localidx3pt_to_globalidx3pt = (int *restrict)malloc(sizeof(int) * ntotchare);
  startlocali = 0;
  startlocalj = 0;
  startlocalk = 0;
  endlocali = Nxx_plus_2NGHOSTS0chare - 1;
  endlocalj = Nxx_plus_2NGHOSTS1chare - 1;
  endlocalk = Nxx_plus_2NGHOSTS2chare - 1;
  for (int localk = startlocalk; localk <= endlocalk; localk++) {
    for (int localj = startlocalj; localj <= endlocalj; localj++) {
      for (int locali = startlocali; locali <= endlocali; locali++) {
        
        int localidx3 =  IDX3GENERAL(locali, localj, localk, Nxx0chare, Nxx1chare);
        
        int globali = mapLocalToGlobalIdx0(thischareindex[0], locali, Nxx0chare);
        int globalj = mapLocalToGlobalIdx0(thischareindex[1], localj, Nxx1chare);
        int globalk = mapLocalToGlobalIdx0(thischareindex[2], localk, Nxx2chare);
        int globalidx3 =  IDX3GENERAL(globali, globalj, globalk, Nxx0, Nxx1);
                      
        charecommstruct->localidx3pt_to_globalidx3pt[localidx3] = globalidx3;
      }
    }
  }
"""
    prefunc = r"""
// The following depend on the chare index
// A global grid point can be part of 2 local grids
static int mapLocalToGlobalIdx0(int chareidx0, int local_idx0, int Nxx0chare) {
    return (chareidx0 * Nxx0chare) + local_idx0;
}
static int mapLocalToGlobalIdx1(int chareidx1, int local_idx1, int Nxx1chare) {
    return (chareidx1 * Nxx1chare) + local_idx1
}
static int mapLocalToGlobalIdx2(int chareidx2, int local_idx2, int Nxx2chare) {
    return (chareidx2 * Nxx2chare) + local_idx2;
}
// Assumes grid point point lies within local grid of chare
static int mapGlobalToLocalIdx0(int chareidx0, int global_idx0, int Nxx0chare) {
    return global_idx0 - (chareidx0 * Nxx0chare);
}
static int mapGlobalToLocalIdx1(int chareidx1, int global_idx1, int Nxx1chare) {
    return global_idx1 - (chareidx1 * Nxx1chare);
}
static int mapGlobalToLocalIdx2(int chareidx2, int global_idx2, int Nxx2chare) {
    return global_idx2 - (chareidx2 * Nxx2chare);
}
"""
    cfc.register_CFunction(
        includes=includes,
        prefunc=prefunc,
        desc=desc,
        c_type=c_type,
        CoordSystem_for_wrapper_func=CoordSystem,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
    )


def chare_comm_register_C_functions(
    list_of_CoordSystems: List[str],    
) -> None:
    """
    :param list_of_CoordSystems: List of coordinate systems to use.    
    :return: None
    """
    for CoordSystem in list_of_CoordSystems:
        # Register C function to set up the boundary condition struct.
        register_CFunction_charecommstruct_set_up(CoordSystem=CoordSystem)

    # Register charecomm_struct's contribution to griddata_struct:
    griddata_commondata.register_griddata_commondata(
        __name__,
        "charecomm_struct charecommstruct",
        "maps that convert between index of a pt in chare's local grid to the index on the global grid, etc",
    )
    
    BHaH_defines_h.register_BHaH_defines(
    __name__,
    r"""#define IDX3_OF_CHARE(i, j, k) ((i) + Nchare0 * ((j) + Nchare1 * ((k))))
#define IDX3GENERAL(i, j, k, Ni, Nj) ((i) + (Ni) * ((j) + (Nj) * (k)))
#define REVERSE_IDX3GENERAL(index, Ni, Nj) \
{ \
	int k = (index) % (Nj); \
	int temp = (index) / (Nj); \
	int j = temp % (Ni); \
	int i = temp / (Ni); \
	i, j, k; \
}

typedef struct __charecomm_struct__ {
  int *restrict globalidx3pt_to_chareidx3;    // which chare is evolving or applying bcs to grid point
  int *restrict globalidx3pt_to_localidx3pt;  // local index of grid point on chare that is evolving or setting bcs for gridpoint
  int *restrict localidx3pt_to_globalidx3pt;  // local to this chare  
} charecomm_struct;    
""",
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
