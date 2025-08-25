"""
Generate horizon_finder.cpp, horizon_finder.h, horizon_finder.ci for the superB infrastructure.

Author: Nishita Jadoo
        njadoo **at** uidaho **dot* edu
"""

from pathlib import Path

from nrpy.helpers.generic import clang_format


def output_horizon_finder_h(
    project_dir: str,
) -> None:
    """
    Generate horizon_finder.h.
    :param project_dir: Directory where the project C code is output
    """
    project_Path = Path(project_dir)
    project_Path.mkdir(parents=True, exist_ok=True)

    file_output_str = """
#ifndef __HORIZON_FINDER_H__
#define __HORIZON_FINDER_H__
#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"
#include "superB/superB_pup_function_prototypes.h"
#include "horizon_finder.decl.h"

class Horizon_finder : public CBase_Horizon_finder {
Horizon_finder_SDAG_CODE

    private :
  /// Member Variables (Object State) ///
  commondata_struct commondata;
  griddata_struct *griddata;
  REAL x_center, y_center, z_center;
  int total_elements;
  int which_horizon;
  REAL *radii;
  REAL (*dst_x0x1x2)[3];
  REAL **dst_data_ptrs;
  bool write_diagnostics_this_step;

  /// Member Functions (private) ///

public:
  /// Constructors ///
  Horizon_finder();
  Horizon_finder(CkMigrateMessage *msg);
  /// Destructor ///
  ~Horizon_finder();
  void process_interpolation_results(InterpBufMsg *msg);

  /// Entry Methods ///
};

#endif //__HORIZON_FINDER_H__


"""
    horizon_finder_h_file = project_Path / "horizon_finder.h"
    with horizon_finder_h_file.open("w", encoding="utf-8") as file:
        file.write(clang_format(file_output_str))


def output_horizon_finder_cpp(
    project_dir: str,
) -> None:
    """
    Generate the horizon_finder.cpp.
    :param project_dir: Directory where the project C code is output
    """
    project_Path = Path(project_dir)
    project_Path.mkdir(parents=True, exist_ok=True)

    file_output_str = r"""
#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"
#include "main.h"
#include "timestepping.h"
#include "interpolator3d.h"
#include "horizon_finder.h"
#include "interp_buf_msg.h"

extern/* readonly */ CProxy_Main mainProxy;
extern/* readonly */ CProxy_Timestepping timesteppingArray;
extern/* readonly */ CProxy_Interpolator3d interpolator3dArray;
extern/* readonly */ CProxy_Horizon_finder horizon_finderProxy;

Horizon_finder::Horizon_finder() {
  CkPrintf("Horizon_finder chare %d created on PE %d\n", thisIndex, CkMyPe());
}

// migration constructor
Horizon_finder::Horizon_finder(CkMigrateMessage *msg) : CBase_Horizon_finder(msg) {}

// destructor
Horizon_finder::~Horizon_finder() {}

void Horizon_finder::process_interpolation_results(InterpBufMsg *msg) {
  char *buf = msg->buf;
  size_t bufSz = msg->len;
  const size_t bytes_per_pt = sizeof(int) + BHAHAHA_NUM_INTERP_GFS * sizeof(REAL);

  // Sanity check
  if (bufSz % bytes_per_pt != 0) {
      CkAbort("Error: Buffer size is not a multiple of bytes_per_pt in report_interpolation_results!");
  }
  int npts = bufSz / bytes_per_pt;
  char *p = buf;
  for (int k = 0; k < npts; ++k) {
      int idx;
      memcpy(&idx, p, sizeof(int));
      p += sizeof(int);
      for (int gf = 0; gf < BHAHAHA_NUM_INTERP_GFS; gf++) {
          memcpy(&dst_data_ptrs[gf][idx], p, sizeof(REAL));
          p += sizeof(REAL);
      }
  }
  delete msg;
}

#include "horizon_finder.def.h"
"""
    horizon_finder_cpp_file = project_Path / "horizon_finder.cpp"
    with horizon_finder_cpp_file.open("w", encoding="utf-8") as file:
        file.write(clang_format(file_output_str))


def output_horizon_finder_ci(
    project_dir: str,
) -> None:
    """
    Generate horizon_finder.ci.

    :param project_dir: Directory where the project C code is output
    """
    project_Path = Path(project_dir)
    project_Path.mkdir(parents=True, exist_ok=True)

    file_output_str = """
module horizon_finder {
  include "BHaH_defines.h";
  include "BHaH_function_prototypes.h";
  include "commondata_object.h";
  include "griddata_object.h";
  include "ckio.h";
  include "pup_stl.h";

  extern message InterpBufMsg {
    int horizon_idx;
    int len;
    char buf[];
  };

  array[1D] Horizon_finder {
    entry Horizon_finder();
    entry void start(CommondataObject & inData1, GriddataObject & inData2) {
      serial {
        commondata = inData1.commondata;
        griddata = inData2.griddata;
      }
      while (commondata.time < commondata.t_final) {
        serial {
          write_diagnostics_this_step = fabs(round(commondata.time / commondata.diagnostics_output_every) * commondata.diagnostics_output_every -
                                             commondata.time) < 0.5 * commondata.dt;
        }
        if (write_diagnostics_this_step) {
          serial {
            which_horizon = thisIndex;
          }
          // In BBH mode, skip horizons that aren't active.
          if (!(commondata.bah_BBH_mode_enable && !commondata.bah_BBH_mode_horizon_active[which_horizon])) {
            serial {
              //bhahaha_find_horizons(&commondata, griddata, &x_center, &y_center, &z_center, &radii, &total_elements, &dst_x0x1x2, &dst_data_ptrs, which_horizon, BHAHAHA_FIND_HORIZONS_SETUP);
            }

            when ready_for_interpolation() {
              serial {
                interpolator3dArray[CkArrayIndex3D(0, 0, 0)].charezero_start_interpolation(thisIndex, BHAHAHA_NUM_INTERP_GFS, total_elements, (REAL*)dst_x0x1x2);
              }
            }

            when report_interpolation_results(InterpBufMsg *msg) {
              serial {
                process_interpolation_results(msg);
              }
            }
            serial {
              //bhahaha_find_horizons(&commondata, griddata, &x_center, &y_center, &z_center, &radii, &total_elements, &dst_x0x1x2, &dst_data_ptrs, which_horizon, BHAHAHA_FIND_HORIZONS_FIND_AND_WRITE_TO_FILE);
              free(radii);
              thisProxy[CkArrayIndex1D(thisIndex)].horizon_finding_complete();
            }

            when horizon_finding_complete() {
              serial{}
            }
          } // END condition: skip horizons that aren't active.
        } // END condition: write_diagnostics_this_step

        serial {
          // Adding dt to commondata.time many times will induce roundoff error,
          //   so here we set time based on the iteration number.
          commondata.time = (REAL)(commondata.nn + 1) * commondata.dt;
          // Finally, increment the timestep n:
          commondata.nn++;
        }
      }
      serial { mainProxy.done(); }
    };
    entry void report_interpolation_results(InterpBufMsg *msg);
    entry void ready_for_interpolation();
    entry void horizon_finding_complete();
   };
};
"""
    horizon_finder_ci_file = project_Path / "horizon_finder.ci"
    with horizon_finder_ci_file.open("w", encoding="utf-8") as file:
        file.write(clang_format(file_output_str))


def output_horizon_finder_h_cpp_ci(
    project_dir: str,
) -> None:
    """
    Generate horizon_finder.h, horizon_finder.cpp and horizon_finder.ci.
    :param project_dir: Directory where the project C code is output.
    """

    output_horizon_finder_h(
        project_dir=project_dir,
    )

    output_horizon_finder_cpp(
        project_dir=project_dir,
    )

    output_horizon_finder_ci(
        project_dir=project_dir,
    )
