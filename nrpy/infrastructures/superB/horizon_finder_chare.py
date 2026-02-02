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
  REAL *x_guess;
  REAL *y_guess;
  REAL *z_guess;
  int total_elements;
  int which_horizon;
  REAL *radii;
  REAL (*dst_x0x1x2)[3];
  REAL **dst_data_ptrs;
  bool do_horizon_find;

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
#include "InterpBufMsg.h"

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
Horizon_finder::~Horizon_finder() {
  if (x_guess) free(x_guess);
  if (y_guess) free(y_guess);
  if (z_guess) free(z_guess);
}

void Horizon_finder::process_interpolation_results(InterpBufMsg *msg) {
  char *buf = msg->buf;
  size_t bufSz = msg->len;
  const int num_gfs = msg->num_gfs;
  // Keep this path quiet during normal runs.
  if (num_gfs != BHAHAHA_NUM_INTERP_GFS) {
      CkAbort("Error: Horizon interpolation received unexpected number of gridfunctions.");
  }
  const size_t bytes_per_pt = sizeof(int) + (size_t)num_gfs * sizeof(REAL);

  // Sanity check
  if (bufSz % bytes_per_pt != 0) {
      CkAbort("Error: Buffer size is not a multiple of bytes_per_pt in report_interpolation_results!");
  }
  if (unpack_interpolation_buffer(num_gfs, buf, bufSz, dst_data_ptrs) != 0) {
      CkAbort("Error: Failed to unpack interpolation buffer in report_interpolation_results!");
  }
  (void)bytes_per_pt;
  delete msg;
}

#include "horizon_finder.def.h"
"""
    horizon_finder_cpp_file = project_Path / "horizon_finder.cpp"
    with horizon_finder_cpp_file.open("w", encoding="utf-8") as file:
        file.write(file_output_str)


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
    int request_type;
    int request_id;
    int num_gfs;
    int len;
    char buf[];
  };

  array[1D] Horizon_finder {
    entry Horizon_finder();
    entry void start(CommondataObject & inData1, GriddataObject & inData2) {
      serial {
        commondata = inData1.commondata;
        griddata = inData2.griddata;
        x_guess = (REAL *)malloc(commondata.bah_max_num_horizons * sizeof(REAL));
        y_guess = (REAL *)malloc(commondata.bah_max_num_horizons * sizeof(REAL));
        z_guess = (REAL *)malloc(commondata.bah_max_num_horizons * sizeof(REAL));
      }
      while (commondata.time < commondata.t_final) {
        serial {
          do_horizon_find = true;  // default: yes, find horizon
          // STEP 1: Check if horizon find is scheduled for the current iteration.
          if (commondata.diagnostics_output_every <= 0 ||
              (commondata.nn % (int)(commondata.diagnostics_output_every / commondata.dt + 0.5)) != 0) {
            int bah_find_every = 1;  // Placeholder: find every iteration. This should be a commondata param.
            if (commondata.diagnostics_output_every > commondata.dt) {
              // A basic way to get find_every from time interval
              bah_find_every = (int)(commondata.diagnostics_output_every / commondata.dt + 0.5);
              if (bah_find_every == 0)
                bah_find_every = 1;
            }
            if (bah_find_every <= 0 || (commondata.nn % bah_find_every) != 0) {
              do_horizon_find = false;  // not scheduled this iteration
            }
          } // END IF: diagnostics_output_every > 0
        }
        if (do_horizon_find) {
          serial {
            which_horizon = thisIndex;
          }
          // In BBH mode, skip horizons that aren't active.
          if (!(commondata.bah_enable_BBH_mode && !commondata.bah_BBH_mode_horizon_active[which_horizon])) {
            serial {
              bhahaha_find_horizons(&commondata, griddata, x_guess, y_guess, z_guess, &radii, &total_elements, &dst_x0x1x2, &dst_data_ptrs, which_horizon, BHAHAHA_FIND_HORIZONS_SETUP);
            }

            when ready_for_interpolation() {
              serial {
                interpolator3dArray.start_interpolation(INTERP_REQUEST_BHAHAHA, thisIndex, BHAHAHA_NUM_INTERP_GFS,
                                                       total_elements, (REAL *)dst_x0x1x2);
              }
            }

            when report_interpolation_results(InterpBufMsg *msg) {
              serial {
                process_interpolation_results(msg);
              }
            }
            serial {
              bhahaha_find_horizons(&commondata, griddata, x_guess, y_guess, z_guess, &radii, &total_elements, &dst_x0x1x2, &dst_data_ptrs, which_horizon, BHAHAHA_FIND_HORIZONS_FIND_AND_WRITE_TO_FILE);
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
