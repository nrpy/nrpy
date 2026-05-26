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
  griddata_struct *griddata = nullptr;
  bool owns_griddata = false;
  REAL *x_guess = nullptr;
  REAL *y_guess = nullptr;
  REAL *z_guess = nullptr;
  int total_elements = 0;
  int which_horizon = 0;
  REAL *radii = nullptr;
  REAL (*dst_x0x1x2)[3] = nullptr;
  REAL **dst_data_ptrs = nullptr;
  bool do_horizon_find = false;

  /// Member Functions (private) ///

public:
  /// Constructors ///
  Horizon_finder();
  Horizon_finder(CkMigrateMessage *msg);
  void pup(PUP::er &p);
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

/**
 * PUP an optional REAL array owned by Horizon_finder.
 *
 * @param[in,out] p Charm++ PUP serializer.
 * @param[in,out] array Pointer to the optional array pointer.
 * @param length Number of REAL entries to serialize.
 * @param[in] name Allocation/error label.
 */
static void horizon_finder_pup_optional_REAL_array(PUP::er &p, REAL **array, const int length, const char *name) {
  int has_array = (*array != nullptr) ? 1 : 0;
  p | has_array;
  if (has_array != 0) {
    if (p.isUnpacking()) {
      *array = (REAL *restrict)malloc(sizeof(REAL) * length);
      if (*array == nullptr && length > 0) {
        CkAbort("%s", name);
      } // END IF: optional REAL array allocation failed
    } // END IF: unpacking optional REAL array
    PUParray(p, *array, length);
  } else if (p.isUnpacking()) {
    *array = nullptr;
  } // END ELSE IF: unpacking absent optional REAL array
} // END FUNCTION: horizon_finder_pup_optional_REAL_array

Horizon_finder::Horizon_finder() {
  CkPrintf("Horizon_finder chare %d created on PE %d\n", thisIndex, CkMyPe());
}

// migration constructor
Horizon_finder::Horizon_finder(CkMigrateMessage *msg) : CBase_Horizon_finder(msg) {}

// PUP routine for class Horizon_finder
void Horizon_finder::pup(PUP::er &p) {
  CBase_Horizon_finder::pup(p);
  __sdag_pup(p);
  pup_commondata_struct(p, commondata);

  int has_griddata = (griddata != nullptr) ? 1 : 0;
  p | has_griddata;
  if (has_griddata != 0) {
    int size_griddata = commondata.NUMGRIDS;
    p | size_griddata;
    if (p.isUnpacking()) {
      owns_griddata = true;
      griddata = (griddata_struct *restrict)malloc(sizeof(griddata_struct) * size_griddata);
      if (griddata == nullptr && size_griddata > 0) {
        CkAbort("Horizon_finder PUP failed to allocate griddata.");
      }
    }
    for (int i = 0; i < size_griddata; i++) {
      pup_griddata(p, griddata[i]);
    }
  } else if (p.isUnpacking()) {
    owns_griddata = false;
    griddata = nullptr;
  }

  p | total_elements;
  p | which_horizon;
  p | do_horizon_find;

  const int num_horizons = commondata.bah_max_num_horizons;
  horizon_finder_pup_optional_REAL_array(p, &x_guess, num_horizons, "Horizon_finder PUP failed to allocate x_guess.");
  horizon_finder_pup_optional_REAL_array(p, &y_guess, num_horizons, "Horizon_finder PUP failed to allocate y_guess.");
  horizon_finder_pup_optional_REAL_array(p, &z_guess, num_horizons, "Horizon_finder PUP failed to allocate z_guess.");

  int has_dst_x0x1x2 = (dst_x0x1x2 != nullptr) ? 1 : 0;
  p | has_dst_x0x1x2;
  if (has_dst_x0x1x2 != 0) {
    if (p.isUnpacking()) {
      dst_x0x1x2 = (REAL(*)[3])malloc(sizeof(REAL) * total_elements * 3);
      if (dst_x0x1x2 == nullptr && total_elements > 0) {
        CkAbort("Horizon_finder PUP failed to allocate dst_x0x1x2.");
      }
    }
    PUParray(p, (REAL *)dst_x0x1x2, total_elements * 3);
  } else if (p.isUnpacking()) {
    dst_x0x1x2 = nullptr;
  }

  int has_dst_data_ptrs = (dst_data_ptrs != nullptr) ? 1 : 0;
  p | has_dst_data_ptrs;
  if (has_dst_data_ptrs != 0) {
    if (p.isUnpacking()) {
      dst_data_ptrs = (REAL **)malloc(sizeof(REAL *) * BHAHAHA_NUM_INTERP_GFS);
      if (dst_data_ptrs == nullptr) {
        CkAbort("Horizon_finder PUP failed to allocate dst_data_ptrs.");
      }
    }
    for (int i = 0; i < BHAHAHA_NUM_INTERP_GFS; i++) {
      if (p.isUnpacking()) {
        dst_data_ptrs[i] = (REAL *restrict)malloc(sizeof(REAL) * total_elements);
        if (dst_data_ptrs[i] == nullptr && total_elements > 0) {
          CkAbort("Horizon_finder PUP failed to allocate dst_data_ptrs[i].");
        }
      }
      PUParray(p, dst_data_ptrs[i], total_elements);
    }
  } else if (p.isUnpacking()) {
    dst_data_ptrs = nullptr;
  }
}

// destructor
Horizon_finder::~Horizon_finder() {
  BHAH_FREE(x_guess);
  BHAH_FREE(y_guess);
  BHAH_FREE(z_guess);
  BHAH_FREE(radii);
  BHAH_FREE(dst_x0x1x2);
  if (dst_data_ptrs != nullptr) {
    for (int i = 0; i < BHAHAHA_NUM_INTERP_GFS; i++) {
      BHAH_FREE(dst_data_ptrs[i]);
    }
    BHAH_FREE(dst_data_ptrs);
  }
  if (owns_griddata && griddata != nullptr) {
    for (int grid = 0; grid < commondata.NUMGRIDS; grid++) {
      BHAH_FREE(griddata[grid].xx[0]);
      BHAH_FREE(griddata[grid].xx[1]);
      BHAH_FREE(griddata[grid].xx[2]);
    }
    BHAH_FREE(griddata);
  }
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
              BHAH_FREE(radii);
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
              BHAH_FREE(dst_x0x1x2);
              if (dst_data_ptrs != nullptr) {
                for (int i = 0; i < BHAHAHA_NUM_INTERP_GFS; i++) {
                  BHAH_FREE(dst_data_ptrs[i]);
                } // END LOOP: for interpolated horizon gridfunctions
                BHAH_FREE(dst_data_ptrs);
              } // END IF: interpolation data pointers allocated
              total_elements = 0;
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
      serial {
        CkCallback cb(CkReductionTarget(Horizon_finder, notify_horizon_finder_completion), thisProxy[CkArrayIndex1D(0)]);
        contribute(cb);
      }
    };
    entry void report_interpolation_results(InterpBufMsg *msg);
    entry void ready_for_interpolation();
    entry void horizon_finding_complete();
    entry [reductiontarget] void notify_horizon_finder_completion() {
      serial {
        mainProxy.horizon_finder_done();
      }
    }
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
