"""
Generate interpolator3d.cpp, interpolator3d.h and interpolator3d.ci for the superB infrastructure.

Author: Nishita Jadoo
        njadoo **at** uidaho **dot* edu
"""

from pathlib import Path

from nrpy.helpers.generic import clang_format


def output_interpolator3d_h(
    project_dir: str,
) -> None:
    """
    Generate interpolator3d.h.
    :param project_dir: Directory where the project C code is output
    """
    project_Path = Path(project_dir)
    project_Path.mkdir(parents=True, exist_ok=True)

    file_output_str = """
#ifndef __INTERPOLATOR3D_H__
#define __INTERPOLATOR3D_H__
#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"
#include "superB/superB_pup_function_prototypes.h"
#include "interpolator3d.decl.h"

class Interpolator3d : public CBase_Interpolator3d {
Interpolator3d_SDAG_CODE

    private :
  /// Member Variables (Object State) ///
  commondata_struct commondata;
  griddata_struct *griddata_chare;
  const int grid = 0;
  REAL (*dst_x0x1x2)[3] = NULL;
  REAL (*dst_x0x1x2_chare)[3];
  REAL *dst_data_ptrs_chare[BHAHAHA_NUM_INTERP_GFS];
  int total_elements_chare;
  const REAL *restrict src_gf_ptrs[BHAHAHA_NUM_INTERP_GFS];
  int *dst_indices_chare = nullptr;
  int iter = 0;
  char **interp_bufs=nullptr;
  int *interp_lens=nullptr;
  int interp_count = 0;
  int interp_total = 0;
  int interp_horizon_idx = 0;

  /// Member Functions (private) ///
  void contribute_interpolation_results(int curr_index_horizonfinder_chare);

public:
  /// Constructors ///
  Interpolator3d();
  Interpolator3d(CkMigrateMessage *msg);
  /// Destructor ///
  ~Interpolator3d();
  void recv_interp_msg(InterpBufMsg *m);
  void send_interp_concat();

  /// Entry Methods ///

};

#endif //__INTERPOLATOR3D_H__
"""
    interpolator3d_h_file = project_Path / "interpolator3d.h"
    with interpolator3d_h_file.open("w", encoding="utf-8") as file:
        file.write(clang_format(file_output_str))


def output_interpolator3d_cpp(
    project_dir: str,
) -> None:
    """
    Generate the interpolator3d.cpp.
    :param project_dir: Directory where the project C code is output
    """
    project_Path = Path(project_dir)
    project_Path.mkdir(parents=True, exist_ok=True)

    file_output_str = r"""
#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"
#include "timestepping.h"
#include "horizon_finder.h"
#include "interpolator3d.h"
#include "InterpBufMsg.h"
#include <cstring>

extern/* readonly */ CProxy_Timestepping timesteppingArray;
extern/* readonly */ CProxy_Horizon_finder horizon_finderProxy;
extern/* readonly */ CProxy_Interpolator3d interpolator3dArray;

Interpolator3d::Interpolator3d() {
  CkPrintf("Interpolator3d chare %d%d%d created on PE %d\n", thisIndex.x, thisIndex.y, thisIndex.z, CkMyPe());
}

// migration constructor
Interpolator3d::Interpolator3d(CkMigrateMessage *msg) : CBase_Interpolator3d(msg) {}

// destructor
Interpolator3d::~Interpolator3d() {}

void Interpolator3d::contribute_interpolation_results(int curr_index_horizonfinder_chare) {
  // We have total_elements_chare number of grid points
  // For each grid point, we send its index in the original dst_x0x1x2 array
  // and BHAHAHA_NUM_INTERP_GFS number of gfs values at that grid point

  // 1) compute how many bytes each point contributes:
  size_t bytes_per_point = sizeof(int) + BHAHAHA_NUM_INTERP_GFS * sizeof(REAL);
  size_t total_bytes = total_elements_chare * bytes_per_point;

  // 2) pack into a raw buffer with memcpy
  char *packbuf = (char*)malloc(total_bytes);
  char *p = packbuf;

  for (int k = 0; k < total_elements_chare; ++k) {
    // copy the integer index
    std::memcpy(p, &dst_indices_chare[k], sizeof(int));
    p += sizeof(int);

    // copy the BHAHAHA_NUM_INTERP_GFS REAL values
    for (int i = 0; i < BHAHAHA_NUM_INTERP_GFS; ++i) {
      std::memcpy(p, &dst_data_ptrs_chare[i][k], sizeof(REAL));
      p += sizeof(REAL);
    }
  }

  {
     InterpBufMsg *m = new (total_bytes) InterpBufMsg();
     m->horizon_idx = curr_index_horizonfinder_chare;
     m->len = (int) total_bytes;
     memcpy(m->buf, packbuf, total_bytes);
     thisProxy[CkArrayIndex3D(0,0,0)].recv_interp_msg(m);
  }
  free(packbuf);
}

void Interpolator3d::recv_interp_msg(InterpBufMsg *m){
  if(interp_count == 0) {
    int Nchare0 = commondata.Nchare0;
    int Nchare1 = commondata.Nchare1;
    int Nchare2 = commondata.Nchare2;
    int Ncharetotal = Nchare0 * Nchare1 * Nchare2;
    interp_total = Ncharetotal;
    interp_horizon_idx = m->horizon_idx;
    interp_bufs = new char*[interp_total];
    interp_lens = new int[interp_total];
  }
  interp_lens[interp_count] = m->len;

  // Allocate a new buffer and copy the data from the message.
  interp_bufs[interp_count] = new char[m->len];
  memcpy(interp_bufs[interp_count], m->buf, m->len);

  interp_count++;
  delete m; // This is now safe because we have our own copy of the data.

  if(interp_count == interp_total){
     send_interp_concat();
   }
}

void Interpolator3d::send_interp_concat(){
  size_t tot=0;
  for(int i=0; i<interp_total; i++){
     tot += interp_lens[i];
  }
  char *agg=(char*)malloc(tot);
  size_t off=0;
  for(int i=0;i<interp_total;i++){
    memcpy(agg+off, interp_bufs[i], interp_lens[i]);
    off += interp_lens[i];
  }
  InterpBufMsg *out = new (tot) InterpBufMsg();
  out->horizon_idx = interp_horizon_idx;
  out->len = (int) tot;
  memcpy(out->buf, agg, tot);
  horizon_finderProxy[CkArrayIndex1D(interp_horizon_idx)].report_interpolation_results(out);

  free(agg);

  // Deallocate each individual buffer that was allocated in recv_interp_msg.
  for (int i = 0; i < interp_total; i++) {
    delete[] interp_bufs[i];
  }
  // Now free the array that held the pointers.
  delete[] interp_bufs;

  delete[] interp_lens;
  interp_bufs=nullptr;
  interp_lens=nullptr;
  interp_count=interp_total=interp_horizon_idx=0;

  thisProxy.interp_concatenation_complete();
}

#include "interpolator3d.def.h"
"""
    interpolator3d_cpp_file = project_Path / "interpolator3d.cpp"
    with interpolator3d_cpp_file.open("w", encoding="utf-8") as file:
        file.write(clang_format(file_output_str))


def output_interpolator3d_ci(
    project_dir: str,
) -> None:
    """
    Generate interpolator3d.ci.

    :param project_dir: Directory where the project C code is output
    """
    project_Path = Path(project_dir)
    project_Path.mkdir(parents=True, exist_ok=True)

    file_output_str = """
module interpolator3d {
  include "BHaH_defines.h";
  include "BHaH_function_prototypes.h";
  include "commondata_object.h";
  include "griddata_object.h";
  include "ckio.h";
  include "pup_stl.h";

  message InterpBufMsg {
    int horizon_idx;
    int len;
    char buf[];
  };

  array[3D] Interpolator3d {
    entry Interpolator3d();
    entry void start(CommondataObject & inData1, GriddataObject & inData2) {
      serial {
        commondata = inData1.commondata;
        griddata_chare = inData2.griddata;
      }
      while (commondata.time < commondata.t_final) {
        when receiv_bhahaha_gfs(int len_tmpBuffer, REAL tmpBuffer[len_tmpBuffer]) {
          serial {
            const int Nxx_plus_2NGHOSTS0 = griddata_chare[grid].params.Nxx_plus_2NGHOSTS0;
            const int Nxx_plus_2NGHOSTS1 = griddata_chare[grid].params.Nxx_plus_2NGHOSTS1;
            const int Nxx_plus_2NGHOSTS2 = griddata_chare[grid].params.Nxx_plus_2NGHOSTS2;
            const int Nxx_plus_2NGHOSTS_tot = Nxx_plus_2NGHOSTS0 * Nxx_plus_2NGHOSTS1 * Nxx_plus_2NGHOSTS2;
            for (int idx = 0; idx < BHAHAHA_NUM_INTERP_GFS; idx++) {
              src_gf_ptrs[idx] = tmpBuffer + idx * Nxx_plus_2NGHOSTS_tot;
            }
          }

          if (thisIndex.x == 0 && thisIndex.y == 0 && thisIndex.z == 0) {
            serial {
              horizon_finderProxy.ready_for_interpolation();
            }
          }
          for (iter = 0; iter < commondata.bah_max_num_horizons; iter++) {
            if (thisIndex.x == 0 && thisIndex.y == 0 && thisIndex.z == 0) {
              when charezero_start_interpolation[iter](int index_horizonfinder_chare, int bhahaha_num_interp_gfs, int total_elements, REAL dst_x0x1x2_linear[3*total_elements]) {
                serial {
                  thisProxy.start_interpolation(index_horizonfinder_chare, bhahaha_num_interp_gfs, total_elements, dst_x0x1x2_linear);
                }
              }
            }

            when start_interpolation(int index_horizonfinder_chare, int bhahaha_num_interp_gfs, int total_elements, REAL dst_x0x1x2_linear[3*total_elements]) {
              serial {
                dst_x0x1x2 = (REAL (*)[3])dst_x0x1x2_linear;
                int count_total_elements_chare = 0;
                for (int i = 0; i < total_elements; i++) {
                  if ((griddata_chare->params.xxmin0 <= dst_x0x1x2[i][0] && dst_x0x1x2[i][0] <= griddata_chare->params.xxmax0) &&
                      (griddata_chare->params.xxmin1 <= dst_x0x1x2[i][1] && dst_x0x1x2[i][1] <= griddata_chare->params.xxmax1) &&
                      (griddata_chare->params.xxmin2 <= dst_x0x1x2[i][2] && dst_x0x1x2[i][2] <= griddata_chare->params.xxmax2)) {

                    count_total_elements_chare++;
                  }
                }
                total_elements_chare = count_total_elements_chare;
                dst_x0x1x2_chare = malloc(total_elements_chare * 3 * sizeof(REAL));
                dst_indices_chare = (int*)malloc(total_elements_chare * sizeof(int));

                count_total_elements_chare = 0;
                for (int i = 0; i < total_elements; i++) {
                  if (griddata_chare->params.xxmin0 <= dst_x0x1x2[i][0] && dst_x0x1x2[i][0] <= griddata_chare->params.xxmax0 &&
                      griddata_chare->params.xxmin1 <= dst_x0x1x2[i][1] && dst_x0x1x2[i][1] <= griddata_chare->params.xxmax1 &&
                      griddata_chare->params.xxmin2 <= dst_x0x1x2[i][2] && dst_x0x1x2[i][2] <= griddata_chare->params.xxmax2) {
                    dst_x0x1x2_chare[count_total_elements_chare][0] = dst_x0x1x2[i][0];
                    dst_x0x1x2_chare[count_total_elements_chare][1] = dst_x0x1x2[i][1];
                    dst_x0x1x2_chare[count_total_elements_chare][2] = dst_x0x1x2[i][2];

                    dst_indices_chare[count_total_elements_chare] = i;

                    count_total_elements_chare++;
                  }
                }
                for (int i = 0; i < BHAHAHA_NUM_INTERP_GFS; i++) {
                  dst_data_ptrs_chare[i] = (REAL *)malloc(total_elements_chare * sizeof(REAL));
                }
                const int Nxx_plus_2NGHOSTS0 = griddata_chare[grid].params.Nxx_plus_2NGHOSTS0;
                const int Nxx_plus_2NGHOSTS1 = griddata_chare[grid].params.Nxx_plus_2NGHOSTS1;
                const int Nxx_plus_2NGHOSTS2 = griddata_chare[grid].params.Nxx_plus_2NGHOSTS2;

                interpolation_3d_general__uniform_src_grid((NGHOSTS - 1), griddata_chare[grid].params.dxx0, griddata_chare[grid].params.dxx1, griddata_chare[grid].params.dxx2, Nxx_plus_2NGHOSTS0,
                                                       Nxx_plus_2NGHOSTS1, Nxx_plus_2NGHOSTS2, BHAHAHA_NUM_INTERP_GFS, griddata_chare[grid].xx, src_gf_ptrs,
                                                       total_elements_chare, dst_x0x1x2_chare, dst_data_ptrs_chare);
                contribute_interpolation_results(index_horizonfinder_chare);
                free(dst_x0x1x2_chare);
                free(dst_indices_chare);
                for (int i = 0; i < BHAHAHA_NUM_INTERP_GFS; i++) {
                    free(dst_data_ptrs_chare[i]);
                }
              }

              // Only continue when custom reduction on interpolator chare 0 is complete
              when interp_concatenation_complete() {
                serial {}
              }
            }
          } // end for (iter = 0; iter < commondata.bah_max_num_horizons; iter++)
        } // end when receiv_bhahaha_gfs
        serial {
          // Adding dt to commondata.time many times will induce roundoff error,
          //   so here we set time based on the iteration number.
          commondata.time = (REAL)(commondata.nn + 1) * commondata.dt;
          // Finally, increment the timestep n:
          commondata.nn++;
        }
      }// end while (commondata.time < commondata.t_final)
    };
    entry void start_interpolation(int index_horizonfinder_chare, int bhahaha_num_interp_gfs, int total_elements, REAL dst_x0x1x2_linear[3*total_elements]);
    entry void charezero_start_interpolation(int index_horizonfinder_chare, int bhahaha_num_interp_gfs, int total_elements, REAL dst_x0x1x2_linear[3*total_elements]);
    entry void receiv_bhahaha_gfs(int len_tmpBuffer, REAL tmpBuffer[len_tmpBuffer]);
    entry void recv_interp_msg(InterpBufMsg *m);
    entry void send_interp_concat();
    entry void interp_concatenation_complete();
   };
};
"""
    interpolator3d_ci_file = project_Path / "interpolator3d.ci"
    with interpolator3d_ci_file.open("w", encoding="utf-8") as file:
        file.write(clang_format(file_output_str))


def output_interpolator3d_h_cpp_ci(
    project_dir: str,
) -> None:
    """
    Generate interpolator3d.h, interpolator3d.cpp and interpolator3d.ci.
    :param project_dir: Directory where the project C code is output.
    """

    output_interpolator3d_h(
        project_dir=project_dir,
    )

    output_interpolator3d_cpp(
        project_dir=project_dir,
    )

    output_interpolator3d_ci(
        project_dir=project_dir,
    )
