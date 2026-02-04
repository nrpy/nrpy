"""
Generate interpolator3d.cpp, interpolator3d.h and interpolator3d.ci for the superB infrastructure.

Author: Nishita Jadoo
        njadoo **at** uidaho **dot* edu
"""

from pathlib import Path

from nrpy.helpers.generic import clang_format


def output_griddata_object_h(
    project_dir: str,
) -> None:
    r"""
    Output header file with definition for class GriddataObject.
    :param project_dir: Directory where the project C code is output
    """
    project_Path = Path(project_dir)
    project_Path.mkdir(parents=True, exist_ok=True)

    file_output_str = """#ifndef __GRIDDATAOBJECT_H__
#define __GRIDDATAOBJECT_H__
#include "BHaH_defines.h"
#include "superB/superB_pup_function_prototypes.h"

class GriddataObject {
public:

  griddata_struct *griddata = nullptr;
  int size_griddata;

  void pup(PUP::er &p) {
    p | size_griddata;

    if (p.isUnpacking())
      griddata = new griddata_struct[size_griddata];

    // pup_griddata only PUPs params and xx in griddata
    for (int i = 0; i < size_griddata; i++)
      pup_griddata(p, griddata[i]);
  }
};
#endif //__GRIDDATAOBJECT_H__
"""
    commondata_object_file = project_Path / "griddata_object.h"
    with commondata_object_file.open("w", encoding="utf-8") as file:
        file.write(clang_format(file_output_str))


def output_interp_buf_msg_h(
    project_dir: str,
) -> None:
    r"""
    Output header file with definition for class InterpBufMsg.
    :param project_dir: Directory where the project C code is output
    """
    project_Path = Path(project_dir)
    project_Path.mkdir(parents=True, exist_ok=True)

    file_output_str = """#ifndef INTERP_BUF_MSG_H
#define INTERP_BUF_MSG_H

class InterpBufMsg : public CMessage_InterpBufMsg {
public:
  int request_type;
  int request_id;
  int num_gfs;
  int len;
  char *buf;
};

#endif // INTERP_BUF_MSG_H
"""
    interp_buf_msg_file = project_Path / "InterpBufMsg.h"
    with interp_buf_msg_file.open("w", encoding="utf-8") as file:
        file.write(clang_format(file_output_str))


def output_interpolator3d_h(
    project_dir: str,
    enable_psi4: bool = False,
) -> None:
    """
    Generate interpolator3d.h.
    :param project_dir: Directory where the project C code is output
    """
    project_Path = Path(project_dir)
    project_Path.mkdir(parents=True, exist_ok=True)

    psi4_member_vars = ""
    if enable_psi4:
        psi4_member_vars = """  psi4_shell_angular_grid_t psi4_shell = {};
  REAL *psi4r_at_R_ext = nullptr;
  REAL *psi4i_at_R_ext = nullptr;
"""

    file_output_str = """
#ifndef __INTERPOLATOR3D_H__
#define __INTERPOLATOR3D_H__
#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"
#include "interpolator3d.decl.h"
#include "superB/superB_pup_function_prototypes.h"

class Interpolator3d : public CBase_Interpolator3d {
Interpolator3d_SDAG_CODE

    private :
    /// Member Variables (Object State) ///
    commondata_struct commondata;
  griddata_struct *griddata_chare;
  bool owns_griddata_chare = false;
  const int grid = 0;
  REAL (*dst_x0x1x2)[3] = NULL;
  REAL (*dst_x0x1x2_chare)[3];
  REAL **dst_data_ptrs_chare = nullptr;
  int total_elements_chare;
  const REAL *restrict *src_gf_ptrs = nullptr;
  int src_gf_ptrs_capacity = 0;
  int *dst_indices_chare = nullptr;
  int iter = 0;
  char **interp_bufs = nullptr;
  int *interp_lens = nullptr;
  int interp_count = 0;
  int interp_total = 0;
  int interp_request_type = -1;
  int interp_request_id = -1;
  int interp_num_gfs = 0;
<<PSI4_MEMBER_VARS>>

  /// Member Functions (private) ///
  void contribute_interpolation_results(int curr_index_horizonfinder_chare);
  void perform_interpolation(int request_type, int request_id, int num_gfs, int total_elements, REAL *dst_x0x1x2_linear);

public:
  /// Constructors ///
  Interpolator3d();
  Interpolator3d(CkMigrateMessage *msg);
  void pup(PUP::er &p);
  /// Destructor ///
  ~Interpolator3d();
  void recv_interp_msg(InterpBufMsg *m);
  void send_interp_concat();

  /// Entry Methods ///
};

#endif //__INTERPOLATOR3D_H__
"""
    file_output_str = file_output_str.replace("<<PSI4_MEMBER_VARS>>", psi4_member_vars)
    interpolator3d_h_file = project_Path / "interpolator3d.h"
    with interpolator3d_h_file.open("w", encoding="utf-8") as file:
        file.write(clang_format(file_output_str))


def output_interpolator3d_cpp(
    project_dir: str,
    enable_psi4: bool = False,
) -> None:
    """
    Generate the interpolator3d.cpp.
    :param project_dir: Directory where the project C code is output
    """
    project_Path = Path(project_dir)
    project_Path.mkdir(parents=True, exist_ok=True)

    psi4_send_concat_log = ""
    psi4_send_concat_block = ""
    if enable_psi4:
        psi4_send_concat_block = r"""
  } else if (interp_request_type == INTERP_REQUEST_PSI4) {
    if (psi4r_at_R_ext == nullptr || psi4i_at_R_ext == nullptr) {
      CkAbort("Error: PSI4 interpolation results received without allocated psi4 buffers.");
    }
    if (interp_num_gfs != 2) {
      CkAbort("Error: PSI4 interpolation expects exactly 2 gridfunctions.");
    }
    REAL *dst_data_ptrs[2] = {psi4r_at_R_ext, psi4i_at_R_ext};
    if (unpack_interpolation_buffer(interp_num_gfs, agg, tot, dst_data_ptrs) != 0) {
      CkAbort("Error: Failed to unpack PSI4 interpolation buffer.");
    }
    const REAL R_ext = commondata.list_of_psi4_extraction_radii[interp_request_id];
    psi4_spinweightm2_decompose_shell(&commondata, &psi4_shell, commondata.time, R_ext, psi4r_at_R_ext, psi4i_at_R_ext);
    free(psi4r_at_R_ext);
    free(psi4i_at_R_ext);
    psi4r_at_R_ext = nullptr;
    psi4i_at_R_ext = nullptr;
    delete out;
"""

    file_output_str = r"""
#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"
#include "timestepping.h"
#include "horizon_finder.h"
#include "interpolator3d.h"
#include "InterpBufMsg.h"
#include <cstring>
#include <cstdlib>

extern/* readonly */ CProxy_Timestepping timesteppingArray;
extern/* readonly */ CProxy_Horizon_finder horizon_finderProxy;
extern/* readonly */ CProxy_Interpolator3d interpolator3dArray;

Interpolator3d::Interpolator3d() {
  CkPrintf("Interpolator3d chare %d,%d,%d created on PE %d\n", thisIndex.x, thisIndex.y, thisIndex.z, CkMyPe());
}

// migration constructor
Interpolator3d::Interpolator3d(CkMigrateMessage *msg) : CBase_Interpolator3d(msg) {}

// PUP routine for class Interpolator3d
void Interpolator3d::pup(PUP::er &p) {
  CBase_Interpolator3d::pup(p);
  __sdag_pup(p);
  pup_commondata_struct(p, commondata);

  int size_griddata = commondata.NUMGRIDS;
  p | size_griddata;
  if (p.isUnpacking()) {
    owns_griddata_chare = true;
    griddata_chare = (griddata_struct *restrict)malloc(sizeof(griddata_struct) * size_griddata);
  }
  for (int i = 0; i < size_griddata; i++) {
    pup_griddata(p, griddata_chare[i]);
  }

  p | iter;
  p | interp_count;
  p | interp_total;
  p | interp_request_type;
  p | interp_request_id;
  p | interp_num_gfs;

  if (p.isUnpacking()) {
    if (interp_total > 0) {
      interp_bufs = new char *[interp_total];
      interp_lens = new int[interp_total];
      for (int i = 0; i < interp_total; i++) {
        interp_bufs[i] = nullptr;
        interp_lens[i] = 0;
      }
    } else {
      interp_bufs = nullptr;
      interp_lens = nullptr;
    }
  }

  for (int i = 0; i < interp_count; i++) {
    p | interp_lens[i];
    if (p.isUnpacking()) {
      interp_bufs[i] = new char[interp_lens[i]];
    }
    p(interp_bufs[i], interp_lens[i]);
  }

  if (p.isUnpacking()) {
    dst_x0x1x2 = nullptr;
    dst_x0x1x2_chare = nullptr;
    dst_data_ptrs_chare = nullptr;
    dst_indices_chare = nullptr;
    src_gf_ptrs = nullptr;
    src_gf_ptrs_capacity = 0;
    total_elements_chare = 0;
  }
}

// destructor
Interpolator3d::~Interpolator3d() {
  if (owns_griddata_chare && griddata_chare != nullptr) {
    free(griddata_chare);
    griddata_chare = nullptr;
  }
  delete[] src_gf_ptrs;
}

void Interpolator3d::perform_interpolation(int request_type, int request_id, int num_gfs, int total_elements, REAL *dst_x0x1x2_linear) {
  dst_x0x1x2 = (REAL(*)[3])dst_x0x1x2_linear;
  int count_total_elements_chare = 0;
  for (int i = 0; i < total_elements; i++) {
    if ((griddata_chare->params.xxmin0 <= dst_x0x1x2[i][0] && dst_x0x1x2[i][0] <= griddata_chare->params.xxmax0) &&
        (griddata_chare->params.xxmin1 <= dst_x0x1x2[i][1] && dst_x0x1x2[i][1] <= griddata_chare->params.xxmax1) &&
        (griddata_chare->params.xxmin2 <= dst_x0x1x2[i][2] && dst_x0x1x2[i][2] <= griddata_chare->params.xxmax2)) {
      count_total_elements_chare++;
    }
  }
  total_elements_chare = count_total_elements_chare;
  dst_x0x1x2_chare = (REAL(*)[3])malloc(total_elements_chare * 3 * sizeof(REAL));
  dst_indices_chare = (int *)malloc(total_elements_chare * sizeof(int));

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

  interp_request_type = request_type;
  interp_request_id = request_id;
  interp_num_gfs = num_gfs;
  dst_data_ptrs_chare = (REAL **)malloc(interp_num_gfs * sizeof(REAL *));
  for (int i = 0; i < interp_num_gfs; i++) {
    dst_data_ptrs_chare[i] = (REAL *)malloc(total_elements_chare * sizeof(REAL));
  }
  const int Nxx_plus_2NGHOSTS0 = griddata_chare[grid].params.Nxx_plus_2NGHOSTS0;
  const int Nxx_plus_2NGHOSTS1 = griddata_chare[grid].params.Nxx_plus_2NGHOSTS1;
  const int Nxx_plus_2NGHOSTS2 = griddata_chare[grid].params.Nxx_plus_2NGHOSTS2;

  const int interp_err = interpolation_3d_general__uniform_src_grid((NGHOSTS), griddata_chare[grid].params.dxx0, griddata_chare[grid].params.dxx1,
                                                                    griddata_chare[grid].params.dxx2, Nxx_plus_2NGHOSTS0, Nxx_plus_2NGHOSTS1, Nxx_plus_2NGHOSTS2,
                                                                    interp_num_gfs, griddata_chare[grid].xx, src_gf_ptrs, total_elements_chare, dst_x0x1x2_chare,
                                                                    dst_data_ptrs_chare);
  if (interp_err != 0) {
    CkAbort("Error: interpolation_3d_general__uniform_src_grid failed in perform_interpolation.");
  }
  contribute_interpolation_results(request_id);
  free(dst_x0x1x2_chare);
  free(dst_indices_chare);
  for (int i = 0; i < interp_num_gfs; i++) {
    free(dst_data_ptrs_chare[i]);
  }
  free(dst_data_ptrs_chare);
  dst_data_ptrs_chare = nullptr;
}

void Interpolator3d::contribute_interpolation_results(int curr_index_horizonfinder_chare) {
  // We have total_elements_chare number of grid points
  // For each grid point, we send its index in the original dst_x0x1x2 array
  // and interp_num_gfs number of gfs values at that grid point

  // 1) compute how many bytes each point contributes:
  size_t bytes_per_point = sizeof(int) + interp_num_gfs * sizeof(REAL);
  size_t total_bytes = total_elements_chare * bytes_per_point;

  // 2) pack into a raw buffer with memcpy
  char *packbuf = (char*)malloc(total_bytes);
  char *p = packbuf;

  for (int k = 0; k < total_elements_chare; ++k) {
    // copy the integer index
    std::memcpy(p, &dst_indices_chare[k], sizeof(int));
    p += sizeof(int);

    // copy the interp_num_gfs REAL values
    for (int i = 0; i < interp_num_gfs; ++i) {
      std::memcpy(p, &dst_data_ptrs_chare[i][k], sizeof(REAL));
      p += sizeof(REAL);
    }
  }

  {
     InterpBufMsg *m = new (total_bytes) InterpBufMsg();
     m->request_type = interp_request_type;
     m->request_id = curr_index_horizonfinder_chare;
     m->num_gfs = interp_num_gfs;
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
    interp_request_type = m->request_type;
    interp_request_id = m->request_id;
    interp_num_gfs = m->num_gfs;
    interp_bufs = new char*[interp_total];
    interp_lens = new int[interp_total];
  } else if (interp_request_type != m->request_type || interp_request_id != m->request_id || interp_num_gfs != m->num_gfs) {
    CkAbort("Error: Interpolator3d received mismatched interpolation messages in recv_interp_msg.");
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
<<PSI4_SEND_CONCAT_LOG>>
  char *agg=(char*)malloc(tot);
  size_t off=0;
  for(int i=0;i<interp_total;i++){
    memcpy(agg+off, interp_bufs[i], interp_lens[i]);
    off += interp_lens[i];
  }
  InterpBufMsg *out = new (tot) InterpBufMsg();
  out->request_type = interp_request_type;
  out->request_id = interp_request_id;
  out->num_gfs = interp_num_gfs;
  out->len = (int) tot;
  memcpy(out->buf, agg, tot);
  if (interp_request_type == INTERP_REQUEST_BHAHAHA) {
    horizon_finderProxy[CkArrayIndex1D(interp_request_id)].report_interpolation_results(out);
<<PSI4_SEND_CONCAT_BLOCK>>
  } else {
    CkAbort("Error: Unknown interpolation request type in send_interp_concat.");
  }

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
  interp_count=interp_total=0;
  interp_request_type = -1;
  interp_request_id = -1;
  interp_num_gfs = 0;

  thisProxy.interp_concatenation_complete();
}

#include "interpolator3d.def.h"
"""
    file_output_str = file_output_str.replace(
        "<<PSI4_SEND_CONCAT_LOG>>", psi4_send_concat_log
    )
    file_output_str = file_output_str.replace(
        "<<PSI4_SEND_CONCAT_BLOCK>>", psi4_send_concat_block
    )
    interpolator3d_cpp_file = project_Path / "interpolator3d.cpp"
    with interpolator3d_cpp_file.open("w", encoding="utf-8") as file:
        file.write(file_output_str)


def output_interpolation_buffer_utils_cpp(
    project_dir: str,
) -> None:
    """
    Generate interpolation_buffer_utils.cpp.
    :param project_dir: Directory where the project C code is output
    """
    project_Path = Path(project_dir)
    project_Path.mkdir(parents=True, exist_ok=True)

    file_output_str = r"""#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"
#include <cstring>
#include <cstdlib>

int unpack_interpolation_buffer(const int num_gfs, const char *buf, const size_t buf_sz, REAL *dst_data_ptrs[]) {
  if (num_gfs <= 0) {
    return -1;
  }
  const size_t bytes_per_pt = sizeof(int) + (size_t)num_gfs * sizeof(REAL);
  if (bytes_per_pt == 0 || (buf_sz % bytes_per_pt) != 0) {
    return -1;
  }

  const int npts = (int)(buf_sz / bytes_per_pt);
  const char *p = buf;
  for (int k = 0; k < npts; ++k) {
    int idx = 0;
    std::memcpy(&idx, p, sizeof(int));
    p += sizeof(int);
    for (int gf = 0; gf < num_gfs; ++gf) {
      std::memcpy(&dst_data_ptrs[gf][idx], p, sizeof(REAL));
      p += sizeof(REAL);
    }
  }
  return 0;
}
"""
    interpolation_buffer_utils_file = project_Path / "interpolation_buffer_utils.cpp"
    with interpolation_buffer_utils_file.open("w", encoding="utf-8") as file:
        file.write(file_output_str)


def output_interpolator3d_ci(
    project_dir: str,
    enable_psi4: bool = False,
) -> None:
    """
    Generate interpolator3d.ci.

    :param project_dir: Directory where the project C code is output
    """
    project_Path = Path(project_dir)
    project_Path.mkdir(parents=True, exist_ok=True)

    psi4_recv_log = ""
    psi4_start_interp_log = ""
    psi4_request_block = "\n          }\n"
    if enable_psi4:
        psi4_recv_log = ""
        psi4_start_interp_log = ""
        psi4_request_block = r"""
          } else if (request_type == INTERP_REQUEST_PSI4) {
            if (thisIndex.x == 0 && thisIndex.y == 0 && thisIndex.z == 0) {
              serial { psi4_spinweightm2_shell_init(&commondata, &psi4_shell); }
            }
            for (iter = 0; iter < commondata.num_psi4_extraction_radii; iter++) {
              if (thisIndex.x == 0 && thisIndex.y == 0 && thisIndex.z == 0) {
                serial {
                  const REAL R_ext = commondata.list_of_psi4_extraction_radii[iter];
                  const int num_pts = psi4_shell.num_pts;
                  REAL(*dst_pts)[3] = (REAL(*)[3])malloc(sizeof(REAL) * num_pts * 3);
                  psi4r_at_R_ext = (REAL *)calloc(num_pts, sizeof(REAL));
                  psi4i_at_R_ext = (REAL *)calloc(num_pts, sizeof(REAL));
                  psi4_spinweightm2_shell_fill_points(&griddata_chare[grid].params, &psi4_shell, R_ext, dst_pts, NULL);
                  thisProxy.start_interpolation(INTERP_REQUEST_PSI4, iter, num_gfs, num_pts, (REAL *)dst_pts);
                  free(dst_pts);
                }
              }
              when start_interpolation(int request_type, int request_id, int num_gfs, int total_elements,
                                       REAL dst_x0x1x2_linear[3 * total_elements]) {
                serial {
                  perform_interpolation(request_type, request_id, num_gfs, total_elements, dst_x0x1x2_linear);
                }

                // Only continue when custom reduction on interpolator chare 0 is complete
                when interp_concatenation_complete() {
                  serial {}
                }
              }
            } // end for (iter = 0; iter < commondata.num_psi4_extraction_radii; iter++)
            if (thisIndex.x == 0 && thisIndex.y == 0 && thisIndex.z == 0) {
              serial { psi4_spinweightm2_shell_free(&psi4_shell); }
            }
"""
        psi4_request_block += "\n          }\n"

    file_output_str = """
module interpolator3d {
  include "BHaH_defines.h";
  include "BHaH_function_prototypes.h";
  include "commondata_object.h";
  include "griddata_object.h";
  include "ckio.h";
  include "pup_stl.h";

  message InterpBufMsg {
    int request_type;
    int request_id;
    int num_gfs;
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
        when receiv_interp_gfs(int request_type, int num_gfs, int nn, int len_tmpBuffer, REAL tmpBuffer[len_tmpBuffer]) {
          serial {
            commondata.nn = nn;
            commondata.time = (REAL)nn * commondata.dt;
            if (request_type == INTERP_REQUEST_PSI4 && thisIndex.x == 0 && thisIndex.y == 0 && thisIndex.z == 0) {
              CkPrintf("Interpolator3d: received PSI4 gfs at nn=%d (t=%e), num_gfs=%d\\n", nn, (double)commondata.time, num_gfs);
            }
            const int Nxx_plus_2NGHOSTS0 = griddata_chare[grid].params.Nxx_plus_2NGHOSTS0;
            const int Nxx_plus_2NGHOSTS1 = griddata_chare[grid].params.Nxx_plus_2NGHOSTS1;
            const int Nxx_plus_2NGHOSTS2 = griddata_chare[grid].params.Nxx_plus_2NGHOSTS2;
            const int Nxx_plus_2NGHOSTS_tot = Nxx_plus_2NGHOSTS0 * Nxx_plus_2NGHOSTS1 * Nxx_plus_2NGHOSTS2;
            interp_num_gfs = num_gfs;
            interp_request_type = request_type;
            if (src_gf_ptrs == nullptr) {
              src_gf_ptrs = new const REAL *[interp_num_gfs];
              src_gf_ptrs_capacity = interp_num_gfs;
            } else if (src_gf_ptrs_capacity < interp_num_gfs) {
              delete[] src_gf_ptrs;
              src_gf_ptrs = new const REAL *[interp_num_gfs];
              src_gf_ptrs_capacity = interp_num_gfs;
            }
            for (int idx = 0; idx < interp_num_gfs; idx++) {
              src_gf_ptrs[idx] = tmpBuffer + idx * Nxx_plus_2NGHOSTS_tot;
            }
          }

          if (request_type == INTERP_REQUEST_BHAHAHA) {
            if (thisIndex.x == 0 && thisIndex.y == 0 && thisIndex.z == 0) {
              serial { horizon_finderProxy.ready_for_interpolation(); }
            }
            for (iter = 0; iter < commondata.bah_max_num_horizons; iter++) {
              when start_interpolation(int request_type, int request_id, int num_gfs, int total_elements,
                                       REAL dst_x0x1x2_linear[3 * total_elements]) {
                serial {
                  if (request_type == INTERP_REQUEST_PSI4 && thisIndex.x == 0 && thisIndex.y == 0 && thisIndex.z == 0) {
                    CkPrintf("Interpolator3d: start_interpolation PSI4 request_id=%d total_elements=%d\\n", request_id, total_elements);
                  }
                  perform_interpolation(request_type, request_id, num_gfs, total_elements, dst_x0x1x2_linear);
                }

                // Only continue when custom reduction on interpolator chare 0 is complete
                when interp_concatenation_complete() {
                  serial {}
                }
              }
            } // end for (iter = 0; iter < commondata.bah_max_num_horizons; iter++)
{psi4_request_block}
        } // end when receiv_interp_gfs
      } // end while (commondata.time < commondata.t_final)
    };
    entry void start_interpolation(int request_type, int request_id, int num_gfs, int total_elements, REAL dst_x0x1x2_linear[3 * total_elements]);
    entry void receiv_interp_gfs(int request_type, int num_gfs, int nn, int len_tmpBuffer, REAL tmpBuffer[len_tmpBuffer]);
    entry void recv_interp_msg(InterpBufMsg * m);
    entry void send_interp_concat();
    entry void interp_concatenation_complete();
  };
};
"""
    file_output_str = file_output_str.replace(
        """            if (request_type == INTERP_REQUEST_PSI4 && thisIndex.x == 0 && thisIndex.y == 0 && thisIndex.z == 0) {
              CkPrintf("Interpolator3d: received PSI4 gfs at nn=%d (t=%e), num_gfs=%d\\n", nn, (double)commondata.time, num_gfs);
            }
""",
        psi4_recv_log,
    )
    file_output_str = file_output_str.replace(
        """                  if (request_type == INTERP_REQUEST_PSI4 && thisIndex.x == 0 && thisIndex.y == 0 && thisIndex.z == 0) {
                    CkPrintf("Interpolator3d: start_interpolation PSI4 request_id=%d total_elements=%d\\n", request_id, total_elements);
                  }
""",
        psi4_start_interp_log,
    )
    file_output_str = file_output_str.replace(
        "{psi4_request_block}", psi4_request_block
    )
    interpolator3d_ci_file = project_Path / "interpolator3d.ci"
    with interpolator3d_ci_file.open("w", encoding="utf-8") as file:
        file.write(file_output_str)


def output_interpolator3d_h_cpp_ci(
    project_dir: str,
    enable_psi4: bool = False,
) -> None:
    """
    Generate interpolator3d.h, interpolator3d.cpp and interpolator3d.ci.
    :param project_dir: Directory where the project C code is output.
    """
    output_griddata_object_h(
        project_dir=project_dir,
    )

    output_interp_buf_msg_h(
        project_dir=project_dir,
    )

    output_interpolator3d_h(
        project_dir=project_dir,
        enable_psi4=enable_psi4,
    )

    output_interpolator3d_cpp(
        project_dir=project_dir,
        enable_psi4=enable_psi4,
    )

    output_interpolation_buffer_utils_cpp(
        project_dir=project_dir,
    )

    output_interpolator3d_ci(
        project_dir=project_dir,
        enable_psi4=enable_psi4,
    )
