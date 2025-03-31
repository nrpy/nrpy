"""
Register the superB PUP routines for Charm++ checkpointing and load balancing.

This script registers the C function "superB_pup_routines" from the file
"superB_pup_routines.cpp", which contains various Pack-Unpack (PUP) routines for
structs used in checkpointing and load balancing in Charm++. These routines handle
the packing and unpacking of data structures such as commondata, params, rfm, boundary
conditions, MoL grid functions, and more.

Author: Nishita Jadoo
        njadoo **at** uidaho **dot* edu
"""

from typing import List

import nrpy.c_function as cfc
import nrpy.params as par
from nrpy.infrastructures.BHaH.MoLtimestepping.MoL_gridfunction_names import (
    generate_gridfunction_names,
)
from nrpy.infrastructures.BHaH.MoLtimestepping.RK_Butcher_Table_Dictionary import (
    generate_Butcher_tables,
)


def register_CFunction_superB_pup_routines(
    MoL_method: str = "RK4",
    enable_psi4_diagnostics: bool = False,
) -> None:
    """
    Register the C function "superB_pup_routines", which is a collection of Pack-Unpack (PUP) routines for various structs.
    These routines are used for checkpointing and load balancing in Charm++.

    :param MoL_method: The Method of Lines (MoL) method to be used (default is "RK4").
    :param enable_psi4_diagnostics: Flag to enable psi4 diagnostics.

    DocTests:
        >>> register_CFunction_superB_pup_routines()
    """
    desc = """This file implements a collection of Pack-Unpack (PUP) routines used in Charm++ for checkpointing,
and load balancing in the superB framework.
It includes routines for serializing and deserializing:
    - commondata_struct and params_struct,
    - rfm_struct with reference metric precomputation and memory allocation,
    - inner and outer boundary condition structures (innerpt_bc_struct, outerpt_bc_struct, bc_info_struct, bc_struct),
    - MoL grid functions (MoL_gridfunctions_struct),
    - chare communication structures (charecomm_struct),
    - diagnostic information (diagnostic_struct),
    - temporary buffers and nonlocal inner boundary conditions (tmpBuffers_struct, nonlocalinnerbc_struct),
    - and grid data structures (griddata_struct, griddata_chare).
This comprehensive set of routines is crucial for efficient data management and communication in high-performance, parallel simulations.
"""
    # prefunc contains most of the C++ source code for the PUP routines.
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    prefunc = """
// PUP routine for struct commondata_struct
void pup_commondata_struct(PUP::er &p, commondata_struct &commondata) {
"""
    struct_list: List[str] = (
        []
    )  # List to store individual struct elements for commondata_struct.
    for parname, CodeParam in par.glb_code_params_dict.items():
        if CodeParam.commondata:
            struct = "commondata"
            CPtype = CodeParam.cparam_type
            comment = f"  // {CodeParam.module}::{parname}"
            if "char" in CPtype and "[" in CPtype and "]" in CPtype:
                chararray_size = CPtype.split("[")[1].replace("]", "")
                c_output = (
                    f"PUParray(p, {struct}.{parname}, {chararray_size});{comment}\n"
                )
            elif "TIMEVAR" in CPtype:
                c_output = f"p|{struct}.{parname}.tv_sec;{comment}\n"
                c_output += f"p|{struct}.{parname}.tv_nsec;{comment}\n"
            else:
                c_output = f"p|{struct}.{parname};{comment}\n"
            struct_list.append(c_output)
    # Sort the lines alphabetically and join them with line breaks.
    prefunc += "// PUP commondata struct\n"
    prefunc += "".join(sorted(struct_list))
    prefunc += """
}"""

    prefunc += """
// PUP routine for struct params_struct
void pup_params_struct(PUP::er &p, params_struct &params) {
"""
    params_struct_list: List[str] = (
        []
    )  # List to store individual struct elements for params_struct.
    for parname, CodeParam in par.glb_code_params_dict.items():
        CPtype = CodeParam.cparam_type
        if not CodeParam.commondata and CPtype != "#define":
            struct = "params"
            comment = f"  // {CodeParam.module}::{parname}"
            if "char" in CPtype and "[" in CPtype and "]" in CPtype:
                chararray_size = CPtype.split("[")[1].replace("]", "")
                c_output = (
                    f"PUParray(p, {struct}.{parname}, {chararray_size});{comment}\n"
                )
            elif "TIMEVAR" in CPtype:
                c_output = f"p|{struct}.{parname}.tv_sec;{comment}\n"
                c_output += f"p|{struct}.{parname}.tv_nsec;{comment}\n"
            else:
                c_output = f"p|{struct}.{parname};{comment}\n"
            params_struct_list.append(c_output)
    # Sort the lines alphabetically and join them with line breaks.
    prefunc += "// PUP params struct\n"
    prefunc += "".join(sorted(params_struct_list))
    prefunc += """
}"""

    prefunc += """
// PUP routine for struct innerpt_bc_struct
void pup_innerpt_bc_struct(PUP::er &p, innerpt_bc_struct &ibc) {
  p | ibc.dstpt;
  p | ibc.srcpt;
  for (int i = 0; i < 10; i++) {
    p | ibc.parity[i];
  }
}

// PUP routine for struct outerpt_bc_struct
void pup_outerpt_bc_struct(PUP::er &p, outerpt_bc_struct &obc) {
  p | obc.i0;
  p | obc.i1;
  p | obc.i2;
  p | obc.FACEX0;
  p | obc.FACEX1;
  p | obc.FACEX2;
}

// PUP routine for struct bc_info_struct
void pup_bc_info_struct(PUP::er &p, bc_info_struct &bci) {
  p | bci.num_inner_boundary_points;
  p | bci.num_inner_boundary_points_nonlocal;
  for (int i = 0; i < NGHOSTS; i++) {
    for (int j = 0; j < 3; j++) {
      p | bci.num_pure_outer_boundary_points[i][j];
    }
  }
  for (int i = 0; i < NGHOSTS; i++) {
    for (int j = 0; j < 6; j++) {
      for (int k = 0; k < 6; k++) {
        p | bci.bc_loop_bounds[i][j][k];
      }
    }
  }
}

// PUP routine for struct bc_struct
void pup_bc_struct(PUP::er &p, bc_struct &bc) {

  pup_bc_info_struct(p, bc.bc_info);

  if (p.isUnpacking()) {

    bc.inner_bc_array = (innerpt_bc_struct *restrict)malloc(sizeof(innerpt_bc_struct) * bc.bc_info.num_inner_boundary_points);
    bc.inner_bc_array_nonlocal = (innerpt_bc_struct *restrict)malloc(sizeof(innerpt_bc_struct) * bc.bc_info.num_inner_boundary_points_nonlocal);

    for (int which_gz = 0; which_gz < NGHOSTS; which_gz++) {
      for (int dirn = 0; dirn < 3; dirn++) {
        bc.pure_outer_bc_array[dirn + (3 * which_gz)] =
            (outerpt_bc_struct *restrict)malloc(bc.bc_info.num_pure_outer_boundary_points[which_gz][dirn] * sizeof(outerpt_bc_struct));
      }
    }
  }

  for (int i = 0; i < bc.bc_info.num_inner_boundary_points; i++) {
    pup_innerpt_bc_struct(p, bc.inner_bc_array[i]);
  }

  for (int i = 0; i < bc.bc_info.num_inner_boundary_points_nonlocal; i++) {
    pup_innerpt_bc_struct(p, bc.inner_bc_array_nonlocal[i]);
  }

  for (int which_gz = 0; which_gz < NGHOSTS; which_gz++) {
    for (int dirn = 0; dirn < 3; dirn++) {
      int num_points = bc.bc_info.num_pure_outer_boundary_points[which_gz][dirn];
      for (int pt = 0; pt < num_points; pt++) {
          pup_outerpt_bc_struct(p, bc.pure_outer_bc_array[dirn + (3 * which_gz)][pt]);
      }
    }
  }
}
"""

    prefunc += """
// PUP routine for struct MoL_gridfunctions_struct
void pup_MoL_gridfunctions_struct(PUP::er &p, MoL_gridfunctions_struct &gridfuncs, const params_struct &params, const commondata_struct &commondata) {"""

    prefunc += r"""
  p | gridfuncs.num_evol_gfs_to_sync;
  p | gridfuncs.num_auxevol_gfs_to_sync;
  p | gridfuncs.num_aux_gfs_to_sync;
  p | gridfuncs.max_sync_gfs;
  PUParray(p, gridfuncs.evol_gfs_to_sync, gridfuncs.num_evol_gfs_to_sync);
  PUParray(p, gridfuncs.auxevol_gfs_to_sync, gridfuncs.num_auxevol_gfs_to_sync);
  PUParray(p, gridfuncs.aux_gfs_to_sync, gridfuncs.num_aux_gfs_to_sync);
"""

    prefunc += """
  const int Nxx_plus_2NGHOSTS_tot = params.Nxx_plus_2NGHOSTS0 * params.Nxx_plus_2NGHOSTS1 * params.Nxx_plus_2NGHOSTS2;
  if (p.isUnpacking()) {
"""
    Butcher_dict = generate_Butcher_tables()
    (
        y_n_gridfunctions,
        non_y_n_gridfunctions_list,
        _,
        diagnostic_gridfunctions2_point_to,
    ) = generate_gridfunction_names(Butcher_dict, MoL_method=MoL_method)
    # Combine y_n_gfs and non_y_n_gfs into a single list.
    gridfunctions_list = [y_n_gridfunctions] + non_y_n_gridfunctions_list
    for gridfunctions in gridfunctions_list:
        num_gfs = (
            "NUM_EVOL_GFS" if gridfunctions != "auxevol_gfs" else "NUM_AUXEVOL_GFS"
        )
        # Do not allocate a zero-sized array.
        if num_gfs == "NUM_AUXEVOL_GFS":
            prefunc += "  if(NUM_AUXEVOL_GFS > 0) "
        prefunc += (
            f"gridfuncs.{gridfunctions} = (REAL *restrict)malloc(sizeof(REAL) * {num_gfs} * "
            "Nxx_plus_2NGHOSTS_tot);\n"
        )
    # In superB, allocate separate memory for diagnostic_output_gfs.
    prefunc += "gridfuncs.diagnostic_output_gfs  = (REAL *restrict)malloc(sizeof(REAL) * NUM_EVOL_GFS * Nxx_plus_2NGHOSTS_tot);\n"
    prefunc += f"gridfuncs.diagnostic_output_gfs2 = gridfuncs.{diagnostic_gridfunctions2_point_to};\n"
    prefunc += """
    initialize_yn_and_non_yn_gfs_to_nan(&commondata, &params, &gridfuncs);
  }"""
    for gridfunctions in gridfunctions_list:
        num_gfs = (
            "NUM_EVOL_GFS" if gridfunctions != "auxevol_gfs" else "NUM_AUXEVOL_GFS"
        )
        if gridfunctions == "y_n_gfs":
            prefunc += (
                f"PUParray(p, gridfuncs.y_n_gfs, {num_gfs} * Nxx_plus_2NGHOSTS_tot);\n"
            )
        elif gridfunctions == "auxevol_gfs":
            prefunc += rf"""
  if (strstr(params.CoordSystemName, "Spherical") != NULL) {{
    if(NUM_AUXEVOL_GFS > 0) {{
      PUParray(p, gridfuncs.auxevol_gfs, {num_gfs} * Nxx_plus_2NGHOSTS_tot);
    }}
  }}"""
    prefunc += """
}
"""

    prefunc += """
// PUP routine for struct charecomm_struct
void pup_charecomm_struct(PUP::er &p, charecomm_struct &cc, const params_struct &params, const params_struct &params_chare) {
  const int ntot = params.Nxx_plus_2NGHOSTS0 * params.Nxx_plus_2NGHOSTS1 * params.Nxx_plus_2NGHOSTS2;
  const int ntotchare = params_chare.Nxx_plus_2NGHOSTS0 * params_chare.Nxx_plus_2NGHOSTS1 * params_chare.Nxx_plus_2NGHOSTS2;
  if (p.isUnpacking()) {
    cc.globalidx3pt_to_chareidx3 = (int *restrict)malloc(sizeof(int) * ntot);
    cc.globalidx3pt_to_localidx3pt = (int *restrict)malloc(sizeof(int) * ntot);
    cc.localidx3pt_to_globalidx3pt = (int *restrict)malloc(sizeof(int) * ntotchare);
  }
  PUParray(p, cc.globalidx3pt_to_chareidx3, ntot);
  PUParray(p, cc.globalidx3pt_to_localidx3pt, ntot);
  PUParray(p, cc.localidx3pt_to_globalidx3pt, ntotchare);
}"""

    prefunc += """
// PUP routine for struct diagnostic_struct
void pup_diagnostic_struct(PUP::er &p, diagnostic_struct &ds, const params_struct &params_chare) {
  p | ds.num_output_quantities;
  p | ds.tot_num_diagnostic_1d_y_pts;
  p | ds.tot_num_diagnostic_1d_z_pts;
  p | ds.tot_num_diagnostic_2d_xy_pts;
  p | ds.tot_num_diagnostic_2d_yz_pts;
  p | ds.num_diagnostic_1d_y_pts;
  p | ds.num_diagnostic_1d_z_pts;
  p | ds.num_diagnostic_2d_xy_pts;
  p | ds.num_diagnostic_2d_yz_pts;

  if (p.isUnpacking()) {
    ds.localidx3_diagnostic_1d_y_pt = (int *restrict)malloc(sizeof(int) * ds.num_diagnostic_1d_y_pts);
    ds.locali0_diagnostic_1d_y_pt = (int *restrict)malloc(sizeof(int) * ds.num_diagnostic_1d_y_pts);
    ds.locali1_diagnostic_1d_y_pt = (int *restrict)malloc(sizeof(int) * ds.num_diagnostic_1d_y_pts);
    ds.locali2_diagnostic_1d_y_pt = (int *restrict)malloc(sizeof(int) * ds.num_diagnostic_1d_y_pts);
    ds.offset_diagnostic_1d_y_pt = (int *restrict)malloc(sizeof(int) * ds.num_diagnostic_1d_y_pts);

    ds.localidx3_diagnostic_1d_z_pt = (int *restrict)malloc(sizeof(int) * ds.num_diagnostic_1d_z_pts);
    ds.locali0_diagnostic_1d_z_pt = (int *restrict)malloc(sizeof(int) * ds.num_diagnostic_1d_z_pts);
    ds.locali1_diagnostic_1d_z_pt = (int *restrict)malloc(sizeof(int) * ds.num_diagnostic_1d_z_pts);
    ds.locali2_diagnostic_1d_z_pt = (int *restrict)malloc(sizeof(int) * ds.num_diagnostic_1d_z_pts);
    ds.offset_diagnostic_1d_z_pt = (int *restrict)malloc(sizeof(int) * ds.num_diagnostic_1d_z_pts);

    ds.localidx3_diagnostic_2d_xy_pt = (int *restrict)malloc(sizeof(int) * ds.num_diagnostic_2d_xy_pts);
    ds.locali0_diagnostic_2d_xy_pt = (int *restrict)malloc(sizeof(int) * ds.num_diagnostic_2d_xy_pts);
    ds.locali1_diagnostic_2d_xy_pt = (int *restrict)malloc(sizeof(int) * ds.num_diagnostic_2d_xy_pts);
    ds.locali2_diagnostic_2d_xy_pt = (int *restrict)malloc(sizeof(int) * ds.num_diagnostic_2d_xy_pts);
    ds.offset_diagnostic_2d_xy_pt = (int *restrict)malloc(sizeof(int) * ds.num_diagnostic_2d_xy_pts);

    ds.localidx3_diagnostic_2d_yz_pt = (int *restrict)malloc(sizeof(int) * ds.num_diagnostic_2d_yz_pts);
    ds.locali0_diagnostic_2d_yz_pt = (int *restrict)malloc(sizeof(int) * ds.num_diagnostic_2d_yz_pts);
    ds.locali1_diagnostic_2d_yz_pt = (int *restrict)malloc(sizeof(int) * ds.num_diagnostic_2d_yz_pts);
    ds.locali2_diagnostic_2d_yz_pt = (int *restrict)malloc(sizeof(int) * ds.num_diagnostic_2d_yz_pts);
    ds.offset_diagnostic_2d_yz_pt = (int *restrict)malloc(sizeof(int) * ds.num_diagnostic_2d_yz_pts);
  }

  PUParray(p, ds.localidx3_diagnostic_1d_y_pt, ds.num_diagnostic_1d_y_pts);
  PUParray(p, ds.locali0_diagnostic_1d_y_pt, ds.num_diagnostic_1d_y_pts);
  PUParray(p, ds.locali1_diagnostic_1d_y_pt, ds.num_diagnostic_1d_y_pts);
  PUParray(p, ds.locali2_diagnostic_1d_y_pt, ds.num_diagnostic_1d_y_pts);
  PUParray(p, ds.offset_diagnostic_1d_y_pt, ds.num_diagnostic_1d_y_pts);

  PUParray(p, ds.localidx3_diagnostic_1d_z_pt, ds.num_diagnostic_1d_z_pts);
  PUParray(p, ds.locali0_diagnostic_1d_z_pt, ds.num_diagnostic_1d_z_pts);
  PUParray(p, ds.locali1_diagnostic_1d_z_pt, ds.num_diagnostic_1d_z_pts);
  PUParray(p, ds.locali2_diagnostic_1d_z_pt, ds.num_diagnostic_1d_z_pts);
  PUParray(p, ds.offset_diagnostic_1d_z_pt, ds.num_diagnostic_1d_z_pts);

  PUParray(p, ds.localidx3_diagnostic_2d_xy_pt, ds.num_diagnostic_2d_xy_pts);
  PUParray(p, ds.locali0_diagnostic_2d_xy_pt, ds.num_diagnostic_2d_xy_pts);
  PUParray(p, ds.locali1_diagnostic_2d_xy_pt, ds.num_diagnostic_2d_xy_pts);
  PUParray(p, ds.locali2_diagnostic_2d_xy_pt, ds.num_diagnostic_2d_xy_pts);
  PUParray(p, ds.offset_diagnostic_2d_xy_pt, ds.num_diagnostic_2d_xy_pts);

  PUParray(p, ds.localidx3_diagnostic_2d_yz_pt, ds.num_diagnostic_2d_yz_pts);
  PUParray(p, ds.locali0_diagnostic_2d_yz_pt, ds.num_diagnostic_2d_yz_pts);
  PUParray(p, ds.locali1_diagnostic_2d_yz_pt, ds.num_diagnostic_2d_yz_pts);
  PUParray(p, ds.locali2_diagnostic_2d_yz_pt, ds.num_diagnostic_2d_yz_pts);
  PUParray(p, ds.offset_diagnostic_2d_yz_pt, ds.num_diagnostic_2d_yz_pts);

  PUParray(p, ds.filename_1d_y, 256);
  PUParray(p, ds.filename_1d_z, 256);
  PUParray(p, ds.filename_2d_xy, 256);
  PUParray(p, ds.filename_2d_yz, 256);
"""

    if enable_psi4_diagnostics:
        prefunc += r"""
  p | ds.num_of_R_exts_chare;
  p | ds.psi4_spinweightm2_sph_harmonics_max_l;
  p | ds.length_localsums_for_psi4_decomp;

  if (p.isUnpacking()) {
    ds.list_of_R_exts_chare = (REAL *restrict)malloc(sizeof(REAL) * ds.num_of_R_exts_chare);
    ds.localsums_for_psi4_decomp = (REAL *restrict)malloc(sizeof(REAL) * ds.length_localsums_for_psi4_decomp);
    ds.globalsums_for_psi4_decomp = (REAL *restrict)malloc(sizeof(REAL) * ds.length_localsums_for_psi4_decomp);
  }
  PUParray(p, ds.list_of_R_exts_chare, ds.num_of_R_exts_chare);

  if (strstr(params_chare.CoordSystemName, "Cylindrical") != NULL) {
    p | ds.tot_N_shell_pts_chare;
    p | ds.dtheta;
    if (p.isUnpacking()) {
      ds.N_shell_pts_chare = (int *restrict)malloc(sizeof(int) * ds.num_of_R_exts_chare);
      ds.N_theta_shell_chare = (int *restrict)malloc(sizeof(int) * ds.num_of_R_exts_chare);
      ds.xx_shell_chare = (REAL ***restrict)malloc(ds.num_of_R_exts_chare * sizeof(REAL **));
      ds.theta_shell_chare = (REAL **restrict)malloc(ds.num_of_R_exts_chare * sizeof(REAL *));
    }
    PUParray(p, ds.N_shell_pts_chare, ds.num_of_R_exts_chare);
    PUParray(p, ds.N_theta_shell_chare, ds.num_of_R_exts_chare);
    if (p.isUnpacking()) {
      ds.theta_shell_chare = (REAL * *restrict)malloc(sizeof(REAL *) * ds.num_of_R_exts_chare);
      for (int i = 0; i < ds.num_of_R_exts_chare; i++) {
        ds.theta_shell_chare[i] = (REAL *restrict)malloc(sizeof(REAL) * ds.N_theta_shell_chare[i]);
      }
      ds.xx_shell_chare = (REAL * **restrict)malloc(sizeof(REAL **) * ds.num_of_R_exts_chare);
      for (int i = 0; i < ds.num_of_R_exts_chare; i++) {
        ds.xx_shell_chare[i] = (REAL * *restrict)malloc(sizeof(REAL *) * ds.N_shell_pts_chare[i]);
        for (int j = 0; j < ds.N_shell_pts_chare[i]; j++) {
          ds.xx_shell_chare[i][j] = (REAL *restrict)malloc(sizeof(REAL) * 3);
        }
      }
    }
    for (int i = 0; i < ds.num_of_R_exts_chare; i++) {
      for (int j = 0; j < ds.N_shell_pts_chare[i]; j++) {
        PUParray(p, ds.xx_shell_chare[i][j], 3);
      }
    }
     for (int i = 0; i < ds.num_of_R_exts_chare; i++) {
      PUParray(p, ds.theta_shell_chare[i], ds.N_theta_shell_chare[i]);
    }
  }"""
    prefunc += r"""
}"""

    prefunc += """
// PUP routine for struct tmpBuffers_struct
void pup_tmpBuffers_struct(PUP::er &p, tmpBuffers_struct &tmpBuffers, const params_struct &params, const nonlocalinnerbc_struct &nonlocalinnerbc, const MoL_gridfunctions_struct &gridfuncs) {
  const int Nxx_plus_2NGHOSTS_face0 = params.Nxx_plus_2NGHOSTS1 * params.Nxx_plus_2NGHOSTS2;
  const int Nxx_plus_2NGHOSTS_face1 = params.Nxx_plus_2NGHOSTS0 * params.Nxx_plus_2NGHOSTS2;
  const int Nxx_plus_2NGHOSTS_face2 = params.Nxx_plus_2NGHOSTS0 * params.Nxx_plus_2NGHOSTS1;
  const int max_sync_gfs = gridfuncs.max_sync_gfs;
  size_t size_EW = static_cast<size_t>(max_sync_gfs) * NGHOSTS * Nxx_plus_2NGHOSTS_face0;
  size_t size_NS = static_cast<size_t>(max_sync_gfs) * NGHOSTS * Nxx_plus_2NGHOSTS_face1;
  size_t size_TB = static_cast<size_t>(max_sync_gfs) * NGHOSTS * Nxx_plus_2NGHOSTS_face2;
  if (p.isUnpacking()) {
    tmpBuffers.tmpBuffer_EW = (REAL *restrict)malloc(sizeof(REAL) * size_EW);
    tmpBuffers.tmpBuffer_NS = (REAL *restrict)malloc(sizeof(REAL) * size_NS);
    tmpBuffers.tmpBuffer_TB = (REAL *restrict)malloc(sizeof(REAL) * size_TB);
  }
  const int tot_num_dst_chares = nonlocalinnerbc.tot_num_dst_chares;
  const int tot_num_src_chares = nonlocalinnerbc.tot_num_src_chares;
  const int *num_srcpts_tosend_each_chare = nonlocalinnerbc.num_srcpts_tosend_each_chare;
  const int *num_srcpts_each_chare = nonlocalinnerbc.num_srcpts_each_chare;
  if (p.isUnpacking()) {
    tmpBuffers.tmpBuffer_innerbc_send = (REAL **)malloc(tot_num_dst_chares * sizeof(REAL *));
    for (int which_chare = 0; which_chare < tot_num_dst_chares; which_chare++) {
      tmpBuffers.tmpBuffer_innerbc_send[which_chare] = (REAL *restrict)malloc(sizeof(REAL) * max_sync_gfs * num_srcpts_tosend_each_chare[which_chare]);
    }
    tmpBuffers.tmpBuffer_innerbc_receiv = (REAL **)malloc(tot_num_src_chares * sizeof(REAL *));
    for (int which_chare = 0; which_chare < tot_num_src_chares; which_chare++) {
      tmpBuffers.tmpBuffer_innerbc_receiv[which_chare] = (REAL *restrict)malloc(sizeof(REAL) * max_sync_gfs * num_srcpts_each_chare[which_chare]);
    }
  }
}

// PUP routine for struct nonlocalinnerbc_struct
void pup_nonlocalinnerbc_struct(PUP::er &p, nonlocalinnerbc_struct &nonlocal, const commondata_struct &commondata) {
  const int Nchare0 = commondata.Nchare0;
  const int Nchare1 = commondata.Nchare1;
  const int Nchare2 = commondata.Nchare2;
  const int tot_num_chares = Nchare0 * Nchare1 * Nchare2;

  p | nonlocal.tot_num_src_chares;
  p | nonlocal.tot_num_dst_chares;
  if (p.isUnpacking()) {
    nonlocal.idx3_of_src_chares = (int *restrict)malloc(sizeof(int) * nonlocal.tot_num_src_chares);
    nonlocal.idx3chare_to_src_chare_id = (int *restrict)malloc(sizeof(int) * tot_num_chares);
    nonlocal.num_srcpts_each_chare = (int *restrict)malloc(sizeof(int) * nonlocal.tot_num_src_chares);

    nonlocal.idx3_of_dst_chares = (int *restrict)malloc(sizeof(int) * nonlocal.tot_num_dst_chares);
    nonlocal.idx3chare_to_dst_chare_id = (int *restrict)malloc(sizeof(int) * nonlocal.tot_num_dst_chares);
    nonlocal.num_srcpts_tosend_each_chare = (int *restrict)malloc(sizeof(int) * nonlocal.tot_num_dst_chares);

    nonlocal.map_srcchare_and_srcpt_id_to_linear_id = (int **)malloc(sizeof(int *) * nonlocal.tot_num_src_chares);
    nonlocal.globalidx3_srcpts = (int **)malloc(sizeof(int *) * nonlocal.tot_num_src_chares);

    nonlocal.globalidx3_srcpts_tosend = (int **)malloc(sizeof(int *) * nonlocal.tot_num_dst_chares);
  }
  PUParray(p, nonlocal.idx3_of_src_chares, nonlocal.tot_num_src_chares);
  PUParray(p, nonlocal.idx3chare_to_src_chare_id, tot_num_chares);
  PUParray(p, nonlocal.num_srcpts_each_chare, nonlocal.tot_num_src_chares);

  PUParray(p, nonlocal.idx3_of_dst_chares, nonlocal.tot_num_dst_chares);
  PUParray(p, nonlocal.idx3chare_to_dst_chare_id, nonlocal.tot_num_dst_chares);
  PUParray(p, nonlocal.num_srcpts_tosend_each_chare, nonlocal.tot_num_dst_chares);

  for (int src_chare = 0; src_chare < nonlocal.tot_num_src_chares; src_chare++) {
    if (p.isUnpacking()) {
      nonlocal.map_srcchare_and_srcpt_id_to_linear_id[src_chare] =
          (int *restrict)malloc(sizeof(int) * nonlocal.num_srcpts_each_chare[src_chare]);
      nonlocal.globalidx3_srcpts[src_chare] =
          (int *restrict)malloc(sizeof(int) * nonlocal.num_srcpts_each_chare[src_chare]);
    }
    PUParray(p, nonlocal.map_srcchare_and_srcpt_id_to_linear_id[src_chare], nonlocal.num_srcpts_each_chare[src_chare]);
    PUParray(p, nonlocal.globalidx3_srcpts[src_chare], nonlocal.num_srcpts_each_chare[src_chare]);
  }
  for (int dst_chare = 0; dst_chare < nonlocal.tot_num_dst_chares; dst_chare++) {
    if (p.isUnpacking()) {
      nonlocal.globalidx3_srcpts_tosend[dst_chare] =
          (int *restrict)malloc(sizeof(int) * nonlocal.num_srcpts_tosend_each_chare[dst_chare]);
    }
    PUParray(p, nonlocal.globalidx3_srcpts_tosend[dst_chare], nonlocal.num_srcpts_tosend_each_chare[dst_chare]);
  }
}

// PUP routine for struct griddata
// During time evolution, need params from griddata which is used to unpack charecomm_struct in griddata_chare and xx for diagnostics.
void pup_griddata(PUP::er &p, griddata_struct &gd) {
  pup_params_struct(p, gd.params);
  if (p.isUnpacking()) {
    gd.xx[0] = (REAL *restrict)malloc(sizeof(REAL) * gd.params.Nxx_plus_2NGHOSTS0);
    gd.xx[1] = (REAL *restrict)malloc(sizeof(REAL) * gd.params.Nxx_plus_2NGHOSTS1);
    gd.xx[2] = (REAL *restrict)malloc(sizeof(REAL) * gd.params.Nxx_plus_2NGHOSTS2);
  }
  PUParray(p, gd.xx[0], gd.params.Nxx_plus_2NGHOSTS0);
  PUParray(p, gd.xx[1], gd.params.Nxx_plus_2NGHOSTS1);
  PUParray(p, gd.xx[2], gd.params.Nxx_plus_2NGHOSTS2);
}

// PUP routine for struct griddata_chare
// For unpacking order is important; unpacked structs are used for unpacking the subsequent structs.
void pup_griddata_chare(PUP::er &p, griddata_struct &gd, const params_struct &params, const commondata_struct &commondata) {

  pup_params_struct(p, gd.params);

  if (p.isUnpacking()) {
    gd.xx[0] = (REAL *restrict)malloc(sizeof(REAL) * gd.params.Nxx_plus_2NGHOSTS0);
    gd.xx[1] = (REAL *restrict)malloc(sizeof(REAL) * gd.params.Nxx_plus_2NGHOSTS1);
    gd.xx[2] = (REAL *restrict)malloc(sizeof(REAL) * gd.params.Nxx_plus_2NGHOSTS2);
  }
  PUParray(p, gd.xx[0], gd.params.Nxx_plus_2NGHOSTS0);
  PUParray(p, gd.xx[1], gd.params.Nxx_plus_2NGHOSTS1);
  PUParray(p, gd.xx[2], gd.params.Nxx_plus_2NGHOSTS2);

  pup_diagnostic_struct(p, gd.diagnosticstruct, gd.params);

  pup_charecomm_struct(p, gd.charecommstruct, params, gd.params);

  pup_bc_struct(p, gd.bcstruct);

  pup_nonlocalinnerbc_struct(p, gd.nonlocalinnerbcstruct, commondata);

  pup_MoL_gridfunctions_struct(p, gd.gridfuncs, gd.params, commondata);

  pup_tmpBuffers_struct(p, gd.tmpBuffers, gd.params, gd.nonlocalinnerbcstruct, gd.gridfuncs);

  if (p.isUnpacking()) {
    gd.rfmstruct = (rfm_struct *)malloc(sizeof(rfm_struct));
    rfm_precompute_malloc(&commondata, &gd.params, gd.rfmstruct);
    rfm_precompute_defines(&commondata, &gd.params, gd.rfmstruct, gd.xx);
  }
}
"""

    name = "superB_pup_routines"
    params = """ """
    body = r"""
  // This space intentionally left blank:
  // NRPy+ requires C files share the same name as the primary function within the C file.
"""
    cfc.register_CFunction(
        subdirectory="superB",
        includes=includes,
        prefunc=prefunc,
        desc=desc,
        name=name,
        params=params,
        body=body,
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
