"""
Construct BHaH_defines.h from data registered to griddata_commondata,
  CodeParameters, and NRPyParameters.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

import sys
from pathlib import Path
from typing import Optional, Dict, List
import nrpy.params as par
import nrpy.grid as gri
from nrpy.infrastructures.BHaH import griddata_commondata
from nrpy.helpers.generic import clang_format


def register_griddata_struct_and_return_griddata_struct_str() -> str:
    """
    Registers contributions to the griddata struct from different modules and constructs the griddata_struct string.

    :return: A string representing the typedef structure for grid data, including contributions from BHaH modules,
             reference_metric, CurviBoundaryConditions, and masking, if applicable.
    """
    # Step 1: Register griddata_struct contributions from all BHaH modules:
    griddata_commondata.register_griddata_commondata(
        "params",
        "params_struct params",
        "BHaH parameters, generated from NRPy+'s CodeParameters",
    )
    griddata_commondata.register_griddata_commondata(
        "MoL", "MoL_gridfunctions_struct gridfuncs", "MoL gridfunctions"
    )

    if any("reference_metric" in key for key in sys.modules):
        griddata_commondata.register_griddata_commondata(
            "reference_metric",
            "char CoordSystemname[100]",
            "the name of the CoordSystem (from reference_metric)",
        )
        griddata_commondata.register_griddata_commondata(
            "reference_metric",
            "char gridname[100]",
            "a user-defined alias for describing the grid",
        )
        griddata_commondata.register_griddata_commondata(
            "reference_metric",
            "rfm_struct rfmstruct",
            "includes e.g., 1D arrays of reference metric quantities",
        )

    if any("CurviBoundaryConditions" in key for key in sys.modules):
        griddata_commondata.register_griddata_commondata(
            "CurviBoundaryConditions",
            "bc_struct bcstruct",
            "all data needed to perform boundary conditions in curvilinear coordinates",
        )

    if any("masking" in key for key in sys.modules):
        griddata_commondata.register_griddata_commondata(
            "masking",
            "int8_t *restrict mask",
            "specifies grid overlaps on multipatch grids",
            is_commondata=True,
        )

    griddata_struct_def = r"""
typedef struct __griddata__ {
  // griddata_struct stores data needed on each grid
  // xx[3] stores the uniform grid coordinates.
  REAL *restrict xx[3];
"""
    for module, item_list in griddata_commondata.glb_griddata_struct_dict.items():
        griddata_struct_def += f"  // NRPy+ MODULE: {module}\n"
        for item in item_list:
            griddata_struct_def += f"  {item.c_declaration};"
            if item.description != "":
                griddata_struct_def += f"// <- {item.description}"
            griddata_struct_def += "\n"
    griddata_struct_def += "} griddata_struct;\n"

    return griddata_struct_def


def register_commondata_struct_and_return_commondata_struct_str() -> str:
    """
    Registers contributions to the commondata struct from different modules and constructs the commondata_struct string.

    :return: A string representing the typedef structure for grid data, including contributions from BHaH modules,
             reference_metric, CurviBoundaryConditions, and masking, if applicable.
    """
    commondata_struct_def = r"""
typedef struct __commondata__ {
  // commondata_struct stores data common to all grids
  int NUMGRIDS;
}
"""
    return commondata_struct_def


def output_BHaH_defines_h(
    project_dir: str,
    REAL_means: str = "double",
    enable_simd: bool = True,
    fin_NGHOSTS_add_one_for_upwinding: bool = False,
    MoL_method: str = "RK4",
    CoordSystem: str = "Cartesian",
    ID_persist_struct_contents_str: str = "",
    include_T4UU: bool = False,
    supplemental_defines_dict: Optional[Dict[str, str]] = None,
    clang_format_options: str = "-style={BasedOnStyle: LLVM, ColumnLimit: 0}",
) -> None:
    """
    Outputs C code header file with macro definitions and other configurations for the project.

    :param project_dir: Directory where the project C code is output
    :param enable_simd: Flag to enable Single Instruction Multiple Data (SIMD) optimizations
    :param fin_NGHOSTS_add_one_for_upwinding: Option to add one extra ghost zone for upwinding
    :param MoL_method: Specific Method of Lines timestepping method, default is "RK4"
    :param supplemental_defines_dict: Additional key-value pairs to be included in the output file
    :param clang_format_options: Options for clang formatting.

    :raises ValueError: If the project directory does not exist

    >>> from nrpy.infrastructures.BHaH.MoLtimestepping import MoL
    >>> import nrpy.finite_difference as fin
    >>> from nrpy.helpers.generic import compress_string_to_base64, decompress_base64_to_string, diff_strings
    >>> project_dir = Path("/tmp", "tmp_BHaH_defines_h")
    >>> project_dir.mkdir(parents=True, exist_ok=True)
    >>> output_BHaH_defines_h(project_dir=str(project_dir))
    >>> expected_string = decompress_base64_to_string("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")
    >>> returned_string = (project_dir / "BHaH_defines.h").read_text()
    >>> if returned_string != expected_string:
    ...    compressed_str = compress_string_to_base64(returned_string)
    ...    error_message = "Trusted BHaH_defines.h string changed!\\n"
    ...    error_message += "Here's the diff:\\n" + diff_strings(expected_string, returned_string) + "\\n"
    ...    raise ValueError(error_message + f"base64-encoded output: {compressed_str}")
    """
    project_Path = Path(project_dir)
    project_Path.mkdir(parents=True, exist_ok=True)

    BHaH_defines_h_dict: Dict[str, str] = {}

    ###############################
    # GENERALLY USEFUL DEFINITIONS
    gen_Nbd_str = """#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <errno.h>
#include <ctype.h>
"""
    if any("progress_indicator" in key for key in sys.modules):
        gen_Nbd_str += "#include <time.h>\n"
    if enable_simd:
        gen_Nbd_str += "// output_BHaH_defines_h(...,enable_simd=True) was called so we #include SIMD intrinsics:\n"
        gen_Nbd_str += """#include "simd/simd_intrinsics.h"\n"""
    gen_Nbd_str += f"#define REAL {REAL_means}\n"
    if any("progress_indicator" in key for key in sys.modules):
        gen_Nbd_str += r"""
#ifdef __linux__
// Timer with nanosecond resolution. Only on Linux.
#define TIMEVAR struct timespec
#define CURRTIME_FUNC(currtime) clock_gettime(CLOCK_REALTIME, currtime)
#define TIME_IN_NS(start, end) (REAL)(1000000000L * (end.tv_sec - start.tv_sec) + end.tv_nsec - start.tv_nsec);
#else
// Low-resolution timer, 1-second resolution. Widely available.
#define TIMEVAR time_t
#define CURRTIME_FUNC(currtime) time(currtime)
#define TIME_IN_NS(start, end) (REAL)(difftime(end, start)*1.0e9 + 1e-6) // Round up to avoid divide-by-zero.
#endif
"""
    gen_Nbd_str += r"""
#define MIN(A, B) ( ((A) < (B)) ? (A) : (B) )
#define MAX(A, B) ( ((A) > (B)) ? (A) : (B) )
"""
    BHaH_defines_h_dict["general"] = gen_Nbd_str

    #####################################
    # PARAMS_STRUCT AND COMMONDATA_STRUCT
    # Generate C code to declare C params_struct & commondata_struct;
    #         output to "declare_Cparameters_struct.h"
    #         We want the elements of this struct to be *sorted*,
    #         to ensure that the struct is consistently ordered
    #         for checkpointing purposes.
    par_Nbd_str = "typedef struct __params_struct__ {\n"
    commondata_Nbd_str = "typedef struct __commondata_struct__ {\n"
    CCodelines_params_struct: List[str] = []
    CCodelines_commondata_struct: List[str] = []

    # Iterate through the global code parameters dictionary
    for CPname, CodeParam in par.glb_code_params_dict.items():
        CPtype: str = CodeParam.c_type_alias
        if CPtype != "#define":
            comment: str = f"  // {CodeParam.module}::{CPname}"
            c_output = f"  {CPtype} {CPname};{comment}\n"
            if "char" in CPtype and "[" in CPtype and "]" in CPtype:
                chararray_size = CPtype.split("[")[1].replace("]", "")
                c_output = f"char {CPname}[{chararray_size}];{comment}\n"

            # Concatenate module and CPname for the comment
            if CodeParam.commondata:
                CCodelines_commondata_struct.append(c_output)
            else:
                CCodelines_params_struct.append(c_output)

    # Sort CCodelines_params_struct and append them to the par_Nbd_str
    for line in sorted(CCodelines_params_struct):
        par_Nbd_str += line
    # Sort CCodelines_commondata_struct and append them to the commondata_Nbd_str
    for line in sorted(CCodelines_commondata_struct):
        commondata_Nbd_str += line

    par_Nbd_str += "} params_struct;\n"
    commondata_Nbd_str += "} commondata_struct;\n"

    BHaH_defines_h_dict["params_struct"] = par_Nbd_str
    BHaH_defines_h_dict["commondata_struct"] = commondata_Nbd_str

    ###############################
    # FINITE DIFFERENCE
    # First register C functions needed by finite_difference

    # Then set up the dictionary entry for finite_difference in BHaH_defines
    if any("finite_difference" in key for key in sys.modules):
        NGHOSTS = int(par.parval_from_str("finite_difference::fd_order") / 2)
        if fin_NGHOSTS_add_one_for_upwinding:
            NGHOSTS += 1
        fin_Nbd_str = f"""
// Set the number of ghost zones
// Note that upwinding in e.g., BSSN requires that NGHOSTS = fd_order/2 + 1 <- Notice the +1.
#define NGHOSTS {NGHOSTS}
"""
        if not enable_simd:
            fin_Nbd_str += """
// When enable_simd = False, this is the UPWIND_ALG() macro:
#define UPWIND_ALG(UpwindVecU) UpwindVecU > 0.0 ? 1.0 : 0.0\n"""
        BHaH_defines_h_dict["finite_difference"] = fin_Nbd_str

    if any("reference_metric" in key for key in sys.modules):
        # pylint: disable=C0415
        import nrpy.infrastructures.BHaH.rfm_precompute as rfmp

        rfm_Nbd_str = rfmp.ReferenceMetricPrecompute(CoordSystem).BHaH_defines_str
        BHaH_defines_h_dict["reference_metric"] = rfm_Nbd_str

    ###############################
    # MoL:
    # if any of the modules contains "MoL" in its name
    if any("MoL" in key for key in sys.modules):
        # pylint: disable=C0415
        from nrpy.infrastructures.BHaH.MoLtimestepping import MoL

        # pylint: disable=C0415
        from nrpy.infrastructures.BHaH.MoLtimestepping.RK_Butcher_Table_Dictionary import (
            generate_Butcher_tables,
        )

        Butcher_dict = generate_Butcher_tables()

        # Generating gridfunction names based on the given MoL method
        (
            y_n_gridfunctions,
            non_y_n_gridfunctions_list,
            _diagnostic_gridfunctions_point_to,
            _diagnostic_gridfunctions2_point_to,
        ) = MoL.generate_gridfunction_names(Butcher_dict, MoL_method=MoL_method)

        # Step 3.b: Create MoL_timestepping struct:
        MoL_Nbd_str = (
            f"typedef struct __MoL_gridfunctions_struct__ {{\n"
            f"REAL *restrict {y_n_gridfunctions};\n"
            + "".join(f"REAL *restrict {gfs};\n" for gfs in non_y_n_gridfunctions_list)
            + r"""REAL *restrict diagnostic_output_gfs;
REAL *restrict diagnostic_output_gfs2;
} MoL_gridfunctions_struct;
#define LOOP_ALL_GFS_GPS(ii) \
_Pragma("omp parallel for") \
  for(int (ii)=0;(ii)<Nxx_plus_2NGHOSTS0*Nxx_plus_2NGHOSTS1*Nxx_plus_2NGHOSTS2*NUM_EVOL_GFS;(ii)++)
"""
        )
        BHaH_defines_h_dict["MoL"] = MoL_Nbd_str

    ###############################
    if any("ADM_Initial_Data_Reader" in key for key in sys.modules):
        # pylint: disable=C0415
        import nrpy.infrastructures.BHaH.general_relativity.ADM_Initial_Data_Reader__BSSN_Converter as ADMID

        ADMID_Nbd_str = ADMID.generate_BHaH_defines_contribution(
            ID_persist_struct_contents_str, include_T4UU
        )
        BHaH_defines_h_dict["ADMID"] = ADMID_Nbd_str

    ###############################
    if any("CurviBoundaryConditions" in key for key in sys.modules):
        # pylint: disable=C0415
        import nrpy.infrastructures.BHaH.CurviBoundaryConditions.CurviBoundaryConditions as CBC

        # First register basic C data structures/macros inside BHaH_defines.h
        CBC_Nbd_str = r"""
// NRPy+ Curvilinear Boundary Conditions: Core data structures
// Documented in: Tutorial-Start_to_Finish-Curvilinear_BCs.ipynb

typedef struct __innerpt_bc_struct__ {
  int dstpt;  // dstpt is the 3D grid index IDX3S(i0,i1,i2) of the inner boundary point (i0,i1,i2)
  int srcpt;  // srcpt is the 3D grid index (a la IDX3S) to which the inner boundary point maps
  int8_t parity[10];  // parity[10] is a calculation of dot products for the 10 independent parity types
} innerpt_bc_struct;

typedef struct __outerpt_bc_struct__ {
  short i0,i1,i2;  // the outer boundary point grid index (i0,i1,i2), on the 3D grid
  int8_t FACEX0,FACEX1,FACEX2;  // 1-byte integers that store
  //                               FACEX0,FACEX1,FACEX2 = +1, 0, 0 if on the i0=i0min face,
  //                               FACEX0,FACEX1,FACEX2 = -1, 0, 0 if on the i0=i0max face,
  //                               FACEX0,FACEX1,FACEX2 =  0,+1, 0 if on the i1=i2min face,
  //                               FACEX0,FACEX1,FACEX2 =  0,-1, 0 if on the i1=i1max face,
  //                               FACEX0,FACEX1,FACEX2 =  0, 0,+1 if on the i2=i2min face, or
  //                               FACEX0,FACEX1,FACEX2 =  0, 0,-1 if on the i2=i2max face,
} outerpt_bc_struct;

typedef struct __bc_info_struct__ {
  int num_inner_boundary_points;  // stores total number of inner boundary points
  int num_pure_outer_boundary_points[NGHOSTS][3];  // stores number of outer boundary points on each
  //                                                  ghostzone level and direction (update min and
  //                                                  max faces simultaneously on multiple cores)
  int bc_loop_bounds[NGHOSTS][6][6];  // stores outer boundary loop bounds. Unused after bcstruct_set_up()
} bc_info_struct;

typedef struct __bc_struct__ {
  innerpt_bc_struct *restrict inner_bc_array;  // information needed for updating each inner boundary point
  outerpt_bc_struct *restrict pure_outer_bc_array[NGHOSTS*3]; // information needed for updating each outer
  //                                                             boundary point
  bc_info_struct bc_info;  // stores number of inner and outer boundary points, needed for setting loop
  //                          bounds and parallelizing over as many boundary points as possible.
} bc_struct;
"""
        CBC_Nbd_str += CBC.BHaH_defines_set_gridfunction_defines_with_parity_types(
            verbose=True
        )
        BHaH_defines_h_dict["CurviBoundaryConditions"] = CBC_Nbd_str

    ###############################
    # GRID, etc.
    # Then set up the dictionary entry for grid in BHaH_defines
    gri_Nbd_str = gri.BHaHGridFunction.gridfunction_defines()
    gri_Nbd_str += (
        r"""
// Declare the IDX4(gf,i,j,k) macro, which enables us to store 4-dimensions of
//   data in a 1D array. In this case, consecutive values of "i"
//   (all other indices held to a fixed value) are consecutive in memory, where
//   consecutive values of "j" (fixing all other indices) are separated by
//   Nxx_plus_2NGHOSTS0 elements in memory. Similarly, consecutive values of
//   "k" are separated by Nxx_plus_2NGHOSTS0*Nxx_plus_2NGHOSTS1 in memory, etc.
#define IDX4(g,i,j,k)                                                  \
  ( (i) + Nxx_plus_2NGHOSTS0 * ( (j) + Nxx_plus_2NGHOSTS1 * ( (k) + Nxx_plus_2NGHOSTS2 * (g) ) ) )
#define IDX4pt(g,idx) ( (idx) + (Nxx_plus_2NGHOSTS0*Nxx_plus_2NGHOSTS1*Nxx_plus_2NGHOSTS2) * (g) )
#define IDX3(i,j,k) ( (i) + Nxx_plus_2NGHOSTS0 * ( (j) + Nxx_plus_2NGHOSTS1 * ( (k) ) ) )
#define LOOP_REGION(i0min,i0max, i1min,i1max, i2min,i2max)              \
  for(int i2=i2min;i2<i2max;i2++) for(int i1=i1min;i1<i1max;i1++) for(int i0=i0min;i0<i0max;i0++)
#define LOOP_OMP(__OMP_PRAGMA__, i0,i0min,i0max, i1,i1min,i1max, i2,i2min,i2max) \
_Pragma(__OMP_PRAGMA__)  \
    for(int (i2)=(i2min);(i2)<(i2max);(i2)++) for(int (i1)=(i1min);(i1)<(i1max);(i1)++) for(int (i0)=(i0min);(i0)<(i0max);(i0)++)
#define LOOP_NOOMP(i0,i0min,i0max, i1,i1min,i1max, i2,i2min,i2max)      \
  for(int (i2)=(i2min);(i2)<(i2max);(i2)++) for(int (i1)=(i1min);(i1)<(i1max);(i1)++) for(int (i0)=(i0min);(i0)<(i0max);(i0)++)
#define LOOP_BREAKOUT(i0,i1,i2, i0max,i1max,i2max) { i0=(i0max); i1=(i1max); i2=(i2max); break; }
#define IS_IN_GRID_INTERIOR(i0i1i2, Nxx_plus_2NGHOSTS0,Nxx_plus_2NGHOSTS1,Nxx_plus_2NGHOSTS2, NG) \
  ( i0i1i2[0] >= (NG) && i0i1i2[0] < (Nxx_plus_2NGHOSTS0)-(NG) &&       \
    i0i1i2[1] >= (NG) && i0i1i2[1] < (Nxx_plus_2NGHOSTS1)-(NG) &&       \
    i0i1i2[2] >= (NG) && i0i1i2[2] < (Nxx_plus_2NGHOSTS2)-(NG) )
"""
        + register_griddata_struct_and_return_griddata_struct_str()
    )
    BHaH_defines_h_dict["grid"] = gri_Nbd_str

    def output_key(key_name: str, item_name: str) -> str:
        return f"""
//********************************************
// Basic definitions for module {key_name}:\n{item_name}
//********************************************"""

    file_output_str = """// BHaH core header file, automatically generated from output_BHaH_defines_h within BHaH_defines_h.py,
//    DO NOT EDIT THIS FILE BY HAND.\n\n"""
    # The ordering here is based largely on data structure dependencies. E.g., griddata_struct contains bc_struct.
    core_modules_list = [
        "general",
        "commondata_struct",
        "params_struct",
        "finite_difference",
        "reference_metric",
        "CurviBoundaryConditions",
        "MoL",
        "interpolate",
        "grid",  # griddata struct depends upon other core modules
    ]
    # Populate BHaH_defines.h with core modules
    for key in core_modules_list:
        if key in BHaH_defines_h_dict:
            file_output_str += output_key(key, BHaH_defines_h_dict[key])

    # Populate BHaH_defines.h with non-core modules
    for key, item in BHaH_defines_h_dict.items():
        if key not in core_modules_list:
            print(f"Outputting non-core modules key = {key} to BHaH_defines.h")
            file_output_str += output_key(key, item)

    # Populate BHaH_defines.h with whatever else is desired.
    if supplemental_defines_dict:
        for key in supplemental_defines_dict:
            file_output_str += output_key(key, supplemental_defines_dict[key])

    bhah_defines_file = project_Path / "BHaH_defines.h"
    with bhah_defines_file.open("w", encoding="utf-8") as file:
        file.write(
            clang_format(file_output_str, clang_format_options=clang_format_options)
        )


if __name__ == "__main__":
    import doctest

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
