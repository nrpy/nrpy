"""
Set up basic functions and loop insertions for precomputed reference metric infrastructure.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""
from typing import List
import sympy as sp

import nrpy.c_function as cfc
import nrpy.reference_metric as refmetric
from nrpy.helpers.generic import superfast_uniq
from nrpy.infrastructures.BHaH.BHaH_defines_h import register_BHaH_defines


class ReferenceMetricPrecompute:
    """
    Base class for reference metric precomputation.

    This class stores contributions to BHaH_defines.h, as well as functions for memory allocation,
    definition, and freeing of rfm precomputation data. It also provides strings for reading rfm
    precompute quantities within loops for both SIMD-ized and non-SIMD loops.
    """

    def __init__(self, CoordSystem: str):
        rfm = refmetric.reference_metric[CoordSystem + "_rfm_precompute"]
        # Step 7: Construct needed C code for declaring rfmstruct, allocating storage for
        #         rfmstruct arrays, defining each element in each array, reading the
        #         rfmstruct data from memory (both with and without SIMD enabled), and
        #         freeing allocated memory for the rfmstruct arrays.

        # BHaH_defines_list: String that contains the body of the rfmstruct struct.
        self.BHaH_defines_list: List[str] = []
        self.rfm_struct__define = ""
        # rfmstruct stores pointers to (so far) 1D arrays. The rfm_struct__malloc string allocates space for the arrays.
        self.rfm_struct__malloc = ""
        self.rfm_struct__freemem = ""

        # readvr_str reads the arrays from memory as needed
        self.readvr_str = ["", "", ""]
        self.readvr_SIMD_outer_str = ["", "", ""]
        self.readvr_SIMD_inner_str = ["", "", ""]

        # Sort freevars_uniq_vals and freevars_uniq_xx_indep, according to alphabetized freevars_uniq_xx_indep.
        #    Without this step, the ordering of elements in rfmstruct would be random, and would change each time
        #    this function was called.
        freevars_uniq_xx_indep: List[sp.Expr] = []
        freevars_uniq_vals: List[sp.Expr] = []
        if len(rfm.freevars_uniq_xx_indep) > 0:
            freevars_uniq_xx_indep, freevars_uniq_vals = (
                list(x)
                for x in zip(
                    *sorted(
                        zip(rfm.freevars_uniq_xx_indep, rfm.freevars_uniq_vals), key=str
                    )
                )
            )

        # Tease out how many variables each function in freevars_uniq_vals
        Nxx_plus_2NGHOSTS = [
            sp.Symbol("Nxx_plus_2NGHOSTS0", real=True),
            sp.Symbol("Nxx_plus_2NGHOSTS1", real=True),
            sp.Symbol("Nxx_plus_2NGHOSTS2", real=True),
        ]
        which_freevar: int = 0
        for expr in freevars_uniq_vals:
            if "_of_xx" in str(freevars_uniq_xx_indep[which_freevar]):
                frees = list(expr.free_symbols)
                frees_uniq = superfast_uniq(frees)
                xx_list: List[sp.Basic] = []
                malloc_size: int = 1
                for i in range(3):
                    if rfm.xx[i] in frees_uniq:
                        xx_list.append(rfm.xx[i])
                        malloc_size *= Nxx_plus_2NGHOSTS[i]

                self.BHaH_defines_list += [
                    f"REAL *restrict {freevars_uniq_xx_indep[which_freevar]};"
                ]
                self.rfm_struct__malloc += f"rfmstruct->{freevars_uniq_xx_indep[which_freevar]} = (REAL *)malloc(sizeof(REAL)*{malloc_size});\n"
                self.rfm_struct__freemem += (
                    f"free(rfmstruct->{freevars_uniq_xx_indep[which_freevar]});\n"
                )
                output_define_and_readvr = False
                for dirn in range(3):
                    if (
                        (rfm.xx[dirn] in frees_uniq)
                        and not (rfm.xx[(dirn + 1) % 3] in frees_uniq)
                        and not (rfm.xx[(dirn + 2) % 3] in frees_uniq)
                    ):
                        self.rfm_struct__define += (
                            f"for(int i{dirn}=0;i{dirn}<Nxx_plus_2NGHOSTS{dirn};i{dirn}++) {{\n"
                            f"  const REAL xx{dirn} = xx[{dirn}][i{dirn}];\n"
                            f"  rfmstruct->{freevars_uniq_xx_indep[which_freevar]}[i{dirn}] = {sp.ccode(freevars_uniq_vals[which_freevar])};\n"
                        )
                        self.rfm_struct__define += "}\n\n"
                        self.readvr_str[
                            dirn
                        ] += f"const REAL {freevars_uniq_xx_indep[which_freevar]} = rfmstruct->{freevars_uniq_xx_indep[which_freevar]}[i{dirn}];\n"
                        self.readvr_SIMD_outer_str[
                            dirn
                        ] += f"const double NOSIMD{freevars_uniq_xx_indep[which_freevar]} = rfmstruct->{freevars_uniq_xx_indep[which_freevar]}[i{dirn}]; "
                        self.readvr_SIMD_outer_str[
                            dirn
                        ] += f"const REAL_SIMD_ARRAY {freevars_uniq_xx_indep[which_freevar]} = ConstSIMD(NOSIMD{freevars_uniq_xx_indep[which_freevar]});\n"
                        self.readvr_SIMD_inner_str[
                            dirn
                        ] += f"const REAL_SIMD_ARRAY {freevars_uniq_xx_indep[which_freevar]} = ReadSIMD(&rfmstruct->{freevars_uniq_xx_indep[which_freevar]}[i{dirn}]);\n"
                        output_define_and_readvr = True

                if (
                    (not output_define_and_readvr)
                    and (rfm.xx[0] in frees_uniq)
                    and (rfm.xx[1] in frees_uniq)
                ):
                    self.rfm_struct__define += f"""
                for(int i1=0;i1<Nxx_plus_2NGHOSTS1;i1++) for(int i0=0;i0<Nxx_plus_2NGHOSTS0;i0++) {{
                  const REAL xx0 = xx[0][i0];
                  const REAL xx1 = xx[1][i1];
                  rfmstruct->{freevars_uniq_xx_indep[which_freevar]}[i0 + Nxx_plus_2NGHOSTS0*i1] = {sp.ccode(freevars_uniq_vals[which_freevar])};
                }}\n\n"""
                    self.readvr_str[
                        0
                    ] += f"const REAL {freevars_uniq_xx_indep[which_freevar]} = rfmstruct->{freevars_uniq_xx_indep[which_freevar]}[i0 + Nxx_plus_2NGHOSTS0*i1];\n"
                    self.readvr_SIMD_outer_str[
                        0
                    ] += f"const double NOSIMD{freevars_uniq_xx_indep[which_freevar]} = rfmstruct->{freevars_uniq_xx_indep[which_freevar]}[i0 + Nxx_plus_2NGHOSTS0*i1]; "
                    self.readvr_SIMD_outer_str[
                        0
                    ] += f"const REAL_SIMD_ARRAY {freevars_uniq_xx_indep[which_freevar]} = ConstSIMD(NOSIMD{freevars_uniq_xx_indep[which_freevar]});\n"
                    self.readvr_SIMD_inner_str[
                        0
                    ] += f"const REAL_SIMD_ARRAY {freevars_uniq_xx_indep[which_freevar]} = ReadSIMD(&rfmstruct->{freevars_uniq_xx_indep[which_freevar]}[i0 + Nxx_plus_2NGHOSTS0*i1]);\n"
                    output_define_and_readvr = True

                if not output_define_and_readvr:
                    raise RuntimeError(
                        f"ERROR: Could not figure out the (xx0,xx1,xx2) dependency within the expression for {freevars_uniq_xx_indep[which_freevar]}: {freevars_uniq_vals[which_freevar]}"
                    )

                if not output_define_and_readvr:
                    raise RuntimeError(
                        f"ERROR: Could not figure out the (xx0,xx1,xx2) dependency within the expression for {freevars_uniq_xx_indep[which_freevar]}: {freevars_uniq_vals[which_freevar]}"
                    )

            which_freevar += 1


def register_CFunctions_rfm_precompute(list_of_CoordSystems: List[str]) -> None:
    """
    Register C functions for reference metric precomputed lookup arrays.

    :param list_of_CoordSystems: List of coordinate systems to register the C functions.
    """
    combined_BHaH_defines_list = []
    for CoordSystem in list_of_CoordSystems:
        rfm_precompute = ReferenceMetricPrecompute(CoordSystem)

        for func in [
            ("malloc", rfm_precompute.rfm_struct__malloc),
            ("defines", rfm_precompute.rfm_struct__define),
            ("free", rfm_precompute.rfm_struct__freemem),
        ]:
            includes = ["BHaH_defines.h"]

            desc = f"rfm_precompute_{func[0]}: reference metric precomputed lookup arrays: {func[0]}"
            c_type = "void"
            name = "rfm_precompute_" + func[0]
            params = "const commondata_struct *restrict commondata, const params_struct *restrict params, rfm_struct *restrict rfmstruct"
            include_CodeParameters_h = True
            if func[0] == "defines":
                params += ", REAL *restrict xx[3]"

            body = " "
            body += func[1]

            combined_BHaH_defines_list.extend(rfm_precompute.BHaH_defines_list)
            cfc.register_CFunction(
                includes=includes,
                desc=desc,
                c_type=c_type,
                CoordSystem_for_wrapper_func=CoordSystem,
                name=name,
                params=params,
                include_CodeParameters_h=include_CodeParameters_h,
                body=body,
            )

    BHaH_defines = "typedef struct __rfmstruct__ {\n"
    BHaH_defines += "\n".join(sorted(superfast_uniq(combined_BHaH_defines_list)))
    BHaH_defines += "\n} rfm_struct;\n"
    register_BHaH_defines("reference_metric", BHaH_defines)
