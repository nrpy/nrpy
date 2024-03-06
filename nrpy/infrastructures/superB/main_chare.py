"""

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
        Nishita Jadoo
        njadoo **at** uidaho **dot* edu
"""

from typing import List, Tuple, Optional, Dict
import nrpy.c_function as cfc
import sys
from pathlib import Path
import nrpy.params as par
import nrpy.grid as gri
from nrpy.infrastructures.BHaH import griddata_commondata
from nrpy.helpers.generic import clang_format
from nrpy.infrastructures.BHaH import BHaH_defines_h

# fmt: off
_ = par.CodeParameter("int", __name__, "Nchare0", 1, add_to_parfile=True, add_to_set_CodeParameters_h=True, commondata=True)
_ = par.CodeParameter("int", __name__, "Nchare1", 1, add_to_parfile=True, add_to_set_CodeParameters_h=True, commondata=True)
_ = par.CodeParameter("int", __name__, "Nchare2", 1, add_to_parfile=True, add_to_set_CodeParameters_h=True, commondata=True)
# fmt: on


def output_commondata_object_h(
    project_dir: str,
    clang_format_options: str = "-style={BasedOnStyle: LLVM, ColumnLimit: 150}",
) -> None:
    r"""
    Output C code header file with macro definitions and other configurations for the project.
    :param project_dir: Directory where the project C code is output
    """
    project_Path = Path(project_dir)
    project_Path.mkdir(parents=True, exist_ok=True)

    file_output_str = """#ifndef __COMMONDATAOBJECT_H__
#define __COMMONDATAOBJECT_H__
#include "BHaH_defines.h"
#include <time.h>
class CommondataObject {
  public:
    commondata_struct commondata;
    void pup(PUP::er &p) {
"""
    struct_list: List[str] = []  # List to store individual struct elements
    for parname, CodeParam in par.glb_code_params_dict.items():
      struct = "commondata"
      CPtype = CodeParam.c_type_alias
      comment = f"  // {CodeParam.module}::{parname}"
      if "char" in CPtype and "[" in CPtype and "]" in CPtype:
        chararray_size = CPtype.split("[")[1].replace("]", "")
        c_output = f'PUParray(p, {struct}.{parname}, {chararray_size});{comment}\n'
      else:
        c_output = f"p|{struct}.{parname};{comment}\n"
      struct_list.append(c_output)
    # Sort the lines alphabetically and join them with line breaks
    file_output_str += "// PUP commondata struct\n"
    file_output_str += "".join(sorted(struct_list))
    file_output_str += "}\n"
    file_output_str += """
};
#endif //__COMMONDATAOBJECT_H__
"""
    commondata_object_file = project_Path / "commondata_object.h"
    with commondata_object_file.open("w", encoding="utf-8") as file:
        file.write(
            clang_format(file_output_str, clang_format_options=clang_format_options)
        )


def output_main_h(
    project_dir: str,
    clang_format_options: str = "-style={BasedOnStyle: LLVM, ColumnLimit: 150}",
) -> None:
    """
    Generate main.h
    :param clang_format_options: Clang formatting options, default is "-style={BasedOnStyle: LLVM, ColumnLimit: 150}".
    """
    project_Path = Path(project_dir)
    project_Path.mkdir(parents=True, exist_ok=True)

    file_output_str = """#ifndef __MAIN_H__
#define __MAIN_H__

#include "main.decl.h"
#include "timetepping.decl.h"

class Main : public CBase_Main {

  private:
    /// Member Variables (Object State) ///
    CProxy_Timetepping_array;

    /// Private Member Functions ///

 public:

  /// Constructors ///
  Main(CkArgMsg* msg);
  Main(CkMigrateMessage* msg);

  /// Entry Methods ///
  void done();
};

#endif //__MAIN_H__
"""
    main_h_file = project_Path / "main.h"
    with main_h_file.open("w", encoding="utf-8") as file:
        file.write(
            clang_format(file_output_str, clang_format_options=clang_format_options)
        )


def output_main_cpp(
    project_dir: str,
    clang_format_options: str = "-style={BasedOnStyle: LLVM, ColumnLimit: 150}",
) -> None:
    """
    Generate the "generic" C main() function for all simulation codes in the BHaH infrastructure.

    :param clang_format_options: Clang formatting options, default is "-style={BasedOnStyle: LLVM, ColumnLimit: 150}".
    """

    # Make sure all required C functions are registered
    missing_functions: List[Tuple[str, str]] = []
    for func_tuple in [
        ("commondata_struct_set_to_default", "CodeParameters.py"),
        ("cmdline_input_and_parfile_parser", "cmdline_input_and_parfiles.py"),
    ]:
        if func_tuple[0] not in cfc.CFunction_dict:
            missing_functions += [func_tuple]
    if missing_functions:
        error_msg = "Error: These functions are required, and are not registered.\n"
        for func_tuple in missing_functions:
            error_msg += (
                f'  {func_tuple[0]}, registered by function within "{func_tuple[1]}"\n'
            )
        raise ValueError(error_msg)

    project_Path = Path(project_dir)
    project_Path.mkdir(parents=True, exist_ok=True)

    file_output_str = """#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"
#include "main.h"
#include "timestepping.decl.h"

/* readonly */ CProxy_Main mainProxy;

/*
 * -={ main() function }=-
 * Step 1.a: Set each commondata CodeParameter to default.
 * Step 1.b: Overwrite default values to parfile values. Then overwrite parfile values with values set at cmd line.
 */

Main::Main(CkArgMsg* msg) {
  mainProxy = thisProxy;

  CommondataObject commondataObj; // commondataObj.commondata contains parameters common to all grids.

  // Step 1.a: Set each commondata CodeParameter to default.
  commondata_struct_set_to_default(&commondataObj.commondata);

  // Step 1.b: Overwrite default values to parfile values. Then overwrite parfile values with values set at cmd line.
  cmdline_input_and_parfile_parser(&commondataObj.commondata, msg->argc, msg->argv);

  timestepping_array = CProxy_Timestepping::ckNew(commondataObj, &commondataObj.commondata.Nchare0, &commondataObj.commondata.Nchare1, &commondataObj.commondata.Nchare2);

  timestepping_array.start();
}

Main::Main(CkMigrateMessage* msg) { }

void Main::done() {
  CkExit();
}

#include "main.def.h"
"""
    main_cpp_file = project_Path / "main.cpp"
    with main_cpp_file.open("w", encoding="utf-8") as file:
        file.write(
            clang_format(file_output_str, clang_format_options=clang_format_options)
        )


def output_main_ci(
    project_dir: str,
    clang_format_options: str = "-style={BasedOnStyle: LLVM, ColumnLimit: 150}",
) -> None:
    """
    Generate main.h
    :param clang_format_options: Clang formatting options, default is "-style={BasedOnStyle: LLVM, ColumnLimit: 150}".
    """
    project_Path = Path(project_dir)
    project_Path.mkdir(parents=True, exist_ok=True)

    file_output_str = """mainmodule main {

  include "commondata_object.h";

  readonly CProxy_Main mainProxy;

  extern module moL_step_forward_in_time;

  mainchare Main {
    entry Main(CkArgMsg* msg);
    entry void done();
  };

}
"""
    main_ci_file = project_Path / "main.ci"
    with main_ci_file.open("w", encoding="utf-8") as file:
        file.write(
            clang_format(file_output_str, clang_format_options=clang_format_options)
        )


def output_commondata_object_h_and_main_h_cpp_ci(
    project_dir: str,
    clang_format_options: str = "-style={BasedOnStyle: LLVM, ColumnLimit: 150}",
) -> None:
    """
    Generate commondata_object.h, main.h, main.cpp and main.ci
    :param clang_format_options: Clang formatting options, default is "-style={BasedOnStyle: LLVM, ColumnLimit: 150}".
    """

    output_commondata_object_h(
      project_dir=project_dir,
    )

    output_main_h(
        project_dir=project_dir,
    )

    output_main_cpp(
        project_dir=project_dir,
    )

    output_main_ci(
        project_dir=project_dir,
    )
