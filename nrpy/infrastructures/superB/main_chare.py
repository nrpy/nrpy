"""
Generate main.cpp, main.h, main.ci and commondata.h for the superB infrastructure.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
        Nishita Jadoo
        njadoo **at** uidaho **dot* edu
"""

from pathlib import Path
from typing import List, Tuple

import nrpy.c_function as cfc
import nrpy.params as par
from nrpy.helpers.generic import clang_format

# fmt: off
_ = par.CodeParameter("int", __name__, "Nchare0", 1, add_to_parfile=True, add_to_set_CodeParameters_h=True, commondata=True)
_ = par.CodeParameter("int", __name__, "Nchare1", 1, add_to_parfile=True, add_to_set_CodeParameters_h=True, commondata=True)
_ = par.CodeParameter("int", __name__, "Nchare2", 1, add_to_parfile=True, add_to_set_CodeParameters_h=True, commondata=True)
# fmt: on


def output_commondata_object_h(
    project_dir: str,
) -> None:
    r"""
    Output header file with definition for class CommondataObject.
    :param project_dir: Directory where the project C code is output
    """
    project_Path = Path(project_dir)
    project_Path.mkdir(parents=True, exist_ok=True)

    file_output_str = """#ifndef __COMMONDATAOBJECT_H__
#define __COMMONDATAOBJECT_H__
#include "BHaH_defines.h"
#include "superB/superB_pup_function_prototypes.h"
#include <time.h>
class CommondataObject {
  public:
    commondata_struct commondata;
    void pup(PUP::er &p) {
    // PUP commondata struct
    pup_commondata_struct(p, commondata);
    }
};
#endif //__COMMONDATAOBJECT_H__
"""
    commondata_object_file = project_Path / "commondata_object.h"
    with commondata_object_file.open("w", encoding="utf-8") as file:
        file.write(clang_format(file_output_str))


def output_main_h(
    project_dir: str,
    enable_charm_checkpointing: bool = False,
    enable_BHaHAHA: bool = False,
) -> None:
    """
    Generate main.h.
    :param project_dir: Directory where the project C code is output
    :param enable_charm_checkpointing: Enable checkpointing using Charm++.
    :param enable_BHaHAHA: If True, include horizon_finder and interpolator3d declarations.
    """
    project_Path = Path(project_dir)
    project_Path.mkdir(parents=True, exist_ok=True)

    file_output_str = """#ifndef __MAIN_H__
#define __MAIN_H__

#include "pup.h"
#include "main.decl.h"
#include "timestepping.decl.h" """
    if enable_BHaHAHA:
        file_output_str += """
#include "horizon_finder.decl.h"
#include "interpolator3d.decl.h"
"""
    file_output_str += """
class Main : public CBase_Main {

  private:
    /// Member Variables (Object State) ///
    REAL start_time;

    /// Private Member Functions ///

 public:

  /// Constructors ///
  Main(CkArgMsg* msg);
  Main(CkMigrateMessage* msg);
"""
    if enable_charm_checkpointing:
        file_output_str += r"""
  void pup(PUP::er &p);
"""
    file_output_str += r"""
  /// Entry Methods ///
  void done();
};

#endif //__MAIN_H__
"""
    main_h_file = project_Path / "main.h"
    with main_h_file.open("w", encoding="utf-8") as file:
        file.write(clang_format(file_output_str))


def output_main_cpp(
    project_dir: str,
    enable_charm_checkpointing: bool = False,
    enable_BHaHAHA: bool = False,
) -> None:
    """
    Generate the "generic" C main() function for all simulation codes in the superB infrastructure.
    :param project_dir: Directory where the project C code is output
    :param enable_charm_checkpointing: Enable checkpointing using Charm++.
    :raises ValueError: Raised if any required function for superB main() is not registered.
    :param enable_BHaHAHA: If True, emit horizon_finder + interpolator3d arrays and the RR mapper.
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

    file_output_str = r"""#include <cstring>
#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"
#include "main.h"
#include "timestepping.decl.h" """

    if enable_BHaHAHA:
        file_output_str += r"""
#include "horizon_finder.decl.h"
#include "interpolator3d.decl.h" """

    file_output_str += r"""
/* readonly */ CProxy_Main mainProxy;
/* readonly */ CProxy_Timestepping timesteppingArray;
"""

    if enable_BHaHAHA:
        file_output_str += r"""
/* readonly */ CProxy_Horizon_finder horizon_finderProxy;
/* readonly */ CProxy_Interpolator3d interpolator3dArray; """

    if enable_BHaHAHA and not enable_charm_checkpointing:
        file_output_str += r"""
//===============================================
// RRMap_with_offset
// - Flatten (x,y,z) -> linear = x + N0*y + N0*N1*z
// - PE = (linear + offset) % CkNumPes()
// - Works for 1D/2D (missing dims treated as 0).
//===============================================
class RRMap_with_offset : public CkArrayMap {
    int offset;
    int Nchare0, Nchare1;

public:
    RRMap_with_offset(int _offset = 0, int _Nchare0 = 1, int _Nchare1 = 1)
        : offset(_offset), Nchare0(_Nchare0), Nchare1(_Nchare1) {}
    RRMap_with_offset(CkMigrateMessage* /*msg*/) {}

    void pup(PUP::er &p) override {
        p | offset;
        p | Nchare0;
        p | Nchare1;
    }

    int procNum(int /*unused*/, const CkArrayIndex& idx) override {
        const int x = idx.data()[0];
        const int y = (idx.nInts > 1) ? idx.data()[1] : 0;
        const int z = (idx.nInts > 2) ? idx.data()[2] : 0;
        const int linear = x + Nchare0 * y + Nchare0 * Nchare1 * z;
        return (linear + offset) % CkNumPes();
    }
};
"""

    file_output_str += r"""
/*
 * -={ main() function }=-
 * Step 1.a: Set each commondata CodeParameter to default.
 * Step 1.b: Overwrite default values to parfile values. Then overwrite parfile values with values set at cmd line.
 */

Main::Main(CkArgMsg* msg) {
  start_time = CkWallTimer();

  mainProxy = thisProxy;

  CommondataObject commondataObj; // commondataObj.commondata contains parameters common to all grids.

  // Step 1.a: Set each commondata CodeParameter to default.
  commondata_struct_set_to_default(&commondataObj.commondata);

  // Step 1.b: Overwrite default values to parfile values. Then overwrite parfile values with values set at cmd line.
  const char** argv_const = new const char*[msg->argc];
  for (int i = 0; i < msg->argc; i++) {
      argv_const[i] = msg->argv[i];
  }
  cmdline_input_and_parfile_parser(&commondataObj.commondata, msg->argc, argv_const);
  delete[] argv_const;
"""

    if not enable_BHaHAHA:
        # Original minimal launch
        file_output_str += r"""
  timesteppingArray = CProxy_Timestepping::ckNew(commondataObj, commondataObj.commondata.Nchare0, commondataObj.commondata.Nchare1, commondataObj.commondata.Nchare2);

  timesteppingArray.start();
}
"""
    else:
        # Horizons + Timestepping + Interpolator
        file_output_str += r"""
  // Dimensions & counts
  const int Nchare0 = commondataObj.commondata.Nchare0;
  const int Nchare1 = commondataObj.commondata.Nchare1;
  const int Nchare2 = commondataObj.commondata.Nchare2;
  const int Ncharetotal = Nchare0 * Nchare1 * Nchare2;
  const int numHorizons = commondataObj.commondata.bah_max_num_horizons;

  CkArrayOptions optsH(numHorizons);
  horizon_finderProxy = CProxy_Horizon_finder::ckNew(optsH);

  CkArrayOptions optsTS(Nchare0, Nchare1, Nchare2);
"""
        if not enable_charm_checkpointing:
            file_output_str += r"""
  CProxy_RRMap_with_offset rrMap = CProxy_RRMap_with_offset::ckNew(0, Nchare0, Nchare1);
  optsTS.setMap(rrMap);
"""
        file_output_str += r"""
  timesteppingArray = CProxy_Timestepping::ckNew(commondataObj, optsTS);

  CkArrayOptions optsI(Nchare0, Nchare1, Nchare2);
"""
        if not enable_charm_checkpointing:
            file_output_str += r"""
  optsI.setMap(rrMap);
"""
        file_output_str += r"""
  interpolator3dArray = CProxy_Interpolator3d::ckNew(optsI);

"""
        file_output_str += r"""
  timesteppingArray.start();
}
"""

    file_output_str += r"""
Main::Main(CkMigrateMessage* msg) {
  mainProxy = thisProxy;
}

void Main::done() {
  CkPrintf("\nTotal wall clock time = %f s.\n", CkWallTimer() - start_time);
  CkExit();
}
"""
    if enable_charm_checkpointing:
        file_output_str += r"""
// PUP routine for checkpointing
void Main::pup(PUP::er &p) {
  CBase_Main::pup(p);
}
"""
    file_output_str += r"""
#include "main.def.h"
"""
    main_cpp_file = project_Path / "main.cpp"
    with main_cpp_file.open("w", encoding="utf-8") as file:
        file.write(clang_format(file_output_str))


def output_main_ci(
    project_dir: str,
    enable_BHaHAHA: bool = False,
    enable_charm_checkpointing: bool = False,
) -> None:
    """
    Generate main.ci.

    :param project_dir: Directory where the project C code is output
    :param enable_BHaHAHA: If True, include horizon_finder and interpolator3d declarations.
    :param enable_charm_checkpointing: Enable checkpointing using Charm++.
    """
    project_Path = Path(project_dir)
    project_Path.mkdir(parents=True, exist_ok=True)

    file_output_str = """mainmodule main {

  include "commondata_object.h";
  include "pup_stl.h";

  readonly CProxy_Main mainProxy;
  readonly CProxy_Timestepping timesteppingArray;
  """
    if enable_BHaHAHA:
        file_output_str += r"""
  readonly CProxy_Horizon_finder horizon_finderProxy;
  readonly CProxy_Interpolator3d interpolator3dArray;
  """
    file_output_str += r"""
  extern module timestepping;
  """
    if enable_BHaHAHA:
        file_output_str += r"""
  extern module horizon_finder;
  extern module interpolator3d;
"""
    file_output_str += r"""
  mainchare Main {
    entry Main(CkArgMsg* msg);
    entry void done();
  };"""
    if enable_BHaHAHA and not enable_charm_checkpointing:
        file_output_str += r"""
  //-----------------------------------------------
  // RRMap_with_offset: round-robin array mapper with a start offset.
  // PE = (flatten(x,y,z; Nchare0,Nchare1) + offset) % CkNumPes()
  //-----------------------------------------------
  group RRMap_with_offset : CkArrayMap {
    entry RRMap_with_offset(int offset, int Nchare0, int Nchare1);
  };"""

    file_output_str += r"""
}
"""
    main_ci_file = project_Path / "main.ci"
    with main_ci_file.open("w", encoding="utf-8") as file:
        file.write(clang_format(file_output_str))


def output_commondata_object_h_and_main_h_cpp_ci(
    project_dir: str,
    enable_charm_checkpointing: bool = False,
    enable_BHaHAHA: bool = False,
) -> None:
    """
    Generate commondata_object.h, main.h, main.cpp and main.ci.
    :param project_dir: Directory where the project C code is output.
    :param enable_charm_checkpointing: Enable checkpointing using Charm++.
    :param enable_BHaHAHA: If True, include horizon_finder and interpolator3d declarations.
    """
    output_commondata_object_h(
        project_dir=project_dir,
    )

    output_main_h(
        project_dir=project_dir,
        enable_charm_checkpointing=enable_charm_checkpointing,
        enable_BHaHAHA=enable_BHaHAHA,
    )

    output_main_cpp(
        project_dir=project_dir,
        enable_charm_checkpointing=enable_charm_checkpointing,
        enable_BHaHAHA=enable_BHaHAHA,
    )

    output_main_ci(
        project_dir=project_dir,
        enable_BHaHAHA=enable_BHaHAHA,
        enable_charm_checkpointing=enable_charm_checkpointing,
    )
