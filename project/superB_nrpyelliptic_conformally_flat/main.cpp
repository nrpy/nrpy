#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"
#include "main.h"
#include "timestepping.decl.h"
/* readonly */ CProxy_Main mainProxy;
/* readonly */ CProxy_Timestepping timesteppingArray;

/*
 * -={ main() function }=-
 * Step 1.a: Set each commondata CodeParameter to default.
 * Step 1.b: Overwrite default values to parfile values. Then overwrite parfile values with values set at cmd line.
 */

Main::Main(CkArgMsg *msg) {
  start_time = CkWallTimer();

  mainProxy = thisProxy;

  CommondataObject commondataObj; // commondataObj.commondata contains parameters common to all grids.

  // Step 1.a: Set each commondata CodeParameter to default.
  commondata_struct_set_to_default(&commondataObj.commondata);

  // Step 1.b: Overwrite default values to parfile values. Then overwrite parfile values with values set at cmd line.
  const char **argv_const = new const char *[msg->argc];
  for (int i = 0; i < msg->argc; i++) {
    argv_const[i] = msg->argv[i];
  }
  cmdline_input_and_parfile_parser(&commondataObj.commondata, msg->argc, argv_const);
  delete[] argv_const;

  timesteppingArray =
      CProxy_Timestepping::ckNew(commondataObj, commondataObj.commondata.Nchare0, commondataObj.commondata.Nchare1, commondataObj.commondata.Nchare2);

  timesteppingArray.start();
}

Main::Main(CkMigrateMessage *msg) { mainProxy = thisProxy; }

void Main::done() {
  CkPrintf("\nTotal wall clock time = %f s.\n", CkWallTimer() - start_time);
  CkExit();
}

#include "main.def.h"
