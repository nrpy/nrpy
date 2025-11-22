#ifndef __MAIN_H__
#define __MAIN_H__

#include "main.decl.h"
#include "pup.h"
#include "timestepping.decl.h"
class Main : public CBase_Main {

private:
  /// Member Variables (Object State) ///
  REAL start_time;

  /// Private Member Functions ///

public:
  /// Constructors ///
  Main(CkArgMsg *msg);
  Main(CkMigrateMessage *msg);

  /// Entry Methods ///
  void done();
};

#endif //__MAIN_H__
