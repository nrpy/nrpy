#ifndef __COMMONDATAOBJECT_H__
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
