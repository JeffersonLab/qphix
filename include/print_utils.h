#ifndef PRINT_UTILS_H
#define PRINT_UTILS_H

#include <cstdio>
#include <cstdarg>
#include "dslash_config_internal.h"
using namespace std;
#ifdef QMP_COMMS
#include <qmp.h>
#endif


namespace CPlusPlusWilsonDslash { 

  void masterPrintf(const char * format, ... );
  void localPrintf(const char * format, ... );

};

#endif
