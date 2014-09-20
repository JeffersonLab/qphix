#ifndef QPHIX_PRINT_UTILS_H
#define QPHIX_PRINT_UTILS_H

#include <cstdio>
#include <cstdarg>
#include "qphix/qphix_config.h"

using namespace std;
#ifdef QPHIX_QMP_COMMS
#include <qmp.h>
#endif


namespace  QPhiX 
{
  void masterPrintf(const char * format, ... );
  void localPrintf(const char * format, ... );
};

#endif
