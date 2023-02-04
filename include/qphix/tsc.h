#pragma once

// Stuff for rdtsc
#include <sys/types.h>
#ifdef QPHIX_HAVE_SYSCONF
// no need for any headers
#elif defined QPHIX_HAVE_SYSCTL
#include <sys/sysctl.h>
#endif

typedef long long TSC_tick;
#define CLOCK_NOW(a) (a) = __rdtsc()
