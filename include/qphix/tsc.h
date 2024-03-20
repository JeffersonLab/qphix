#pragma once

#ifdef TIMING_CG
// Stuff for rdtsc
#include <sys/types.h>
#include <sys/sysctl.h>

#define CLOCK_NOW(a) (a) = __rdtsc()
#endif
typedef long long TSC_tick;
