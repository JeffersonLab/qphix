#pragma once

#ifdef TIMING_CG
// Stuff for rdtsc
#include <sys/types.h>
#include <sys/sysctl.h>

typedef long long TSC_tick;
#define CLOCK_NOW(a) (a) = __rdtsc()
#endif
