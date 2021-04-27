#pragma once

#include <cstdio>
#include <cstdarg>

#include "qphix/qphix_config.h"

namespace QPhiX
{
void masterPrintf(const char *format, ...);
void localPrintf(const char *format, ...);
}
