#pragma once

namespace QPhiX
{
void setThreadAffinity(int nCores, int threadsPerCore);
void reportAffinity(void);
}
