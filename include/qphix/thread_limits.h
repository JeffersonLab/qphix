#ifndef THREAD_LIMITS_H
#define THREAD_LIMITS_H

namespace QPhiX 
{

#define MAX_CORES       128
#define MAX_THREADS_PER_CORE    4
#define MAX_THREADS     (MAX_CORES * MAX_THREADS_PER_CORE)

};
#endif
