#ifndef BLAS_H
#define BLAS_H

#include "dslash_config.h"
#ifdef QMP_COMMS
#include "qmp.h"
#endif

// Generic OpenMP templated
#include "blas_c.h"
#include "blas_new_c.h"

// MIC Specializations
#ifdef MIC_SOURCE
#include "blas_mic.h"
#endif

// SSE Specializations
#ifdef SSE_SOURCE
#endif

// AVX Specialization
#ifdef AVX_SOURCE
#endif

#endif
