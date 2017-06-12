#ifndef QPHIX_BLAS_H
#define QPHIX_BLAS_H

// Generic OpenMP templated
#include "qphix/blas_c.h"
#include "qphix/blas_new_c.h"

// MIC Specializations
#ifdef QPHIX_MIC_SOURCE
#include "qphix/blas_mic.h"
#endif

// SSE Specializations
#ifdef QPHIX_SSE_SOURCE
#endif

// AVX Specialization
#ifdef QPHIX_AVX_SOURCE
#endif

#ifdef QPHIX_AVX2_SOURCE
#endif
#endif
