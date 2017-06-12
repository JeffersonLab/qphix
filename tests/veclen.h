#pragma once

#ifndef QPHIX_SOALEN
#error "QPHIX_SOALEN is not defined"
#endif

#if defined(QPHIX_MIC_SOURCE) || defined(QPHIX_AVX512_SOURCE)

#define VECLEN_SP 16
#define VECLEN_HP 16
#define VECLEN_DP 8

#elif defined(QPHIX_AVX_SOURCE) || defined(QPHIX_AVX2_SOURCE)

#define VECLEN_SP 8
#define VECLEN_HP 8
#define VECLEN_DP 4

#elif defined(QPHIX_SSE_SOURCE)

#define VECLEN_SP 4
#define VECLEN_DP 2

#elif defined(QPHIX_SCALAR_SOURCE)
#define VECLEN_DP 1
#define VECLEN_SP 1
#define QPHIX_SOALEN 1

#elif defined(QPHIX_QPX_SOURCE)
#define VECLEN_DP 4
#define QPHIX_SOALEN 4

#endif
