#ifndef QPHIX_DSLASH_GENERATED_H
#define QPHIX_DSLASH_GENERATED_H

#if defined (QPHIX_MIC_SOURCE)
#include "qphix/mic/dslash_mic_complete_specialization.h"

#elif defined(QPHIX_AVX512_SOURCE)
#include "qphix/avx512/dslash_avx512_complete_specialization.h"

#elif defined(QPHIX_AVX_SOURCE)
#warning "including dslash_avx_complete_specializations.h"
#include "qphix/avx/dslash_avx_complete_specialization.h"

#elif defined(QPHIX_AVX2_SOURCE)
#warning "including dslash_avx2_complete_specializations.h"
#include "qphix/avx2/dslash_avx2_complete_specialization.h"

#elif defined(QPHIX_SCALAR_SOURCE)
#include "qphix/scalar/dslash_scalar_complete_specialization.h"

#elif defined(QPHIX_QPX_SOURCE)
#include "qphix/qpx/dslash_qpx_complete_specialization.h"

#elif defined(QPHIX_SSE_SOURCE)
#include "qphix/sse/dslash_sse_complete_specialization.h"

#endif


#endif
