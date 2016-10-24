#ifndef QPHIX_TM_CLOV_DSLASH_C_GENERATED
#define QPHIX_TM_CLOV_DSLASH_C_GENERATED


#if defined (QPHIX_MIC_SOURCE)
#include "qphix/mic/tm_clov_dslash_mic_complete_specialization.h"

#elif defined (QPHIX_AVX512_SOURCE)
#include "qphix/avx512/tm_clov_dslash_avx512_complete_specialization.h"

#elif defined(QPHIX_AVX_SOURCE)
#include "qphix/avx/tm_clov_dslash_avx_complete_specialization.h"

#elif defined(QPHIX_AVX2_SOURCE)
#include "qphix/avx2/tm_clov_dslash_avx2_complete_specialization.h"

#elif defined(QPHIX_SCALAR_SOURCE)
#include "qphix/scalar/tm_clov_dslash_scalar_complete_specialization.h"

#elif defined(QPHIX_QPX_SOURCE)
#include "qphix/qpx/tm_clov_dslash_qpx_complete_specialization.h"

#elif defined(QPHIX_SSE_SOURCE)
#include "qphix/sse/tm_clov_dslash_sse_complete_specialization.h"

#endif // SOURCES (ARCHs)


#endif // Include guard
