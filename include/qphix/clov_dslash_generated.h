// This file has been automatically generated. Do not change it manually, rather look
// for the template in qphix-codegen.

#pragma once

#include "qphix/diagnostics.h"

#if defined(QPHIX_AVX_SOURCE)
QPHIX_MESSAGE("including qphix/avx/clov_dslash_avx_complete_specialization.h")
#include "qphix/avx/clov_dslash_avx_complete_specialization.h"
#endif

#if defined(QPHIX_AVX2_SOURCE)
QPHIX_MESSAGE("including qphix/avx2/clov_dslash_avx2_complete_specialization.h")
#include "qphix/avx2/clov_dslash_avx2_complete_specialization.h"
#endif

#if defined(QPHIX_AVX512_SOURCE)
QPHIX_MESSAGE("including qphix/avx512/clov_dslash_avx512_complete_specialization.h")
#include "qphix/avx512/clov_dslash_avx512_complete_specialization.h"
#endif

#if defined(QPHIX_MIC_SOURCE)
QPHIX_MESSAGE("including qphix/mic/clov_dslash_mic_complete_specialization.h")
#include "qphix/mic/clov_dslash_mic_complete_specialization.h"
#endif

#if defined(QPHIX_QPX_SOURCE)
QPHIX_MESSAGE("including qphix/qpx/clov_dslash_qpx_complete_specialization.h")
#include "qphix/qpx/clov_dslash_qpx_complete_specialization.h"
#endif

#if defined(QPHIX_SCALAR_SOURCE)
QPHIX_MESSAGE("including qphix/scalar/clov_dslash_scalar_complete_specialization.h")
#include "qphix/scalar/clov_dslash_scalar_complete_specialization.h"
#endif

#if defined(QPHIX_SSE_SOURCE)
QPHIX_MESSAGE("including qphix/sse/clov_dslash_sse_complete_specialization.h")
#include "qphix/sse/clov_dslash_sse_complete_specialization.h"
#endif
