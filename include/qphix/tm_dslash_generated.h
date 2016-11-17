#ifndef QPHIX_TM_DSLASH_GENERATED_H
#define QPHIX_TM_DSLASH_GENERATED_H


#if defined (QPHIX_MIC_SOURCE)
#warning "including tm_dslash_mic_complete_specializations.h"
#include "qphix/mic/tm_dslash_mic_complete_specialization.h"

#elif defined (QPHIX_AVX512_SOURCE)
#warning "including tm_dslash_avx512_complete_specializations.h"
#include "qphix/avx512/tm_dslash_avx512_complete_specialization.h"

#elif defined(QPHIX_AVX_SOURCE)
#warning "including tm_dslash_avx_complete_specializations.h"
#include "qphix/avx/tm_dslash_avx_complete_specialization.h"

#elif defined(QPHIX_AVX2_SOURCE)
#warning "including tm_dslash_avx2_complete_specializations.h"
#include "qphix/avx2/tm_dslash_avx2_complete_specialization.h"

#elif defined(QPHIX_SCALAR_SOURCE)
#warning "including tm_dslash_scalar_complete_specializations.h"
#include "qphix/scalar/tm_dslash_scalar_complete_specialization.h"


#endif // SOURCES (ARCHs)
#endif // Include guard
