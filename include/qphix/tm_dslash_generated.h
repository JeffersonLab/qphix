#ifndef QPHIX_TM_DSLASH_GENERATED_H
#define QPHIX_TM_DSLASH_GENERATED_H

#if defined (QPHIX_MIC_SOURCE)
#include "qphix/mic/tm_dslash_mic_complete_specialization.h"
#elif defined(QPHIX_AVX_SOURCE)
#warning "including tm_dslash_avx_complete_specializations.h"
#include "qphix/avx/tm_dslash_avx_complete_specialization.h"
#elif defined(QPHIX_SCALAR_SOURCE)
#warning "Scalar source is not implemented, including tm_dslash_avx_complete_specializations.h"
#include "qphix/avx/tm_dslash_avx_complete_specialization.h"
#endif


#endif
