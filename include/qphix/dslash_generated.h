#ifndef QPHIX_DSLASH_GENERATED_H
#define QPHIX_DSLASH_GENERATED_H

#if defined (QPHIX_MIC_SOURCE)
#include "qphix/mic/dslash_mic_complete_specialization.h"
#elif defined(QPHIX_AVX_SOURCE)
#warning including dslash_avx_complete_specializations.h"
#include "qphix/avx/dslash_avx_complete_specialization.h"
#endif


#endif
