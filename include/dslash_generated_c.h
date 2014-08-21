#ifndef DSLASH_C_GENERATED
#define DSLASH_C_GENERATED

#if defined (MIC_SOURCE)
#include "mic/dslash_mic_complete_specialization.h"
// #include "mic/dslash_mic_face_complete_specialization.h"
#elif defined(AVX_SOURCE)
#include "avx/dslash_avx_complete_specialization.h"
// #include "avx/dslash_avx_face_complete_specialization.h"
#endif


#endif
