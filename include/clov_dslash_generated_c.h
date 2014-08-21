#ifndef CLOV_DSLASH_C_GENERATED
#define CLOV_DSLASH_C_GENERATED

#if defined (MIC_SOURCE)
#include "mic/clov_dslash_mic_complete_specialization.h"
// #include "mic/dslash_mic_face_complete_specialization.h"

#elif defined(AVX_SOURCE)
#include "avx/clov_dslash_avx_complete_specialization.h"
//#include "avx/dslash_avx_face_complete_specialization.h"
#endif


#endif
