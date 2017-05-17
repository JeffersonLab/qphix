// This file has been automatically generated. Do not change it manually, rather look for the template in qphix-codegen.

#include "qphix/dslash_def.h"

#include "qphix/dslash_body.h"

#ifdef QPHIX_DO_COMMS
// Disable comms for nokw
#include "qphix/clov_face.h"
#endif

namespace QPhiX {


#if defined (QPHIX_AVX_SOURCE)
template class Dslash<double, 4, 2, true>;
template class Dslash<double, 4, 2, false>;
template class Dslash<double, 4, 4, true>;
template class Dslash<double, 4, 4, false>;
template class Dslash<float, 8, 4, true>;
template class Dslash<float, 8, 4, false>;
template class Dslash<float, 8, 8, true>;
template class Dslash<float, 8, 8, false>;
#endif

#if defined (QPHIX_AVX2_SOURCE)
template class Dslash<double, 4, 2, true>;
template class Dslash<double, 4, 2, false>;
template class Dslash<double, 4, 4, true>;
template class Dslash<double, 4, 4, false>;
template class Dslash<float, 8, 4, true>;
template class Dslash<float, 8, 4, false>;
template class Dslash<float, 8, 8, true>;
template class Dslash<float, 8, 8, false>;
template class Dslash<half, 8, 4, true>;
template class Dslash<half, 8, 4, false>;
template class Dslash<half, 8, 8, true>;
template class Dslash<half, 8, 8, false>;
#endif

#if defined (QPHIX_AVX512_SOURCE)
template class Dslash<double, 8, 4, true>;
template class Dslash<double, 8, 4, false>;
template class Dslash<double, 8, 8, true>;
template class Dslash<double, 8, 8, false>;
template class Dslash<float, 16, 4, true>;
template class Dslash<float, 16, 4, false>;
template class Dslash<float, 16, 8, true>;
template class Dslash<float, 16, 8, false>;
template class Dslash<float, 16, 16, true>;
template class Dslash<float, 16, 16, false>;
template class Dslash<half, 16, 4, true>;
template class Dslash<half, 16, 4, false>;
template class Dslash<half, 16, 8, true>;
template class Dslash<half, 16, 8, false>;
template class Dslash<half, 16, 16, true>;
template class Dslash<half, 16, 16, false>;
#endif

#if defined (QPHIX_MIC_SOURCE)
template class Dslash<double, 8, 4, true>;
template class Dslash<double, 8, 4, false>;
template class Dslash<double, 8, 8, true>;
template class Dslash<double, 8, 8, false>;
template class Dslash<float, 16, 4, true>;
template class Dslash<float, 16, 4, false>;
template class Dslash<float, 16, 8, true>;
template class Dslash<float, 16, 8, false>;
template class Dslash<float, 16, 16, true>;
template class Dslash<float, 16, 16, false>;
template class Dslash<half, 16, 4, true>;
template class Dslash<half, 16, 4, false>;
template class Dslash<half, 16, 8, true>;
template class Dslash<half, 16, 8, false>;
template class Dslash<half, 16, 16, true>;
template class Dslash<half, 16, 16, false>;
#endif

#if defined (QPHIX_QPX_SOURCE)
template class Dslash<double, 4, 4, true>;
template class Dslash<double, 4, 4, false>;
#endif

#if defined (QPHIX_SCALAR_SOURCE)
template class Dslash<double, 1, 1, true>;
template class Dslash<double, 1, 1, false>;
template class Dslash<float, 1, 1, true>;
template class Dslash<float, 1, 1, false>;
#endif

#if defined (QPHIX_SSE_SOURCE)
template class Dslash<double, 2, 2, true>;
template class Dslash<double, 2, 2, false>;
template class Dslash<float, 4, 4, true>;
template class Dslash<float, 4, 4, false>;
#endif

}