#include "qphix/tm_dslash_def.h"

#include "qphix/tm_dslash_body.h"

#ifdef QPHIX_QMP_COMMS
#include "qphix/tm_dslash_face.h"
#endif

namespace QPhiX {
template class TMDslash<double, 4, 1, false>;
template class TMDslash<double, 4, 1, true>;
template class TMDslash<double, 4, 2, false>;
template class TMDslash<double, 4, 2, true>;
template class TMDslash<double, 4, 4, false>;
template class TMDslash<double, 4, 4, true>;
template class TMDslash<float, 8, 1, false>;
template class TMDslash<float, 8, 1, true>;
template class TMDslash<float, 8, 2, false>;
template class TMDslash<float, 8, 2, true>;
template class TMDslash<float, 8, 4, false>;
template class TMDslash<float, 8, 4, true>;
template class TMDslash<float, 8, 8, false>;
template class TMDslash<float, 8, 8, true>;
}
