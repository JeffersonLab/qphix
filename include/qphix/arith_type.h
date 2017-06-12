#pragma once

#include "qphix/geometry.h"

namespace QPhiX
{
// Typically the arithmetic type is the same as the type
template <typename FT>
struct ArithType {
  typedef FT Type;
};

// But sometimes it is not.
template <>
struct ArithType<half> {
  typedef float Type;
};
}
