#pragma once

#include <qphix_codegen/decl_common.h>

#include "qdp.h"

template <typename T>
struct tolerance {
    static const QDP::Double small;
    static const double value;
};

template <typename T>
struct rsdTarget {
  static const double value;
};
