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

template <>
const QDP::Double tolerance<QPhiX::half>::small = QDP::Double(5.0e-3);
template <>
const double tolerance<QPhiX::half>::value = 5.0e-3;

template <>
const QDP::Double tolerance<float>::small = QDP::Double(1.0e-6);
template <>
const double tolerance<float>::value = 1.0e-6;

template <>
const QDP::Double tolerance<double>::small = QDP::Double(1.0e-7);
template <>
const double tolerance<double>::value = 1.0e-7;

template <>
const double rsdTarget<QPhiX::half>::value = (double)(1.0e-3);

template <>
const double rsdTarget<float>::value = (double)(1.0e-7);

template <>
const double rsdTarget<double>::value = (double)(1.0e-12);
