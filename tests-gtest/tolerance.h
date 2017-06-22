#pragma once

#include <qphix/geometry.h>

template <typename T>
struct tolerance {
  static const double small; // Always fail
};

template <>
const double tolerance<QPhiX::half>::small = double(5.0e-3);

template <>
const double tolerance<float>::small = double(1.0e-6);

template <>
const double tolerance<double>::small = double(1.0e-13);

template <typename T>
struct rsdTarget {
  static const double value;
};

template <>
const double rsdTarget<QPhiX::half>::value = (double)(1.0e-3);

template <>
const double rsdTarget<float>::value = (double)(1.0e-7);

template <>
const double rsdTarget<double>::value = (double)(1.0e-12);
