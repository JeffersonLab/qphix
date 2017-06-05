#pragma once

#include "qdp.h"

template <typename T>
struct tolerance {
  static const QDP::Double small; // Always fail
};

template <>
const QDP::Double tolerance<half>::small = QDP::Double(5.0e-3);

template <>
const QDP::Double tolerance<float>::small = QDP::Double(1.0e-6);

template <>
const QDP::Double tolerance<double>::small = QDP::Double(1.0e-7);

template <typename T>
struct rsdTarget {
  static const double value;
};

template <>
const double rsdTarget<half>::value = (double)(1.0e-3);

template <>
const double rsdTarget<float>::value = (double)(1.0e-7);

template <>
const double rsdTarget<double>::value = (double)(1.0e-12);
