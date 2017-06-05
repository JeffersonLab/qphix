#pragma once

#include "qdp.h"

using namespace QDP;

template <typename T>
struct tolerance {
  static const Double small; // Always fail
};

template <>
const Double tolerance<half>::small = Double(5.0e-3);

template <>
const Double tolerance<float>::small = Double(1.0e-6);

template <>
const Double tolerance<double>::small = Double(1.0e-7);

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
