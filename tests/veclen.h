#pragma once

#ifndef QPHIX_SOALEN
#define QPHIX_SOALEN 4
#endif

#if defined(QPHIX_MIC_SOURCE) || defined(QPHIX_AVX512_SOURCE)

#define VECLEN_SP 16
#define VECLEN_HP 16
#define VECLEN_DP 8

#elif defined(QPHIX_AVX_SOURCE) || defined(QPHIX_AVX2_SOURCE)

#define VECLEN_SP 8
#define VECLEN_HP 8
#define VECLEN_DP 4

#elif defined(QPHIX_SSE_SOURCE)

#define VECLEN_SP 4
#define VECLEN_DP 2

#elif defined(QPHIX_SCALAR_SOURCE)
#define VECLEN_DP 1
#define VECLEN_SP 1
#ifdef QPHIX_SOALEN
#undef QPHIX_SOALEN
#endif
#define QPHIX_SOALEN 1
#elif defined(QPHIX_QPX_SOURCE)
#define VECLEN_DP 4
#ifdef QPHIX_SOALEN
#undef QPHIX_SOALEN
#endif 
#define QPHIX_SOALEN 4

#endif

#include "qphix/geometry.h"

#include "qdp.h"

template <typename FT>
constexpr int get_veclen()
{
  return 0;
};

#if defined(QPHIX_MIC_SOURCE) || defined(QPHIX_AVX512_SOURCE)
template <>
constexpr int get_veclen<QPhiX::half>()
{
  return 16;
}
template <>
constexpr int get_veclen<float>()
{
  return 16;
}
template <>
constexpr int get_veclen<double>()
{
  return 8;
}

#elif defined(QPHIX_AVX_SOURCE)
template <>
constexpr int get_veclen<QPhiX::half>()
{
  return 8;
}
template <>
constexpr int get_veclen<float>()
{
  return 8;
}
template <>
constexpr int get_veclen<double>()
{
  return 4;
}

#elif defined(QPHIX_AVX2_SOURCE)
template <>
constexpr int get_veclen<float>()
{
  return 8;
}
template <>
constexpr int get_veclen<double>()
{
  return 4;
}

#elif defined(QPHIX_SSE_SOURCE)
template <>
constexpr int get_veclen<float>()
{
  return 4;
}
template <>
constexpr int get_veclen<double>()
{
  return 2;
}

#elif defined(QPHIX_SCALAR_SOURCE)
template <>
constexpr int get_veclen<float>()
{
  return 1;
}
template <>
constexpr int get_veclen<double>()
{
  return 1;
}

#elif defined(QPHIX_QPX_SOURCE)
template <>
constexpr int get_veclen<double>()
{
  return 4;
}

#endif
