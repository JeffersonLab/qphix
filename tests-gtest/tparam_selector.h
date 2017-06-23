#pragma once

#include "veclen.h"
#include "prec.h"

#include <qphix/print_utils.h>

#include <qdp.h>
#include <qphix_codegen/decl_common.h>

using namespace QDP;

template <typename FT>
struct QdpTypes {
};

template <>
struct QdpTypes<QPhiX::half> {
  typedef LatticeColorMatrixF QdpGauge;
  typedef LatticeDiracFermionF QdpSpinor;
};

template <>
struct QdpTypes<float> {
  typedef LatticeColorMatrixF QdpGauge;
  typedef LatticeDiracFermionF QdpSpinor;
};

template <>
struct QdpTypes<double> {
  typedef LatticeColorMatrixD QdpGauge;
  typedef LatticeDiracFermionD QdpSpinor;
};

template <class TestClass>
static void call(TestClass &instance, Prec prec, int soalen, bool compress12);

template <class TestClass, typename FT>
static void call_1(TestClass &instance, int soalen, bool compress12);

template <int soalen>
struct Call2 {
  template <class TestClass, typename FT>
  static void call_2(TestClass &instance, bool compress12);
};

template <class TestClass>
void call(TestClass &instance, Prec prec, int soalen, bool compress12)
{
  if (prec == HALF_PREC) {
#if defined(QPHIX_MIC_SOURCE) || defined(QPHIX_AVX2_SOURCE) ||                      \
    defined(QPHIX_AVX512_SOURCE)
    call_1<TestClass, QPhiX::half>(instance, soalen, compress12);
#else
    QPhiX::masterPrintf("This architecture does not support half-precision floating "
                        "point operations.\n");
    std::abort();
#endif
  } else if (prec == FLOAT_PREC) {
    call_1<TestClass, float>(instance, soalen, compress12);
  } else {
    call_1<TestClass, double>(instance, soalen, compress12);
  }
}

template <class TestClass, typename FT>
void call_1(TestClass &instance, int soalen, bool compress12)
{
  constexpr int veclen = get_veclen<FT>();
  constexpr int veclen_half = veclen / 2;

  if (soalen == veclen) {
    Call2<veclen>::template call_2<TestClass, FT>(instance, compress12);
  } else if (soalen == veclen_half) {
    Call2<veclen_half>::template call_2<TestClass, FT>(instance, compress12);
  }
  // On AVX512, we also have a `soalen == 4` variant for `half` and `double`,
  // where `veclen == 16`. Enable this as well.
#if defined(QPHIX_MIC_SOURCE) || defined(QPHIX_AVX512_SOURCE)
  else if (veclen == 16 && soalen == 4) {
    Call2<4>::template call_2<TestClass, FT>(instance, compress12);
  }
#endif
  else {
    QPhiX::masterPrintf(
        "soalen %d is not implemented for a floating point type of size %d.\n",
        soalen,
        sizeof(FT));
    std::abort();
  }
}

template <int soalen>
template <class TestClass, typename FT>
void Call2<soalen>::call_2(TestClass &instance, bool compress12)
{
  constexpr int veclen = get_veclen<FT>();

  if (compress12) {
    instance.template operator()<FT,
                                 veclen,
                                 soalen,
                                 true,
                                 typename QdpTypes<FT>::QdpGauge,
                                 typename QdpTypes<FT>::QdpSpinor>();
  } else {
    instance.template operator()<FT,
                                 veclen,
                                 soalen,
                                 false,
                                 typename QdpTypes<FT>::QdpGauge,
                                 typename QdpTypes<FT>::QdpSpinor>();
  }
}

template <>
template <class TestClass, typename FT>
void Call2<0>::call_2(TestClass &instance, bool compress12){};
