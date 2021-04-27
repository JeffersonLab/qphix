#pragma once

#include "qphix/blas_c.h"
#include "qphix/qphix_config.h"
#include "qphix/blas_utils.h"
#include "qphix/print_utils.h"
#include "qphix/arith_type.h"
#include <array>

namespace QPhiX
{

template <typename FT, int V, int S, bool compress>
class CopyFunctor
{
 public:
  CopyFunctor(typename Geometry<FT, V, S, compress>::FourSpinorBlock *res_,
              const typename Geometry<FT, V, S, compress>::FourSpinorBlock *src_)
      : res(res_), src(src_)
  {
  }

  CopyFunctor(const CopyFunctor<FT, V, S, compress> &rhs)
      : res(rhs.res), src(rhs.src)
  {
  }

  ~CopyFunctor() {}

  inline void func(int block) const
  {
    int nvec_in_spinor = (3 * 4 * 2 * S) / V;
    const FT *srcbase = &src[block][0][0][0][0];
    FT *resbase = &res[block][0][0][0][0];

#if defined(QPHIX_MIC_SOURCE)
    // Intel MIC
    const int prefdist1 = 12;
    const char *pref1base = (const char *)srcbase + prefdist1 * 64;

    const int prefdist2 = 64;
    const char *pref2base = (const char *)srcbase + prefdist2 * 64;

    for (int numv = 0; numv < nvec_in_spinor; numv++) {
      _mm_prefetch(&pref1base[numv * V * sizeof(FT)], _MM_HINT_T0);
      _mm_prefetch(&pref2base[numv * V * sizeof(FT)], _MM_HINT_T1);

#pragma omp simd aligned(resbase, srcbase : V)
      for (int s = 0; s < V; s++) {
        resbase[numv * V + s] = srcbase[numv * V + s];
      }
    }
#else
    // Generic
    for (int numv = 0; numv < nvec_in_spinor; numv++) {
#pragma omp simd aligned(resbase, srcbase : V)
      for (int s = 0; s < V; s++) {
        resbase[numv * V + s] = srcbase[numv * V + s];
      }
    }
#endif
  }

 private:
  typename Geometry<FT, V, S, compress>::FourSpinorBlock *res;
  const typename Geometry<FT, V, S, compress>::FourSpinorBlock *src;
};

template <typename FT, int V, int S, bool compress>
class ZeroFunctor
{
 public:
  typedef typename ArithType<FT>::Type AT;

  ZeroFunctor(typename Geometry<FT, V, S, compress>::FourSpinorBlock *res_)
      : res(res_)
  {
  }
  ~ZeroFunctor() {}

  inline void func(int block) const
  {
    int nvec_in_spinor = (3 * 4 * 2 * S) / V;
    FT *resbase = &res[block][0][0][0][0];

#if defined(__GNUG__) && !defined(__INTEL_COMPILER)
    typename Geometry<AT, V, S, compress>::FourSpinorBlock res_spinor
        __attribute__((aligned(QPHIX_LLC_CACHE_ALIGN)));
#else
    __declspec(align(QPHIX_LLC_CACHE_ALIGN))
        typename Geometry<AT, V, S, compress>::FourSpinorBlock res_spinor;
#endif

    // Now we are hopefully both in L1 and in the right layout so
    for (int col = 0; col < 3; col++) {
      for (int spin = 0; spin < 4; spin++) {
        for (int reim = 0; reim < 2; reim++) {
#pragma omp simd aligned(res_spinor : QPHIX_LLC_CACHE_ALIGN)
          for (int i = 0; i < S; i++) {
            res_spinor[col][spin][reim][i] = rep<AT, double>((double)0);
          }
        }
      }
    }
    BLASUtils::streamOutSpinor<FT, V>(
        resbase, (const AT *)res_spinor, nvec_in_spinor);
  }

 private:
  typename Geometry<FT, V, S, compress>::FourSpinorBlock *res;
};

template <typename FT, int V, int S, bool compress>
class AYPXFunctor
{
 public:
  typedef typename ArithType<FT>::Type AT;

  AYPXFunctor(double a_,
              const typename Geometry<FT, V, S, compress>::FourSpinorBlock *x_,
              typename Geometry<FT, V, S, compress>::FourSpinorBlock *y_)
      : a(rep<AT, double>(a_)), x(x_), y(y_)
  {
  }

  AYPXFunctor(const AYPXFunctor<FT, V, S, compress> &rhs)
      : a(rhs.a), x(rhs.x), y(rhs.y)
  {
  }

  ~AYPXFunctor() {}

  inline void func(int block) const
  {
    int nvec_in_spinor = (3 * 4 * 2 * S) / V;
    const FT *xbase = &x[block][0][0][0][0];
    FT *ybase = &y[block][0][0][0][0];

// Temporary storage to stream into and out of
#if defined(__GNUG__) && !defined(__INTEL_COMPILER)
    typename Geometry<AT, V, S, compress>::FourSpinorBlock x_spinor
        __attribute__((aligned(QPHIX_LLC_CACHE_ALIGN)));
    typename Geometry<AT, V, S, compress>::FourSpinorBlock y_spinor
        __attribute__((aligned(QPHIX_LLC_CACHE_ALIGN)));
#else
    __declspec(align(QPHIX_LLC_CACHE_ALIGN))
        typename Geometry<AT, V, S, compress>::FourSpinorBlock x_spinor;
    __declspec(align(QPHIX_LLC_CACHE_ALIGN))
        typename Geometry<AT, V, S, compress>::FourSpinorBlock y_spinor;
#endif

    BLASUtils::streamInSpinor<FT, V>((AT *)x_spinor, xbase, nvec_in_spinor);
    BLASUtils::streamInSpinor<FT, V>((AT *)y_spinor, ybase, nvec_in_spinor);

    // Now we are hopefully both in L1 and in the right layout so
    for (int col = 0; col < 3; col++) {
      for (int spin = 0; spin < 4; spin++) {
        for (int reim = 0; reim < 2; reim++) {
#pragma omp simd aligned(y_spinor, x_spinor : QPHIX_LLC_CACHE_ALIGN)
          for (int i = 0; i < S; i++) {
            y_spinor[col][spin][reim][i] =
                a * y_spinor[col][spin][reim][i] + x_spinor[col][spin][reim][i];
          }
        }
      }
    }

    BLASUtils::streamOutSpinor<FT, V>(ybase, (const AT *)y_spinor, nvec_in_spinor);
  }

 private:
  AT a;
  const typename Geometry<FT, V, S, compress>::FourSpinorBlock *x;
  typename Geometry<FT, V, S, compress>::FourSpinorBlock *y;
};

template <typename FT, int V, int S, bool compress>
class YPEQXFunctor
{
 public:
  typedef typename ArithType<FT>::Type AT;

  YPEQXFunctor(const typename Geometry<FT, V, S, compress>::FourSpinorBlock *x_,
              typename Geometry<FT, V, S, compress>::FourSpinorBlock *y_)
      :  x(x_), y(y_)
  {
  }

  YPEQXFunctor(const YPEQXFunctor<FT, V, S, compress> &rhs)
      : x(rhs.x), y(rhs.y)
  {
  }

  ~YPEQXFunctor() {}

  inline void func(int block) const
  {
    int nvec_in_spinor = (3 * 4 * 2 * S) / V;
    const FT *xbase = &x[block][0][0][0][0];
    FT *ybase = &y[block][0][0][0][0];

// Temporary storage to stream into and out of
#if defined(__GNUG__) && !defined(__INTEL_COMPILER)
    typename Geometry<AT, V, S, compress>::FourSpinorBlock x_spinor
        __attribute__((aligned(QPHIX_LLC_CACHE_ALIGN)));
    typename Geometry<AT, V, S, compress>::FourSpinorBlock y_spinor
        __attribute__((aligned(QPHIX_LLC_CACHE_ALIGN)));
#else
    __declspec(align(QPHIX_LLC_CACHE_ALIGN))
        typename Geometry<AT, V, S, compress>::FourSpinorBlock x_spinor;
    __declspec(align(QPHIX_LLC_CACHE_ALIGN))
        typename Geometry<AT, V, S, compress>::FourSpinorBlock y_spinor;
#endif

    BLASUtils::streamInSpinor<FT, V>((AT *)x_spinor, xbase, nvec_in_spinor);
    BLASUtils::streamInSpinor<FT, V>((AT *)y_spinor, ybase, nvec_in_spinor);

    // Now we are hopefully both in L1 and in the right layout so
    for (int col = 0; col < 3; col++) {
      for (int spin = 0; spin < 4; spin++) {
        for (int reim = 0; reim < 2; reim++) {
#pragma omp simd aligned(y_spinor, x_spinor : QPHIX_LLC_CACHE_ALIGN)
          for (int i = 0; i < S; i++) {
            y_spinor[col][spin][reim][i] += x_spinor[col][spin][reim][i];
          }
        }
      }
    }

    BLASUtils::streamOutSpinor<FT, V>(ybase, (const AT *)y_spinor, nvec_in_spinor);
  }

 private:
  const typename Geometry<FT, V, S, compress>::FourSpinorBlock *x;
  typename Geometry<FT, V, S, compress>::FourSpinorBlock *y;
};

template <typename FT, int V, int S, bool compress>
class YMEQXFunctor
{
 public:
  typedef typename ArithType<FT>::Type AT;

  YMEQXFunctor(const typename Geometry<FT, V, S, compress>::FourSpinorBlock *x_,
              typename Geometry<FT, V, S, compress>::FourSpinorBlock *y_)
      :  x(x_), y(y_)
  {
  }

  YMEQXFunctor(const YPEQXFunctor<FT, V, S, compress> &rhs)
      : x(rhs.x), y(rhs.y)
  {
  }

  ~YMEQXFunctor() {}

  inline void func(int block) const
  {
    int nvec_in_spinor = (3 * 4 * 2 * S) / V;
    const FT *xbase = &x[block][0][0][0][0];
    FT *ybase = &y[block][0][0][0][0];

// Temporary storage to stream into and out of
#if defined(__GNUG__) && !defined(__INTEL_COMPILER)
    typename Geometry<AT, V, S, compress>::FourSpinorBlock x_spinor
        __attribute__((aligned(QPHIX_LLC_CACHE_ALIGN)));
    typename Geometry<AT, V, S, compress>::FourSpinorBlock y_spinor
        __attribute__((aligned(QPHIX_LLC_CACHE_ALIGN)));
#else
    __declspec(align(QPHIX_LLC_CACHE_ALIGN))
        typename Geometry<AT, V, S, compress>::FourSpinorBlock x_spinor;
    __declspec(align(QPHIX_LLC_CACHE_ALIGN))
        typename Geometry<AT, V, S, compress>::FourSpinorBlock y_spinor;
#endif

    BLASUtils::streamInSpinor<FT, V>((AT *)x_spinor, xbase, nvec_in_spinor);
    BLASUtils::streamInSpinor<FT, V>((AT *)y_spinor, ybase, nvec_in_spinor);

    // Now we are hopefully both in L1 and in the right layout so
    for (int col = 0; col < 3; col++) {
      for (int spin = 0; spin < 4; spin++) {
        for (int reim = 0; reim < 2; reim++) {
#pragma omp simd aligned(y_spinor, x_spinor : QPHIX_LLC_CACHE_ALIGN)
          for (int i = 0; i < S; i++) {
            y_spinor[col][spin][reim][i] -= x_spinor[col][spin][reim][i];
          }
        }
      }
    }

    BLASUtils::streamOutSpinor<FT, V>(ybase, (const AT *)y_spinor, nvec_in_spinor);
  }

 private:
  const typename Geometry<FT, V, S, compress>::FourSpinorBlock *x;
  typename Geometry<FT, V, S, compress>::FourSpinorBlock *y;
};


template <typename FT, int V, int S, bool compress>
class AXPYFunctor
{
 public:
  typedef typename ArithType<FT>::Type AT;

  AXPYFunctor(double a_,
              const typename Geometry<FT, V, S, compress>::FourSpinorBlock *x_,
              typename Geometry<FT, V, S, compress>::FourSpinorBlock *y_)
      : a(rep<AT, double>(a_)), x(x_), y(y_)
  {
  }

  AXPYFunctor(const AXPYFunctor<FT, V, S, compress> &rhs)
      : a(rhs.a), x(rhs.x), y(rhs.y)
  {
  }

  ~AXPYFunctor() {}

  inline void func(int block) const
  {
    int nvec_in_spinor = (3 * 4 * 2 * S) / V;
    const FT *xbase = &x[block][0][0][0][0];
    FT *ybase = &y[block][0][0][0][0];

// Temporary storage to stream into and out of
#if defined(__GNUG__) && !defined(__INTEL_COMPILER)
    typename Geometry<AT, V, S, compress>::FourSpinorBlock x_spinor
        __attribute__((aligned(QPHIX_LLC_CACHE_ALIGN)));
    typename Geometry<AT, V, S, compress>::FourSpinorBlock y_spinor
        __attribute__((aligned(QPHIX_LLC_CACHE_ALIGN)));
#else
    __declspec(align(QPHIX_LLC_CACHE_ALIGN))
        typename Geometry<AT, V, S, compress>::FourSpinorBlock x_spinor;
    __declspec(align(QPHIX_LLC_CACHE_ALIGN))
        typename Geometry<AT, V, S, compress>::FourSpinorBlock y_spinor;
#endif

    BLASUtils::streamInSpinor<FT, V>((AT *)x_spinor, xbase, nvec_in_spinor);
    BLASUtils::streamInSpinor<FT, V>((AT *)y_spinor, ybase, nvec_in_spinor);

    // Now we are hopefully both in L1 and in the right layout so
    for (int col = 0; col < 3; col++) {
      for (int spin = 0; spin < 4; spin++) {
        for (int reim = 0; reim < 2; reim++) {
#pragma omp simd aligned(y_spinor, x_spinor : QPHIX_LLC_CACHE_ALIGN)
          for (int i = 0; i < S; i++) {
            y_spinor[col][spin][reim][i] =
                a * x_spinor[col][spin][reim][i] + y_spinor[col][spin][reim][i];
          }
        }
      }
    }

    BLASUtils::streamOutSpinor<FT, V>(ybase, (const AT *)y_spinor, nvec_in_spinor);
  }

 private:
  AT a;
  const typename Geometry<FT, V, S, compress>::FourSpinorBlock *x;
  typename Geometry<FT, V, S, compress>::FourSpinorBlock *y;
};

template <typename FT, int V, int S, bool compress>
class AXPBYFunctor
{
 public:
  typedef typename ArithType<FT>::Type AT;

  AXPBYFunctor(double a_,
               const typename Geometry<FT, V, S, compress>::FourSpinorBlock *x_,
               double b_,
               typename Geometry<FT, V, S, compress>::FourSpinorBlock *y_)
      : a(rep<AT, double>(a_)), b(rep<AT, double>(b_)), x(x_), y(y_)
  {
  }
  AXPBYFunctor(const AXPBYFunctor<FT, V, S, compress> &rhs)
      : a(rhs.a), b(rhs.b), x(rhs.x), y(rhs.y)
  {
  }
  ~AXPBYFunctor() {}

  inline void func(int block) const
  {
    int nvec_in_spinor = (3 * 4 * 2 * S) / V;
    const FT *xbase = &x[block][0][0][0][0];
    FT *ybase = &y[block][0][0][0][0];

// Temporary storage to stream into and out of
#if defined(__GNUG__) && !defined(__INTEL_COMPILER)
    typename Geometry<AT, V, S, compress>::FourSpinorBlock x_spinor
        __attribute__((aligned(QPHIX_LLC_CACHE_ALIGN)));
    typename Geometry<AT, V, S, compress>::FourSpinorBlock y_spinor
        __attribute__((aligned(QPHIX_LLC_CACHE_ALIGN)));
#else
    __declspec(align(QPHIX_LLC_CACHE_ALIGN))
        typename Geometry<AT, V, S, compress>::FourSpinorBlock x_spinor;
    __declspec(align(QPHIX_LLC_CACHE_ALIGN))
        typename Geometry<AT, V, S, compress>::FourSpinorBlock y_spinor;
#endif

    BLASUtils::streamInSpinor<FT, V>((AT *)x_spinor, xbase, nvec_in_spinor);
    BLASUtils::streamInSpinor<FT, V>((AT *)y_spinor, ybase, nvec_in_spinor);

    // Now we are hopefully both in L1 and in the right layout so
    for (int col = 0; col < 3; col++) {
      for (int spin = 0; spin < 4; spin++) {
        for (int reim = 0; reim < 2; reim++) {
#pragma omp simd aligned(y_spinor, x_spinor : QPHIX_LLC_CACHE_ALIGN)
          for (int i = 0; i < S; i++) {
            y_spinor[col][spin][reim][i] =
                a * x_spinor[col][spin][reim][i] + b * y_spinor[col][spin][reim][i];
          }
        }
      }
    }
    BLASUtils::streamOutSpinor<FT, V>(ybase, (const AT *)y_spinor, nvec_in_spinor);
  }

 private:
  AT a;
  AT b;
  const typename Geometry<FT, V, S, compress>::FourSpinorBlock *x;
  typename Geometry<FT, V, S, compress>::FourSpinorBlock *y;
};

template <typename FT, int V, int S, bool compress>
class AXYFunctor
{
 public:
  typedef typename ArithType<FT>::Type AT;

  AXYFunctor(double a_,
             const typename Geometry<FT, V, S, compress>::FourSpinorBlock *x_,
             typename Geometry<FT, V, S, compress>::FourSpinorBlock *y_)
      : a(rep<AT, double>(a_)), x(x_), y(y_)
  {
  }
  ~AXYFunctor() {}

  inline void func(int block) const
  {
    int nvec_in_spinor = (3 * 4 * 2 * S) / V;
    const FT *xbase = &x[block][0][0][0][0];
    FT *ybase = &y[block][0][0][0][0];

// Temporary storage to stream into and out of
#if defined(__GNUG__) && !defined(__INTEL_COMPILER)
    typename Geometry<AT, V, S, compress>::FourSpinorBlock x_spinor
        __attribute__((aligned(QPHIX_LLC_CACHE_ALIGN)));
    typename Geometry<AT, V, S, compress>::FourSpinorBlock y_spinor
        __attribute__((aligned(QPHIX_LLC_CACHE_ALIGN)));
#else
    __declspec(align(QPHIX_LLC_CACHE_ALIGN))
        typename Geometry<AT, V, S, compress>::FourSpinorBlock x_spinor;
    __declspec(align(QPHIX_LLC_CACHE_ALIGN))
        typename Geometry<AT, V, S, compress>::FourSpinorBlock y_spinor;
#endif

    BLASUtils::streamInSpinor<FT, V>((AT *)x_spinor, xbase, nvec_in_spinor);

    // Now we are hopefully both in L1 and in the right layout so
    for (int col = 0; col < 3; col++) {
      for (int spin = 0; spin < 4; spin++) {
        for (int reim = 0; reim < 2; reim++) {
#pragma omp simd aligned(y_spinor, x_spinor : QPHIX_LLC_CACHE_ALIGN)
          for (int i = 0; i < S; i++) {
            y_spinor[col][spin][reim][i] = a * x_spinor[col][spin][reim][i];
          }
        }
      }
    }

    BLASUtils::streamOutSpinor<FT, V>(ybase, (const AT *)y_spinor, nvec_in_spinor);
  }

 private:
  AT a;
  const typename Geometry<FT, V, S, compress>::FourSpinorBlock *x;
  typename Geometry<FT, V, S, compress>::FourSpinorBlock *y;
};

template <typename FT, int V, int S, bool compress>
class AXFunctor
{
 public:
  typedef typename ArithType<FT>::Type AT;

  AXFunctor(double a_,
            typename Geometry<FT, V, S, compress>::FourSpinorBlock *x_)
      : a(rep<AT, double>(a_)), x(x_)
  {
  }
  ~AXFunctor() {}

  inline void func(int block) const
  {
    int nvec_in_spinor = (3 * 4 * 2 * S) / V;
    FT *xbase = &x[block][0][0][0][0];

// Temporary storage to stream into and out of
#if defined(__GNUG__) && !defined(__INTEL_COMPILER)
    typename Geometry<AT, V, S, compress>::FourSpinorBlock x_spinor
        __attribute__((aligned(QPHIX_LLC_CACHE_ALIGN)));

#else
    __declspec(align(QPHIX_LLC_CACHE_ALIGN))
        typename Geometry<AT, V, S, compress>::FourSpinorBlock x_spinor;

#endif

    BLASUtils::streamInSpinor<FT, V>((AT *)x_spinor, xbase, nvec_in_spinor);

    // Now we are hopefully both in L1 and in the right layout so
    for (int col = 0; col < 3; col++) {
      for (int spin = 0; spin < 4; spin++) {
        for (int reim = 0; reim < 2; reim++) {
#pragma omp simd aligned(x_spinor : QPHIX_LLC_CACHE_ALIGN)
          for (int i = 0; i < S; i++) {
            x_spinor[col][spin][reim][i] *= a;
          }
        }
      }
    }

    BLASUtils::streamOutSpinor<FT, V>(xbase, (const AT *)x_spinor, nvec_in_spinor);
  }

 private:
  AT a;
  typename Geometry<FT, V, S, compress>::FourSpinorBlock *x;
};


template <typename FT, int V, int S, bool compress>
class Norm2Functor
{
 public:
  typedef typename ArithType<FT>::Type AT;

  Norm2Functor(const typename Geometry<FT, V, S, compress>::FourSpinorBlock *x_)
      : x(x_)
  {
  }
  ~Norm2Functor() {}

  inline void func(int block, std::array<double, S> &reduction) const
  {
    int nvec_in_spinor = (3 * 4 * 2 * S) / V;
    const FT *xbase = &x[block][0][0][0][0];

#if defined(__GNUG__) && !defined(__INTEL_COMPILER)
    typename Geometry<AT, V, S, compress>::FourSpinorBlock x_spinor
        __attribute__((aligned(QPHIX_LLC_CACHE_ALIGN)));
#else
    __declspec(align(QPHIX_LLC_CACHE_ALIGN))
        typename Geometry<AT, V, S, compress>::FourSpinorBlock x_spinor;
#endif

    BLASUtils::streamInSpinor<FT, V>((AT *)x_spinor, xbase, nvec_in_spinor);

    for (int col = 0; col < 3; col++) {
      for (int spin = 0; spin < 4; spin++) {
        for (int reim = 0; reim < 2; reim++) {
          for (int s = 0; s < S; s++) {

            double xfoo = rep<double, AT>(x_spinor[col][spin][reim][s]);
            reduction[s] += xfoo * xfoo;
          }
        }
      }
    }
  }

 private:
  const typename Geometry<FT, V, S, compress>::FourSpinorBlock *x;
};

template <typename FT, int V, int S, bool compress>
class XMYNorm2Functor
{
 public:
  typedef typename ArithType<FT>::Type AT;

  XMYNorm2Functor(typename Geometry<FT, V, S, compress>::FourSpinorBlock *res_,
                  const typename Geometry<FT, V, S, compress>::FourSpinorBlock *x_,
                  const typename Geometry<FT, V, S, compress>::FourSpinorBlock *y_)
      : res(res_), x(x_), y(y_)
  {
  }

  ~XMYNorm2Functor() {}

  inline void func(int block, std::array<double, S> &red) const
  {
    int nvec_in_spinor = (3 * 4 * 2 * S) / V;
    const FT *xbase = &x[block][0][0][0][0];
    const FT *ybase = &y[block][0][0][0][0];
    FT *resbase = &res[block][0][0][0][0];

// Temporary storage to stream into and out of
#if defined(__GNUG__) && !defined(__INTEL_COMPILER)
    typename Geometry<AT, V, S, compress>::FourSpinorBlock x_spinor
        __attribute__((aligned(QPHIX_LLC_CACHE_ALIGN)));
    typename Geometry<AT, V, S, compress>::FourSpinorBlock y_spinor
        __attribute__((aligned(QPHIX_LLC_CACHE_ALIGN)));
    typename Geometry<AT, V, S, compress>::FourSpinorBlock res_spinor
        __attribute__((aligned(QPHIX_LLC_CACHE_ALIGN)));
#else
    __declspec(align(QPHIX_LLC_CACHE_ALIGN))
        typename Geometry<AT, V, S, compress>::FourSpinorBlock x_spinor;
    __declspec(align(QPHIX_LLC_CACHE_ALIGN))
        typename Geometry<AT, V, S, compress>::FourSpinorBlock y_spinor;
    __declspec(align(QPHIX_LLC_CACHE_ALIGN))
        typename Geometry<AT, V, S, compress>::FourSpinorBlock res_spinor;
#endif

    BLASUtils::streamInSpinor<FT, V>((AT *)x_spinor, xbase, nvec_in_spinor);
    BLASUtils::streamInSpinor<FT, V>((AT *)y_spinor, ybase, nvec_in_spinor);

    // Now we are hopefully both in L1 and in the right layout so
    for (int col = 0; col < 3; col++) {
      for (int spin = 0; spin < 4; spin++) {
        for (int reim = 0; reim < 2; reim++) {
#pragma omp simd aligned(res_spinor, x_spinor, y_spinor : QPHIX_LLC_CACHE_ALIGN)
          for (int i = 0; i < S; i++) {
            res_spinor[col][spin][reim][i] =
                x_spinor[col][spin][reim][i] - y_spinor[col][spin][reim][i];
          }
        }
      }
    }

#if 1
    for (int col = 0; col < 3; col++) {
      for (int spin = 0; spin < 4; spin++) {
        for (int reim = 0; reim < 2; reim++) {
          for (int i = 0; i < S; i++) {
            red[i] += (double)res_spinor[col][spin][reim][i] *
                      (double)res_spinor[col][spin][reim][i];
          }
        }
      }
    }
#endif
    BLASUtils::streamOutSpinor<FT, V>(
        resbase, (const AT *)res_spinor, nvec_in_spinor);
  }

 private:
  typename Geometry<FT, V, S, compress>::FourSpinorBlock *res;
  const typename Geometry<FT, V, S, compress>::FourSpinorBlock *x;
  const typename Geometry<FT, V, S, compress>::FourSpinorBlock *y;
};

template <typename FT, int V, int S, bool compress>
class XMY2Norm2Functor
{
 public:
  typedef typename ArithType<FT>::Type AT;

  XMY2Norm2Functor(const typename Geometry<FT, V, S, compress>::FourSpinorBlock *x_,
                  typename Geometry<FT, V, S, compress>::FourSpinorBlock *y_)
      :  x(x_), y(y_)
  {
  }

  ~XMY2Norm2Functor() {}

  inline void func(int block, std::array<double, S> &red) const
  {
    int nvec_in_spinor = (3 * 4 * 2 * S) / V;
    const FT *xbase = &x[block][0][0][0][0];
    FT *ybase = &y[block][0][0][0][0];

// Temporary storage to stream into and out of
#if defined(__GNUG__) && !defined(__INTEL_COMPILER)
    typename Geometry<AT, V, S, compress>::FourSpinorBlock x_spinor
        __attribute__((aligned(QPHIX_LLC_CACHE_ALIGN)));
    typename Geometry<AT, V, S, compress>::FourSpinorBlock y_spinor
        __attribute__((aligned(QPHIX_LLC_CACHE_ALIGN)));

#else
    __declspec(align(QPHIX_LLC_CACHE_ALIGN))
        typename Geometry<AT, V, S, compress>::FourSpinorBlock x_spinor;
    __declspec(align(QPHIX_LLC_CACHE_ALIGN))
        typename Geometry<AT, V, S, compress>::FourSpinorBlock y_spinor;

#endif

    BLASUtils::streamInSpinor<FT, V>((AT *)x_spinor, xbase, nvec_in_spinor);
    BLASUtils::streamInSpinor<FT, V>((AT *)y_spinor, ybase, nvec_in_spinor);

    // Now we are hopefully both in L1 and in the right layout so
    for (int col = 0; col < 3; col++) {
      for (int spin = 0; spin < 4; spin++) {
        for (int reim = 0; reim < 2; reim++) {
#pragma omp simd aligned(x_spinor, y_spinor : QPHIX_LLC_CACHE_ALIGN)
          for (int i = 0; i < S; i++) {
            y_spinor[col][spin][reim][i] -=
                x_spinor[col][spin][reim][i];
          }
        }
      }
    }

#if 1
    for (int col = 0; col < 3; col++) {
      for (int spin = 0; spin < 4; spin++) {
        for (int reim = 0; reim < 2; reim++) {
#pragma omp simd safelen(S)
          for (int i = 0; i < S; i++) {
            red[i] += (double)y_spinor[col][spin][reim][i] *
                      (double)y_spinor[col][spin][reim][i];
          }
        }
      }
    }
#endif
    BLASUtils::streamOutSpinor<FT, V>(
        ybase, (const AT *)y_spinor, nvec_in_spinor);
  }

 private:

  const typename Geometry<FT, V, S, compress>::FourSpinorBlock *x;
  typename Geometry<FT, V, S, compress>::FourSpinorBlock *y;
};

template <typename FT, int V, int S, bool compress>
class AXPYNorm2Functor
{
 public:
  typedef typename ArithType<FT>::Type AT;

  AXPYNorm2Functor(double a_,
                   const typename Geometry<FT, V, S, compress>::FourSpinorBlock *x_,
                   typename Geometry<FT, V, S, compress>::FourSpinorBlock *y_)
      : a(rep<AT, double>(a_)), x(x_), y(y_)
  {
  }
  ~AXPYNorm2Functor() {}

  inline void func(int block, std::array<double, S> &reduction) const
  {
    int nvec_in_spinor = (3 * 4 * 2 * S) / V;
    const FT *xbase = &x[block][0][0][0][0];
    FT *ybase = &y[block][0][0][0][0];

// Temporary storage to stream into and out of
#if defined(__GNUG__) && !defined(__INTEL_COMPILER)
    typename Geometry<AT, V, S, compress>::FourSpinorBlock x_spinor
        __attribute__((aligned(QPHIX_LLC_CACHE_ALIGN)));
    typename Geometry<AT, V, S, compress>::FourSpinorBlock y_spinor
        __attribute__((aligned(QPHIX_LLC_CACHE_ALIGN)));
#else
    __declspec(align(QPHIX_LLC_CACHE_ALIGN))
        typename Geometry<AT, V, S, compress>::FourSpinorBlock x_spinor;
    __declspec(align(QPHIX_LLC_CACHE_ALIGN))
        typename Geometry<AT, V, S, compress>::FourSpinorBlock y_spinor;
#endif

    BLASUtils::streamInSpinor<FT, V>((AT *)x_spinor, xbase, nvec_in_spinor);
    BLASUtils::streamInSpinor<FT, V>((AT *)y_spinor, ybase, nvec_in_spinor);

    // Now we are hopefully both in L1 and in the right layout so
    for (int col = 0; col < 3; col++) {
      for (int spin = 0; spin < 4; spin++) {
        for (int reim = 0; reim < 2; reim++) {
#pragma omp simd aligned(y_spinor, x_spinor : QPHIX_LLC_CACHE_ALIGN)
          for (int i = 0; i < S; i++) {
            y_spinor[col][spin][reim][i] =
                a * x_spinor[col][spin][reim][i] + y_spinor[col][spin][reim][i];
          }
        }
      }
    }

#if 1
    for (int col = 0; col < 3; col++) {
      for (int spin = 0; spin < 4; spin++) {
        for (int reim = 0; reim < 2; reim++) {
          for (int i = 0; i < S; i++) {
            reduction[i] += rep<double, AT>(y_spinor[col][spin][reim][i]) *
                            rep<double, AT>(y_spinor[col][spin][reim][i]);
          }
        }
      }
    }
#endif

    BLASUtils::streamOutSpinor<FT, V>(ybase, (const AT *)y_spinor, nvec_in_spinor);
  }

 private:
  AT a;
  const typename Geometry<FT, V, S, compress>::FourSpinorBlock *x;
  typename Geometry<FT, V, S, compress>::FourSpinorBlock *y;
};

template <typename FT, int V, int S, bool compress>
class XMYFunctor
{
 public:
  typedef typename ArithType<FT>::Type AT;

  XMYFunctor(const typename Geometry<FT, V, S, compress>::FourSpinorBlock *x_,
             typename Geometry<FT, V, S, compress>::FourSpinorBlock *y_)
      : x(x_), y(y_)
  {
  }

  ~XMYFunctor() {}

  inline void func(int block) const
  {
    int nvec_in_spinor = (3 * 4 * 2 * S) / V;
    const FT *xbase = &x[block][0][0][0][0];
    FT *ybase = &y[block][0][0][0][0];

#if defined(__GNUG__) && !defined(__INTEL_COMPILER)
    typename Geometry<AT, V, S, compress>::FourSpinorBlock x_spinor
        __attribute__((aligned(QPHIX_LLC_CACHE_ALIGN)));
    typename Geometry<AT, V, S, compress>::FourSpinorBlock y_spinor
        __attribute__((aligned(QPHIX_LLC_CACHE_ALIGN)));
#else
    __declspec(align(QPHIX_LLC_CACHE_ALIGN))
        typename Geometry<AT, V, S, compress>::FourSpinorBlock x_spinor;
    __declspec(align(QPHIX_LLC_CACHE_ALIGN))
        typename Geometry<AT, V, S, compress>::FourSpinorBlock y_spinor;
#endif

    BLASUtils::streamInSpinor<FT, V>((AT *)x_spinor, xbase, nvec_in_spinor);
    BLASUtils::streamInSpinor<FT, V>((AT *)y_spinor, ybase, nvec_in_spinor);
    // Now we are hopefully both in L1 and in the right layout so
    for (int col = 0; col < 3; col++) {
      for (int spin = 0; spin < 4; spin++) {
        for (int reim = 0; reim < 2; reim++) {
#pragma omp simd aligned(y_spinor, x_spinor : QPHIX_LLC_CACHE_ALIGN)
          for (int i = 0; i < S; i++) {
            y_spinor[col][spin][reim][i] =
                x_spinor[col][spin][reim][i] - y_spinor[col][spin][reim][i];
          }
        }
      }
    }
    BLASUtils::streamOutSpinor<FT, V>(ybase, (const AT *)y_spinor, nvec_in_spinor);
  }

 private:
  const typename Geometry<FT, V, S, compress>::FourSpinorBlock *x;
  typename Geometry<FT, V, S, compress>::FourSpinorBlock *y;
};

template <typename FT, int V, int S, bool compress>
class RmammpNorm2rxpapFunctor
{
 public:
  typedef typename ArithType<FT>::Type AT;

  RmammpNorm2rxpapFunctor(
      double a_,
      typename Geometry<FT, V, S, compress>::FourSpinorBlock *r_,
      const typename Geometry<FT, V, S, compress>::FourSpinorBlock *mmp_,
      typename Geometry<FT, V, S, compress>::FourSpinorBlock *x_,
      const typename Geometry<FT, V, S, compress>::FourSpinorBlock *p_)
      : a(rep<AT, double>(a_)), r(r_), mmp(mmp_), x(x_), p(p_)
  {
  }

  RmammpNorm2rxpapFunctor(const RmammpNorm2rxpapFunctor<FT, V, S, compress> &rhs)
      : a(rhs.a), r(rhs.r), mmp(rhs.mmp), x(rhs.x), p(rhs.p)
  {
  }

  ~RmammpNorm2rxpapFunctor() {}

  inline void func(int block, std::array<double, S> &reduction) const
  {
    int nvec_in_spinor = (3 * 4 * 2 * S) / V;
    FT *rbase = &r[block][0][0][0][0];
    const FT *mmpbase = &mmp[block][0][0][0][0];
    FT *xbase = &x[block][0][0][0][0];
    const FT *pbase = &p[block][0][0][0][0];

#if defined(__GNUG__) && !defined(__INTEL_COMPILER)
    typename Geometry<AT, V, S, compress>::FourSpinorBlock r_spinor
        __attribute__((aligned(QPHIX_LLC_CACHE_ALIGN)));
    typename Geometry<AT, V, S, compress>::FourSpinorBlock mmp_spinor
        __attribute__((aligned(QPHIX_LLC_CACHE_ALIGN)));
    typename Geometry<AT, V, S, compress>::FourSpinorBlock x_spinor
        __attribute__((aligned(QPHIX_LLC_CACHE_ALIGN)));
    typename Geometry<AT, V, S, compress>::FourSpinorBlock p_spinor
        __attribute__((aligned(QPHIX_LLC_CACHE_ALIGN)));
#else
    __declspec(align(QPHIX_LLC_CACHE_ALIGN))
        typename Geometry<AT, V, S, compress>::FourSpinorBlock r_spinor;
    __declspec(align(QPHIX_LLC_CACHE_ALIGN))
        typename Geometry<AT, V, S, compress>::FourSpinorBlock mmp_spinor;
    __declspec(align(QPHIX_LLC_CACHE_ALIGN))
        typename Geometry<AT, V, S, compress>::FourSpinorBlock x_spinor;
    __declspec(align(QPHIX_LLC_CACHE_ALIGN))
        typename Geometry<AT, V, S, compress>::FourSpinorBlock p_spinor;
#endif

    BLASUtils::streamInSpinor<FT, V>((AT *)r_spinor, rbase, nvec_in_spinor);
    BLASUtils::streamInSpinor<FT, V>((AT *)mmp_spinor, mmpbase, nvec_in_spinor);
    BLASUtils::streamInSpinor<FT, V>((AT *)x_spinor, xbase, nvec_in_spinor);
    BLASUtils::streamInSpinor<FT, V>((AT *)p_spinor, pbase, nvec_in_spinor);

    // Now we are hopefully both in L1 and in the right layout so
    for (int col = 0; col < 3; col++) {
      for (int spin = 0; spin < 4; spin++) {
        for (int reim = 0; reim < 2; reim++) {
#pragma omp simd aligned(                                                           \
    r_spinor, mmp_spinor, x_spinor, p_spinor : QPHIX_LLC_CACHE_ALIGN)
          for (int i = 0; i < S; i++) {
            r_spinor[col][spin][reim][i] =
                r_spinor[col][spin][reim][i] - a * mmp_spinor[col][spin][reim][i];
            x_spinor[col][spin][reim][i] =
                x_spinor[col][spin][reim][i] + a * p_spinor[col][spin][reim][i];
          }
        }
      }
    }

    for (int col = 0; col < 3; col++) {
      for (int spin = 0; spin < 4; spin++) {
        for (int reim = 0; reim < 2; reim++) {
          for (int i = 0; i < S; i++) {
            reduction[i] += (double)r_spinor[col][spin][reim][i] *
                            (double)r_spinor[col][spin][reim][i];
          }
        }
      }
    }

    BLASUtils::streamOutSpinor<FT, V>(rbase, (const AT *)r_spinor, nvec_in_spinor);
    BLASUtils::streamOutSpinor<FT, V>(xbase, (const AT *)x_spinor, nvec_in_spinor);
  }

 private:
  AT a;
  typename Geometry<FT, V, S, compress>::FourSpinorBlock *r;
  const typename Geometry<FT, V, S, compress>::FourSpinorBlock *mmp;
  typename Geometry<FT, V, S, compress>::FourSpinorBlock *x;
  const typename Geometry<FT, V, S, compress>::FourSpinorBlock *p;
};

template <typename FT, int V, int S, bool compress>
class RichardsonRXUpdateNormRFunctor
{
 public:
  typedef typename ArithType<FT>::Type AT;

  RichardsonRXUpdateNormRFunctor(
      typename Geometry<FT, V, S, compress>::FourSpinorBlock *x_,
      typename Geometry<FT, V, S, compress>::FourSpinorBlock *r_,

      const typename Geometry<FT, V, S, compress>::FourSpinorBlock *delta_x_,
      const typename Geometry<FT, V, S, compress>::FourSpinorBlock *delta_r_)
      : x(x_), r(r_), delta_x(delta_x_), delta_r(delta_r_)
  {
  }

  ~RichardsonRXUpdateNormRFunctor() {}

  inline void func(int block, std::array<double, S> &reduction) const
  {
    FT *xbase = &x[block][0][0][0][0];
    FT *rbase = &r[block][0][0][0][0];
    const FT *delta_xbase = &delta_x[block][0][0][0][0];
    const FT *delta_rbase = &delta_r[block][0][0][0][0];
    int nvec_in_spinor = (3 * 4 * 2 * S) / V;

#if defined(__GNUG__) && !defined(__INTEL_COMPILER)
    typename Geometry<AT, V, S, compress>::FourSpinorBlock x_spinor
        __attribute__((aligned(QPHIX_LLC_CACHE_ALIGN)));
    typename Geometry<AT, V, S, compress>::FourSpinorBlock r_spinor
        __attribute__((aligned(QPHIX_LLC_CACHE_ALIGN)));
    typename Geometry<AT, V, S, compress>::FourSpinorBlock delta_x_spinor
        __attribute__((aligned(QPHIX_LLC_CACHE_ALIGN)));
    typename Geometry<AT, V, S, compress>::FourSpinorBlock delta_r_spinor
        __attribute__((aligned(QPHIX_LLC_CACHE_ALIGN)));
#else
    __declspec(align(QPHIX_LLC_CACHE_ALIGN))
        typename Geometry<AT, V, S, compress>::FourSpinorBlock x_spinor;
    __declspec(align(QPHIX_LLC_CACHE_ALIGN))
        typename Geometry<AT, V, S, compress>::FourSpinorBlock r_spinor;
    __declspec(align(QPHIX_LLC_CACHE_ALIGN))
        typename Geometry<AT, V, S, compress>::FourSpinorBlock delta_x_spinor;
    __declspec(align(QPHIX_LLC_CACHE_ALIGN))
        typename Geometry<AT, V, S, compress>::FourSpinorBlock delta_r_spinor;
#endif

    BLASUtils::streamInSpinor<FT, V>((AT *)x_spinor, xbase, nvec_in_spinor);
    BLASUtils::streamInSpinor<FT, V>((AT *)r_spinor, rbase, nvec_in_spinor);
    BLASUtils::streamInSpinor<FT, V>(
        (AT *)delta_x_spinor, delta_xbase, nvec_in_spinor);
    BLASUtils::streamInSpinor<FT, V>(
        (AT *)delta_r_spinor, delta_rbase, nvec_in_spinor);

    // Now we are hopefully both in L1 and in the right layout so
    for (int col = 0; col < 3; col++) {
      for (int spin = 0; spin < 4; spin++) {
        for (int reim = 0; reim < 2; reim++) {
#pragma omp simd aligned(                                                           \
    x_spinor, r_spinor, delta_x_spinor, delta_r_spinor : QPHIX_LLC_CACHE_ALIGN)
          for (int i = 0; i < S; i++) {
            x_spinor[col][spin][reim][i] += delta_x_spinor[col][spin][reim][i];
            r_spinor[col][spin][reim][i] -= delta_r_spinor[col][spin][reim][i];
          }
        }
      }
    }

    for (int col = 0; col < 3; col++) {
      for (int spin = 0; spin < 4; spin++) {
        for (int reim = 0; reim < 2; reim++) {
          for (int i = 0; i < S; i++) {
            reduction[i] += (double)r_spinor[col][spin][reim][i] *
                            (double)r_spinor[col][spin][reim][i];
          }
        }
      }
    }

    BLASUtils::writeSpinor<FT, V>(xbase, (const AT *)x_spinor, nvec_in_spinor);
    BLASUtils::writeSpinor<FT, V>(rbase, (const AT *)r_spinor, nvec_in_spinor);
  }

 private:
  typename Geometry<FT, V, S, compress>::FourSpinorBlock *x;
  typename Geometry<FT, V, S, compress>::FourSpinorBlock *r;
  const typename Geometry<FT, V, S, compress>::FourSpinorBlock *delta_x;
  const typename Geometry<FT, V, S, compress>::FourSpinorBlock *delta_r;
};

}; // Namespace
