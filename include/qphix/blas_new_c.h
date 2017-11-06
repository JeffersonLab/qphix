#pragma once

#include "qphix/geometry.h"
#include "qphix/comm.h"
#include "qphix/qphix_config.h"

#include "qphix/site_loops.h"
#include "qphix/real_functors.h"
#include "qphix/complex_functors.h"
#include "qphix/arith_type.h"

#include <omp.h>

namespace QPhiX
{

template <typename FT, int V, int S, bool compress>
void copySpinor(typename Geometry<FT, V, S, compress>::FourSpinorBlock *res,
                const typename Geometry<FT, V, S, compress>::FourSpinorBlock *src,
                const Geometry<FT, V, S, compress> &geom,
                int n_blas_simt)
{
  CopyFunctor<FT, V, S, compress> f(res, src);
  siteLoopNoReduction<FT, V, S, compress, CopyFunctor<FT, V, S, compress>>(
      f, geom, n_blas_simt);
}

/**
  \see The article \ref intel-cpp-compiler-workaround contains an explanation
  of the `typename Spinor1`, `enable_if`, and `is_same` constructs.
  */
template <typename FT, int V, int S, bool compress, int num_flav, typename Spinor1>
typename std::enable_if<
    std::is_same<const typename Geometry<FT, V, S, compress>::FourSpinorBlock,
                 const Spinor1>::value,
    void>::type
copySpinor(
    typename Geometry<FT, V, S, compress>::FourSpinorBlock *const *res,
    Spinor1 *const *src,
    const Geometry<FT, V, S, compress> &geom,
    int n_blas_simt)
{
  for (uint8_t f = 0; f < num_flav; ++f) {
    copySpinor(res[f], src[f], geom, n_blas_simt);
  }
}

template <typename FT, int V, int S, bool compress>
void zeroSpinor(typename Geometry<FT, V, S, compress>::FourSpinorBlock *res,
                const Geometry<FT, V, S, compress> &geom,
                int n_blas_simt)
{
  ZeroFunctor<FT, V, S, compress> f(res);
  siteLoopNoReduction<FT, V, S, compress, ZeroFunctor<FT, V, S, compress>>(
      f, geom, n_blas_simt);
}

template <typename FT, int V, int S, bool compress, int num_flav>
void zeroSpinor(
    typename Geometry<FT, V, S, compress>::FourSpinorBlock *const res[num_flav],
    const Geometry<FT, V, S, compress> &geom,
    int n_blas_simt)
{
  for (uint8_t f = 0; f < num_flav; ++f) {
    zeroSpinor(res[f], geom, n_blas_simt);
  }
}

template <typename FT, int V, int S, bool compress>
void axy(const double alpha,
         const typename Geometry<FT, V, S, compress>::FourSpinorBlock *x,
         typename Geometry<FT, V, S, compress>::FourSpinorBlock *y,
         const Geometry<FT, V, S, compress> &geom,
         int n_blas_simt)
{
  AXYFunctor<FT, V, S, compress> f(alpha, x, y);
  siteLoopNoReduction<FT, V, S, compress, AXYFunctor<FT, V, S, compress>>(
      f, geom, n_blas_simt);
}

/**
  \see The article \ref intel-cpp-compiler-workaround contains an explanation
  of the `typename Spinor1`, `enable_if`, and `is_same` constructs.
  */
template <typename FT, int V, int S, bool compress, int num_flav, typename Spinor1>
typename std::enable_if<
    std::is_same<const typename Geometry<FT, V, S, compress>::FourSpinorBlock,
                 const Spinor1>::value,
    void>::type
axy(const double alpha,
    Spinor1 *const x[num_flav],
    typename Geometry<FT, V, S, compress>::FourSpinorBlock *const y[num_flav],
    const Geometry<FT, V, S, compress> &geom,
    int n_blas_simt)
{
  for (int f = 0; f < num_flav; ++f) {
    axy(alpha, x[f], y[f], geom, n_blas_simt);
  }
}

template <typename FT, int V, int S, bool compress>
void ax(const double alpha,
    typename Geometry<FT, V, S, compress>::FourSpinorBlock *x,
    const Geometry<FT, V, S, compress> &geom,
    int n_blas_simt)
{

  AXFunctor<FT, V, S, compress> f(alpha, x);
  siteLoopNoReduction<FT, V, S, compress, AXFunctor<FT, V, S, compress>>(
      f, geom, n_blas_simt);
}

template <typename FT, int V, int S, bool compress>
void ypeqx(typename Geometry<FT, V, S, compress>::FourSpinorBlock *x,
    typename Geometry<FT, V, S, compress>::FourSpinorBlock *y,
    const Geometry<FT, V, S, compress> &geom,
    int n_blas_simt)
{

  YPEQXFunctor<FT, V, S, compress> f(x,y);
  siteLoopNoReduction<FT, V, S, compress, YPEQXFunctor<FT, V, S, compress>>(
      f, geom, n_blas_simt);
}

template <typename FT, int V, int S, bool compress>
void ymeqx(typename Geometry<FT, V, S, compress>::FourSpinorBlock *x,
    typename Geometry<FT, V, S, compress>::FourSpinorBlock *y,
    const Geometry<FT, V, S, compress> &geom,
    int n_blas_simt)
{

  YMEQXFunctor<FT, V, S, compress> f(x,y);
  siteLoopNoReduction<FT, V, S, compress, YMEQXFunctor<FT, V, S, compress>>(
      f, geom, n_blas_simt);
}


template <typename FT, int V, int S, bool compress>
void aypx(const double alpha,
          const typename Geometry<FT, V, S, compress>::FourSpinorBlock *x,
          typename Geometry<FT, V, S, compress>::FourSpinorBlock *y,
          const Geometry<FT, V, S, compress> &geom,
          int n_blas_simt)
{
  AYPXFunctor<FT, V, S, compress> f(alpha, x, y);
  siteLoopNoReduction<FT, V, S, compress, AYPXFunctor<FT, V, S, compress>>(
      f, geom, n_blas_simt);
}

/**
  Multiplies the spinor with \f$ \alpha + \mathrm i \mu \gamma_5 \f$.

  \param[in] apimu Array of length 2. Real and imaginary part of \f$ \alpha +
  \mathrm i \mu \f$, the first array element is α, the second is μ.
  */
template <typename FT, int V, int S, bool compress>
void twisted_mass(const double apimu[2],
                  const typename Geometry<FT, V, S, compress>::FourSpinorBlock *x,
                  typename Geometry<FT, V, S, compress>::FourSpinorBlock *y,
                  const Geometry<FT, V, S, compress> &geom,
                  int n_blas_simt)
{
  TwistedMassFunctor<FT, V, S, compress> f(apimu, x, y);
  siteLoopNoReduction<FT, V, S, compress, TwistedMassFunctor<FT, V, S, compress>>(
     f, geom, n_blas_simt);
}

/**
  Multiplies the spinor with (α + iμγ⁵τ³ + ετ¹).
  */
template <typename FT, int V, int S, bool compress, typename Spinor1>
typename std::enable_if<
    std::is_same<const typename Geometry<FT, V, S, compress>::FourSpinorBlock,
                 const Spinor1>::value,
    void>::type
two_flav_twisted_mass(
    const double apimu[2],
    const double epsilon,
    Spinor1 *const *x,
    typename Geometry<FT, V, S, compress>::FourSpinorBlock *const *y,
    const Geometry<FT, V, S, compress> &geom,
    int n_blas_simt)
{
  TwoFlavTwistedMassFunctor<FT, V, S, compress> f(apimu, epsilon, x, y);
  siteLoopNoReduction<FT,
                      V,
                      S,
                      compress,
                      TwoFlavTwistedMassFunctor<FT, V, S, compress>>(
      f, geom, n_blas_simt);
}

/**
  \see The article \ref intel-cpp-compiler-workaround contains an explanation
  of the `typename Spinor1`, `enable_if`, and `is_same` constructs.
  */
template <typename FT, int V, int S, bool compress, int num_flav, typename Spinor1>
typename std::enable_if<
    std::is_same<const typename Geometry<FT, V, S, compress>::FourSpinorBlock,
                 const Spinor1>::value,
    void>::type
aypx(const double alpha,
     Spinor1 *const x[num_flav],
     typename Geometry<FT, V, S, compress>::FourSpinorBlock *const y[num_flav],
     const Geometry<FT, V, S, compress> &geom,
     int n_blas_simt)
{
  for (int f = 0; f < num_flav; ++f) {
    aypx(alpha, x[f], y[f], geom, n_blas_simt);
  }
}

template <typename FT, int V, int S, bool compress>
void axpy(const double alpha,
          const typename Geometry<FT, V, S, compress>::FourSpinorBlock *x,
          typename Geometry<FT, V, S, compress>::FourSpinorBlock *y,
          const Geometry<FT, V, S, compress> &geom,
          int n_blas_simt)
{

  AXPYFunctor<FT, V, S, compress> f(alpha, x, y);
  siteLoopNoReduction<FT, V, S, compress, AXPYFunctor<FT, V, S, compress>>(
      f, geom, n_blas_simt);
}

/**
  \see The article \ref intel-cpp-compiler-workaround contains an explanation
  of the `typename Spinor1`, `enable_if`, and `is_same` constructs.
  */
template <typename FT, int V, int S, bool compress, int num_flav, typename Spinor1>
typename std::enable_if<
    std::is_same<const typename Geometry<FT, V, S, compress>::FourSpinorBlock,
                 const Spinor1>::value,
    void>::type
axpy(const double alpha,
     Spinor1 *const x[num_flav],
     typename Geometry<FT, V, S, compress>::FourSpinorBlock *y[num_flav],
     const Geometry<FT, V, S, compress> &geom,
     int n_blas_simt)
{
  for (int f = 0; f < num_flav; ++f) {
    axpy(alpha, x[f], y[f], geom, n_blas_simt);
  }
}

template <typename FT, int V, int S, bool compress>
void axpby(const double alpha,
           const typename Geometry<FT, V, S, compress>::FourSpinorBlock *x,
           const double beta,
           typename Geometry<FT, V, S, compress>::FourSpinorBlock *y,
           const Geometry<FT, V, S, compress> &geom,
           int n_blas_simt)
{

  AXPBYFunctor<FT, V, S, compress> f(alpha, x, beta, y);
  siteLoopNoReduction<FT, V, S, compress, AXPBYFunctor<FT, V, S, compress>>(
      f, geom, n_blas_simt);
}

template <typename FT, int V, int S, bool compress, int num_flav, typename Spinor1>
typename std::enable_if<
    std::is_same<const typename Geometry<FT, V, S, compress>::FourSpinorBlock,
                 const Spinor1>::value,
    void>::type
axpby(const double alpha,
      Spinor1 *const x[num_flav],
      const double beta,
      typename Geometry<FT, V, S, compress>::FourSpinorBlock *const y[num_flav],
      const Geometry<FT, V, S, compress> &geom,
      int n_blas_simt)
{
  for (int f = 0; f < num_flav; ++f) {
    axpby(alpha, x[f], beta, y[f], geom, n_blas_simt);
  }
}

template <typename FT, int V, int S, bool compress>
void norm2Spinor(double &n2,
                 const typename Geometry<FT, V, S, compress>::FourSpinorBlock *x,
                 const Geometry<FT, V, S, compress> &geom,
                 int n_blas_simt,
                 bool globalSum=true)
{
  Norm2Functor<FT, V, S, compress> f(x);
  siteLoop1Reduction<FT, V, S, compress, Norm2Functor<FT, V, S, compress>>(
      f, n2, geom, n_blas_simt,globalSum);
}

/**
  \see The article \ref intel-cpp-compiler-workaround contains an explanation
  of the `typename Spinor1`, `enable_if`, and `is_same` constructs.
  */
template <typename FT, int V, int S, bool compress, int num_flav, typename Spinor1>
typename std::enable_if<
    std::is_same<const typename Geometry<FT, V, S, compress>::FourSpinorBlock,
                 const Spinor1>::value,
    void>::type
norm2Spinor(double &n2,
            Spinor1 *const x[num_flav],
            Geometry<FT, V, S, compress> &geom,
            int n_blas_simt)
{
  n2 = 0.0;
  for (uint8_t f = 0; f < num_flav; ++f) {
    double local_n2;
    norm2Spinor(local_n2, x[f], geom, n_blas_simt);
    n2 += local_n2;
  }
}

template <typename FT, int V, int S, bool compress>
void axpyNorm2(const double alpha,
               const typename Geometry<FT, V, S, compress>::FourSpinorBlock *x,
               typename Geometry<FT, V, S, compress>::FourSpinorBlock *y,
               double &norm2y,
               const Geometry<FT, V, S, compress> &geom,
               int n_blas_simt)
{
  AXPYNorm2Functor<FT, V, S, compress> f(alpha, x, y);
  siteLoop1Reduction<FT, V, S, compress, AXPYNorm2Functor<FT, V, S, compress>>(
      f, norm2y, geom, n_blas_simt);
} // End of Function.

template <typename FT, int V, int S, bool compress, int num_flav, typename Spinor1>
typename std::enable_if<
    std::is_same<const typename Geometry<FT, V, S, compress>::FourSpinorBlock,
                 const Spinor1>::value,
    void>::type
axpyNorm2(const double alpha,
          Spinor1 *const x[num_flav],
          typename Geometry<FT, V, S, compress>::FourSpinorBlock *const y[num_flav],
          double &norm2y,
          const Geometry<FT, V, S, compress> &geom,
          int n_blas_simt)
{
  norm2y = 0.0;
  for (uint8_t f = 0; f < num_flav; ++f) {
    double local_norm2y;
    axpyNorm2(alpha, x[f], y[f], local_norm2y, geom, n_blas_simt);
    norm2y += local_norm2y;
  }
}

template <typename FT, int V, int S, bool compress>
void xmyNorm2Spinor(typename Geometry<FT, V, S, compress>::FourSpinorBlock *res,
                    const typename Geometry<FT, V, S, compress>::FourSpinorBlock *x,
                    typename Geometry<FT, V, S, compress>::FourSpinorBlock *y,
                    double &n2res,
                    const Geometry<FT, V, S, compress> &geom,
                    int n_blas_simt,
                    bool globalSum=true)
{
  XMYNorm2Functor<FT, V, S, compress> f(res, x, y);
  siteLoop1Reduction<FT, V, S, compress, XMYNorm2Functor<FT, V, S, compress>>(
      f, n2res, geom, n_blas_simt, globalSum);
} // End of Function.


template <typename FT, int V, int S, bool compress>
void xmy2Norm2Spinor(const typename Geometry<FT, V, S, compress>::FourSpinorBlock *x,
                      typename Geometry<FT, V, S, compress>::FourSpinorBlock *y,
                      double &n2res,
                      const Geometry<FT, V, S, compress> &geom,
                      int n_blas_simt,
                      bool globalSum=true)
{
  XMY2Norm2Functor<FT, V, S, compress> f(x, y);
  siteLoop1Reduction<FT, V, S, compress, XMY2Norm2Functor<FT, V, S, compress>>(
      f, n2res, geom, n_blas_simt, globalSum);
} // End of Function.
/**
  \see The article \ref intel-cpp-compiler-workaround contains an explanation
  of the `typename Spinor1`, `enable_if`, and `is_same` constructs.
  */
template <typename FT, int V, int S, bool compress, int num_flav, typename Spinor1>
typename std::enable_if<
    std::is_same<const typename Geometry<FT, V, S, compress>::FourSpinorBlock,
                 const Spinor1>::value,
    void>::type
xmyNorm2Spinor(
    typename Geometry<FT, V, S, compress>::FourSpinorBlock *const res[num_flav],
    Spinor1 *const x[num_flav],
    typename Geometry<FT, V, S, compress>::FourSpinorBlock *const y[num_flav],
    double &n2res,
    const Geometry<FT, V, S, compress> &geom,
    int n_blas_simt)
{
  n2res = 0;
  for (uint8_t f = 0; f < num_flav; ++f) {
    double local_n2res;
    xmyNorm2Spinor(res[f], x[f], y[f], local_n2res, geom, n_blas_simt);
    n2res += local_n2res;
  }
}

template <typename FT, int V, int S, bool compress>
void rmammpNorm2rxpap(typename Geometry<FT, V, S, compress>::FourSpinorBlock *r,
                      const double &ar,
                      typename Geometry<FT, V, S, compress>::FourSpinorBlock *mmp,
                      double &cp,
                      typename Geometry<FT, V, S, compress>::FourSpinorBlock *x,
                      typename Geometry<FT, V, S, compress>::FourSpinorBlock *p,
                      const Geometry<FT, V, S, compress> &geom,
                      int n_blas_simt)
{
  RmammpNorm2rxpapFunctor<FT, V, S, compress> f(ar, r, mmp, x, p);
  siteLoop1Reduction<FT,
                     V,
                     S,
                     compress,
                     RmammpNorm2rxpapFunctor<FT, V, S, compress>>(
      f, cp, geom, n_blas_simt);
} // End of Function.

template <typename FT, int V, int S, bool compress, int num_flav>
void rmammpNorm2rxpap(
    typename Geometry<FT, V, S, compress>::FourSpinorBlock *const r[num_flav],
    const double &ar,
    typename Geometry<FT, V, S, compress>::FourSpinorBlock *const mmp[num_flav],
    double &cp,
    typename Geometry<FT, V, S, compress>::FourSpinorBlock *const x[num_flav],
    typename Geometry<FT, V, S, compress>::FourSpinorBlock *const p[num_flav],
    const Geometry<FT, V, S, compress> &geom,
    int n_blas_simt)
{
  cp = 0;
  for (uint8_t f = 0; f < num_flav; ++f) {
    double local_cp;
    rmammpNorm2rxpap(r[f], ar, mmp[f], local_cp, x[f], p[f], geom, n_blas_simt);
    cp += local_cp;
  }
}

template <typename FT, int V, int S, bool compress>
void richardson_rxupdateNormR(
    typename Geometry<FT, V, S, compress>::FourSpinorBlock *x,
    typename Geometry<FT, V, S, compress>::FourSpinorBlock *r,
    const typename Geometry<FT, V, S, compress>::FourSpinorBlock *delta_x,
    const typename Geometry<FT, V, S, compress>::FourSpinorBlock *delta_r,
    double &cp,
    const Geometry<FT, V, S, compress> &geom,
    int n_blas_simt)
{
  RichardsonRXUpdateNormRFunctor<FT, V, S, compress> f(x, r, delta_x, delta_r);
  siteLoop1Reduction<FT,
                     V,
                     S,
                     compress,
                     RichardsonRXUpdateNormRFunctor<FT, V, S, compress>>(
      f, cp, geom, n_blas_simt);
} // End of Function.

/**
  \see The article \ref intel-cpp-compiler-workaround contains an explanation
  of the `typename Spinor1`, `enable_if`, and `is_same` constructs.
  */
template <typename FT,
          int V,
          int S,
          bool compress,
          int num_flav,
          typename Spinor1,
          typename Spinor2>
typename std::enable_if<
    std::is_same<const typename Geometry<FT, V, S, compress>::FourSpinorBlock,
                 const Spinor1>::value &&
        std::is_same<const typename Geometry<FT, V, S, compress>::FourSpinorBlock,
                     const Spinor2>::value,
    void>::type
richardson_rxupdateNormR(
    typename Geometry<FT, V, S, compress>::FourSpinorBlock *const x[num_flav],
    typename Geometry<FT, V, S, compress>::FourSpinorBlock *const r[num_flav],
    Spinor1 *const delta_x[num_flav],
    Spinor2 *const delta_r[num_flav],
    double &cp,
    const Geometry<FT, V, S, compress> &geom,
    int n_blas_simt)
{
  cp = 0;
  for (uint8_t f = 0; f < num_flav; ++f) {
    double local_cp;
    richardson_rxupdateNormR(
        x[f], r[f], delta_x[f], delta_r[f], local_cp, geom, n_blas_simt);
    cp += local_cp;
  }
}

template <typename FT, int V, int S, bool compress>
void bicgstab_xmy(const typename Geometry<FT, V, S, compress>::FourSpinorBlock *x,
                  typename Geometry<FT, V, S, compress>::FourSpinorBlock *y,
                  const Geometry<FT, V, S, compress> &geom,
                  int n_blas_simt)
{
  XMYFunctor<FT, V, S, compress> f(x, y);
  siteLoopNoReduction<FT, V, S, compress, XMYFunctor<FT, V, S, compress>>(
      f, geom, n_blas_simt);
} // End of Function.

/**
  \see The article \ref intel-cpp-compiler-workaround contains an explanation
  of the `typename Spinor1`, `enable_if`, and `is_same` constructs.
  */
template <typename FT, int V, int S, bool compress, int num_flav, typename Spinor1>
typename std::enable_if<
    std::is_same<const typename Geometry<FT, V, S, compress>::FourSpinorBlock,
                 const Spinor1>::value,
    void>::type
bicgstab_xmy(
    Spinor1 *const x[num_flav],
    typename Geometry<FT, V, S, compress>::FourSpinorBlock *const y[num_flav],
    const Geometry<FT, V, S, compress> &geom,
    int n_blas_simt)
{
  for (uint8_t f = 0; f < num_flav; ++f) {
    bicgstab_xmy(x[f], y[f], geom, n_blas_simt);
  }
}

template <typename FT, int V, int S, bool compress>
void innerProduct(double results[2],
                  const typename Geometry<FT, V, S, compress>::FourSpinorBlock *x,
                  const typename Geometry<FT, V, S, compress>::FourSpinorBlock *y,
                  const Geometry<FT, V, S, compress> &geom,
                  int n_blas_simt,
                  bool globalSum=true)
{
  results[0] = (double) 0;
  results[1] = (double) 0;

  InnerProductFunctor<FT, V, S, compress> f(x, y);
  siteLoop2Reductions<FT, V, S, compress, InnerProductFunctor<FT, V, S, compress>>(
      f, results, geom, n_blas_simt, globalSum);
}



/**
  \see The article \ref intel-cpp-compiler-workaround contains an explanation
  of the `typename Spinor1`, `enable_if`, and `is_same` constructs.
  */
template <typename FT,
          int V,
          int S,
          bool compress,
          int num_flav,
          typename Spinor1,
          typename Spinor2>
typename std::enable_if<
    std::is_same<const typename Geometry<FT, V, S, compress>::FourSpinorBlock,
                 const Spinor1>::value &&
        std::is_same<const typename Geometry<FT, V, S, compress>::FourSpinorBlock,
                     const Spinor2>::value,
    void>::type
innerProduct(double results[2],
             Spinor1 *const x[num_flav],
             Spinor2 *const y[num_flav],
             const Geometry<FT, V, S, compress> &geom,
             int n_blas_simt)
{
  results[0] = 0;
  results[1] = 0;

  for (uint8_t f = 0; f < num_flav; ++f) {
    double local_results[2];
    innerProduct(local_results, x[f], y[f], geom, n_blas_simt);

    results[0] += local_results[0];
    results[1] += local_results[1];
  }
}

template <typename FT, int V, int S, bool compress>
void innerProductNorm(double results[3],
                  const typename Geometry<FT, V, S, compress>::FourSpinorBlock *x,
                  const typename Geometry<FT, V, S, compress>::FourSpinorBlock *y,
                  const Geometry<FT, V, S, compress> &geom,
                  int n_blas_simt)
{
  results[0] = (double) 0;
  results[1] = (double) 0;
  results[2] = (double) 0;

  InnerProductNormFunctor<FT, V, S, compress> f(x, y);
  siteLoop3Reductions<FT, V, S, compress, InnerProductNormFunctor<FT, V, S, compress>>(
      f, results, geom, n_blas_simt);
}

/**
  \see The article \ref intel-cpp-compiler-workaround contains an explanation
  of the `typename Spinor1`, `enable_if`, and `is_same` constructs.
  */
template <typename FT,
          int V,
          int S,
          bool compress,
          int num_flav,
          typename Spinor1,
          typename Spinor2>
typename std::enable_if<
    std::is_same<const typename Geometry<FT, V, S, compress>::FourSpinorBlock,
                 const Spinor1>::value &&
        std::is_same<const typename Geometry<FT, V, S, compress>::FourSpinorBlock,
                     const Spinor2>::value,
    void>::type
innerProductNorm(double results[3],
             Spinor1 *const x[num_flav],
             Spinor2 *const y[num_flav],
             const Geometry<FT, V, S, compress> &geom,
             int n_blas_simt)
{
  results[0] = 0;
  results[1] = 0;
  results[2] = 0;

  for (uint8_t f = 0; f < num_flav; ++f) {
    double local_results[3];
    innerProductNorm(local_results, x[f], y[f], geom, n_blas_simt);

    results[0] += local_results[0];
    results[1] += local_results[1];
    results[2] += local_results[2];
  }
}


template <typename FT, int V, int S, bool compress>
void bicgstab_p_update(
    const typename Geometry<FT, V, S, compress>::FourSpinorBlock *r,
    typename Geometry<FT, V, S, compress>::FourSpinorBlock *p,
    const typename Geometry<FT, V, S, compress>::FourSpinorBlock *v,
    double beta[2],
    double omega[2],
    const Geometry<FT, V, S, compress> &geom,
    int n_blas_simt)
{
  BiCGStabPUpdateFunctor<FT, V, S, compress> f(r, p, v, beta, omega);
  siteLoopNoReduction<FT,
                      V,
                      S,
                      compress,
                      BiCGStabPUpdateFunctor<FT, V, S, compress>>(
      f, geom, n_blas_simt);
}

/**
  \see The article \ref intel-cpp-compiler-workaround contains an explanation
  of the `typename Spinor1`, `enable_if`, and `is_same` constructs.
  */
template <typename FT,
          int V,
          int S,
          bool compress,
          int num_flav,
          typename Spinor1,
          typename Spinor2>
typename std::enable_if<
    std::is_same<const typename Geometry<FT, V, S, compress>::FourSpinorBlock,
                 const Spinor1>::value &&
        std::is_same<const typename Geometry<FT, V, S, compress>::FourSpinorBlock,
                     const Spinor2>::value,
    void>::type
bicgstab_p_update(
    Spinor1 *const r[num_flav],
    typename Geometry<FT, V, S, compress>::FourSpinorBlock *const p[num_flav],
    Spinor2 *const v[num_flav],
    double beta[2],
    double omega[2],
    const Geometry<FT, V, S, compress> &geom,
    int n_blas_simt)
{
  // TODO (Martin Ueding): Check whether this is the correct generalization to
  // multiple flavors.
  for (uint8_t f = 0; f < num_flav; ++f) {
    bicgstab_p_update(r[f], p[f], v[f], beta, omega, geom, n_blas_simt);
  }
}

template <typename FT, int V, int S, bool compress>
void bicgstab_s_update(
    double alpha[2],
    typename Geometry<FT, V, S, compress>::FourSpinorBlock *s,
    const typename Geometry<FT, V, S, compress>::FourSpinorBlock *v,
    const Geometry<FT, V, S, compress> &geom,
    int n_blas_simt)
{
  BiCGStabSUpdateFunctor<FT, V, S, compress> f(alpha, s, v);
  siteLoopNoReduction<FT,
                      V,
                      S,
                      compress,
                      BiCGStabSUpdateFunctor<FT, V, S, compress>>(
      f, geom, n_blas_simt);
}

template <typename FT, int V, int S, bool compress>
void caxpySpinor(
    double alpha[2],
    const typename Geometry<FT, V, S, compress>::FourSpinorBlock *x,
    typename Geometry<FT, V, S, compress>::FourSpinorBlock *y,
    const Geometry<FT, V, S, compress> &geom,
    int n_blas_simt)
{
  CAXPYFunctor<FT, V, S, compress> f(alpha, x, y);
  siteLoopNoReduction<FT,
                      V,
                      S,
                      compress,
                      CAXPYFunctor<FT, V, S, compress>>(
      f, geom, n_blas_simt);
}

/**
  \see The article \ref intel-cpp-compiler-workaround contains an explanation
  of the `typename Spinor1`, `enable_if`, and `is_same` constructs.
  */
template <typename FT, int V, int S, bool compress, int num_flav, typename Spinor1>
typename std::enable_if<
    std::is_same<const typename Geometry<FT, V, S, compress>::FourSpinorBlock,
                 const Spinor1>::value,
    void>::type
bicgstab_s_update(
    double alpha[2],
    typename Geometry<FT, V, S, compress>::FourSpinorBlock *const s[num_flav],
    Spinor1 *const v[num_flav],
    const Geometry<FT, V, S, compress> &geom,
    int n_blas_simt)
{
  // TODO (Martin Ueding): Check whether this is the correct generalization to
  // multiple flavors.
  for (uint8_t f = 0; f < num_flav; ++f) {
    bicgstab_s_update(alpha, s[f], v[f], geom, n_blas_simt);
  }
}

template <typename FT, int V, int S, bool compress>
void bicgstab_rxupdate(
    typename Geometry<FT, V, S, compress>::FourSpinorBlock *x,
    typename Geometry<FT, V, S, compress>::FourSpinorBlock *r,
    const typename Geometry<FT, V, S, compress>::FourSpinorBlock *t,
    const typename Geometry<FT, V, S, compress>::FourSpinorBlock *p,
    double omega[2],
    double alpha[2],
    double &r_norm,
    const Geometry<FT, V, S, compress> &geom,
    int n_blas_simt)
{
  BiCGStabRXUpdateFunctor<FT, V, S, compress> f(x, r, t, p, omega, alpha);
  siteLoop1Reduction<FT,
                     V,
                     S,
                     compress,
                     BiCGStabRXUpdateFunctor<FT, V, S, compress>>(
      f, r_norm, geom, n_blas_simt);
}

/**
  \see The article \ref intel-cpp-compiler-workaround contains an explanation
  of the `typename Spinor1`, `enable_if`, and `is_same` constructs.
  */
template <typename FT,
          int V,
          int S,
          bool compress,
          int num_flav,
          typename Spinor1,
          typename Spinor2>
typename std::enable_if<
    std::is_same<const typename Geometry<FT, V, S, compress>::FourSpinorBlock,
                 const Spinor1>::value &&
        std::is_same<const typename Geometry<FT, V, S, compress>::FourSpinorBlock,
                     const Spinor2>::value,
    void>::type
bicgstab_rxupdate(
    typename Geometry<FT, V, S, compress>::FourSpinorBlock *const x[num_flav],
    typename Geometry<FT, V, S, compress>::FourSpinorBlock *const r[num_flav],
    Spinor1 *const t[num_flav],
    Spinor2 *const p[num_flav],
    double omega[2],
    double alpha[2],
    double &r_norm,
    const Geometry<FT, V, S, compress> &geom,
    int n_blas_simt)
{
  // TODO (Martin Ueding): Check whether this is the correct generalization to
  // multiple flavors.
  r_norm = 0;
  for (uint8_t f = 0; f < num_flav; ++f) {
    double local_r_norm;
    bicgstab_rxupdate(
        x[f], r[f], t[f], p[f], omega, alpha, local_r_norm, geom, n_blas_simt);
    r_norm += local_r_norm;
  }
}

template <typename FT, int V, int S, bool compress>
void mr_rxupdate(
    typename Geometry<FT, V, S, compress>::FourSpinorBlock *x,
    typename Geometry<FT, V, S, compress>::FourSpinorBlock *r,
    const typename Geometry<FT, V, S, compress>::FourSpinorBlock *Mr,
    double a[2],
    const Geometry<FT, V, S, compress> &geom,
    int n_blas_simt)
{
  MRRXUpdateFunctor<FT, V, S, compress> f(x, r, Mr, a);
  siteLoopNoReduction<FT,
                     V,
                     S,
                     compress,
                     MRRXUpdateFunctor<FT, V, S, compress>>(
      f, geom, n_blas_simt);
}



/**
  \see The article \ref intel-cpp-compiler-workaround contains an explanation
  of the `typename Spinor1`, `enable_if`, and `is_same` constructs.
  */
template <typename FT,
          int V,
          int S,
          bool compress,
          int num_flav,
          typename Spinor1>
typename std::enable_if<
    std::is_same<const typename Geometry<FT, V, S, compress>::FourSpinorBlock,
                 const Spinor1>::value,void>::type
mr_rxupdate(
        typename Geometry<FT, V, S, compress>::FourSpinorBlock *const  x[num_flav],
        typename Geometry<FT, V, S, compress>::FourSpinorBlock *const  r[num_flav],
    Spinor1 *const Mr[num_flav],
    double a[2],
    const Geometry<FT, V, S, compress> &geom,
    int n_blas_simt)
{
  // TODO Check whether this is the correct generalization to multiple
  // flavors.

  for (uint8_t f = 0; f < num_flav; ++f) {
    mr_rxupdate(
		x[f], r[f], Mr[f], a, geom, n_blas_simt);
  }
}


template <typename FTOut,
          int VOut,
          int SOut,
          bool CompressOut,
          typename FTIn,
          int VIn,
          int SIn,
          bool CompressIn>
void convert(
    typename Geometry<FTOut, VOut, SOut, CompressOut>::FourSpinorBlock *spinor_out,
    double scale_factor,
    const typename Geometry<FTIn, VIn, SIn, CompressIn>::FourSpinorBlock *spinor_in,
    const Geometry<FTOut, VOut, SOut, CompressOut> &geom_out,
    const Geometry<FTIn, VIn, SIn, CompressIn> &geom_in,
    int n_blas_threads)
{
  typedef typename ArithType<FTOut>::Type AT_out;

  // Get the subgrid latt size.
  int Nt = geom_out.Nt();
  int Nz = geom_out.Nz();
  int Ny = geom_out.Ny();
  int Nxh = geom_out.Nxh();
  int nvecs_out = geom_out.nVecs();
  int Pxy_out = geom_out.getPxy();
  int Pxyz_out = geom_out.getPxyz();

  int nvecs_in = geom_in.nVecs();
  int Pxy_in = geom_in.getPxy();
  int Pxyz_in = geom_in.getPxyz();

#pragma omp parallel for collapse(4)
  for (int t = 0; t < Nt; t++) {
    for (int z = 0; z < Nz; z++) {
      for (int y = 0; y < Ny; y++) {
        for (int s = 0; s < nvecs_out; s++) {
          for (int col = 0; col < 3; col++) {
            for (int spin = 0; spin < 4; spin++) {
              for (int x = 0; x < SOut; x++) {
                // Index of the soavector from the output.
                int const ind_out = t * Pxyz_out + z * Pxy_out + y * nvecs_out + s;

                // Actual X coordinate is given by the number of elements in each
                // soavector `SOut` and the position within the soavector `x`.
                int const x_coord = s * SOut + x;

                // The soavector on the current line of X elements is computed from the
                // actual X coordinate `x_coord` and the input soalen `SIn`.
                int const s_in = x_coord / SIn;

                // The position within the input soavector is given by subtracting the X
                // coordinate of the head of the input soavector from the actual X
                // coordinate.
                int const x_in = x_coord - SIn * s_in;

                // The index of the input soavector is computed analogous to the output
                // soavector.
                int const ind_in = t * Pxyz_in + z * Pxy_in + y * nvecs_in + s_in;

                // Iterate over real and imaginary part.
                for (int const reim : {0, 1}) {
                    // Extract source and target from the data structures.
                    auto const source = spinor_in[ind_in][col][spin][reim][x_in];
                    auto &target = spinor_out[ind_out][col][spin][reim][x];

                    // Convert the input numbers into the arithmetic type.
                    auto const scale_factor_rep = rep<AT_out, double>(scale_factor);
                    auto const source_rep = rep<AT_out, FTIn>(source);

                    // Perform the multiplication within the arithmetic types.
                    auto const result_rep = scale_factor_rep * source_rep;

                    // Convert the result back into the storage type.
                    auto const result = rep<FTOut, AT_out>(result_rep);

                    // Assign the result into the target array element.
                    target = result;
                }
              }
            }
          }
        }
      }
    }
  }
}

/**
  \see The article \ref intel-cpp-compiler-workaround contains an explanation
  of the `typename Spinor1`, `enable_if`, and `is_same` constructs.
  */
template <typename FTOut,
          int VOut,
          int SOut,
          bool CompressOut,
          typename FTIn,
          int VIn,
          int SIn,
          bool CompressIn,
          int num_flav,
          typename Spinor1>
typename std::enable_if<
    std::is_same<
        const typename Geometry<FTIn, VIn, SIn, CompressIn>::FourSpinorBlock,
        const Spinor1>::value,
    void>::type
convert(typename Geometry<FTOut, VOut, SOut, CompressOut>::FourSpinorBlock
            *const spinor_out[num_flav],
        double scale_factor,
        Spinor1 *const spinor_in[num_flav],
        const Geometry<FTOut, VOut, SOut, CompressOut> &geom_out,
        const Geometry<FTIn, VIn, SIn, CompressIn> &geom_in,
        int n_blas_threads)
{
  for (uint8_t f = 0; f < num_flav; ++f) {
    convert(spinor_out[f],
            scale_factor,
            spinor_in[f],
            geom_out,
            geom_in,
            n_blas_threads);
  }
}

}; // Namespace
