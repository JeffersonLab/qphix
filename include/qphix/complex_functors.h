#pragma once

#include "qphix/blas_utils.h"
#include "qphix/arith_type.h"
#include <array>

namespace QPhiX
{

template <typename FT, int V, int S, bool compress>
class InnerProductFunctor
{
 public:
  typedef typename ArithType<FT>::Type AT;

  InnerProductFunctor(const typename Geometry<FT, V, S, compress>::FourSpinorBlock *l_,
                      const typename Geometry<FT, V, S, compress>::FourSpinorBlock *r_)
      : l(l_), r(r_)
  {
  }

  ~InnerProductFunctor() {}

  inline void func(int block, std::array<double, S> &reduction_re, std::array<double, S> &reduction_im) const
  {
    int nvec_in_spinor = (3 * 4 * 2 * S) / V;
    const FT *lbase = &l[block][0][0][0][0];
    const FT *rbase = &r[block][0][0][0][0];

#if defined(__GNUG__) && !defined(__INTEL_COMPILER)
    typename Geometry<AT, V, S, compress>::FourSpinorBlock l_spinor
        __attribute__((aligned(QPHIX_LLC_CACHE_ALIGN)));
    typename Geometry<AT, V, S, compress>::FourSpinorBlock r_spinor
        __attribute__((aligned(QPHIX_LLC_CACHE_ALIGN)));
#else
    __declspec(align(QPHIX_LLC_CACHE_ALIGN))
        typename Geometry<AT, V, S, compress>::FourSpinorBlock l_spinor;
    __declspec(align(QPHIX_LLC_CACHE_ALIGN))
        typename Geometry<AT, V, S, compress>::FourSpinorBlock r_spinor;
#endif

    BLASUtils::streamInSpinor<FT, V>((AT *)l_spinor, lbase, nvec_in_spinor);
    BLASUtils::streamInSpinor<FT, V>((AT *)r_spinor, rbase, nvec_in_spinor);

    // Accumulate the inner product from this spinor

    for (int col = 0; col < 3; col++) {
      for (int spin = 0; spin < 4; spin++) {
#pragma omp simd aligned(l_spinor, r_spinor : V)
        for (int s = 0; s < S; s++) {
          reduction_re[s] += l_spinor[col][spin][0][s] * r_spinor[col][spin][0][s] +
                             l_spinor[col][spin][1][s] * r_spinor[col][spin][1][s];
          reduction_im[s] += l_spinor[col][spin][0][s] * r_spinor[col][spin][1][s] -
                             l_spinor[col][spin][1][s] * r_spinor[col][spin][0][s];
        }
      }
    }
  }

 private:
  const typename Geometry<FT, V, S, compress>::FourSpinorBlock *l;
  const typename Geometry<FT, V, S, compress>::FourSpinorBlock *r;
};

template <typename FT, int V, int S, bool compress>
class InnerProductNormFunctor
{
 public:
  typedef typename ArithType<FT>::Type AT;

  InnerProductNormFunctor(
      const typename Geometry<FT, V, S, compress>::FourSpinorBlock *l_,
      const typename Geometry<FT, V, S, compress>::FourSpinorBlock *r_)
      : l(l_), r(r_)
  {
  }

  ~InnerProductNormFunctor() {}

  inline void func(int block, std::array<double, S> &reduction_re, std::array<double, S> &reduction_im, std::array<double, S> &norm_l) const
  {
    int nvec_in_spinor = (3 * 4 * 2 * S) / V;
    const FT *lbase = &l[block][0][0][0][0];
    const FT *rbase = &r[block][0][0][0][0];

#if defined(__GNUG__) && !defined(__INTEL_COMPILER)
    typename Geometry<AT, V, S, compress>::FourSpinorBlock l_spinor
        __attribute__((aligned(QPHIX_LLC_CACHE_ALIGN)));
    typename Geometry<AT, V, S, compress>::FourSpinorBlock r_spinor
        __attribute__((aligned(QPHIX_LLC_CACHE_ALIGN)));
#else
    __declspec(align(QPHIX_LLC_CACHE_ALIGN))
        typename Geometry<AT, V, S, compress>::FourSpinorBlock l_spinor;
    __declspec(align(QPHIX_LLC_CACHE_ALIGN))
        typename Geometry<AT, V, S, compress>::FourSpinorBlock r_spinor;
#endif

    BLASUtils::streamInSpinor<FT, V>((AT *)l_spinor, lbase, nvec_in_spinor);
    BLASUtils::streamInSpinor<FT, V>((AT *)r_spinor, rbase, nvec_in_spinor);

    // Accumulate the inner product from this spinor

    for (int col = 0; col < 3; col++) {
      for (int spin = 0; spin < 4; spin++) {
#pragma omp simd aligned(l_spinor, r_spinor : V)
        for (int s = 0; s < S; s++) {
          reduction_re[s] += l_spinor[col][spin][0][s] * r_spinor[col][spin][0][s] +
                             l_spinor[col][spin][1][s] * r_spinor[col][spin][1][s];
          reduction_im[s] += l_spinor[col][spin][0][s] * r_spinor[col][spin][1][s] -
                             l_spinor[col][spin][1][s] * r_spinor[col][spin][0][s];

	  norm_l[s] += l_spinor[col][spin][0][s] * l_spinor[col][spin][0][s] +
	    l_spinor[col][spin][1][s] * l_spinor[col][spin][1][s];

        }
      }
    }
  }

 private:
  const typename Geometry<FT, V, S, compress>::FourSpinorBlock *l;
  const typename Geometry<FT, V, S, compress>::FourSpinorBlock *r;
};


template <typename FT, int V, int S, bool compress>
class BiCGStabPUpdateFunctor
{
 public:
  typedef typename ArithType<FT>::Type AT;

  BiCGStabPUpdateFunctor(const typename Geometry<FT, V, S, compress>::FourSpinorBlock *r_,
                         typename Geometry<FT, V, S, compress>::FourSpinorBlock *p_,
                         const typename Geometry<FT, V, S, compress>::FourSpinorBlock *v_,
                         double beta_[2],
                         double omega_[2])
      : r(r_), p(p_), v(v_)
  {
    beta[0] = rep<AT, double>(beta_[0]);
    beta[1] = rep<AT, double>(beta_[1]);
    omega[0] = rep<AT, double>(omega_[0]);
    omega[1] = rep<AT, double>(omega_[1]);
  }

  ~BiCGStabPUpdateFunctor() {}

  inline void func(int block) const
  {
    int nvec_in_spinor = (3 * 4 * 2 * S) / V;
    const FT *rbase = &r[block][0][0][0][0];
    FT *pbase = &p[block][0][0][0][0];
    const FT *vbase = &v[block][0][0][0][0];

#if defined(__GNUG__) && !defined(__INTEL_COMPILER)
    typename Geometry<AT, V, S, compress>::FourSpinorBlock r_spinor
        __attribute__((aligned(QPHIX_LLC_CACHE_ALIGN)));
    typename Geometry<AT, V, S, compress>::FourSpinorBlock p_spinor
        __attribute__((aligned(QPHIX_LLC_CACHE_ALIGN)));
    typename Geometry<AT, V, S, compress>::FourSpinorBlock v_spinor
        __attribute__((aligned(QPHIX_LLC_CACHE_ALIGN)));
#else
    __declspec(align(QPHIX_LLC_CACHE_ALIGN))
        typename Geometry<AT, V, S, compress>::FourSpinorBlock r_spinor;
    __declspec(align(QPHIX_LLC_CACHE_ALIGN))
        typename Geometry<AT, V, S, compress>::FourSpinorBlock p_spinor;
    __declspec(align(QPHIX_LLC_CACHE_ALIGN))
        typename Geometry<AT, V, S, compress>::FourSpinorBlock v_spinor;
#endif

    BLASUtils::streamInSpinor<FT, V>((AT *)r_spinor, rbase, nvec_in_spinor);
    BLASUtils::streamInSpinor<FT, V>((AT *)p_spinor, pbase, nvec_in_spinor);
    BLASUtils::streamInSpinor<FT, V>((AT *)v_spinor, vbase, nvec_in_spinor);

    for (int col = 0; col < 3; col++) {
      for (int spin = 0; spin < 4; spin++) {

#if defined(__GNUG__) && !defined(__INTEL_COMPILER)
        AT tmp_cmpx[2][S] __attribute__((aligned(QPHIX_LLC_CACHE_ALIGN)));
#else
        __declspec(align(QPHIX_LLC_CACHE_ALIGN)) AT tmp_cmpx[2][S];
#endif

        // Evaluate: p = r + beta (p - omega v)
        //           p = r + beta tmp
        // with tmp = p - omega

        // tmp = -omega v + p = p - omega v
        BLASUtils::cnmadd<AT, S>(
            tmp_cmpx, omega, v_spinor[col][spin], p_spinor[col][spin]);

        // p = r + beta tmp
        BLASUtils::cmadd<AT, S>(p_spinor[col][spin], beta, tmp_cmpx, r_spinor[col][spin]);
      }
    }

    // Should this be a stream out spinor?
    BLASUtils::writeSpinor<FT, V>(pbase, (const AT *)p_spinor, nvec_in_spinor);
  }

 private:
  const typename Geometry<FT, V, S, compress>::FourSpinorBlock *r;
  typename Geometry<FT, V, S, compress>::FourSpinorBlock *p;
  const typename Geometry<FT, V, S, compress>::FourSpinorBlock *v;
  AT beta[2];
  AT omega[2];
};

template <typename FT, int V, int S, bool compress>
class BiCGStabSUpdateFunctor
{
 public:
  typedef typename ArithType<FT>::Type AT;

  BiCGStabSUpdateFunctor(double alpha_[2],
                         typename Geometry<FT, V, S, compress>::FourSpinorBlock *s_,
                         const typename Geometry<FT, V, S, compress>::FourSpinorBlock *v_)
      : s(s_), v(v_)
  {
    alpha[0] = rep<AT, double>(alpha_[0]);
    alpha[1] = rep<AT, double>(alpha_[1]);
  }

  ~BiCGStabSUpdateFunctor() {}

  inline void func(int block) const
  {
    int nvec_in_spinor = (3 * 4 * 2 * S) / V;
    FT *sbase = &s[block][0][0][0][0];
    const FT *vbase = &v[block][0][0][0][0];

#if defined(__GNUG__) && !defined(__INTEL_COMPILER)
    typename Geometry<AT, V, S, compress>::FourSpinorBlock s_spinor
        __attribute__((aligned(QPHIX_LLC_CACHE_ALIGN)));
    typename Geometry<AT, V, S, compress>::FourSpinorBlock v_spinor
        __attribute__((aligned(QPHIX_LLC_CACHE_ALIGN)));
#else
    __declspec(align(QPHIX_LLC_CACHE_ALIGN))
        typename Geometry<AT, V, S, compress>::FourSpinorBlock s_spinor;
    __declspec(align(QPHIX_LLC_CACHE_ALIGN))
        typename Geometry<AT, V, S, compress>::FourSpinorBlock v_spinor;
#endif

    BLASUtils::streamInSpinor<FT, V>((AT *)s_spinor, sbase, nvec_in_spinor);
    BLASUtils::streamInSpinor<FT, V>((AT *)v_spinor, vbase, nvec_in_spinor);

    for (int col = 0; col < 3; col++) {
      for (int spin = 0; spin < 4; spin++) {
        // s = s - alpha v
        //   = -alpha v + s
        //
        // 's' is not involved in the complex mul, so its components
        // are not mixed, so I will leave it as an alias to the cnmadd
        BLASUtils::cnmadd<AT, S>(
            s_spinor[col][spin], alpha, v_spinor[col][spin], s_spinor[col][spin]);
      }
    }
    BLASUtils::writeSpinor<FT, V>(sbase, (const AT *)s_spinor, nvec_in_spinor);
  }

 private:
  typename Geometry<FT, V, S, compress>::FourSpinorBlock *s;
  const typename Geometry<FT, V, S, compress>::FourSpinorBlock *v;
  AT alpha[2];
};

template <typename FT, int V, int S, bool compress>
class CAXPYFunctor
{
 public:
  typedef typename ArithType<FT>::Type AT;

  CAXPYFunctor(
      double alpha_[2],
      const typename Geometry<FT, V, S, compress>::FourSpinorBlock *x_,
      typename Geometry<FT, V, S, compress>::FourSpinorBlock *y_)
      : x(x_), y(y_)
  {
    alpha[0] = rep<AT, double>(alpha_[0]);
    alpha[1] = rep<AT, double>(alpha_[1]);
  }

  ~CAXPYFunctor() {}

  inline void func(int block) const
  {
    int nvec_in_spinor = (3 * 4 * 2 * S) / V;
    FT *ybase = &y[block][0][0][0][0];
    const FT *xbase = &x[block][0][0][0][0];

#if defined(__GNUG__) && !defined(__INTEL_COMPILER)
    typename Geometry<AT, V, S, compress>::FourSpinorBlock y_spinor
        __attribute__((aligned(QPHIX_LLC_CACHE_ALIGN)));
    typename Geometry<AT, V, S, compress>::FourSpinorBlock x_spinor
        __attribute__((aligned(QPHIX_LLC_CACHE_ALIGN)));
#else
    __declspec(align(QPHIX_LLC_CACHE_ALIGN))
        typename Geometry<AT, V, S, compress>::FourSpinorBlock y_spinor;
    __declspec(align(QPHIX_LLC_CACHE_ALIGN))
        typename Geometry<AT, V, S, compress>::FourSpinorBlock x_spinor;
#endif

    BLASUtils::streamInSpinor<FT, V>((AT *)y_spinor, ybase, nvec_in_spinor);
    BLASUtils::streamInSpinor<FT, V>((AT *)x_spinor, xbase, nvec_in_spinor);

    for (int col = 0; col < 3; col++) {
      for (int spin = 0; spin < 4; spin++) {
        // y =  alpha x + y
        // are not mixed, so I will leave it as an alias to the cnmadd
        AT tmp_cmpx[2][S];
        BLASUtils::cmadd<AT, S>(
            y_spinor[col][spin], alpha, x_spinor[col][spin], y_spinor[col][spin]);
      }
    }

    BLASUtils::writeSpinor<FT, V>(ybase, (const AT *)y_spinor, nvec_in_spinor);
  }

 private:

  const typename Geometry<FT, V, S, compress>::FourSpinorBlock *x;
  typename Geometry<FT, V, S, compress>::FourSpinorBlock *y;
  AT alpha[2];
};








template <typename FT, int V, int S, bool compress>
class BiCGStabRXUpdateFunctor
{
 public:
  typedef typename ArithType<FT>::Type AT;
  BiCGStabRXUpdateFunctor(
      typename Geometry<FT, V, S, compress>::FourSpinorBlock *x_,
      typename Geometry<FT, V, S, compress>::FourSpinorBlock *r_,
      const typename Geometry<FT, V, S, compress>::FourSpinorBlock *t_,
      const typename Geometry<FT, V, S, compress>::FourSpinorBlock *p_,
      double omega_[2],
      double alpha_[2])
      : x(x_), r(r_), t(t_), p(p_)
  {
    omega[0] = rep<AT, double>(omega_[0]);
    omega[1] = rep<AT, double>(omega_[1]);

    alpha[0] = rep<AT, double>(alpha_[0]);
    alpha[1] = rep<AT, double>(alpha_[1]);
  }

  ~BiCGStabRXUpdateFunctor() {}

  inline void func(int block, std::array<double, S> &reduction) const
  {
    int nvec_in_spinor = (3 * 4 * 2 * S) / V;
    FT *xbase = &x[block][0][0][0][0];
    FT *rbase = &r[block][0][0][0][0];
    const FT *tbase = &t[block][0][0][0][0];
    const FT *pbase = &p[block][0][0][0][0];

#if defined(__GNUG__) && !defined(__INTEL_COMPILER)
    typename Geometry<AT, V, S, compress>::FourSpinorBlock x_spinor
        __attribute__((aligned(QPHIX_LLC_CACHE_ALIGN)));
    typename Geometry<AT, V, S, compress>::FourSpinorBlock r_spinor
        __attribute__((aligned(QPHIX_LLC_CACHE_ALIGN)));
    typename Geometry<AT, V, S, compress>::FourSpinorBlock t_spinor
        __attribute__((aligned(QPHIX_LLC_CACHE_ALIGN)));
    typename Geometry<AT, V, S, compress>::FourSpinorBlock p_spinor
        __attribute__((aligned(QPHIX_LLC_CACHE_ALIGN)));
#else
    __declspec(align(QPHIX_LLC_CACHE_ALIGN))
        typename Geometry<AT, V, S, compress>::FourSpinorBlock x_spinor;
    __declspec(align(QPHIX_LLC_CACHE_ALIGN))
        typename Geometry<AT, V, S, compress>::FourSpinorBlock r_spinor;
    __declspec(align(QPHIX_LLC_CACHE_ALIGN))
        typename Geometry<AT, V, S, compress>::FourSpinorBlock t_spinor;
    __declspec(align(QPHIX_LLC_CACHE_ALIGN))
        typename Geometry<AT, V, S, compress>::FourSpinorBlock p_spinor;
#endif

    BLASUtils::streamInSpinor<FT, V>((AT *)x_spinor, xbase, nvec_in_spinor);
    BLASUtils::streamInSpinor<FT, V>((AT *)r_spinor, rbase, nvec_in_spinor);
    BLASUtils::streamInSpinor<FT, V>((AT *)t_spinor, tbase, nvec_in_spinor);
    BLASUtils::streamInSpinor<FT, V>((AT *)p_spinor, pbase, nvec_in_spinor);

    for (int col = 0; col < 3; col++) {
      for (int spin = 0; spin < 4; spin++) {
        AT tmp_cmpx[2][S];

        /* tmp = alpha p + x */
        BLASUtils::cmadd<AT, S>(
            tmp_cmpx, alpha, p_spinor[col][spin], x_spinor[col][spin]);
        /* x = omega r + tmp = omega r + alpha p + x */
        BLASUtils::cmadd<AT, S>(
            x_spinor[col][spin], omega, r_spinor[col][spin], tmp_cmpx);

        /* r = -omega t + r = r - omega t */
        BLASUtils::cnmadd<AT, S>(
            r_spinor[col][spin], omega, t_spinor[col][spin], r_spinor[col][spin]);

        /* accumulate new r_norm into reduction */
        for (int cmpx = 0; cmpx < 2; cmpx++) {
          for (int s = 0; s < S; s++) {
            reduction[s] += (double)r_spinor[col][spin][cmpx][s] *
                            (double)r_spinor[col][spin][cmpx][s];
          }
        }
      }
    }
    /* Write back x and r */
    /* Should these be streamouts */
    BLASUtils::writeSpinor<FT, V>(xbase, (const AT *)x_spinor, nvec_in_spinor);
    BLASUtils::writeSpinor<FT, V>(rbase, (const AT *)r_spinor, nvec_in_spinor);
  }

 private:
  typename Geometry<FT, V, S, compress>::FourSpinorBlock *x;
  typename Geometry<FT, V, S, compress>::FourSpinorBlock *r;
  const typename Geometry<FT, V, S, compress>::FourSpinorBlock *t;
  const typename Geometry<FT, V, S, compress>::FourSpinorBlock *p;
  AT alpha[2];
  AT omega[2];
};

template <typename FT, int V, int S, bool compress>
class TwistedMassFunctor
{
 public:
  typedef typename ArithType<FT>::Type AT;

  TwistedMassFunctor(double const apimu_[2],
                     const typename Geometry<FT, V, S, compress>::FourSpinorBlock *x_,
                     typename Geometry<FT, V, S, compress>::FourSpinorBlock *y_)
      : apimu{rep<AT, double>(apimu_[0]), rep<AT, double>(apimu_[1])}, x(x_), y(y_)
  {
  }

  ~TwistedMassFunctor() {}

  inline void func(int block) const
  {
    BLASUtils::StreamInSpinor<FT, AT, V, S, compress> x_spinor(&x[block][0][0][0][0]);

    BLASUtils::StreamSpinor<FT, AT, V, S, compress> y_spinor(
        false, BLASUtils::StreamOut::stream, &y[block][0][0][0][0]);

    // Now we are hopefully both in L1 and in the right layout so
    for (int col = 0; col < 3; col++) {
      for (int spin = 0; spin < 4; spin++) {
        (spin < 2 ? BLASUtils::cm<AT, S>(
                        y_spinor.get()[col][spin], apimu, x_spinor.get()[col][spin])
                  : BLASUtils::cconjm<AT, S>(
                        y_spinor.get()[col][spin], apimu, x_spinor.get()[col][spin]));
      }
    }
  }

 private:
  AT apimu[2];
  const typename Geometry<FT, V, S, compress>::FourSpinorBlock *x;
  typename Geometry<FT, V, S, compress>::FourSpinorBlock *y;
};

template <typename FT, int V, int S, bool compress>
class TwoFlavTwistedMassFunctor
{
 public:
  typedef typename ArithType<FT>::Type AT;

  TwoFlavTwistedMassFunctor(
      double const *const apimu_,
      double const epsilon_,
      typename Geometry<FT, V, S, compress>::FourSpinorBlock const *const *x_,
      typename Geometry<FT, V, S, compress>::FourSpinorBlock *const *y_)
      : apimu{rep<AT, double>(apimu_[0]), rep<AT, double>(apimu_[1])},
        epsilon(rep<AT, double>(epsilon_)), x{x_[0], x_[1]}, y{y_[0], y_[1]}
  {
  }

#ifdef __INTEL_COMPILER
  TwoFlavTwistedMassFunctor(
      double const * apimu_,
      double const epsilon_,
      typename Geometry<FT, V, S, compress>::FourSpinorBlock *const *x_,
      typename Geometry<FT, V, S, compress>::FourSpinorBlock *const *y_)
      : apimu{rep<AT, double>(apimu_[0]), rep<AT, double>(apimu_[1])},
        epsilon(rep<AT, double>(epsilon_)), x{x_[0], x_[1]}, y{y_[0], y_[1]}
  {
  }
#endif

  inline void func(int block) const
  {
    int nvec_in_spinor = (3 * 4 * 2 * S) / V;
    const FT *x_up_base = &x[0][block][0][0][0][0];
    const FT *x_dn_base = &x[1][block][0][0][0][0];
    FT *y_up_base = &y[0][block][0][0][0][0];
    FT *y_dn_base = &y[1][block][0][0][0][0];

// Temporary storage to stream into and out of
#if defined(__GNUG__) && !defined(__INTEL_COMPILER)
    typename Geometry<AT, V, S, compress>::FourSpinorBlock x_up_spinor
        __attribute__((aligned(QPHIX_LLC_CACHE_ALIGN)));
    typename Geometry<AT, V, S, compress>::FourSpinorBlock x_dn_spinor
        __attribute__((aligned(QPHIX_LLC_CACHE_ALIGN)));
    typename Geometry<AT, V, S, compress>::FourSpinorBlock y_up_spinor
        __attribute__((aligned(QPHIX_LLC_CACHE_ALIGN)));
    typename Geometry<AT, V, S, compress>::FourSpinorBlock y_dn_spinor
        __attribute__((aligned(QPHIX_LLC_CACHE_ALIGN)));
#else
    __declspec(align(QPHIX_LLC_CACHE_ALIGN))
        typename Geometry<AT, V, S, compress>::FourSpinorBlock x_up_spinor;
    __declspec(align(QPHIX_LLC_CACHE_ALIGN))
        typename Geometry<AT, V, S, compress>::FourSpinorBlock y_up_spinor;
    __declspec(align(QPHIX_LLC_CACHE_ALIGN))
        typename Geometry<AT, V, S, compress>::FourSpinorBlock x_dn_spinor;
    __declspec(align(QPHIX_LLC_CACHE_ALIGN))
        typename Geometry<AT, V, S, compress>::FourSpinorBlock y_dn_spinor;
#endif

    BLASUtils::streamInSpinor<FT, V>((AT *)x_up_spinor, x_up_base, nvec_in_spinor);
    BLASUtils::streamInSpinor<FT, V>((AT *)x_dn_spinor, x_dn_base, nvec_in_spinor);

    // Now we are hopefully both in L1 and in the right layout so
    for (int col = 0; col < 3; col++) {
      for (int spin = 0; spin < 4; spin++) {
        // (a + i mu gamma_5 tau_3 + epsilon tau_1) \psi
        (spin < 2 ? BLASUtils::tau3cm_tau1_scaleadd<AT, S>(y_up_spinor[col][spin],
                                                           y_dn_spinor[col][spin],
                                                           apimu,
                                                           epsilon,
                                                           x_up_spinor[col][spin],
                                                           x_dn_spinor[col][spin])
                  : BLASUtils::tau3cconjm_tau1_scaleadd<AT, S>(y_up_spinor[col][spin],
                                                               y_dn_spinor[col][spin],
                                                               apimu,
                                                               epsilon,
                                                               x_up_spinor[col][spin],
                                                               x_dn_spinor[col][spin]));
      }
    }

    BLASUtils::streamOutSpinor<FT, V>(y_up_base, (const AT *)y_up_spinor, nvec_in_spinor);
    BLASUtils::streamOutSpinor<FT, V>(y_dn_base, (const AT *)y_dn_spinor, nvec_in_spinor);
  }

 private:
  AT apimu[2];
  AT epsilon;
  const typename Geometry<FT, V, S, compress>::FourSpinorBlock *x[2];
  typename Geometry<FT, V, S, compress>::FourSpinorBlock *y[2];
};

template <typename FT, int V, int S, bool compress>
class MRRXUpdateFunctor
{
 public:
  typedef typename ArithType<FT>::Type AT;
  MRRXUpdateFunctor(
      typename Geometry<FT, V, S, compress>::FourSpinorBlock *x_,
      typename Geometry<FT, V, S, compress>::FourSpinorBlock *r_,
      const typename Geometry<FT, V, S, compress>::FourSpinorBlock *Mr_,
      double a_[2])
      : x(x_), r(r_), Mr(Mr_)
  {
    a[0] = rep<AT, double>(a_[0]);
    a[1] = rep<AT, double>(a_[1]);
  }

  ~MRRXUpdateFunctor() {}

  inline void func(int block) const
  {
    int nvec_in_spinor = (3 * 4 * 2 * S) / V;
    FT *xbase = &x[block][0][0][0][0];
    FT *rbase = &r[block][0][0][0][0];
    const FT *Mrbase = &Mr[block][0][0][0][0];

#if defined(__GNUG__) && !defined(__INTEL_COMPILER)
    typename Geometry<AT, V, S, compress>::FourSpinorBlock x_spinor
        __attribute__((aligned(QPHIX_LLC_CACHE_ALIGN)));
    typename Geometry<AT, V, S, compress>::FourSpinorBlock r_spinor
        __attribute__((aligned(QPHIX_LLC_CACHE_ALIGN)));
    typename Geometry<AT, V, S, compress>::FourSpinorBlock Mr_spinor
        __attribute__((aligned(QPHIX_LLC_CACHE_ALIGN)));
#else
    __declspec(align(QPHIX_LLC_CACHE_ALIGN))
        typename Geometry<AT, V, S, compress>::FourSpinorBlock x_spinor;
    __declspec(align(QPHIX_LLC_CACHE_ALIGN))
        typename Geometry<AT, V, S, compress>::FourSpinorBlock r_spinor;
    __declspec(align(QPHIX_LLC_CACHE_ALIGN))
        typename Geometry<AT, V, S, compress>::FourSpinorBlock Mr_spinor;
#endif

    BLASUtils::streamInSpinor<FT, V>((AT *)x_spinor, xbase, nvec_in_spinor);
    BLASUtils::streamInSpinor<FT, V>((AT *)r_spinor, rbase, nvec_in_spinor);
    BLASUtils::streamInSpinor<FT, V>((AT *)Mr_spinor,Mrbase, nvec_in_spinor);


    for (int col = 0; col < 3; col++) {
      for (int spin = 0; spin < 4; spin++) {

        /* x += a  r  */
        BLASUtils::cmadd<AT, S>(
            x_spinor[col][spin], a, r_spinor[col][spin], x_spinor[col][spin]);

        /* r  += -a Mr   */
        BLASUtils::cnmadd<AT, S>(
            r_spinor[col][spin], a, Mr_spinor[col][spin], r_spinor[col][spin]);

               
      }
    }
    /* Write back x and r */
    /* Should these be streamouts */
    BLASUtils::writeSpinor<FT, V>(xbase, (const AT *)x_spinor, nvec_in_spinor);
    BLASUtils::writeSpinor<FT, V>(rbase, (const AT *)r_spinor, nvec_in_spinor);
  }

 private:
  typename Geometry<FT, V, S, compress>::FourSpinorBlock *x;
  typename Geometry<FT, V, S, compress>::FourSpinorBlock *r;
  const typename Geometry<FT, V, S, compress>::FourSpinorBlock *Mr;
  AT a[2];
};

}; // Namespace
