#pragma once

#include "qphix/arith_type.h"

namespace QPhiX
{
namespace BLASUtils
{

template <typename FT, int S>
inline void
cm(FT res[2][S], const FT alpha[2], const FT x[2][S])
{
#pragma omp simd aligned(res, x : S)
  for( int s = 0; s < S; s++){
    res[0][s] = alpha[0] * x[0][s] - alpha[1] * x[1][s];
    res[1][s] = alpha[0] * x[1][s] + alpha[1] * x[0][s];
  }
}

template <typename FT, int S>
inline void
cconjm(FT res[2][S], const FT alpha[2], const FT x[2][S])
{
#pragma omp simd aligned(res, x : S)
  for( int s = 0; s < S; s++){
    res[0][s] = alpha[0] * x[0][s] + alpha[1] * x[1][s];
    res[1][s] = alpha[0] * x[1][s] - alpha[1] * x[0][s];
  }
}

//  res = alpha x + y
//  alpha is complex
template <typename FT, int S>
inline void
cmadd(FT res[2][S], const FT alpha[2], const FT x[2][S], const FT y[2][S])
{
//  (a[RE] x[RE] - a[IM] y[IM])  + res[RE]
//  (a[RE] y[IM] + a[IM] y[RE])  + res[IM]
#pragma omp simd aligned(res, x, y : S)
  for (int s = 0; s < S; s++) {
    res[0][s] = alpha[0] * x[0][s] - alpha[1] * x[1][s] + y[0][s];
    res[1][s] = alpha[0] * x[1][s] + alpha[1] * x[0][s] + y[1][s];
  }
}

// res = -alpha x + y
// res = y - alpha x
// alpha is complex
template <typename FT, int S>
inline void
cnmadd(FT res[2][S], const FT alpha[2], const FT x[2][S], const FT y[2][S])
{
//  res[RE] -(a[RE] x[RE] - a[IM] y[IM])
// =res[RE] - a[RE] x[RE] + a[IM] y[IM]

//  res[IM] -(a[RE] y[IM] + a[IM] y[RE])
// =res[IM] -a[RE]y[IM] - a[IM] y[RE]
#pragma omp simd aligned(res, x, y : S)
  for (int s = 0; s < S; s++) {
    res[0][s] = y[0][s] - alpha[0] * x[0][s] + alpha[1] * x[1][s];
    res[1][s] = y[1][s] - alpha[0] * x[1][s] - alpha[1] * x[0][s];
  }
}

/**
  Computes (α + iμτ³ + βτ¹) x.

  \param[out] res_up Array with complex index and then \c soalen index.
  \param[out] res_dn Similar, down flavor
  \param[in] alpha Complex number α
  \param[in] beta Real number β
  \param[in] x_up Input up flavor spinor
  \param[in] x_dn Similar, down flavor
  */
template <typename AT, int soalen>
inline void tau3cm_tau1_scaleadd(AT res_up[2][soalen],
                                 AT res_dn[2][soalen],
                                 const AT alpha[2],
                                 const AT beta,
                                 const AT x_up[2][soalen],
                                 const AT x_dn[2][soalen])
{
#pragma omp simd aligned(res_up, res_dn, x_up, x_dn : soalen)
  for (int s = 0; s < soalen; s++) {
    res_up[RE][s] = alpha[RE] * x_up[RE][s] - alpha[IM] * x_up[IM][s] + beta*x_dn[RE][s];
    res_up[IM][s] = alpha[RE] * x_up[IM][s] + alpha[IM] * x_up[RE][s] + beta*x_dn[IM][s];

    res_dn[RE][s] = alpha[RE] * x_dn[RE][s] + alpha[IM] * x_dn[IM][s] + beta*x_up[RE][s];
    res_dn[IM][s] = alpha[RE] * x_dn[IM][s] - alpha[IM] * x_dn[RE][s] + beta*x_up[IM][s];
  }
}

/**
  Like \ref tau3cm_tau1_scaleadd but with iμ complex conjugated in the process.
  */
template <typename AT, int S>
inline void tau3cconjm_tau1_scaleadd(AT res_up[2][S],
                                     AT res_dn[2][S],
                                     const AT alpha[2],
                                     const AT beta,
                                     const AT x_up[2][S],
                                     const AT x_dn[2][S])
{
#pragma omp simd aligned(res_up, res_dn, x_up, x_dn : S)
  for (int s = 0; s < S; s++) {
    res_up[RE][s] = alpha[RE] * x_up[RE][s] + alpha[IM] * x_up[IM][s] + beta*x_dn[RE][s];
    res_up[IM][s] = alpha[RE] * x_up[IM][s] - alpha[IM] * x_up[RE][s] + beta*x_dn[IM][s];

    res_dn[RE][s] = alpha[RE] * x_dn[RE][s] - alpha[IM] * x_dn[IM][s] + beta*x_up[RE][s];
    res_dn[IM][s] = alpha[RE] * x_dn[IM][s] + alpha[IM] * x_dn[RE][s] + beta*x_up[IM][s];
  }
}

// Generic stream in
template <typename FT, int V>
inline void streamInSpinor(FT *dst, const FT *src, int numvec)
{

#if defined(__MIC__)
  // Intel MIC
  const int prefdist1 = 12;
  const int prefdist2 = 64;

  const char *prefl1base = (const char *)src + prefdist1 * 64;
  const char *prefl2base = (const char *)src + prefdist2 * 64;

  for (int v = 0; v < numvec; v++) {
    _mm_prefetch(&prefl1base[v * V * sizeof(FT)], _MM_HINT_T0);
    _mm_prefetch(&prefl2base[v * V * sizeof(FT)], _MM_HINT_T1);
#pragma omp simd aligned(dst, src : V)
    for (int s = 0; s < V; s++) {
      dst[v * V + s] = src[v * V + s];
    }
  }
#else
  // Generic
  for (int v = 0; v < numvec; v++) {
#pragma omp simd aligned(dst, src : V)
    for (int s = 0; s < V; s++) {
      dst[v * V + s] = src[v * V + s];
    }
  }
#endif
}

// Generic write  out
template <typename FT, int V>
inline void writeSpinor(FT *dst, const FT *src, int numvec)
{
  for (int v = 0; v < numvec; v++) {
#pragma omp simd aligned(dst, src : V)
    for (int s = 0; s < V; s++) {
      dst[v * V + s] = src[v * V + s];
    }
  }
}

// Generic stream out
template <typename FT, int V>
inline void streamOutSpinor(FT *dst, const FT *src, int numvec)
{

  for (int v = 0; v < numvec; v++) {
#pragma vector temporal(dst)
#pragma omp simd aligned(dst, src : V)
    for (int s = 0; s < V; s++) {
      dst[v * V + s] = src[v * V + s];
    }
  }
}

// Stream In to a different type
template <typename FT, int V>
inline void
streamInSpinor(typename ArithType<FT>::Type *dst, const FT *src, int numvec)
{

#if defined(__MIC__)
  // Intel MIC
  const int prefdist1 = 12;
  const int prefdist2 = 64;

  const char *prefl1base = (const char *)src + prefdist1 * 64;
  const char *prefl2base = (const char *)src + prefdist2 * 64;

  for (int v = 0; v < numvec; v++) {
    _mm_prefetch(&prefl1base[v * V * sizeof(FT)], _MM_HINT_T0);
    _mm_prefetch(&prefl2base[v * V * sizeof(FT)], _MM_HINT_T1);

#pragma omp simd
    for (int s = 0; s < V; s++) {
      dst[v * V + s] = src[v * V + s];
    }
  }
#else
  // Generic
  for (int v = 0; v < numvec; v++) {
#pragma vector temporal(dst)
#pragma omp simd aligned(dst, src : V)
    for (int s = 0; s < V; s++) {
      dst[v * V + s] = src[v * V + s];
    }
  }
#endif
}

// Write out to a different type
template <typename FT, int V>
inline void writeSpinor(FT *dst, const typename ArithType<FT>::Type *src, int numvec)
{

  for (int v = 0; v < numvec; v++) {
#pragma omp simd aligned(dst, src : V)
    for (int s = 0; s < V; s++) {
      dst[v * V + s] = src[v * V + s];
    }
  }
}

// Stream out to a different type
template <typename FT, int V>
inline void
streamOutSpinor(FT *dst, const typename ArithType<FT>::Type *src, int numvec)
{

  for (int v = 0; v < numvec; v++) {
#pragma omp simd aligned(dst, src : V)
    for (int s = 0; s < V; s++) {
      dst[v * V + s] = src[v * V + s];
    }
  }
}

#if defined(__MIC__)
#include <immintrin.h>

// Half prec specicialize
template <>
inline void streamInSpinor<half, 16>(typename ArithType<half>::Type *dst,
                                     const half *src,
                                     int numvec)
{

  const int prefdist1 = 12;
  const int prefdist2 = 64;

  const char *prefl1base = (const char *)src + prefdist1 * 64;
  const char *prefl2base = (const char *)src + prefdist2 * 64;

  for (int v = 0; v < numvec; v++) {
    _mm_prefetch(&prefl1base[v * 16 * sizeof(half)], _MM_HINT_T0);
    _mm_prefetch(&prefl2base[v * 16 * sizeof(half)], _MM_HINT_T1);

    __m512 r = _mm512_extload_ps((void *)&src[v * 16],
                                 _MM_UPCONV_PS_FLOAT16,
                                 _MM_BROADCAST32_NONE,
                                 _MM_HINT_T0);
    _mm512_store_ps((void *)&dst[v * 16], r);
  }
}

template <>
inline void writeSpinor<half, 16>(half *dst,
                                  const typename ArithType<half>::Type *src,
                                  int numvec)
{
  const int prefdist1 = 12;
  const int prefdist2 = 64;

  const char *prefl1base = (const char *)src + prefdist1 * 64;
  const char *prefl2base = (const char *)src + prefdist2 * 64;

  const char *prefl1baseo = (const char *)dst + prefdist1 * 64;
  const char *prefl2baseo = (const char *)dst + prefdist2 * 64;

  for (int v = 0; v < numvec; v++) {
    _mm_prefetch(&prefl1baseo[v * 16 * sizeof(half)],
                 _MM_HINT_T0); // Prefetch for write
    _mm_prefetch(&prefl2baseo[v * 16 * sizeof(half)],
                 _MM_HINT_T1); // Prefetch for write

    __m512 r = _mm512_load_ps((void *)&src[v * 16]);
    _mm512_extstore_ps(
        (void *)&dst[v * 16], r, _MM_DOWNCONV_PS_FLOAT16, _MM_HINT_T0);
  }
}

template <>
inline void streamOutSpinor<half, 16>(half *dst,
                                      const typename ArithType<half>::Type *src,
                                      int numvec)
{
  const int prefdist1 = 12;
  const int prefdist2 = 64;

  const char *prefl1base = (const char *)src + prefdist1 * 64;
  const char *prefl2base = (const char *)src + prefdist2 * 64;

  for (int v = 0; v < numvec; v++) {
    _mm_prefetch(&prefl1base[v * 16 * sizeof(half)], _MM_HINT_T0);
    _mm_prefetch(&prefl2base[v * 16 * sizeof(half)], _MM_HINT_T1);

    __m512 r = _mm512_load_ps((void *)&src[v * 16]);
    _mm512_extstore_ps(
        (void *)&dst[v * 16], r, _MM_DOWNCONV_PS_FLOAT16, _MM_HINT_NT);
  }
}

#endif // defined MIC

enum class StreamOut { none, stream, write };

/**
  Holds a spinor that is streamed in or out or both.

  The BLAS functors are implemented without intrinsics. This means that they
  cannot do any calculations with half-precision data. The streaming will
  convert from the underlying type (`FT`) into an arithmetic type (`AT`) such
  that normal `float` operations can be done. This will probably not yield the
  performance of the generated kernels using the half-precision intrinsics. The
  BLAS routines will become easier to write this way, one just needs `#pragma
  omp simd`.

  \tparam FT Type of the spinor data structure for the whole volume.
  \tparam AT Type of the temporary site that is used within the BLAS function.
  */
template <typename FT, typename AT, int veclen, int soalen, bool compress12>
class StreamInSpinor
{
 public:
#if defined(__GNUG__) && !defined(__INTEL_COMPILER)
  typedef typename Geometry<AT, veclen, soalen, compress12>::FourSpinorBlock Spinor
      __attribute__((aligned(QPHIX_LLC_CACHE_ALIGN)));
#else
  typedef __declspec(align(QPHIX_LLC_CACHE_ALIGN))
      typename Geometry<AT, veclen, soalen, compress12>::FourSpinorBlock Spinor;
#endif

  StreamInSpinor(FT const *const base)
      : data_base_(base),
        tmp_base_(static_cast<AT *>(&tmp_[0][0][0][0]))
  {
      streamInSpinor<FT, veclen>(tmp_base_, data_base_, nvec_in_spinor_);
  }

  Spinor &get() { return tmp_; }

 private :
  static int constexpr nvec_in_spinor_ = (3 * 4 * 2 * soalen) / veclen;

  Spinor tmp_;
  FT const *const data_base_;
  AT *const tmp_base_;
};

template <typename FT, typename AT, int veclen, int soalen, bool compress12>
class StreamSpinor
{
 public:
#if defined(__GNUG__) && !defined(__INTEL_COMPILER)
  typedef typename Geometry<AT, veclen, soalen, compress12>::FourSpinorBlock Spinor
      __attribute__((aligned(QPHIX_LLC_CACHE_ALIGN)));
#else
  typedef __declspec(align(QPHIX_LLC_CACHE_ALIGN))
      typename Geometry<AT, veclen, soalen, compress12>::FourSpinorBlock Spinor;
#endif

  StreamSpinor(bool const stream_in, StreamOut const stream_out, FT *const base)
      : stream_in_(stream_in), stream_out_(stream_out), data_base_(base),
        tmp_base_(static_cast<AT *>(&tmp_[0][0][0][0]))
  {
    if (stream_in_) {
      streamInSpinor<FT, veclen>(tmp_base_, data_base_, nvec_in_spinor_);
    }
  }

  /**
    XXX (Martin Ueding): One must not throw any exceptions in the destructor
    because cleanup must always work. Therefore it might be nicer to have a
    `finalize` method here that one has to call manually. However, the
    streaming or writing out will just do some array element assignment. So
    that should work just fine without any exceptions.
    */
  ~StreamSpinor()
  {
    if (stream_out_ == StreamOut::stream) {
      streamOutSpinor<FT, veclen>(data_base_, tmp_base_, nvec_in_spinor_);
    } else if (stream_out_ == StreamOut::write) {
      writeSpinor<FT, veclen>(data_base_, tmp_base_, nvec_in_spinor_);
    }
  }

  Spinor &get() { return tmp_; }

 private :
  static int constexpr nvec_in_spinor_ = (3 * 4 * 2 * soalen) / veclen;

  bool const stream_in_;
  StreamOut const stream_out_;

  Spinor tmp_;
  FT *const data_base_;
  AT *const tmp_base_;
};

} // Namespace BLASUtils
} // Namespace QPhiX
