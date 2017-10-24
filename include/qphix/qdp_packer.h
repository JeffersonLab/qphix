#pragma once

#include "qdp.h"
#include "qphix/geometry.h"
#include "qphix/full_spinor.h"
#include "qphix/dslash_def.h"
#include "qphix/qphix_config.h"

#ifdef QPHIX_BUILD_QDPJIT
#include "qphix/qdp_packer_qdpjit.h"
#else
#include "qphix/qdp_packer_parscalar.h"
#endif

using namespace QDP;

namespace QPhiX
{

template <typename FT, int veclen, int soalen, bool compress, typename QDPSpinor>
void qdp_pack_spinor(
    const QDPSpinor &psi_in,
    typename Geometry<FT, veclen, soalen, compress>::FourSpinorBlock *psi_even,
    typename Geometry<FT, veclen, soalen, compress>::FourSpinorBlock *psi_odd,
    const Geometry<FT, veclen, soalen, compress> &s)
{
  qdp_pack_cb_spinor(psi_in, psi_even, s, 0);
  qdp_pack_cb_spinor(psi_in, psi_odd, s, 1);
}

template <typename FT, int veclen, int soalen, bool compress, typename QDPSpinor>
void qdp_pack_spinor(
    const QDPSpinor &psi_in,
    FullSpinor<FT, veclen, soalen, compress>& psi,
    const Geometry<FT, veclen, soalen, compress> &s)
{
  qdp_pack_cb_spinor(psi_in, psi.getCBData(0), s, 0);
  qdp_pack_cb_spinor(psi_in, psi.getCBData(1), s, 1);
}

template <typename FT, int veclen, int soalen, bool compress, typename QDPSpinor>
void qdp_unpack_spinor(
    typename Geometry<FT, veclen, soalen, compress>::FourSpinorBlock *chi_even,
    typename Geometry<FT, veclen, soalen, compress>::FourSpinorBlock *chi_odd,
    QDPSpinor &chi,
    const Geometry<FT, veclen, soalen, compress> &s)
{
  qdp_unpack_cb_spinor(chi_even, chi, s, 0);
  qdp_unpack_cb_spinor(chi_odd, chi, s, 1);
}

template <typename FT, int veclen, int soalen, bool compress, typename QDPSpinor>
void qdp_unpack_spinor(
    const FullSpinor<FT,veclen,soalen,compress>& chi_qphix,
    QDPSpinor &chi,
    const Geometry<FT, veclen, soalen, compress> &s)
{
  qdp_unpack_cb_spinor(chi_qphix.getCBData(0), chi, s, 0);
  qdp_unpack_cb_spinor(chi_qphix.getCBData(1), chi, s, 1);
}
#if defined(QPHIX_MIC_SOURCE) || defined(QPHIX_AVX512_SOURCE)

// Downconvert an array of float-vecs to an array of float 16 vecs
inline void downconvert_array(const float *from, half *to, const unsigned int nvecs)
{
#pragma omp parallel for shared(from, to)
  for (int i = 0; i < nvecs; i++) {
    __m512 in = _mm512_load_ps((from + 16 * i));
#if defined(QPHIX_MIC_SOURCE)
    _mm512_extstore_ps((to + 16 * i), in, _MM_DOWNCONV_PS_FLOAT16, _MM_HINT_NT);
#elif defined(QPHIX_AVX512_SOURCE)
    _mm256_store_si256((__m256i *)(to + 16 * i),
                       _mm512_cvtps_ph(in, _MM_FROUND_TO_NEAREST_INT));
#endif
  }
}

// Upconvert an array of float16 vecs to an array of float 32
inline void upconvert_array(const half *from, float *to, const unsigned int nvecs)
{
#pragma omp parallel for shared(from, to)
  for (int i = 0; i < nvecs; i++) {
#if defined(QPHIX_MIC_SOURCE)
    __m512 in = _mm512_extload_ps(
        (from + 16 * i), _MM_UPCONV_PS_FLOAT16, _MM_BROADCAST32_NONE, _MM_HINT_T0);
    _mm512_storenrngo_ps((to + 16 * i), in);
#elif defined(QPHIX_AVX512_SOURCE)
    __m512 in = _mm512_cvtph_ps(_mm256_load_si256((__m256i *)(from + 16 * i)));
    _mm512_stream_ps((to + 16 * i), in);
#endif
  }
}

// Half precision packers...
template <int soalen, bool compress, typename QDPGauge>
void qdp_pack_gauge(
    const QDPGauge &u,
    typename Geometry<half, 16, soalen, compress>::SU3MatrixBlock *u_cb0,
    typename Geometry<half, 16, soalen, compress>::SU3MatrixBlock *u_cb1,
    const Geometry<half, 16, soalen, compress> &s)
{

  typedef typename Geometry<float, 16, soalen, compress>::SU3MatrixBlock GaugeF;

  int latt_size[4];
  int By, Bz, NCores, Sy, Sz, PadXY, PatXYZ, MinCt;
  latt_size[0] = s.Nx();
  latt_size[1] = s.Ny();
  latt_size[2] = s.Nz();
  latt_size[3] = s.Nt();
  // Make a geometry for allocating
  Geometry<float, 16, soalen, compress> g_float(latt_size,
                                                s.getBy(),
                                                s.getBz(),
                                                s.getNumCores(),
                                                s.getSy(),
                                                s.getSz(),
                                                s.getPadXY(),
                                                s.getPadXYZ(),
                                                s.getMinCt());

  GaugeF *tmp_cb0 = (GaugeF *)g_float.allocCBGauge();
  GaugeF *tmp_cb1 = (GaugeF *)g_float.allocCBGauge();
  // CHECK THESE ALLOCATIONS

  // OK Now pack the float
  qdp_pack_gauge<float, 16, soalen, compress, QDPGauge>(
      u, tmp_cb0, tmp_cb1, g_float);

  // This is copied out of the allocation routine.
  // Basically we take all that we have allocated
  // and divide by veclen*sizeof(float) to get the number
  // of vectors to downconvert
  unsigned int n_f_vecs = (((s.getPxyz() * s.Nt() * soalen) / 16) * sizeof(GaugeF)) /
                          (16 * sizeof(float));

  downconvert_array((float *)tmp_cb0, (half *)u_cb0, n_f_vecs);
  downconvert_array((float *)tmp_cb1, (half *)u_cb1, n_f_vecs);

  g_float.free(tmp_cb0);
  g_float.free(tmp_cb1);
}

template <int soalen, bool compress, typename QDPSpinor>
void qdp_pack_spinor(
    const QDPSpinor &psi_in,
    typename Geometry<half, 16, soalen, compress>::FourSpinorBlock *psi_even,
    typename Geometry<half, 16, soalen, compress>::FourSpinorBlock *psi_odd,
    const Geometry<half, 16, soalen, compress> &s)
{
  typedef typename Geometry<float, 16, soalen, compress>::FourSpinorBlock SpinorF;

  int latt_size[4];
  int By, Bz, NCores, Sy, Sz, PadXY, PatXYZ, MinCt;
  latt_size[0] = s.Nx();
  latt_size[1] = s.Ny();
  latt_size[2] = s.Nz();
  latt_size[3] = s.Nt();
  // Make a geometry for allocating
  Geometry<float, 16, soalen, compress> g_float(latt_size,
                                                s.getBy(),
                                                s.getBz(),
                                                s.getNumCores(),
                                                s.getSy(),
                                                s.getSz(),
                                                s.getPadXY(),
                                                s.getPadXYZ(),
                                                s.getMinCt());

  SpinorF *tmp_cb0_alloc = (SpinorF *)g_float.allocCBFourSpinor();
  SpinorF *tmp_cb1_alloc = (SpinorF *)g_float.allocCBFourSpinor();

  SpinorF *tmp_cb0 = tmp_cb0_alloc;
  SpinorF *tmp_cb1 = tmp_cb1_alloc;

  // CHECK THESE ALLOCATIONS

  // OK Now pack the float
  qdp_pack_spinor<float, 16, soalen, compress, QDPSpinor>(
      psi_in, tmp_cb0, tmp_cb1, g_float);

  // This is copied out of the allocation routine.
  // Basically we take all that we have allocated
  // and divide by veclen*sizeof(float) to get the number
  // of vectors to downconvert

  unsigned int n_f_vecs =
      ((s.getPxyz() * s.Nt()) * sizeof(SpinorF)) / (16 * sizeof(float));

  downconvert_array((float *)tmp_cb0, (half *)psi_even, n_f_vecs);
  downconvert_array((float *)tmp_cb1, (half *)psi_odd, n_f_vecs);

  g_float.free(tmp_cb0_alloc);
  g_float.free(tmp_cb1_alloc);
}

template <int soalen, bool compress, typename QDPSpinor>
void qdp_unpack_spinor(
    typename Geometry<half, 16, soalen, compress>::FourSpinorBlock *chi_even,
    typename Geometry<half, 16, soalen, compress>::FourSpinorBlock *chi_odd,
    QDPSpinor &chi,
    const Geometry<half, 16, soalen, compress> &s)
{

  typedef typename Geometry<float, 16, soalen, compress>::FourSpinorBlock SpinorF;

  int latt_size[4];
  int By, Bz, NCores, Sy, Sz, PadXY, PatXYZ, MinCt;
  latt_size[0] = s.Nx();
  latt_size[1] = s.Ny();
  latt_size[2] = s.Nz();
  latt_size[3] = s.Nt();
  // Make a geometry for allocating
  Geometry<float, 16, soalen, compress> g_float(latt_size,
                                                s.getBy(),
                                                s.getBz(),
                                                s.getNumCores(),
                                                s.getSy(),
                                                s.getSz(),
                                                s.getPadXY(),
                                                s.getPadXYZ(),
                                                s.getMinCt());

  // FIXME: CHECK THESE ALLOCATIONS
  SpinorF *tmp_cb0_alloc = (SpinorF *)g_float.allocCBFourSpinor();
  SpinorF *tmp_cb1_alloc = (SpinorF *)g_float.allocCBFourSpinor();

  SpinorF *tmp_cb0 = tmp_cb0_alloc;
  SpinorF *tmp_cb1 = tmp_cb1_alloc;

  // This is copied out of the allocation routine.
  // Basically we take all that we have allocated
  // and divide by veclen*sizeof(float) to get the number
  // of vectors to downconvert

  unsigned int n_f_vecs =
      ((s.getPxyz() * s.Nt()) * sizeof(SpinorF)) / (16 * sizeof(float));

  upconvert_array((half *)chi_even, (float *)tmp_cb0, n_f_vecs);
  upconvert_array((half *)chi_odd, (float *)tmp_cb1, n_f_vecs);

  // OK Now pack the float
  qdp_unpack_spinor<float, 16, soalen, compress, QDPSpinor>(
      tmp_cb0, tmp_cb1, chi, g_float);

  // std::cout << "spinor up free in" << std::endl;
  g_float.free(tmp_cb0_alloc);
  g_float.free(tmp_cb1_alloc);
  // std::cout << "spinor up free out" << std::endl;
}

#if defined(QPHIX_BUILD_CLOVER) || defined(QPHIX_BUILD_TWISTED_MASS_WITH_CLOVER)
template <int soalen, bool compress, typename ClovTerm>
void qdp_pack_clover(
    const ClovTerm &qdp_clov_in,
    typename Geometry<half, 16, soalen, compress>::CloverBlock *cl_out,
    const Geometry<half, 16, soalen, compress> &s,
    int cb)
{
  typedef typename Geometry<float, 16, soalen, compress>::CloverBlock ClovF;

  int latt_size[4];
  int By, Bz, NCores, Sy, Sz, PadXY, PatXYZ, MinCt;
  latt_size[0] = s.Nx();
  latt_size[1] = s.Ny();
  latt_size[2] = s.Nz();
  latt_size[3] = s.Nt();
  // Make a geometry for allocating
  Geometry<float, 16, soalen, compress> g_float(latt_size,
                                                s.getBy(),
                                                s.getBz(),
                                                s.getNumCores(),
                                                s.getSy(),
                                                s.getSz(),
                                                s.getPadXY(),
                                                s.getPadXYZ(),
                                                s.getMinCt());

  ClovF *tmp_clov = (ClovF *)g_float.allocCBClov();

  // This is copied out of the allocation routine.
  // Basically we take all that we have allocated
  // and divide by veclen*sizeof(float) to get the number
  // of vectors to downconvert

  unsigned int n_f_vecs = (((s.getPxyz() * s.Nt() * soalen) / 16) * sizeof(ClovF)) /
                          (16 * sizeof(float));

  // OK Now pack the float
  qdp_pack_clover<float, 16, soalen, compress, ClovTerm>(
      qdp_clov_in, tmp_clov, g_float, cb);

  downconvert_array((float *)tmp_clov, (half *)cl_out, n_f_vecs);
  // std::cout << "clov down free in" << std::endl;
  g_float.free(tmp_clov);
  // std::cout << "clov down free out" << std::endl;
}
#endif // Build Clover

#endif // if defined(QPHIX_MIC_SOURCE) || defined(QPHIX_AVX512_SOURCE)
};
