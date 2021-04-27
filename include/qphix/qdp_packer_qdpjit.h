#pragma once

#include "qphix/diagnostics.h"

QPHIX_MESSAGE("building with QDPJIT packers")

#include "qdp.h"
#include "qphix/geometry.h"

#include "qphix/dslash_def.h"
#include "qphix/qphix_config.h"

#if defined(QPHIX_MIC_SOURCE) || defined(QPHIX_AVX512_SOURCE)
#include <immintrin.h>
#endif

#if defined(QPHIX_BUILD_CLOVER) || defined(QPHIX_BUILD_TWISTED_MASS_WITH_CLOVER)
#include "qphix/clover_dslash_def.h"
#endif

#undef DEBUG_PACKER
using namespace QDP;

#define MANUAL_COLLAPSE 1

namespace QPhiX
{

template <typename FT, int veclen, int soalen, bool compress, typename QDPGauge>
void qdp_pack_gauge(
    const QDPGauge &u,
    typename Geometry<FT, veclen, soalen, compress>::SU3MatrixBlock *u_cb0,
    typename Geometry<FT, veclen, soalen, compress>::SU3MatrixBlock *u_cb1,
    const Geometry<FT, veclen, soalen, compress> &s)
{
  // Get the subgrid latt size.
  int64_t Nt = (int64_t)s.Nt();
  int64_t Nz = (int64_t)s.Nz();
  int64_t Ny = (int64_t)s.Ny();
  int64_t nvecs = s.nVecs();
  int64_t nyg = s.nGY();
  int64_t Pxy = s.getPxy();
  int64_t Pxyz = s.getPxyz();
  int64_t qdp_inner_dim = (int64_t)getDataLayoutInnerSize();
  const int64_t n_complex_dim = 2;

#if 1
  int64_t volcb = rb[0].numSiteTable();
#endif

  // Shift the lattice to get U(x-mu)
  QDPGauge u_minus(4);
  for (int mu = 0; mu < 4; mu++) {
    u_minus[mu] = shift(u[mu], BACKWARD, mu);
  }

  // This ought to be the underlying precision type (float/double) in use
  // by QDP-JIT and may well be different from FT
  // So casting will be needed
  using WT = typename WordType<typename QDPGauge::InnerType_t>::Type_t;

#ifdef DEBUG_PACKER
  QDPIO::cout << "QDP GaugePacker: sizeof(FT)=" << sizeof(FT)
              << " sizeof(WT)=" << sizeof(WT) << " qdp_inner=" << qdp_inner_dim
              << " SOA=" << soalen << " veclen=" << veclen << std::endl;
#endif

#ifndef MANUAL_COLLAPSE
#pragma omp parallel for collapse(4)
  for (int64_t t = 0; t < Nt; t++) {
    for (int64_t z = 0; z < Nz; z++) {
      for (int64_t y = 0; y < Ny; y++) {
        for (int64_t s = 0; s < nvecs; s++) {
#else

  // cindex = combined index = s + nvecs*(y + Ny*(z + Nz*t ) )
  // in case collapse(4) doesnt work in OpenMP
  int max_cindex = Nt * Nz * Ny * nvecs;

#pragma omp parallel for
  for (int cindex = 0; cindex < max_cindex; cindex++) {
    int64_t t1 = cindex / nvecs;
    int64_t s = cindex - nvecs * t1;

    int64_t t2 = t1 / Ny;
    int64_t y = t1 - Ny * t2;

    int64_t t = t2 / Nz;
    int64_t z = t2 - Nz * t;

#endif
          for (int64_t mu = 0; mu < 4; mu++) {

            const WT *u_ptr = (const WT *)((u[mu]).getFjit());
            const WT *u_minus_ptr = (const WT *)((u_minus[mu]).getFjit());

            int64_t outer_c = 3;
            if (compress) {
              outer_c = 2;
            }
            for (int64_t c = 0; c < outer_c; c++) {
              for (int64_t c2 = 0; c2 < 3; c2++) {
                for (int64_t x = 0; x < soalen; x++) {

                  int64_t block = (t * Pxyz + z * Pxy) / nyg + (y / nyg) * nvecs + s;
                  int64_t xx = (y % nyg) * soalen + x;

                  // QDP-JIT in OCSRI layout:
                  //  Outer x Color x Spin x Real x Inner
                  //
                  // However, it hardly matters maybe because
                  // There is no spin on the gauge fields

                  int64_t qdpsite =
                      x + soalen * (s + nvecs * (y + Ny * (z + Nz * t)));

#if 0
										int64_t qdpsite_rb0 = qdpsite+rb[0].start();
#else
            int64_t qdpsite_rb0 = qdpsite;
#endif
                  int64_t qdp_outer_rb0 = qdpsite_rb0 / qdp_inner_dim;
                  int64_t qdp_inner_rb0 = qdpsite_rb0 % qdp_inner_dim;
                  int64_t rb0_offset =
                      qdp_inner_rb0 +
                      qdp_inner_dim * n_complex_dim * Nc * Nc * qdp_outer_rb0;

                  u_cb0[block][2 * mu][c][c2][RE][xx] =
                      (FT)u_minus_ptr[rb0_offset +
                                      qdp_inner_dim *
                                          (RE + n_complex_dim * (c + Nc * c2))];
                  u_cb0[block][2 * mu][c][c2][IM][xx] =
                      (FT)u_minus_ptr[rb0_offset +
                                      qdp_inner_dim *
                                          (IM + n_complex_dim * (c + Nc * c2))];
                  u_cb0[block][2 * mu + 1][c][c2][RE][xx] = (FT)
                      u_ptr[rb0_offset +
                            qdp_inner_dim * (RE + n_complex_dim * (c + Nc * c2))];
                  u_cb0[block][2 * mu + 1][c][c2][IM][xx] = (FT)
                      u_ptr[rb0_offset +
                            qdp_inner_dim * (IM + n_complex_dim * (c + Nc * c2))];

#if 0
										int64_t qdpsite_rb1 = qdpsite+rb[1].start();
#else
            int64_t qdpsite_rb1 = qdpsite + rb[1].numSiteTable();
#endif
                  int64_t qdp_outer_rb1 = qdpsite_rb1 / qdp_inner_dim;
                  int64_t qdp_inner_rb1 = qdpsite_rb1 % qdp_inner_dim;
                  int64_t rb1_offset =
                      qdp_inner_rb1 +
                      qdp_inner_dim * n_complex_dim * Nc * Nc * qdp_outer_rb1;

                  u_cb1[block][2 * mu][c][c2][RE][xx] =
                      (FT)u_minus_ptr[rb1_offset +
                                      qdp_inner_dim *
                                          (RE + n_complex_dim * (c + Nc * c2))];
                  u_cb1[block][2 * mu][c][c2][IM][xx] =
                      (FT)u_minus_ptr[rb1_offset +
                                      qdp_inner_dim *
                                          (IM + n_complex_dim * (c + Nc * c2))];
                  u_cb1[block][2 * mu + 1][c][c2][RE][xx] = (FT)
                      u_ptr[rb1_offset +
                            qdp_inner_dim * (RE + n_complex_dim * (c + Nc * c2))];
                  u_cb1[block][2 * mu + 1][c][c2][IM][xx] = (FT)
                      u_ptr[rb1_offset +
                            qdp_inner_dim * (IM + n_complex_dim * (c + Nc * c2))];

                } // x
              } // c2
            } // c
          } // mu

#ifndef MANUAL_COLLAPSE
        } // collapsed s
      } // collapsed y
    } // collapsed z
  } // collapsed t
#else
  } // Hanc Collapsed OMP
#endif
}

#if defined(QPHIX_BUILD_CLOVER) || defined(QPHIX_BUILD_TWISTED_MASS_WITH_CLOVER)

// This accesses the Internals of the LLVMCloverTerm
template <typename FT, int veclen, int soalen, bool compress, typename ClovTerm>
void qdp_pack_clover(
    const ClovTerm &qdp_clov_in,
    typename ClovDslash<FT, veclen, soalen, compress>::CloverBlock *cl_out,
    const Geometry<FT, veclen, soalen, compress> &s,
    int cb)
{
  // Get the subgrid latt size.
  int64_t Nt = s.Nt();
  int64_t Nz = s.Nz();
  int64_t Ny = s.Ny();
  int64_t nvecs = s.nVecs();
  int64_t nyg = s.nGY();
  int64_t Pxy = s.getPxy();
  int64_t Pxyz = s.getPxyz();

  // Sanity Check
  // QDP Type is
  // Outer x 2 chiral blocks x 6 floats x inner sites
  // JIT-LLVM CLOVER will need to expose DiagType and OffDiagType
  // It does in the testcase here but make sure chroma version does it too.

  using DiagType = typename ClovTerm::DiagType;
  using OffDiagType = typename ClovTerm::OffDiagType;

  // This is the base type used here. Either double or float
  using WT = typename WordType<DiagType>::Type_t;

  const DiagType &diag_term = qdp_clov_in.getDiagBuffer();
  const OffDiagType &off_diag_term = qdp_clov_in.getOffDiagBuffer();

  const WT *diag_buf = (const WT *)diag_term.getFjit();
  const WT *off_diag_buf = (const WT *)off_diag_term.getFjit();

  int64_t qdp_inner_dim = (int64_t)getDataLayoutInnerSize();
  const int64_t n_complex_dim = 2;
  const int64_t diag_dim = 6;
  const int64_t offdiag_dim = 15;
  const int64_t chiral_dim = 2;
  const int64_t diag_offset = qdp_inner_dim * diag_dim * chiral_dim;
  const int64_t offdiag_offset =
      qdp_inner_dim * offdiag_dim * n_complex_dim * chiral_dim;

#ifdef DEBUG_PACKER
  QDPIO::cout << "QDP CloverPacker: sizeof(FT)=" << sizeof(FT)
              << " sizeof(WT)=" << sizeof(WT) << " qdp_inner=" << qdp_inner_dim
              << " SOA=" << soalen << " veclen=" << veclen << std::endl;
#endif
// No elem calls in parallel region

#ifndef MANUAL_COLLAPSE
#pragma omp parallel for collapse(4)
  for (int64_t t = 0; t < Nt; t++) {
    for (int64_t z = 0; z < Nz; z++) {
      for (int64_t y = 0; y < Ny; y++) {
        for (int64_t s = 0; s < nvecs; s++) {
#else
  // cindex = combined index = s + nvecs*(y + Ny*(z + Nz*t ) )
  // in case collapse(4) doesnt work in OpenMP
  int max_cindex = Nt * Nz * Ny * nvecs;

#pragma omp parallel for
  for (int cindex = 0; cindex < max_cindex; cindex++) {
    int64_t t1 = cindex / nvecs;
    int64_t s = cindex - nvecs * t1;

    int64_t t2 = t1 / Ny;
    int64_t y = t1 - Ny * t2;

    int64_t t = t2 / Nz;
    int64_t z = t2 - Nz * t;

#endif
          for (int64_t x = 0; x < soalen; x++) {

            int64_t block = (t * Pxyz + z * Pxy) / nyg + (y / nyg) * nvecs + s;
            int64_t xx = (y % nyg) * soalen + x;
            int64_t qdpsite = x + soalen * (s + nvecs * (y + Ny * (z + Nz * t))) +
                              cb * rb[cb].numSiteTable();
            int64_t qdp_outer_idx = qdpsite / qdp_inner_dim;
            int64_t qdp_inner_idx = qdpsite % qdp_inner_dim;

            const WT *diag_base = &diag_buf[diag_offset * qdp_outer_idx];
            const WT *offdiag_base = &off_diag_buf[offdiag_offset * qdp_outer_idx];

            // WARNING: THIS WORKS ONLY IN OCSRI Layout in QDP_JIT
            // For OSCRI Layout: chiral component and diagonal must be transpose
            //                   chiral component and off diagonal must also be
            //                   transposed
            // Logical Outer x Spin x Color x Inner => Outer x Comp x Diag x Inner
            // But we are using OCSRI => Outer Diag Comp Inner

            for (int64_t d = 0; d < 6; d++) {
              cl_out[block].diag1[d][xx] = (FT)(
                  diag_base[qdp_inner_idx + qdp_inner_dim * (0 + chiral_dim * d)]);
            }
            for (int64_t od = 0; od < 15; od++) {
              cl_out[block].off_diag1[od][RE][xx] = (FT)(
                  offdiag_base[qdp_inner_idx +
                               qdp_inner_dim *
                                   (RE + n_complex_dim * (0 + chiral_dim * od))]);

              cl_out[block].off_diag1[od][IM][xx] = (FT)(
                  offdiag_base[qdp_inner_idx +
                               qdp_inner_dim *
                                   (IM + n_complex_dim * (0 + chiral_dim * od))]);
            }
            for (int64_t d = 0; d < 6; d++) {
              cl_out[block].diag2[d][xx] = (FT)(
                  diag_base[qdp_inner_idx + qdp_inner_dim * (1 + chiral_dim * d)]);
            }
            for (int64_t od = 0; od < 15; od++) {
              cl_out[block].off_diag2[od][RE][xx] = (FT)(
                  offdiag_base[qdp_inner_idx +
                               qdp_inner_dim *
                                   (RE + n_complex_dim * (1 + chiral_dim * od))]);

              cl_out[block].off_diag2[od][IM][xx] = (FT)(
                  offdiag_base[qdp_inner_idx +
                               qdp_inner_dim *
                                   (IM + n_complex_dim * (1 + chiral_dim * od))]);
            }
          } // x
#ifndef MANUAL_COLLAPSE
        } // collapsed s
      } // collapsed y
    } // collapsed z
  } // collapsed t
#else
  } // Hand collapsed cindex
#endif
} // function
#endif // IFDEF BUILD CLOVER

template <typename FT, int veclen, int soalen, bool compress, typename QDPSpinor>
void qdp_pack_cb_spinor(
    const QDPSpinor &psi_in,
    typename Geometry<FT, veclen, soalen, compress>::FourSpinorBlock *psi,
    const Geometry<FT, veclen, soalen, compress> &s,
    int cb)
{
  // Get the subgrid latt size.
  int64_t Nt = s.Nt();
  int64_t Nz = s.Nz();
  int64_t Ny = s.Ny();
  int64_t Nxh = s.Nxh();
  int64_t nvecs = s.nVecs();
  int64_t Pxy = s.getPxy();
  int64_t Pxyz = s.getPxyz();

  int64_t qdp_inner_dim = (int64_t)getDataLayoutInnerSize();
  const int64_t n_complex_dim = 2;
  using WT = typename WordType<QDPSpinor>::Type_t;

  const int64_t outer_block_size = qdp_inner_dim * n_complex_dim * Nc * Ns;
#ifdef DEBUG_PACKER
  QDPIO::cout << "QDP CBSpinorPacker: sizeof(FT)=" << sizeof(FT)
              << " sizeof(WT)=" << sizeof(WT) << " qdp_inner=" << qdp_inner_dim
              << " SOA=" << soalen << " veclen=" << veclen << std::endl;
#endif

#ifndef MANUAL_COLLAPSE
#pragma omp parallel for collapse(4)
  for (int64_t t = 0; t < Nt; t++) {
    for (int64_t z = 0; z < Nz; z++) {
      for (int64_t y = 0; y < Ny; y++) {
        for (int64_t s = 0; s < nvecs; s++) {
#else
  // cindex = combined index = s + nvecs*(y + Ny*(z + Nz*t ) )
  // in case collapse(4) doesnt work in OpenMP
  int max_cindex = Nt * Nz * Ny * nvecs;

#pragma omp parallel for
  for (int cindex = 0; cindex < max_cindex; cindex++) {
    int64_t t1 = cindex / nvecs;
    int64_t s = cindex - nvecs * t1;

    int64_t t2 = t1 / Ny;
    int64_t y = t1 - Ny * t2;

    int64_t t = t2 / Nz;
    int64_t z = t2 - Nz * t;
#endif

          for (int64_t col = 0; col < 3; col++) {
            for (int64_t spin = 0; spin < 4; spin++) {
              for (int64_t x = 0; x < soalen; x++) {

                int64_t ind =
                    t * Pxyz + z * Pxy + y * nvecs + s; //((t*Nz+z)*Ny+y)*nvecs+s;
                int64_t x_coord = s * soalen + x;
#if 0

									int64_t qdp_index = ((t*Nz + z)*Ny + y)*Nxh + x_coord + rb[cb].start();
#else
          int64_t qdp_index =
              ((t * Nz + z) * Ny + y) * Nxh + x_coord + cb * rb[cb].numSiteTable();
#endif

                int64_t qdp_inner_idx = qdp_index % qdp_inner_dim;
                int64_t qdp_outer_idx = qdp_index / qdp_inner_dim;

                WT *psiptr = (WT *)(psi_in.getFjit());

                // ASSUME OCSRI order. So the fixed points our
                // qdp_outer*block-size + qdp_inner
                // and the iteration is over the spinor as Inner x Complex x Spin x
                // Color

                int64_t offset = qdp_inner_idx + outer_block_size * qdp_outer_idx;
                psi[ind][col][spin][0][x] = (FT)
                    psiptr[offset +
                           qdp_inner_dim * (RE + n_complex_dim * (spin + Ns * col))];

                psi[ind][col][spin][1][x] = (FT)
                    psiptr[offset +
                           qdp_inner_dim * (IM + n_complex_dim * (spin + Ns * col))];
              } // x
            } // spin
          } // col
#ifndef MANUAL_COLLAPSE
        } // s
      } // y
    } // z
  } // t
#else
  } // cindex
#endif
}

template <typename FT, int veclen, int soalen, bool compress, typename QDPSpinor>
void qdp_unpack_cb_spinor(
    typename Geometry<FT, veclen, soalen, compress>::FourSpinorBlock *chi_packed,
    QDPSpinor &chi,
    const Geometry<FT, veclen, soalen, compress> &s,
    int cb)
{
  int64_t Nt = s.Nt();
  int64_t Nz = s.Nz();
  int64_t Ny = s.Ny();
  int64_t Nxh = s.Nxh();
  int64_t nvecs = s.nVecs();
  int64_t Pxy = s.getPxy();
  int64_t Pxyz = s.getPxyz();

  int64_t qdp_inner_dim = (int64_t)getDataLayoutInnerSize();
  const int64_t n_complex_dim = 2;
  using WT = typename WordType<QDPSpinor>::Type_t;

  const int64_t outer_block_size = qdp_inner_dim * n_complex_dim * Nc * Ns;

#ifdef DEBUG_PACKER
  QDPIO::cout << "QDP CBSpinorUnPacker: sizeof(FT)=" << sizeof(FT)
              << " sizeof(WT)=" << sizeof(WT) << " qdp_inner=" << qdp_inner_dim
              << " SOA=" << soalen << " veclen=" << veclen << std::endl;
#endif

#ifndef MANUAL_COLLAPSE
#pragma omp parallel for collapse(4)
  for (int64_t t = 0; t < Nt; t++) {
    for (int64_t z = 0; z < Nz; z++) {
      for (int64_t y = 0; y < Ny; y++) {
        for (int64_t s = 0; s < nvecs; s++) {
#else
  // cindex = combined index = s + nvecs*(y + Ny*(z + Nz*t ) )
  // in case collapse(4) doesnt work in OpenMP
  int max_cindex = Nt * Nz * Ny * nvecs;

#pragma omp parallel for
  for (int cindex = 0; cindex < max_cindex; cindex++) {
    int64_t t1 = cindex / nvecs;
    int64_t s = cindex - nvecs * t1;

    int64_t t2 = t1 / Ny;
    int64_t y = t1 - Ny * t2;

    int64_t t = t2 / Nz;
    int64_t z = t2 - Nz * t;
#endif
          for (int64_t spin = 0; spin < 4; spin++) {
            for (int64_t col = 0; col < 3; col++) {
              for (int64_t x = 0; x < soalen; x++) {

                int64_t ind =
                    t * Pxyz + z * Pxy + y * nvecs + s; //((t*Nz+z)*Ny+y)*nvecs+s;
                int64_t x_coord = s * soalen + x;

#if 0
									int64_t qdp_index = ((t*Nz + z)*Ny + y)*Nxh + x_coord + rb[cb].start();
#else
          int64_t qdp_index =
              ((t * Nz + z) * Ny + y) * Nxh + x_coord + cb * rb[cb].numSiteTable();
#endif
                int64_t qdp_inner_idx = qdp_index % qdp_inner_dim;
                int64_t qdp_outer_idx = qdp_index / qdp_inner_dim;

                WT *chiptr = (WT *)(chi.getFjit());
                // ASSUME OCSRI order. So the fixed points our
                // qdp_outer*block-size + qdp_inner
                // and the iteration is over the spinor as Inner x Complex x Spin x
                // Color

                int64_t offset = qdp_inner_idx + outer_block_size * qdp_outer_idx;

                chiptr[offset +
                       qdp_inner_dim * (RE + n_complex_dim * (spin + Ns * col))] =
                    (WT)chi_packed[ind][col][spin][0][x];
                chiptr[offset +
                       qdp_inner_dim * (IM + n_complex_dim * (spin + Ns * col))] =
                    (WT)chi_packed[ind][col][spin][1][x];

              } // x
            } // spin
          } // col
#ifndef MANUAL_COLLAPSE
        } // s
      } // y
    } // z
  } // t
#else
  } // cindex
#endif

#ifdef DEBUG_PACKER
  QDPIO::cout << "QDP CBSpinorUnPacker: Done" << std::endl;
#endif
} // function
};
