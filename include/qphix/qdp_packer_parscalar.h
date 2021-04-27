#pragma once

#include "qphix/diagnostics.h"

QPHIX_MESSAGE("building with parscalar packers")

#include "qdp.h"
#include "qphix/geometry.h"

#include "qphix/dslash_def.h"
#include "qphix/qphix_config.h"

#if defined(QPHIX_MIC_SOURCE) || defined(QPHIX_AVX512_SOURCE)
#include <immintrin.h>
#endif

#ifdef QPHIX_BUILD_CLOVER
#include "qphix/clover_dslash_def.h"
#endif

#ifdef QPHIX_BUILD_TWISTED_MASS_WITH_CLOVER
#include "qphix/tm_clov_dslash_def.h"
#endif

#define MANUAL_COLLAPSE 1

using namespace QDP;

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
  int Nt = s.Nt();
  int Nz = s.Nz();
  int Ny = s.Ny();
  int nvecs = s.nVecs();
  int nyg = s.nGY();
  int Pxy = s.getPxy();
  int Pxyz = s.getPxyz();

  // Shift the lattice to get U(x-mu)
  QDPGauge u_minus(4);
  for (int mu = 0; mu < 4; mu++) {
    u_minus[mu] = shift(u[mu], BACKWARD, mu);
  }

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

          for (int mu = 0; mu < 4; mu++) {
            int outer_c = 3;
            if (compress) {
              outer_c = 2;
            }
            for (int c = 0; c < outer_c; c++) {
              for (int c2 = 0; c2 < 3; c2++) {
                for (int x = 0; x < soalen; x++) {

                  //#ifndef USE_PACKED_GAUGES
                  // int xx = x;
                  // int block = ((t*Nz+z)*Ny+y)*nvecs+s;

                  //#endif
                  //#else // USE_PACKED_GAUGES
                  int block = (t * Pxyz + z * Pxy) / nyg + (y / nyg) * nvecs + s;
                  int xx = (y % nyg) * soalen + x;
                  // #endif // USE_PACKED_GAUGES

                  int qdpsite = x + soalen * (s + nvecs * (y + Ny * (z + Nz * t)));
                  u_cb0[block][2 * mu][c][c2][0][xx] =
                      u_minus[mu]
                          .elem(rb[0].start() + qdpsite)
                          .elem()
                          .elem(c2, c)
                          .real();
                  u_cb0[block][2 * mu][c][c2][1][xx] =
                      u_minus[mu]
                          .elem(rb[0].start() + qdpsite)
                          .elem()
                          .elem(c2, c)
                          .imag();
                  u_cb0[block][2 * mu + 1][c][c2][0][xx] =
                      u[mu].elem(rb[0].start() + qdpsite).elem().elem(c2, c).real();
                  u_cb0[block][2 * mu + 1][c][c2][1][xx] =
                      u[mu].elem(rb[0].start() + qdpsite).elem().elem(c2, c).imag();

                  u_cb1[block][2 * mu][c][c2][0][xx] =
                      u_minus[mu]
                          .elem(rb[1].start() + qdpsite)
                          .elem()
                          .elem(c2, c)
                          .real();
                  u_cb1[block][2 * mu][c][c2][1][xx] =
                      u_minus[mu]
                          .elem(rb[1].start() + qdpsite)
                          .elem()
                          .elem(c2, c)
                          .imag();
                  u_cb1[block][2 * mu + 1][c][c2][0][xx] =
                      u[mu].elem(rb[1].start() + qdpsite).elem().elem(c2, c).real();
                  u_cb1[block][2 * mu + 1][c][c2][1][xx] =
                      u[mu].elem(rb[1].start() + qdpsite).elem().elem(c2, c).imag();
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
  } // cindex
#endif

} // end of function

#ifdef QPHIX_BUILD_CLOVER
template <typename FT, int veclen, int soalen, bool compress, typename ClovTerm>
void qdp_pack_clover(
    const ClovTerm &qdp_clov_in,
    typename ClovDslash<FT, veclen, soalen, compress>::CloverBlock *cl_out,
    const Geometry<FT, veclen, soalen, compress> &s,
    int cb)
{
  // Get the subgrid latt size.
  int Nt = s.Nt();
  int Nz = s.Nz();
  int Ny = s.Ny();
  int nvecs = s.nVecs();
  int nyg = s.nGY();
  int Pxy = s.getPxy();
  int Pxyz = s.getPxyz();
  const auto &qdp_clov_in_buf = qdp_clov_in.getTriBuffer();

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

            int block = (t * Pxyz + z * Pxy) / nyg + (y / nyg) * nvecs + s;
            int xx = (y % nyg) * soalen + x;
            int qdpsite =
                x + soalen * (s + nvecs * (y + Ny * (z + Nz * t))) + rb[cb].start();

            for (int d = 0; d < 6; d++) {
              cl_out[block].diag1[d][xx] =
                  qdp_clov_in_buf[qdpsite].diag[0][d].elem();
            }
            for (int od = 0; od < 15; od++) {
              cl_out[block].off_diag1[od][RE][xx] =
                  qdp_clov_in_buf[qdpsite].offd[0][od].real();
              cl_out[block].off_diag1[od][IM][xx] =
                  qdp_clov_in_buf[qdpsite].offd[0][od].imag();
            }

            for (int d = 0; d < 6; d++) {
              cl_out[block].diag2[d][xx] =
                  qdp_clov_in_buf[qdpsite].diag[1][d].elem();
            }
            for (int od = 0; od < 15; od++) {
              cl_out[block].off_diag2[od][RE][xx] =
                  qdp_clov_in_buf[qdpsite].offd[1][od].real();
              cl_out[block].off_diag2[od][IM][xx] =
                  qdp_clov_in_buf[qdpsite].offd[1][od].imag();
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
}
#endif // IFDEF BUILD CLOVER

#ifdef QPHIX_BUILD_TWISTED_MASS_WITH_CLOVER
template <typename FT, int veclen, int soalen, bool compress, typename ClovTerm>
void qdp_pack_full_clover(
    const ClovTerm &qdp_clov_in,
    typename TMClovDslash<FT, veclen, soalen, compress>::FullCloverBlock *cl_out,
    const Geometry<FT, veclen, soalen, compress> &s,
    int cb)
{
  int index[2][15] = {{1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 5},
                      {0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4}};

  // Get the subgrid latt size.
  int Nt = s.Nt();
  int Nz = s.Nz();
  int Ny = s.Ny();
  int nvecs = s.nVecs();
  int nyg = s.nGY();
  int Pxy = s.getPxy();
  int Pxyz = s.getPxyz();
  const auto &qdp_clov_in_buf = qdp_clov_in.getTriBuffer();

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

            int block = (t * Pxyz + z * Pxy) / nyg + (y / nyg) * nvecs + s;
            int xx = (y % nyg) * soalen + x;
            int qdpsite =
                x + soalen * (s + nvecs * (y + Ny * (z + Nz * t))) + rb[cb].start();

            for (int ii = 0; ii < 6; ii++) { // block1, diagonal
              cl_out[block].block1[ii][ii][RE][xx] =
                  qdp_clov_in_buf[qdpsite].diag[0][ii].elem();
              cl_out[block].block1[ii][ii][IM][xx] = 0;
            }
            for (int od = 0; od < 15; od++) { // block1, off-diagonal
              int ii = index[0][od];
              int jj = index[1][od];
              cl_out[block].block1[ii][jj][RE][xx] =
                  qdp_clov_in_buf[qdpsite].offd[0][od].real();
              cl_out[block].block1[ii][jj][IM][xx] =
                  qdp_clov_in_buf[qdpsite].offd[0][od].imag();
              cl_out[block].block1[jj][ii][RE][xx] =
                  qdp_clov_in_buf[qdpsite].offd[0][od].real();
              cl_out[block].block1[jj][ii][IM][xx] =
                  -qdp_clov_in_buf[qdpsite].offd[0][od].imag();
            }

            for (int ii = 0; ii < 6; ii++) { // block2, diagonal
              cl_out[block].block2[ii][ii][RE][xx] =
                  qdp_clov_in_buf[qdpsite].diag[1][ii].elem();
              cl_out[block].block2[ii][ii][IM][xx] = 0;
            }
            for (int od = 0; od < 15; od++) { // block2, off-diagonal
              int ii = index[0][od];
              int jj = index[1][od];
              cl_out[block].block2[ii][jj][RE][xx] =
                  qdp_clov_in_buf[qdpsite].offd[1][od].real();
              cl_out[block].block2[ii][jj][IM][xx] =
                  qdp_clov_in_buf[qdpsite].offd[1][od].imag();
              cl_out[block].block2[jj][ii][RE][xx] =
                  qdp_clov_in_buf[qdpsite].offd[1][od].real();
              cl_out[block].block2[jj][ii][IM][xx] =
                  -qdp_clov_in_buf[qdpsite].offd[1][od].imag();
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
}
#endif // IFDEF BUILD TWISTED-MASS CLOVER

template <typename FT, int veclen, int soalen, bool compress, typename QDPSpinor>
void qdp_pack_cb_spinor(
    const QDPSpinor &psi_in,
    typename Geometry<FT, veclen, soalen, compress>::FourSpinorBlock *psi,
    const Geometry<FT, veclen, soalen, compress> &s,
    int cb)
{
  // Get the subgrid latt size.
  int Nt = s.Nt();
  int Nz = s.Nz();
  int Ny = s.Ny();
  int Nxh = s.Nxh();
  int nvecs = s.nVecs();
  int Pxy = s.getPxy();
  int Pxyz = s.getPxyz();

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

          for (int col = 0; col < 3; col++) {
            for (int spin = 0; spin < 4; spin++) {
              for (int x = 0; x < soalen; x++) {

                int ind =
                    t * Pxyz + z * Pxy + y * nvecs + s; //((t*Nz+z)*Ny+y)*nvecs+s;
                int x_coord = s * soalen + x;

                int qdp_ind = ((t * Nz + z) * Ny + y) * Nxh + x_coord;
                psi[ind][col][spin][0][x] = psi_in.elem(rb[cb].start() + qdp_ind)
                                                .elem(spin)
                                                .elem(col)
                                                .real();
                psi[ind][col][spin][1][x] = psi_in.elem(rb[cb].start() + qdp_ind)
                                                .elem(spin)
                                                .elem(col)
                                                .imag();
              }
            }
          }
#ifndef MANUAL_COLLAPSE
        } // s
      } // y
    } // z
  } // t
#else
  }
#endif
}

template <typename FT, int veclen, int soalen, bool compress, typename QDPSpinor>
void qdp_unpack_cb_spinor(
    typename Geometry<FT, veclen, soalen, compress>::FourSpinorBlock *chi_packed,
    QDPSpinor &chi,
    const Geometry<FT, veclen, soalen, compress> &s,
    int cb)
{
  int Nt = s.Nt();
  int Nz = s.Nz();
  int Ny = s.Ny();
  int Nxh = s.Nxh();
  int nvecs = s.nVecs();
  int Pxy = s.getPxy();
  int Pxyz = s.getPxyz();

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

          for (int spin = 0; spin < 4; spin++) {
            for (int col = 0; col < 3; col++) {
              for (int x = 0; x < soalen; x++) {

                int ind =
                    t * Pxyz + z * Pxy + y * nvecs + s; //((t*Nz+z)*Ny+y)*nvecs+s;
                int x_coord = s * soalen + x;

                int qdp_ind = ((t * Nz + z) * Ny + y) * Nxh + x_coord;

                chi.elem(rb[cb].start() + qdp_ind).elem(spin).elem(col).real() =
                    chi_packed[ind][col][spin][0][x];
                chi.elem(rb[cb].start() + qdp_ind).elem(spin).elem(col).imag() =
                    chi_packed[ind][col][spin][1][x];
              }
            }
          }
#ifndef MANUAL_COLLAPSE
        }
      }
    }
  }
#else
  }
#endif
}
};
