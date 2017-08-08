#include "timeTWMClover.h"
#include <omp.h>

#include "qphix/tm_clov_dslash_def.h"
#include "qphix/tm_clov_dslash_body.h"
#include "qphix/twisted_mass_clover.h"
#include "qphix/blas.h"
#include "qphix/invcg.h"
#include "qphix/invbicgstab.h"

#include <cstdlib>

using namespace std;
using namespace QPhiX;

#include "veclen.h"
#include "tolerance.h"
#include "tparam_selector.h"

void TimeTMClover::run()
{
  call(*this, args_.prec, args_.soalen, args_.compress12);
}

template <typename FT, int V, int S, bool compress, typename QdpGauge, typename QdpSpinor>
void TimeTMClover::operator()()
{
    runTest<FT, V, S, compress>();
}

template <typename FT, int V, int S, bool compress>
void TimeTMClover::runTest()
{

  typedef typename Geometry<FT, V, S, compress>::FourSpinorBlock Spinor;
  typedef typename Geometry<FT, V, S, compress>::SU3MatrixBlock Gauge;
  typedef typename Geometry<FT, V, S, compress>::CloverBlock CloverBlock;
  typedef typename Geometry<FT, V, S, compress>::FullCloverBlock FullCloverBlock;

  bool verbose = false;
  // Work out local lattice size
  int subLattSize[4];
  for (int mu = 0; mu < 4; mu++) {
    subLattSize[mu] = args_.nrow_in[mu] / args_.qmp_geometry[mu];
  }

  // Work out the size of checkerboarded X-dimension
  int X1h = args_.nrow_in[0] / 2;
  int Nx = args_.nrow_in[0];
  int Ny = args_.nrow_in[1];
  int Nz = args_.nrow_in[2];
  int Nt = args_.nrow_in[3];

  int lX1h = subLattSize[0] / 2;
  int lY = subLattSize[1];
  int lZ = subLattSize[2];
  int lT = subLattSize[3];

  masterPrintf("Initializing Dslash\n");

  double t_boundary = (double)(1);
  double coeff_s = (double)(1);
  double coeff_t = (double)(1);

  // Create Scalar Dslash Class
  Geometry<FT, V, S, compress> geom(subLattSize,
                                    args_.By,
                                    args_.Bz,
                                    args_.NCores,
                                    args_.Sy,
                                    args_.Sz,
                                    args_.PadXY,
                                    args_.PadXYZ,
                                    args_.MinCt,
                                    true);
  TMClovDslash<FT, V, S, compress> D32(&geom, t_boundary, coeff_s, coeff_t);

  // Allocate data for the gauges
  Gauge *packed_gauge_cb0 = (Gauge *)geom.allocCBGauge();
  Gauge *packed_gauge_cb1 = (Gauge *)geom.allocCBGauge();

  Gauge *u_packed[2];
  u_packed[0] = packed_gauge_cb0;
  u_packed[1] = packed_gauge_cb1;

  double factor = 0.08;
  masterPrintf("Initializing Fake Gauge Field: ");

  int nvecs = geom.nVecs();
  int nyg = geom.nGY();
  int Pxy = geom.getPxy();
  int Pxyz = geom.getPxyz();

  double start = omp_get_wtime();

#pragma omp parallel for collapse(4)
  for (int t = 0; t < lT; t++) {
    for (int z = 0; z < lZ; z++) {
      for (int y = 0; y < lY; y++) {
        for (int s = 0; s < nvecs; s++) {
          for (int mu = 0; mu < 8; mu++) {
            for (int c = 0; c < (compress ? 2 : 3); c++) {
              for (int c2 = 0; c2 < 3; c2++) {
                for (int x = 0; x < S; x++) {

                  int block = (t * Pxyz + z * Pxy) / nyg + (y / nyg) * nvecs + s;

                  // This will work out to be between 0 and
                  // veclen
                  int xx = (y % nyg) * S + x;

                  double d1 = factor * (drand48() - 0.5);
                  double d2 = factor * (drand48() - 0.5);
                  double d3 = factor * (drand48() - 0.5);
                  double d4 = factor * (drand48() - 0.5);

                  if (c == c2) {
                    u_packed[0][block][mu][c][c2][RE][xx] =
                        rep<FT, double>((double)1 + d1);
                    u_packed[1][block][mu][c][c2][RE][xx] =
                        rep<FT, double>((double)1 + d3);
                  } else {
                    u_packed[0][block][mu][c][c2][RE][xx] = rep<FT, double>(d1);
                    u_packed[1][block][mu][c][c2][RE][xx] = rep<FT, double>(d3);
                  }

                  u_packed[0][block][mu][c][c2][IM][xx] = rep<FT, double>(d2);
                  u_packed[1][block][mu][c][c2][IM][xx] = rep<FT, double>(d4);
                }
              }
            } // row
          }
        }
      }
    }
  }

  if (!compress) {
#pragma omp parallel for collapse(4)
    for (int t = 0; t < lT; t++) {
      for (int z = 0; z < lZ; z++) {
        for (int y = 0; y < lY; y++) {
          for (int s = 0; s < nvecs; s++) {

            int block = (t * Pxyz + z * Pxy) / nyg + (y / nyg) * nvecs + s;

            for (int mu = 0; mu < 8; mu++) {
              for (int row = 0; row < 2; row++) {

                double norm_row_cb0[V];
                double norm_row_cb1[V];

                for (int x = 0; x < V; x++) {
                  norm_row_cb0[x] = 0;
                  norm_row_cb1[x] = 0;
                }

                // This will work out to be between 0 and veclen
                // Accumulate the norms
                for (int col = 0; col < 3; col++) {
                  for (int x = 0; x < S; x++) {
                    int xx = (y % nyg) * S + x;

                    double u0_re =
                        rep<double, FT>(u_packed[0][block][mu][row][col][RE][xx]);
                    double u1_re =
                        rep<double, FT>(u_packed[1][block][mu][row][col][RE][xx]);

                    double u0_im =
                        rep<double, FT>(u_packed[0][block][mu][row][col][IM][xx]);
                    double u1_im =
                        rep<double, FT>(u_packed[1][block][mu][row][col][IM][xx]);

                    norm_row_cb0[xx] += (u0_re * u0_re) + (u0_im * u0_im);

                    norm_row_cb1[xx] += (u1_re * u1_re) + (u1_im * u1_im);

                  } // x
                } // col

                for (int x = 0; x < V; x++) {
                  norm_row_cb0[x] = sqrt(norm_row_cb0[x]);
                  norm_row_cb1[x] = sqrt(norm_row_cb1[x]);
                }

                // Normalize each component.

                for (int col = 0; col < 3; col++) {
                  for (int x = 0; x < S; x++) {
                    int xx = (y % nyg) * S + x;

                    double u0_re =
                        rep<double, FT>(u_packed[0][block][mu][row][col][RE][xx]) /
                        norm_row_cb0[xx];
                    double u1_re =
                        rep<double, FT>(u_packed[1][block][mu][row][col][RE][xx]) /
                        norm_row_cb1[xx];

                    double u0_im =
                        rep<double, FT>(u_packed[0][block][mu][row][col][IM][xx]) /
                        norm_row_cb0[xx];
                    double u1_im =
                        rep<double, FT>(u_packed[1][block][mu][row][col][IM][xx]) /
                        norm_row_cb1[xx];

                    u_packed[0][block][mu][row][col][RE][xx] =
                        rep<FT, double>(u0_re);
                    u_packed[0][block][mu][row][col][IM][xx] =
                        rep<FT, double>(u0_im);

                    u_packed[1][block][mu][row][col][RE][xx] =
                        rep<FT, double>(u1_re);
                    u_packed[1][block][mu][row][col][IM][xx] =
                        rep<FT, double>(u1_im);
                  } // x
                } // col
              } // row

              {
                for (int x = 0; x < S; x++) {
                  // 3rd row reconstruction.
                  int xx = (y % nyg) * S + x;
                  double ar = rep<double, FT>(u_packed[0][block][mu][0][0][RE][xx]);
                  double ai = rep<double, FT>(u_packed[0][block][mu][0][0][IM][xx]);

                  double br = rep<double, FT>(u_packed[0][block][mu][0][1][RE][xx]);
                  double bi = rep<double, FT>(u_packed[0][block][mu][0][1][IM][xx]);

                  double cr = rep<double, FT>(u_packed[0][block][mu][0][2][RE][xx]);
                  double ci = rep<double, FT>(u_packed[0][block][mu][0][2][IM][xx]);

                  double dr = rep<double, FT>(u_packed[0][block][mu][1][0][RE][xx]);
                  double di = rep<double, FT>(u_packed[0][block][mu][1][0][IM][xx]);

                  double er = rep<double, FT>(u_packed[0][block][mu][1][1][RE][xx]);
                  double ei = rep<double, FT>(u_packed[0][block][mu][1][1][IM][xx]);

                  double fr = rep<double, FT>(u_packed[0][block][mu][1][2][RE][xx]);
                  double fi = rep<double, FT>(u_packed[0][block][mu][1][2][IM][xx]);

                  u_packed[0][block][mu][2][0][RE][xx] =
                      rep<FT, double>(br * fr - bi * fi - er * cr + ei * ci);
                  u_packed[0][block][mu][2][0][IM][xx] =
                      rep<FT, double>(er * ci + ei * cr - br * fi - bi * fr);
                  u_packed[0][block][mu][2][1][RE][xx] =
                      rep<FT, double>(dr * cr - di * ci - ar * fr + ai * fi);
                  u_packed[0][block][mu][2][1][IM][xx] =
                      rep<FT, double>(ar * fi + ai * fr - dr * ci - di * cr);
                  u_packed[0][block][mu][2][2][RE][xx] =
                      rep<FT, double>(ar * er - ai * ei - dr * br + di * bi);
                  u_packed[0][block][mu][2][2][IM][xx] =
                      rep<FT, double>(dr * bi + di * br - ar * ei - ai * er);
                }
              }

              {
                for (int x = 0; x < S; x++) {
                  int xx = (y % nyg) * S + x;
                  // 3rd row reconstruction.
                  double ar = rep<double, FT>(u_packed[1][block][mu][0][0][RE][xx]);
                  double ai = rep<double, FT>(u_packed[1][block][mu][0][0][IM][xx]);

                  double br = rep<double, FT>(u_packed[1][block][mu][0][1][RE][xx]);
                  double bi = rep<double, FT>(u_packed[1][block][mu][0][1][IM][xx]);

                  double cr = rep<double, FT>(u_packed[1][block][mu][0][2][RE][xx]);
                  double ci = rep<double, FT>(u_packed[1][block][mu][0][2][IM][xx]);

                  double dr = rep<double, FT>(u_packed[1][block][mu][1][0][RE][xx]);
                  double di = rep<double, FT>(u_packed[1][block][mu][1][0][IM][xx]);

                  double er = rep<double, FT>(u_packed[1][block][mu][1][1][RE][xx]);
                  double ei = rep<double, FT>(u_packed[1][block][mu][1][1][IM][xx]);

                  double fr = rep<double, FT>(u_packed[1][block][mu][1][2][RE][xx]);
                  double fi = rep<double, FT>(u_packed[1][block][mu][1][2][IM][xx]);

                  u_packed[1][block][mu][2][0][RE][xx] =
                      rep<FT, double>(br * fr - bi * fi - er * cr + ei * ci);
                  u_packed[1][block][mu][2][0][IM][xx] =
                      rep<FT, double>(er * ci + ei * cr - br * fi - bi * fr);
                  u_packed[1][block][mu][2][1][RE][xx] =
                      rep<FT, double>(dr * cr - di * ci - ar * fr + ai * fi);
                  u_packed[1][block][mu][2][1][IM][xx] =
                      rep<FT, double>(ar * fi + ai * fr - dr * ci - di * cr);
                  u_packed[1][block][mu][2][2][RE][xx] =
                      rep<FT, double>(ar * er - ai * ei - dr * br + di * bi);
                  u_packed[1][block][mu][2][2][IM][xx] =
                      rep<FT, double>(dr * bi + di * br - ar * ei - ai * er);
                } // x
              }

            } // mu
          } // s
        } // y
      } // z
    } // t

  } // end if ! compress

  double end = omp_get_wtime();
  masterPrintf(" %g sec\n", end - start);

  // Allocate data for the spinors
  Spinor *p_even = (Spinor *)geom.allocCBFourSpinor();
  Spinor *p_odd = (Spinor *)geom.allocCBFourSpinor();
  Spinor *c_even = (Spinor *)geom.allocCBFourSpinor();
  Spinor *c_odd = (Spinor *)geom.allocCBFourSpinor();

  // Point to the second block of the array. Now there is padding on both
  // ends.
  Spinor *psi_s[2] = {p_even, p_odd};
  Spinor *chi_s[2] = {c_even, c_odd};

  masterPrintf("Filling Input spinor: ");

  start = omp_get_wtime();

#pragma omp parallel for collapse(4)
  for (int t = 0; t < lT; t++) {
    for (int z = 0; z < lZ; z++) {
      for (int y = 0; y < lY; y++) {
        for (int s = 0; s < nvecs; s++) {
          for (int spin = 0; spin < 4; spin++) {
            for (int col = 0; col < 3; col++) {
              for (int x = 0; x < S; x++) {
                double d1 = drand48() - 0.5;
                double d2 = drand48() - 0.5;
                double d3 = drand48() - 0.5;
                double d4 = drand48() - 0.5;

                int ind =
                    t * Pxyz + z * Pxy + y * nvecs + s; //((t*Nz+z)*Ny+y)*nvecs+s;
                psi_s[0][ind][col][spin][0][x] = rep<FT, double>(d1);
                psi_s[0][ind][col][spin][1][x] = rep<FT, double>(d2);
                psi_s[1][ind][col][spin][0][x] = rep<FT, double>(d3);
                psi_s[1][ind][col][spin][1][x] = rep<FT, double>(d4);
              }
            }
          }
        }
      }
    }
  }

  end = omp_get_wtime();

  masterPrintf(" %g sec\n", end - start);
  masterPrintf("Zeroing output spinor: ");

  start = omp_get_wtime();

#pragma omp parallel for collapse(4)
  for (int t = 0; t < lT; t++) {
    for (int z = 0; z < lZ; z++) {
      for (int y = 0; y < lY; y++) {
        for (int s = 0; s < nvecs; s++) {
          for (int spin = 0; spin < 4; spin++) {
            for (int col = 0; col < 3; col++) {
              for (int x = 0; x < S; x++) {

                int ind =
                    t * Pxyz + z * Pxy + y * nvecs + s; //((t*Nz+z)*Ny+y)*nvecs+s;
                double d = 0;
                chi_s[0][ind][col][spin][0][x] = rep<FT, double>(d);
                chi_s[0][ind][col][spin][1][x] = rep<FT, double>(d);
                chi_s[1][ind][col][spin][0][x] = rep<FT, double>(d);
                chi_s[1][ind][col][spin][1][x] = rep<FT, double>(d);
              }
            }
          }
        }
      }
    }
  }
  end = omp_get_wtime();

  masterPrintf(" %g sec\n", end - start);

  masterPrintf("Allocating TWO Clover Terms and TWO Clover Inverses\n");
  CloverBlock *A_cb0_plus = geom.allocCBClov();
  CloverBlock *A_cb0_minus = geom.allocCBClov();
  FullCloverBlock *A_inv_cb1_plus = geom.allocCBFullClov();
  FullCloverBlock *A_inv_cb1_minus = geom.allocCBFullClov();

  CloverBlock *A_cb0[2] = {A_cb0_plus, A_cb0_minus};
  FullCloverBlock *A_inv_cb1[2] = {A_inv_cb1_plus, A_inv_cb1_minus};

  FullCloverBlock *A_inv_multi_cb0[2] = {A_inv_cb1[0], A_inv_cb1[1]};
  FullCloverBlock *A_inv_multi_cb1[2] = {A_inv_cb1[0], A_inv_cb1[1]};
  FullCloverBlock **A_inv_multi[2] = {A_inv_multi_cb0, A_inv_multi_cb1};

  masterPrintf("Done\n");

  masterPrintf("Filling Clover Terms with junk\n");

#pragma omp parallel for collapse(4)
  for (int t = 0; t < lT; t++) {
    for (int z = 0; z < lZ; z++) {
      for (int y = 0; y < lY; y++) {
        for (int s = 0; s < nvecs; s++) {
          for (int x = 0; x < S; x++) {

            int block = (t * Pxyz + z * Pxy) / nyg + (y / nyg) * nvecs + s;

            // This will work out to be between 0 and veclen
            int xx = (y % nyg) * S + x;

            for (int pm : {0, 1}) {
              for (int i = 0; i < 6; i++) {
                A_cb0[pm][block].diag1[i][xx] = rep<FT, double>((double)4.1);
                A_cb0[pm][block].diag2[i][xx] = rep<FT, double>((double)4.1);
              }

              for (int i = 0; i < 15; i++) {
                for (int reim = 0; reim < 2; reim++) {
                  A_cb0[pm][block].off_diag1[i][reim][xx] = rep<FT, double>(0);
                  A_cb0[pm][block].off_diag2[i][reim][xx] = rep<FT, double>(0);
                  ;
                }
              }
            }

            // Zero out both clover terms & inverses for both blocks
            // respectively
            for (int pm = 0; pm < 2; ++pm) {
              for (int i = 0; i < 6; ++i) {
                for (int j = 0; j < 6; ++j) {
                  for (int reim = 0; reim < 2; ++reim) {

                    A_inv_cb1[pm][block].block1[i][j][reim][xx] =
                        rep<FT, double>((double)0);
                    A_inv_cb1[pm][block].block2[i][j][reim][xx] =
                        rep<FT, double>((double)0);

                  } // Real / imaginary part
                } // rows
              } // columns
            } // plus / minus

            // Set the diagonal and its inverse
            for (int pm = 0; pm < 2; ++pm) {
              for (int i = 0; i < 6; ++i) {

                A_inv_cb1[pm][block].block1[i][i][0][xx] =
                    rep<FT, double>((double)1 / double(4.1));
                A_inv_cb1[pm][block].block1[i][i][1][xx] =
                    rep<FT, double>((double)0);
                A_inv_cb1[pm][block].block2[i][i][0][xx] =
                    rep<FT, double>((double)1 / double(4.1));
                A_inv_cb1[pm][block].block2[i][i][1][xx] =
                    rep<FT, double>((double)0);

              } // rows & columns
            } // plus / minus
          }
        }
      }
    }
  } // Filling clover terms

  masterPrintf("Done\n\n");
  masterPrintf("Starting Timings\n");

  masterPrintf("Creating Even/Odd Twisted Mass Clover Operator\n");
  EvenOddTMCloverOperator<FT, V, S, compress> M(
      u_packed, A_cb0[0], A_inv_cb1, &geom, t_boundary, coeff_s, coeff_t);

  if (args_.do_dslash) {

    for (int isign = 1; isign >= -1; isign -= 2) {
      for (int cb = 0; cb <= 1; ++cb) {

        int source_cb = 1 - cb;
        int target_cb = cb;

        masterPrintf("Timing on cb = %d isign = %d\n", cb, isign);
        masterPrintf("=================================\n");

        for (int repeat = 0; repeat < 3; repeat++) {

          double start = omp_get_wtime();

          for (int i = 0; i < args_.iters; i++) {
            // Apply Optimized Dslash
            D32.dslash(chi_s[target_cb],
                       psi_s[source_cb],
                       u_packed[target_cb],
                       (const FullCloverBlock **)A_inv_cb1,
                       isign,
                       target_cb);
          }

          double end = omp_get_wtime();
          double time = end - start;
          CommsUtils::sumDouble(&time);
          time /= (double)CommsUtils::numNodes();

          masterPrintf("\t timing %d of 3\n", repeat + 1);
          masterPrintf("\t %d iterations in %e seconds\n", args_.iters, time);
          masterPrintf("\t %e usec/iteration\n", 1.0e6 * time / (double)args_.iters);
          double Gflops =
              1824.0f * (double)(args_.iters) * (double)(X1h * Ny * Nz * Nt) / 1.0e9;
          double perf = Gflops / time;
          masterPrintf("\t Performance: %g GFLOPS total\n", perf);
          masterPrintf("\t              %g GFLOPS / node\n",
                       perf / (double)CommsUtils::numNodes());

        } // repeats

      } // isign
    } // cb

  } // do dslash

  if (args_.do_m) {

    for (int isign = 1; isign >= -1; isign -= 2) {

      masterPrintf("Timing M: isign = %d\n", isign);
      masterPrintf("===============================\n");

      for (int repeat = 0; repeat < 3; repeat++) {

        double start = omp_get_wtime();

        for (int i = 0; i < args_.iters; i++) {
          // Apply Full Fermion Matrix (w/ Optimized Dslash)
          M(chi_s[0], psi_s[0], isign);
        }

        double end = omp_get_wtime();
        double time = end - start;

        CommsUtils::sumDouble(&time);
        time /= (double)CommsUtils::numNodes();

        masterPrintf("\t timing %d of 3\n", repeat + 1);
        masterPrintf("\t %d iterations in %e seconds\n", args_.iters, time);
        masterPrintf("\t %e usec/iteration\n", 1.0e6 * time / (double)args_.iters);
        double flops_per_iter = 3696.0f;
        double Gflops =
            flops_per_iter * (double)(args_.iters) * (double)(X1h * Ny * Nz * Nt) / 1.0e9;
        double perf = Gflops / time;
        masterPrintf("\t Performance: %g GFLOPS total\n", perf);
        masterPrintf("\t              %g GFLOPS / node\n",
                     perf / (double)CommsUtils::numNodes());
      }
    }

  } // Full Fermion Matrix M

  double rsd_target = rsdTarget<FT>::value;
  int max_iters = 250;
  int niters;
  double rsd_final;

  if (args_.do_cg) {

    InvCG<FT, V, S, compress> solver(M, max_iters);

    for (int solves = 0; solves < 5; solves++) {

      masterPrintf("Zeroing solution\n");
      masterPrintf("Starting solver\n");

      unsigned long site_flops;
      unsigned long mv_apps;

#pragma omp parallel for collapse(4)
      for (int t = 0; t < lT; t++) {
        for (int z = 0; z < lZ; z++) {
          for (int y = 0; y < lY; y++) {
            for (int s = 0; s < nvecs; s++) {
              for (int spin = 0; spin < 4; spin++) {
                for (int col = 0; col < 3; col++) {
                  for (int x = 0; x < S; x++) {

                    int ind = t * Pxyz + z * Pxy + y * nvecs +
                              s; //((t*Nz+z)*Ny+y)*nvecs+s;
                    int x_coord = s * S + x;
                    double d1 = drand48() - 0.5;
                    double d2 = drand48() - 0.5;
                    double d3 = drand48() - 0.5;
                    double d4 = drand48() - 0.5;
                    psi_s[0][ind][col][spin][0][x] = rep<FT, double>(d1);
                    psi_s[0][ind][col][spin][1][x] = rep<FT, double>(d2);
                    psi_s[1][ind][col][spin][0][x] = rep<FT, double>(d3);
                    psi_s[1][ind][col][spin][1][x] = rep<FT, double>(d4);
                    chi_s[0][ind][col][spin][0][x] = rep<FT, double>(0);
                    chi_s[0][ind][col][spin][1][x] = rep<FT, double>(0);
                    chi_s[1][ind][col][spin][0][x] = rep<FT, double>(0);
                    chi_s[1][ind][col][spin][1][x] = rep<FT, double>(0);
                  }
                }
              }
            }
          }
        }
      }

      int isign = 1;
      start = omp_get_wtime();
      solver(chi_s[0],
             psi_s[0],
             rsd_target,
             niters,
             rsd_final,
             site_flops,
             mv_apps,
             isign,
             verbose);
      end = omp_get_wtime();

      unsigned long num_cb_sites = X1h * Ny * Nz * Nt;
      unsigned long total_flops = (site_flops + (3696) * mv_apps) * num_cb_sites;

      masterPrintf("Solver iterations = %d\n", niters);
      masterPrintf("Solver time = %g(s)\n", (end - start));
      masterPrintf("CG GFLOPS = %g\n",
                   1.0e-9 * (double)(total_flops) / (end - start));
    }

  } // CG Solver

  if (args_.do_bicgstab) {

    InvBiCGStab<FT, V, S, compress> solver2(M, max_iters);

    for (int solver = 0; solver < 5; solver++) {

      unsigned long site_flops;
      unsigned long mv_apps;

#pragma omp parallel for collapse(4)
      for (int t = 0; t < lT; t++) {
        for (int z = 0; z < lZ; z++) {
          for (int y = 0; y < lY; y++) {
            for (int s = 0; s < nvecs; s++) {
              for (int spin = 0; spin < 4; spin++) {
                for (int col = 0; col < 3; col++) {
                  for (int x = 0; x < S; x++) {

                    int ind = t * Pxyz + z * Pxy + y * nvecs +
                              s; //((t*Nz+z)*Ny+y)*nvecs+s;
                    int x_coord = s * S + x;

                    double d1 = drand48() - 0.5;
                    double d2 = drand48() - 0.5;
                    double d3 = drand48() - 0.5;
                    double d4 = drand48() - 0.5;

                    psi_s[0][ind][col][spin][0][x] = rep<FT, double>(d1);
                    psi_s[0][ind][col][spin][1][x] = rep<FT, double>(d2);
                    psi_s[1][ind][col][spin][0][x] = rep<FT, double>(d3);
                    psi_s[1][ind][col][spin][1][x] = rep<FT, double>(d4);

                    chi_s[0][ind][col][spin][0][x] = rep<FT, double>(0);
                    chi_s[0][ind][col][spin][1][x] = rep<FT, double>(0);
                    chi_s[1][ind][col][spin][0][x] = rep<FT, double>(0);
                    chi_s[1][ind][col][spin][1][x] = rep<FT, double>(0);
                  }
                }
              }
            }
          }
        }
      }

      int isign = 1;

      start = omp_get_wtime();
      solver2(chi_s[0],
              psi_s[0],
              rsd_target,
              niters,
              rsd_final,
              site_flops,
              mv_apps,
              isign,
              verbose);
      end = omp_get_wtime();

      unsigned long num_cb_sites = X1h * Ny * Nz * Nt;
      unsigned long total_flops = (site_flops + (3696) * mv_apps) * num_cb_sites;

      masterPrintf("Solver iterations = %d\n", niters);
      masterPrintf("Solver time = %g(s)\n", (end - start));
      masterPrintf("BICGSTAB GFLOPS = %g\n",
                   1.0e-9 * (double)(total_flops) / (end - start));
    }

  } // BiCGSTAB Solver

  masterPrintf("Cleaning up\n");

  geom.free(packed_gauge_cb0);
  geom.free(packed_gauge_cb1);
  geom.free(p_even);
  geom.free(p_odd);
  geom.free(c_even);
  geom.free(c_odd);
  geom.free(A_cb0_plus);
  geom.free(A_cb0_minus);
  geom.free(A_inv_cb1_plus);
  geom.free(A_inv_cb1_minus);
}
