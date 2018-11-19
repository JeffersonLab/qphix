
#include "timeTWMNoQDP.h"

#include <omp.h>
#include "qphix/twisted_mass.h"
#if 1
#include "qphix/blas.h"
#include "qphix/invcg.h"
#include "qphix/invbicgstab.h"
#include "qphix/print_utils.h"
#endif

#include <cstdlib>

using namespace std;
using namespace QPhiX;

#include "veclen.h"
#include "tolerance.h"
#include "tparam_selector.h"

void TimeTMDslash::run()
{
  call(*this, args_.prec, args_.soalen, args_.compress12);
}

template <typename FT, int V, int S, bool compress, typename QdpGauge, typename QdpSpinor>
void TimeTMDslash::operator()()
{
    runTest<FT, V, S, compress>();
}


template <typename FT, int V, int S, bool compress>
void TimeTMDslash::runTest()
{

  typedef typename Geometry<FT, V, S, compress>::SU3MatrixBlock Gauge;
  typedef typename Geometry<FT, V, S, compress>::FourSpinorBlock Spinor;

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

  // Diagnostic information:
  masterPrintf("VECLEN=%d SOALEN=%d\n", V, S);
  masterPrintf("Global Lattice Size = ");
  for (int mu = 0; mu < 4; mu++) {
    masterPrintf(" %d", args_.nrow_in[mu]);
  }
  masterPrintf("\n");

  masterPrintf("Local Lattice Size = ");
  for (int mu = 0; mu < 4; mu++) {
    masterPrintf(" %d", subLattSize[mu]);
  }
  masterPrintf("\n");

  masterPrintf("Block Sizes: args_.By= %d args_.Bz=%d\n", args_.By, args_.Bz);
  masterPrintf("Cores = %d\n", args_.NCores);
  masterPrintf("SMT Grid: args_.Sy=%d args_.Sz=%d\n", args_.Sy, args_.Sz);
  masterPrintf("Pad Factors: args_.PadXY=%d args_.PadXYZ=%d\n", args_.PadXY, args_.PadXYZ);

  masterPrintf("Initializing TMDslash\n");

  const double t_boundary = 1.0;
  const double coeff_s = 1.0;
  const double coeff_t = 1.0;

  const double Mass = 0.1;
  const double TwistedMass = 0.1;

  const double alpha = 4.0 + Mass;
  const double beta = 0.25;
  const double Mu = TwistedMass / alpha;
  const double MuInv = alpha / (alpha * alpha + TwistedMass * TwistedMass);

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
  TMDslash<FT, V, S, compress> D32(
      &geom, t_boundary, coeff_s, coeff_t, Mass, TwistedMass);

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
                double d = 0;
                int ind =
                    t * Pxyz + z * Pxy + y * nvecs + s; //((t*Nz+z)*Ny+y)*nvecs+s;
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

#if 1
  // Go through the test cases -- apply SSE dslash versus, QDP Dslash
  for (int isign = 1; isign >= -1; isign -= 2) {
    for (int cb = 0; cb < 2; cb++) {
      int source_cb = 1 - cb;
      int target_cb = cb;
      masterPrintf("Timing on cb=%d isign=%d\n", cb, isign);
      masterPrintf("=============================\n");

      for (int repeat = 0; repeat < 3; repeat++) {
        double start = omp_get_wtime();

        for (int i = 0; i < args_.iters; i++) {
          // Apply Optimized Dslash
          D32.dslash(chi_s[target_cb],
                     psi_s[source_cb],
                     u_packed[target_cb],
                     isign,
                     target_cb);
        }

        double end = omp_get_wtime();
        double time = end - start;
        CommsUtils::sumDouble(&time);
        time /= (double)CommsUtils::numNodes();

        masterPrintf("\t timing %d of 3\n", repeat);
        masterPrintf("\t %d iterations in %e seconds\n", args_.iters, time);
        masterPrintf("\t %e usec/iteration\n", 1.0e6 * time / (double)args_.iters);
        double Gflops = (1320.0f + 72.0f) * (double)(args_.iters) *
                        (double)(X1h * Ny * Nz * Nt) / 1.0e9;
        double perf = Gflops / time;
        masterPrintf("\t Performance: %g GFLOPS total\n", perf);
      }
    }
  }
#endif

#if 1
  masterPrintf("Creating Twisted Mass Wilson Op\n");
  //  double Mass=0.1;
  EvenOddTMWilsonOperator<FT, V, S, compress> M(
      Mass, TwistedMass, u_packed, &geom, t_boundary, coeff_s, coeff_t);

  // Go through the test cases -- apply SSE dslash versus, QDP Dslash
  for (int isign = 1; isign >= -1; isign -= 2) {

    masterPrintf("Timing M: isign=%d\n", isign);
    masterPrintf("=============================\n");

    for (int repeat = 0; repeat < 3; repeat++) {
      double start = omp_get_wtime();

      for (int i = 0; i < args_.iters; i++) {
        // Apply Optimized Dslash
        M(chi_s[0], psi_s[0], isign);
      }

      double end = omp_get_wtime();
      double time = end - start;

      CommsUtils::sumDouble(&time);
      time /= (double)CommsUtils::numNodes();

      masterPrintf("\t timing %d of 3\n", repeat);
      masterPrintf("\t %d iterations in %e seconds\n", args_.iters, time);
      masterPrintf("\t %e usec/iteration\n", 1.0e6 * time / (double)args_.iters);
      double flops_per_iter = (1320.0f + 72.0f) * 2.0 + 24.0 * 3.0;
      double Gflops =
          flops_per_iter * (double)(args_.iters) * (double)(X1h * Ny * Nz * Nt) / 1.0e9;
      double perf = Gflops / time;
      masterPrintf("\t Performance: %g GFLOPS total\n", perf);
      masterPrintf("\t              %g GFLOPS / node\n",
                   perf / (double)CommsUtils::numNodes());
    }
  }
#endif

#if 1
  double rsd_target = rsdTarget<FT>::value;
  int max_iters = 100;
  int niters;
  double rsd_final;
  int len = (geom.getPxyz() * geom.Nt() * sizeof(Spinor)) / sizeof(FT);
  FT *c_s0 = (FT *)chi_s[0];

  {
    masterPrintf("Creating Solver\n");
    InvCG<FT, V, S, compress> solver(M, max_iters);

    masterPrintf("Tuning Solver\n");

    for (int solve = 0; solve < 1; solve++) {
      masterPrintf("Starting solver\n");
      unsigned long site_flops = 0;
      unsigned long mv_apps = 0;

      FT *psi_0 = (FT *)psi_s[0];
#pragma omp parallel for simd
      for (int i = 0; i < len; i++) {
        c_s0[i] = rep<FT, double>(0);
        psi_0[i] = rep<FT, double>(0);
      }

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
                  }
                }
              }
            }
          }
        }
      }

      start = omp_get_wtime();
      solver(chi_s[0],
             psi_s[0],
             rsd_target,
             niters,
             rsd_final,
             site_flops,
             mv_apps,
             1,
             verbose);
      end = omp_get_wtime();

      unsigned long num_cb_sites = X1h * Ny * Nz * Nt;
      unsigned long total_flops =
          (site_flops + (72 + 2 * (1320 + 72 /*twist*/)) * mv_apps) * num_cb_sites;
      masterPrintf("Solver Time=%g(s)\n", (end - start));
      masterPrintf("CG GFLOPS=%g\n", 1.0e-9 * (double)(total_flops) / (end - start));
    }
  } // Solver
#endif

#if 1
  {
    masterPrintf("Creating BiCGStab Solver\n");
    InvBiCGStab<FT, V, S, compress> solver2(M, max_iters);
    masterPrintf("Tuning BiCGStab Solver\n");

    for (int solve = 0; solve < 5; solve++) {
      unsigned long site_flops;
      unsigned long mv_apps;

      FT *psi_0 = (FT *)psi_s[0];
#pragma omp parallel for simd
      for (int i = 0; i < len; i++) {
        c_s0[i] = rep<FT, double>(0);
        psi_0[i] = rep<FT, double>(0);
      }

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
                  }
                }
              }
            }
          }
        }
      }

      start = omp_get_wtime();
      solver2(chi_s[0],
              psi_s[0],
              rsd_target,
              niters,
              rsd_final,
              site_flops,
              mv_apps,
              1,
              verbose,
              1);
      end = omp_get_wtime();

      unsigned long num_cb_sites = X1h * Ny * Nz * Nt;
      unsigned long total_flops =
          (site_flops + (72 + 2 * (1320 + 72 /*twist*/)) * mv_apps) * num_cb_sites;

      masterPrintf("Solver Time=%g(s)\n", (end - start));
      masterPrintf("BICGSTAB GFLOPS=%g\n",
                   1.0e-9 * (double)(total_flops) / (end - start));
    }
  }
#endif

  masterPrintf("Cleaning up\n");

  geom.free(packed_gauge_cb0);
  geom.free(packed_gauge_cb1);
  geom.free(p_even);
  geom.free(p_odd);
  geom.free(c_even);
  geom.free(c_odd);
}
