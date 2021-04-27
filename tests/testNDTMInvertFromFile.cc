#include "testNDTMInvertFromFile.h"

#include "unittest.h"

#include "qdp.h"
using namespace QDP;

#include "dslashm_w.h"

#include "reunit.h"

#include "qphix/qdp_packer.h"
#include "clover_term.h"

#if 1
#include "qphix/ndtm_reuse_operator_clover.h"
#include "qphix/twisted_mass_clover.h"
#include "qphix/wilson.h"
#include "qphix/invcg.h"
#include "qphix/invbicgstab.h"
#endif

#include <omp.h>

using namespace Assertions;
using namespace std;
using namespace QPhiX;

#ifndef QPHIX_SOALEN
#define QPHIX_SOALEN 4
#endif

#if defined(QPHIX_MIC_SOURCE) || defined(QPHIX_AVX512_SOURCE)

#define VECLEN_SP 16
#define VECLEN_HP 16
#define VECLEN_DP 8

#elif defined(QPHIX_AVX_SOURCE) || defined(QPHIX_AVX2_SOURCE)

#define VECLEN_SP 8
#define VECLEN_DP 4

#elif defined(QPHIX_SSE_SOURCE)
#define VECLEN_SP 4
#define VECLEN_DP 2

#elif defined(QPHIX_SCALAR_SOURCE)
#warning SCALAR_SOURCE
#define VECLEN_DP 1
#define VECLEN_SP 1

#elif defined(QPHIX_QPX_SOURCE)
#define VECLEN_DP 4
#define VECLEN_SP 4

#endif

#include "qphix/blas_new_c.h"
#include "packed_spinor.h"

#include "tolerance.h"

void MesPlq(const multi1d<LatticeColorMatrixD> &u, Double &w_plaq, Double &link)
{
  w_plaq = link = Double(0);

  // Compute the average plaquettes
  for (int mu = 1; mu < Nd; ++mu) {
    for (int nu = 0; nu < mu; ++nu) {
      /* tmp_0 = u(x+mu,nu)*u_dag(x+nu,mu) */
      LatticeColorMatrix tmp_0 =
          shift(u[nu], FORWARD, mu) * adj(shift(u[mu], FORWARD, nu));

      /* tmp_1 = tmp_0*u_dag(x,nu)=u(x+mu,nu)*u_dag(x+nu,mu)*u_dag(x,nu)
       */
      LatticeColorMatrix tmp_1 = tmp_0 * adj(u[nu]);

      /* tmp =
       * sum(tr(u(x,mu)*tmp_1=u(x,mu)*u(x+mu,nu)*u_dag(x+nu,mu)*u_dag(x,nu)))
       */
      Double tmp = sum(real(trace(u[mu] * tmp_1)));

      w_plaq += tmp;
    }
  }
  // Normalize
  w_plaq *= 2.0 / double(Layout::vol() * Nd * (Nd - 1) * Nc);

  // Compute the average link
  for (int mu = 0; mu < Nd; ++mu)
    link += sum(real(trace(u[mu])));

  link /= double(12 * Layout::vol());
}

/**
  Complex fused-multiply-accumulate.

  Computes \f$ (r + \mathrm i i) += (a + \mathrm i b) * (c + \mathrm i d) \f$.
  */
template <typename FT>
void cplx_fmacc(
    FT const &a, FT const &b, FT const &c, FT const &d, FT &r_out, FT &i_out)
{
  r_out += a * c - b * d;
  i_out += a * d + b * c;
}

template <typename FT, int V, int S, bool compress>
void inv_clover_spinor_multiplication(
    typename Geometry<FT, V, S, compress>::FullCloverBlock *A_inv[2][2][2],
    typename Geometry<FT, V, S, compress>::FourSpinorBlock *in[2],
    typename Geometry<FT, V, S, compress>::FourSpinorBlock *out[2],
    Geometry<FT, V, S, compress> &geom,
    int isign = 1)
{
  zeroSpinor<FT, V, S, compress, 2>(out, geom, 1);

  // Iterate through all the block.
  auto const num_blocks = geom.get_num_blocks();
  for (auto block = 0u; block < num_blocks; ++block) {
    // Input flavor.
    for (auto f_in : {0, 1}) {
      // Output flavor.
      for (auto f_out : {0, 1}) {
        // The clover term is block-diagonal in spin. Therefore we need
        // to iterate over the two blocks of spin.
        for (auto s_block : {0, 1}) {
          // Extract the spin block as a handy alias.
          auto const &block_in = s_block == 0
                                     ? A_inv[f_out][f_in][isign][block].block1
                                     : A_inv[f_out][f_in][isign][block].block2;
          // Input two-spinor component.
          for (auto s_in : {0, 1}) {
            // Reconstruct four spinor index.
            auto const four_s_in = 2 * s_block + s_in;
            // Output two-spinor component.
            for (auto s_out : {0, 1}) {
              // Reconstruct four spinor index.
              auto const four_s_out = 2 * s_block + s_out;
              // Input color.
              for (auto c_in : {0, 1, 2}) {
                // Spin-color index (0, ..., 5).
                auto const sc_in = 3 * s_in + c_in;
                // Output color.
                for (auto c_out : {0, 1, 2}) {
                  // Spin-color index (0, ..., 5).
                  auto const sc_out = 3 * s_out + c_out;
                  // SIMD vector.
                  for (auto v = 0; v < V; ++v) {
                    cplx_fmacc(block_in[sc_out][sc_in][0][v],
                               block_in[sc_out][sc_in][1][v],
                               in[f_in][block][c_in][four_s_in][0][v],
                               in[f_in][block][c_in][four_s_in][1][v],
                               out[f_out][block][c_out][four_s_out][0][v],
                               out[f_out][block][c_out][four_s_out][1][v]);
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}

template <typename FT, int V, int S, bool compress, typename U, typename Phi>
void testNDTMInvertFromFile::runTest(double mass,
                                     double clov_coeff,
                                     const std::string &filename)
{

  typedef typename Geometry<FT, V, S, compress>::FourSpinorBlock Spinor;
  typedef typename Geometry<FT, V, S, compress>::SU3MatrixBlock Gauge;
  typedef typename Geometry<FT, V, S, compress>::FullCloverBlock FullCloverBlock;

  bool verbose = true;
  QDPIO::cout << endl << "ENTERING CLOVER DSLASH TEST" << endl;

  // Diagnostic information:
  const multi1d<int> &lattSize = Layout::subgridLattSize();
  int Nx = lattSize[0];
  int Ny = lattSize[1];
  int Nz = lattSize[2];
  int Nt = lattSize[3];

  QDPIO::cout << "VECLEN=" << V << "  SOALEN=" << S << endl;
  QDPIO::cout << "Lattice Size: ";
  for (int mu = 0; mu < lattSize.size(); mu++) {
    QDPIO::cout << " " << lattSize[mu];
  }
  QDPIO::cout << endl;

  QDPIO::cout << "Block Sizes: By=" << By << " Bz=" << Bz << endl;
  QDPIO::cout << "N Cores" << NCores << endl;
  QDPIO::cout << "SMT Grid: Sy=" << Sy << " Sz=" << Sz << endl;
  QDPIO::cout << "Pad Factors: PadXY=" << PadXY << " PadXYZ=" << PadXYZ << endl;
  QDPIO::cout << "MinCt=" << MinCt << endl;
  QDPIO::cout << "Threads_per_core = " << N_simt << endl;

  multi1d<LatticeDiracFermionF> control_solution(1);

  QDPIO::cout << "Reading inverted from file" << endl;
  XMLReader control_file_xml, control_record_xml;
  QDPFileReader control_from(
      control_file_xml, "../data/test8.tmlqcd.lime", QDPIO_SERIAL);
  read(control_from, control_record_xml, control_solution);
  close(control_from);

  // What we consider to be small enough...
  QDPIO::cout << "Inititalizing QDP++ gauge field" << endl;

  multi1d<U> u(4);

  QDPIO::cout << "Reading field from file " << endl;
  XMLReader file_xml, record_xml;
  QDPFileReader from(file_xml, "../data/test8.gauge.lime", QDPIO_SERIAL);
  read(from, record_xml, u);
  close(from);

  Double w_plaq, link;
  MesPlq(u, w_plaq, link);
  QDPIO::cout << "Plaquette on reading: w=" << w_plaq << " link=" << link << endl;
  QDPIO::cout << "Reunitarizing links\n";
  for (int mu = 0; mu < 4; mu++) {
    reunit(u[mu]);
  }
  Double w_plaq2, link2;
  MesPlq(u, w_plaq2, link2);
  QDPIO::cout << "Plaquette after reunit: w=" << w_plaq2 << " link=" << link2
              << endl;

  QDPIO::cout << "Delta w_plaq" << w_plaq2 - w_plaq << endl;

  double t_boundary = double(-1);

  QDPIO::cout << "TODO" << endl;
  Geometry<FT, V, S, compress> geom(Layout::subgridLattSize().slice(),
                                    By,
                                    Bz,
                                    NCores,
                                    Sy,
                                    Sz,
                                    PadXY,
                                    PadXYZ,
                                    MinCt);

  QDPIO::cout << "Allocating packged gauge fields" << endl;
  Gauge *packed_gauge_cb0 = (Gauge *)geom.allocCBGauge();
  Gauge *packed_gauge_cb1 = (Gauge *)geom.allocCBGauge();

  QDPIO::cout << "Allocating clover blocks" << endl;
  FullCloverBlock *A[2][2];
  FullCloverBlock *A_inv[2][2][2];

  for (auto f1 : {0, 1}) {
    for (auto isign : {0, 1}) {
      A[f1][isign] = geom.allocCBFullClov();
      for (auto f2 : {0, 1}) {
        A_inv[f1][f2][isign] = geom.allocCBFullClov();
      }
    }
  }

  FullCloverBlock const *const const_A_inv[2][2][2] = {
      {{A_inv[0][0][0], A_inv[0][0][1]}, {A_inv[0][1][0], A_inv[0][1][1]}},
      {{A_inv[1][0][0], A_inv[1][0][1]}, {A_inv[1][1][0], A_inv[1][1][1]}},
  };

  double const kappa = 0.16065;
  double const mu = 0.000000;
  double const epsilon = 0.000;

  // double const kappa = 0.1632550;
  // double const mu = 0.150;
  // double const epsilon = 0.190;

  double const alpha = 1 / (2 * kappa);
  mass = alpha - 4;

  QDPIO::cout << "kappa = " << kappa << endl;
  QDPIO::cout << "mu = " << mu << endl;
  QDPIO::cout << "epsilon = " << epsilon << endl;
  QDPIO::cout << "alpha = " << alpha << endl;
  QDPIO::cout << "mass = " << mass << endl;

  double const denominator = alpha * alpha + mu * mu - epsilon * epsilon;
  double const tau_0_factor = alpha / denominator;
  double const tau_3_factor = mu / denominator;
  double const tau_1_factor = epsilon / denominator;

  // Populate clover term and the inverse.
  QDPIO::cout << "Populate clover terms" << endl;
  auto const num_blocks = geom.get_num_blocks();
  for (auto block = 0u; block < num_blocks; ++block) {
    // QDPIO::cout << "A Block " << block << "/" << num_blocks << endl;
    for (auto f : {0, 1}) {
      for (auto isign_idx : {0, 1}) {
        for (auto s : {0, 1}) {
          for (auto c : {0, 1, 2}) {
            auto const sc = 3 * s + c;
            for (auto v = 0; v < V; ++v) {
              A[f][isign_idx][block].block1[sc][sc][0][v] = alpha;
              A[f][isign_idx][block].block1[sc][sc][1][v] =
                  mu * (isign_idx == 0 ? 1 : -1) * (f == 0 ? 1 : -1);
              A[f][isign_idx][block].block2[sc][sc][0][v] = alpha;
              A[f][isign_idx][block].block2[sc][sc][1][v] =
                  -mu * (isign_idx == 0 ? 1 : -1) * (f == 0 ? 1 : -1);
            }
          }
        }
      }
    }
  }

  QDPIO::cout << "Populate clover inverse" << endl;
  for (auto block = 0u; block < num_blocks; ++block) {
    // QDPIO::cout << "A_inv Block " << block << "/" << num_blocks << endl;
    for (auto f1 : {0, 1}) {
      for (auto f2 : {0, 1}) {
        for (auto isign_idx : {0, 1}) {
          for (auto s : {0, 1}) {
            for (auto c : {0, 1, 2}) {
              auto const sc = 3 * s + c;
              for (auto v = 0; v < V; ++v) {
                A_inv[f1][f2][isign_idx][block].block1[sc][sc][0][v] = tau_0_factor;
                A_inv[f1][f2][isign_idx][block].block2[sc][sc][0][v] = tau_0_factor;

                if (f1 == f2) {
                  A_inv[f1][f2][isign_idx][block].block1[sc][sc][1][v] =
                      tau_3_factor * (isign_idx == 0 ? 1 : -1) * (f1 == 0 ? 1 : -1);
                  A_inv[f1][f2][isign_idx][block].block2[sc][sc][1][v] =
                      -tau_3_factor * (isign_idx == 0 ? 1 : -1) * (f1 == 0 ? 1 : -1);
                } else {
                  A_inv[f1][f2][isign_idx][block].block1[sc][sc][1][v] =
                      tau_1_factor;
                  A_inv[f1][f2][isign_idx][block].block2[sc][sc][1][v] =
                      tau_1_factor;
                }
              }
            }
          }
        }
      }
    }
  }

  QDPIO::cout << "Packing gauge field" << endl;
  qdp_pack_gauge<>(u, packed_gauge_cb0, packed_gauge_cb1, geom);

  Gauge *u_packed[2];
  u_packed[0] = packed_gauge_cb0;
  u_packed[1] = packed_gauge_cb1;

  // Modify u (antiperiodic BC's)
  QDPIO::cout << "Implement antiperiodic boundary conditions" << endl;
  u[3] *= where(Layout::latticeCoordinate(3) == (Layout::lattSize()[3] - 1),
                Real(t_boundary),
                Real(1));

  int max_iters = 5000;

  QDPIO::cout << "Create the operator M" << endl;

  // TODO What do those anisotropy factors really mean?
  EvenOddNDTMCloverReuseOperator<FT, V, S, compress> M_ndtm(
      u_packed, A, A_inv, epsilon, &geom, t_boundary, 1.0, 1.0);

  EvenOddTMCloverOperator<FT, V, S, compress> M_tm_clov(
      u_packed, A[0], A_inv[0][0], &geom, t_boundary, 1.0, 1.0);

  EvenOddWilsonOperator<FT, V, S, compress> M_wilson(
      mass, u_packed, &geom, t_boundary, 1.0, 1.0);

  TMClovDslash<FT, V, S, compress> Dslash(&geom, t_boundary, 1.0, 1.0);

  QDPIO::cout << "Create the solver" << endl;

  InvCG<FT, V, S, compress, TwoFlavEvenOddLinearOperator<FT, V, S, compress>>
      solver_ndtm(M_ndtm, max_iters);
  InvCG<FT, V, S, compress> solver_tm_clov(M_tm_clov, max_iters);
  InvCG<FT, V, S, compress> solver_wilson(M_wilson, max_iters);

  QDPIO::cout << "Tune the solver" << endl;
  solver_ndtm.tune();
  solver_tm_clov.tune();
  solver_wilson.tune();

  double rsd_target = 1.0e-10;
  Real betaFactor = Real(0.25);

  int niters;
  double rsd_final;
  unsigned long site_flops;
  unsigned long mv_apps;

  QDPIO::cout << "Create input and output spinors" << endl;
  PackedSpinor<FT, V, S, compress, Phi> spinor_out[2](geom);
  PackedSpinor<FT, V, S, compress, Phi> spinor_in[2](geom);

  PackedSpinor<FT, V, S, compress, Phi> spinor_b_o_tilde[2](geom);

  PackedSpinor<FT, V, S, compress, Phi> spinor_tmp1[2](geom);
  PackedSpinor<FT, V, S, compress, Phi> spinor_tmp2[2](geom);

  QDPIO::cout << "Zero the input spinors" << endl;
  for (auto f : {0, 1}) {
    spinor_in[f].all = zero;
    // bo_tilde[f].all = zero;
  }

  QDPIO::cout << "Create a point source" << endl;
  spinor_in[0] // First flavor
      .all // Take the whole QDP++ spinor
      .elem(0) // Lattice
      .elem(0) // Spin
      .elem(0) // Color
      .real() = 1.0;
  spinor_in[1] // Second flavor
      .all // Take the whole QDP++ spinor
      .elem(0) // Lattice
      .elem(0) // Spin
      .elem(0) // Color
      .real() = 1.0;

  QDPIO::cout << "Pack the input spinors" << endl;
  for (auto f : {0, 1}) {
    spinor_in[f].pack();
  }

  Spinor *b_e[2] = {spinor_in[0][0], spinor_in[1][0]};
  Spinor *b_o[2] = {spinor_in[0][1], spinor_in[1][1]};

  Spinor *tmp1[2] = {spinor_tmp1[0][0], spinor_tmp1[1][0]};
  Spinor *tmp2[2] = {spinor_tmp2[0][1], spinor_tmp2[1][1]};

  Spinor *b_o_tilde[2] = {spinor_b_o_tilde[0][1], spinor_b_o_tilde[1][1]};
  Spinor *x_o[2] = {spinor_out[0][1], spinor_out[1][1]};
  Spinor *x_e[2] = {spinor_out[0][0], spinor_out[1][0]};

  double norm;

  norm2Spinor<FT, V, S, compress, 2>(norm, b_e, geom, 1);
  QDPIO::cout << "|| b_e ||^2 = " << norm << endl;

  norm2Spinor<FT, V, S, compress, 2>(norm, b_o, geom, 1);
  QDPIO::cout << "|| b_o ||^2 = " << norm << endl;

  QDPIO::cout << "Prepare the source" << endl;
  // \tilde b_o = - M_oe M_ee^{-1} b_e + b_o
  inv_clover_spinor_multiplication(A_inv, b_e, tmp1, geom);
  // FIXME Check the checkerboard indices.
  Dslash.two_flav_dslash(tmp2, tmp1, u_packed[0], const_A_inv, 1, 1);
  aypx<FT, V, S, compress, 2>(-1.0, b_o, tmp2, geom, 1);

  norm2Spinor<FT, V, S, compress, 2>(norm, tmp2, geom, 1);
  QDPIO::cout << "|| \\tilde b_o ||^2 = " << norm << endl;

  norm2Spinor<FT, V, S, compress, 2>(norm, x_e, geom, 1);
  QDPIO::cout << "|| x_e ||^2 = " << norm << endl;

  norm2Spinor<FT, V, S, compress, 2>(norm, x_o, geom, 1);
  QDPIO::cout << "|| x_o ||^2 = " << norm << endl;

  int isign = 1;

  auto cb = 1;

  QDPIO::cout << "Apply M^\\dagger" << endl;
  Dslash.two_flav_dslash(b_o_tilde, tmp2, u_packed[0], const_A_inv, -1, 1);

  norm2Spinor<FT, V, S, compress, 2>(norm, b_o_tilde, geom, 1);
  QDPIO::cout << "|| \\tilde b_o ||^2 = " << norm << endl;

  double start = omp_get_wtime();
  QDPIO::cout << "Call the solver" << endl;

  bool const use_two_flav = false;
  if (use_two_flav) {
    solver_ndtm(x_o,
                b_o_tilde,
                rsd_target,
                niters,
                rsd_final,
                site_flops,
                mv_apps,
                isign,
                verbose);
  } else {
    solver_wilson(x_o[0],
                  b_o_tilde[0],
                  rsd_target,
                  niters,
                  rsd_final,
                  site_flops,
                  mv_apps,
                  isign,
                  verbose);
  }

  double end = omp_get_wtime();

  QDPIO::cout << "Solve took " << (end - start) << " seconds" << endl;

  // Reconstruct the solution.
  // M_ee^{-1} (b_e - M_eo x_o)
  // FIXME Check the checkerboard indices.
  Dslash.two_flav_dslash(tmp1, x_o, u_packed[1], const_A_inv, 1, 0);
  aypx<FT, V, S, compress, 2>(-1.0, b_e, tmp1, geom, 1);
  inv_clover_spinor_multiplication(A_inv, tmp1, x_e, geom);

  QDPIO::cout << "Unpack the output spinors" << endl;
  for (auto f : {0, 1}) {
    spinor_out[f].unpack();
  }

  // TODO Write out the propagators.

  // TODO Adjust the values of the masses. The actual values are encoded in
  // the filename.

  XMLBufferWriter out_file_xml, out_record_xml;
  QDPFileWriter file_writer(
      out_file_xml, "ndtm_inverted.lime", QDPIO_MULTIFILE, QDPIO_SERIAL);
  for (auto f : {0, 1}) {
    QDPIO::cout << "Writing spinor" << endl;
    write(file_writer, out_record_xml, spinor_out[f].all);
  }
  close(file_writer);

  geom.free(packed_gauge_cb0);
  geom.free(packed_gauge_cb1);
}

void testNDTMInvertFromFile::run(void)
{
  typedef LatticeColorMatrixF UF;
  typedef LatticeDiracFermionF PhiF;

  typedef LatticeColorMatrixD UD;
  typedef LatticeDiracFermionD PhiD;

  // std::string
  // filename("/home/mu/Dokumente/Studium/Master_Science_Physik/A100.24_L24_T48_beta190_mul0100_musig150_mudel190_kappa1632550/conf.2999");
  std::string filename(
      "/home/mu/Dokumente/Studium/Master_Science_Physik/test8x8x8x16/conf.0000");
  double mass = -0.3550;
  double clov_coeff = 1.90497469553511;
#if defined(QPHIX_MIC_SOURCE) || defined(QPHIX_AVX512_SOURCE)
  runTest<double, VECLEN_DP, 8, true, UD, PhiD>(mass, clov_coeff, filename);
#else
  runTest<double, VECLEN_DP, 4, true, UD, PhiD>(mass, clov_coeff, filename);
#endif
}
