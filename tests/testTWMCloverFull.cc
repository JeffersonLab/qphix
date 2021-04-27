#include "unittest.h"
#include "testTWMCloverFull.h"

#include "qphix/tm_clov_dslash_def.h"
#include "qphix/tm_clov_dslash_body.h"

#include "qdp.h"
using namespace QDP;

#include "dslashm_w.h"

#include "reunit.h"

#include "qphix/qdp_packer.h"
#include "./clover_term.h"

#include "qphix/twisted_mass_clover.h"
#include "qphix/invcg.h"
#include "qphix/invbicgstab.h"

#include <omp.h>

using namespace Assertions;
using namespace std;
using namespace QPhiX;

#include "tolerance.h"
#include "veclen.h"
#include "compare_qdp_spinors_custom.h"
#include "cli_args.h"
#include "tparam_selector.h"
#include "RandomGauge.h"

void TestTMClover::run()
{
  call(*this, args_.prec, args_.soalen, args_.compress12);
}

namespace
{
/**
  Computes \f$ \psi_\text{out} = A^{-1} D \psi_\text{in} \f$.
  */
template <typename QdpClover, typename QdpGauge, typename QdpSpinor>
void qdp_dslash(QdpSpinor &out,
                QdpSpinor const &in,
                QDP::multi1d<QdpGauge> const &u_aniso,
                QdpClover const &inv_clov,
                int const isign,
                int const target_cb)
{
  QdpSpinor D_psi = zero;
  dslash(D_psi, u_aniso, in, isign, target_cb);
  inv_clov.apply(out, D_psi, isign, target_cb);
}

/**
  Computes \f$ \psi_\text{out} = A \chi_\text{chi} - 0.25 D \psi_\text{psi} \f$.
  */
template <typename QdpClover, typename QdpGauge, typename QdpSpinor>
void qdp_achimbdpsi(QdpSpinor &out,
                    QdpSpinor const &chi,
                    QdpSpinor const &psi,
                    QDP::multi1d<QdpGauge> const &u_aniso,
                    QdpClover const &clov,
                    int const isign,
                    int const target_cb)
{
  int const other_cb = 1 - target_cb;

  QdpSpinor D_psi = zero;
  dslash(D_psi, u_aniso, psi, isign, target_cb);
  clov.apply(out, chi, isign, target_cb);
  out[rb[target_cb]] -= 0.25 * D_psi;
}

/**
  Applies the full operator, \f$ A - 0.25 D A^{-1} D \f$.
  */
template <typename QdpClover, typename QdpGauge, typename QdpSpinor>
void qdp_apply_operator(QdpSpinor &out,
                        QdpSpinor const &in,
                        QDP::multi1d<QdpGauge> const &u_aniso,
                        QdpClover const &clov,
                        QdpClover const &inv_clov,
                        int const isign,
                        int const target_cb)
{
  int const other_cb = 1 - target_cb;

  QdpSpinor A_inv_D_psi = zero;
  qdp_dslash(A_inv_D_psi, in, u_aniso, inv_clov, isign, other_cb);

  qdp_achimbdpsi(out, in, A_inv_D_psi, u_aniso, clov, isign, target_cb);
}
}

template <typename FT, int V, int S, bool compress, typename U, typename Phi>
void TestTMClover::operator()()
{

  bool const verbose = false;
  QDPIO::cout << endl << "ENTERING TWISTED MASS CLOVER TEST" << endl;
  QDPIO::cout << "(testing tm clover with mu = 0 here...)" << endl << endl;

  // Diagnostic information:
  multi1d<int> const &lattSize = Layout::subgridLattSize();
  int const Nx = lattSize[0];
  int const Ny = lattSize[1];
  int const Nz = lattSize[2];
  int const Nt = lattSize[3];

  QDPIO::cout << "Inititalizing QDP++ gauge field." << endl;

  Geometry<FT, V, S, compress> geom(Layout::subgridLattSize().slice(),
                                    args_.By,
                                    args_.Bz,
                                    args_.NCores,
                                    args_.Sy,
                                    args_.Sz,
                                    args_.PadXY,
                                    args_.PadXYZ,
                                    args_.MinCt,
                                    true);

  RandomGauge<FT, V, S, compress> gauge(geom, 1, 0.1);

  // The TM clover operator expects an array of inverse clover terms such that
  // it is the normal and hermitian conjugate. Since there is no twisted mass
  // in this test case (because QDP++ would not support it), these two terms
  // will be the same. Here we build an array such that we can index with
  // checkerboard index and obtain the desired array.
  typename Geometry<FT, V, S, compress>::FullCloverBlock *gauge_inv_clov_array_cb0[2] = {
      gauge.invfullclov_packed[0], gauge.invfullclov_packed[0]};
  typename Geometry<FT, V, S, compress>::FullCloverBlock *gauge_inv_clov_array_cb1[2] = {
      gauge.invfullclov_packed[1], gauge.invfullclov_packed[1]};
  typename Geometry<FT, V, S, compress>::FullCloverBlock **gauge_inv_clov_array[2] = {
      gauge_inv_clov_array_cb0, gauge_inv_clov_array_cb1};

  HybridSpinor<FT, V, S, compress> hs_source(geom), hs_qphix1(geom),
      hs_qphix2(geom), hs_qdp1(geom), hs_qdp2(geom);

  gaussian(hs_source.qdp());
  //make_point_source(hs_source.qdp());

  hs_source.pack();

  TMClovDslash<FT, V, S, compress> D32(
      &geom, gauge.t_boundary, gauge.aniso_fac_s, gauge.aniso_fac_t);

// This tests dslash
#if 1
  QDPIO::cout << std::endl << "TESTING TWISTED-MASS CLOVER DSLASH:" << std::endl;

  int const cbsize_in_blocks = rb[0].numSiteTable() / S;

  for (int const isign : {1, -1}) {
    for (int const target_cb : {0, 1}) {
      int const source_cb = 1 - target_cb;

      // (a) Apply Optimized Dslash, Pack & Unpack
      D32.dslash(hs_qphix1[target_cb],
                 hs_source[source_cb],
                 gauge.u_packed[target_cb],
                 gauge_inv_clov_array[target_cb],
                 isign,
                 target_cb);
      hs_qphix1.unpack();

      // (b) Apply QDP Dslash
      qdp_dslash(hs_qdp1.qdp(),
                 hs_source.qdp(),
                 gauge.u_aniso,
                 gauge.invclov_qdp,
                 isign,
                 target_cb);

      expect_near(hs_qdp1.qdp(),
                  hs_qphix1.qdp(),
                  1e-6,
                  geom,
                  target_cb,
                  "Twisted-mass clover dslash");
    } // cb
  } // isign
#endif

// This tests achimbdpsi
#if 1
  QDPIO::cout << std::endl << "TESTING A CHI - b D PSI" << std::endl;

  double const beta = (double)(0.25); // Always 0.25

  for (int const isign : {1, -1}) {
    for (int const target_cb : {0, 1}) {
      int const source_cb = 1 - target_cb;

      // (a) Apply Optimized Dslash, Pack & Unpack
      D32.dslashAChiMinusBDPsi(hs_qphix1[target_cb],
                               hs_source[source_cb],
                               hs_source[target_cb],
                               gauge.u_packed[target_cb],
                               gauge.clov_packed[target_cb],
                               beta,
                               isign,
                               target_cb);
      hs_qphix1.unpack();

      // (b) Apply QDP Dslash
      qdp_achimbdpsi(hs_qdp1.qdp(),
                     hs_source.qdp(),
                     hs_source.qdp(),
                     gauge.u_aniso,
                     gauge.clov_qdp,
                     isign,
                     target_cb);

      expect_near(hs_qdp1.qdp(), hs_qphix1.qdp(), 1e-6, geom, target_cb, "Twisted-mass clover achimbdpsi");
    } // cb
  } // isign
#endif

  for (int const cb : {1, 0}) {
    // This tests the full even/odd operators
    QDPIO::cout << std::endl
                << std::endl
                << "TESTING TWISTED MASS EVEN/ODD (DEGENERATE) DSLASH CLOVER OPERATOR"
                << std::endl;

    int const other_cb = 1 - cb;

    EvenOddTMCloverOperator<FT, V, S, compress> M(gauge.u_packed,
                                                  gauge.clov_packed[cb],
                                                  gauge_inv_clov_array[other_cb],
                                                  &geom,
                                                  gauge.t_boundary,
                                                  gauge.aniso_fac_s,
                                                  gauge.aniso_fac_t);
    Phi ltmp = zero;
    Real betaFactor = Real(0.25);

    for (int const isign : {1, -1}) {
      masterPrintf("isign: %d, cb: %d\n", isign, cb);

      // (a) Apply QPhiX Operator
      M(hs_qphix1[cb], hs_source[cb], isign, cb);
      hs_qphix1.unpack();

      // (b) Apply QDP Dslash
      qdp_apply_operator(hs_qdp1.qdp(),
                         hs_source.qdp(),
                         gauge.u_aniso,
                         gauge.clov_qdp,
                         gauge.invclov_qdp,
                         isign,
                         cb);

      expect_near(hs_qdp1.qdp(), hs_qphix1.qdp(), 1e-7, geom, cb, "TM clover linop");

    } // isign

// This tests the CG
#if 1
  {
    QDPIO::cout << std::endl
                << std::endl
                << "TESTING CG SOLVER W/ TWISTED MASS CLOVER OPERATOR" << std::endl;

    double rsd_target = rsdTarget<FT>::value;
    int max_iters = 500;
    int niters;
    double rsd_final;
    unsigned long site_flops;
    unsigned long mv_apps;

    InvCG<FT, V, S, compress> solver(M, max_iters);

    double start = omp_get_wtime();
    solver(hs_qphix1[cb],
           hs_source[cb],
           rsd_target,
           niters,
           rsd_final,
           site_flops,
           mv_apps,
           1,
           verbose,
           cb);
    double end = omp_get_wtime();

    hs_qphix1.unpack();

    // Multiply back
    // chi2 = M chi
    qdp_apply_operator(hs_qphix2.qdp(),
                       hs_qphix1.qdp(),
                       gauge.u_aniso,
                       gauge.clov_qdp,
                       gauge.invclov_qdp,
                       1,
                       cb);

    // chi3 = M^\dagger chi2
    qdp_apply_operator(hs_qdp1.qdp(),
                       hs_qphix2.qdp(),
                       gauge.u_aniso,
                       gauge.clov_qdp,
                       gauge.invclov_qdp,
                       -1,
                       cb);

    expect_near(hs_source.qdp(), hs_qdp1.qdp(), 1e-7, geom, cb, "TM Clover CG");

    int Nxh = Nx / 2;
    unsigned long num_cb_sites = Nxh * Ny * Nz * Nt;

    unsigned long total_flops =
        (site_flops + (1320 + 504 + 1320 + 504 + 48) * mv_apps) * num_cb_sites;

    masterPrintf("GFLOPS=%e\n", 1.0e-9 * (double)(total_flops) / (end - start));
  }
#endif

// This tests the BiCG
#if 1
  {
    QDPIO::cout << std::endl
                << std::endl
                << "TESTING BICGSTAB SOLVER W/ TWISTED MASS CLOVER OPERATOR"
                << std::endl;

    double rsd_target = rsdTarget<FT>::value;
    int max_iters = 500;
    int niters;
    double rsd_final;
    unsigned long site_flops;
    unsigned long mv_apps;

    InvBiCGStab<FT, V, S, compress> solver(M, max_iters);
    const int isign = 1;
    double start = omp_get_wtime();
    solver(hs_qphix1[cb],
           hs_source[cb],
           rsd_target,
           niters,
           rsd_final,
           site_flops,
           mv_apps,
           isign,
           verbose,
           cb);
    double end = omp_get_wtime();
    hs_qphix1.unpack();

    // Multiply back
    // chi2 = M chi
    qdp_apply_operator(hs_qphix2.qdp(),
                       hs_qphix1.qdp(),
                       gauge.u_aniso,
                       gauge.clov_qdp,
                       gauge.invclov_qdp,
                       1,
                       cb);

    expect_near(hs_source.qdp(), hs_qphix2.qdp(), 1e-7, geom, cb, "TM Clover BiCGStab");

    int Nxh = Nx / 2;
    unsigned long num_cb_sites = Nxh * Ny * Nz * Nt;

    unsigned long total_flops =
        (site_flops + (1320 + 504 + 1320 + 504 + 48) * mv_apps) * num_cb_sites;
    masterPrintf("GFLOPS=%e\n", 1.0e-9 * (double)(total_flops) / (end - start));
  }
#endif
  }

  QDPIO::cout << std::endl;
}
