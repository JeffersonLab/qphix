#include "unittest.h"
#include "testClovDslashFull.h"

#include "qphix/clover_dslash_def.h"
#include "qphix/clover_dslash_body.h"

#include "qdp.h"
using namespace QDP;

#include "dslashm_w.h"
#include "reunit.h"

#include "qphix/qdp_packer.h"
#include "./clover_term.h"

#if 1
#include "qphix/clover.h"
#include "qphix/invcg.h"
#include "qphix/invbicgstab.h"
#endif

#include <omp.h>

using namespace Assertions;
using namespace std;
using namespace QPhiX;

#include "RandomGauge.h"
#include "compare_qdp_spinors_custom.h"
#include "veclen.h"
#include "tolerance.h"

void TestClover::run()
{
  call(*this, args_.prec, args_.soalen, args_.compress12);
}

template <typename FT,
          int V,
          int S,
          bool compress,
          typename QdpGauge,
          typename QdpSpinor>
void TestClover::operator()()
{

  typedef typename Geometry<FT, V, S, compress>::FourSpinorBlock Spinor;
  typedef typename Geometry<FT, V, S, compress>::SU3MatrixBlock Gauge;
  typedef typename Geometry<FT, V, S, compress>::CloverBlock Clover;

  bool verbose = true;
  QDPIO::cout << endl << "ENTERING CLOVER DSLASH TEST" << endl;

  int Nx = args_.nrow_in[0];
  int Ny = args_.nrow_in[1];
  int Nz = args_.nrow_in[2];
  int Nt = args_.nrow_in[3];

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

  RandomGauge<FT, V, S, compress, QdpGauge, QdpSpinor> gauge(geom, 1.0);

  ClovDslash<FT, V, S, compress> D32(
      &geom, gauge.t_boundary, gauge.aniso_fac_s, gauge.aniso_fac_t);

  HybridSpinor<FT, V, S, compress, QdpSpinor> hs_source(geom), hs_qphix1(geom),
      hs_qphix2(geom), hs_qdp1(geom), hs_qdp2(geom);

  bool const point_source = false;
  if (point_source) {
    hs_source.zero();
    hs_source.qdp()
        .elem(0) // Lattice
        .elem(0) // Spin
        .elem(0) // Color
        .real() = 1.0;
    hs_source.qdp()
        .elem(0) // Lattice
        .elem(0) // Spin
        .elem(0) // Color
        .imag() = 0.0;
  } else {
    gaussian(hs_source.qdp());
  }
  hs_source.pack();

  {
    double norm;
    norm2Spinor(norm, hs_source[0], geom, 1);
    masterPrintf("Norm-sq of source: %g\n", norm);
    assert(norm > 0.0 && "Norm-sq of source must be positive");
  }

#if 1
  // Test only Dslash operator.
  // For clover this will be: A^{-1}_(1-cb,1-cb) D_(1-cb, cb)  psi_cb
  QDPIO::cout << "Testing Dslash \n" << endl;

  for (int isign = 1; isign >= -1; isign -= 2) {
    for (int cb = 0; cb < 2; cb++) {
      int source_cb = 1 - cb;
      int target_cb = cb;

#if 1
      // Apply Optimized Dslash
      hs_qphix1.zero();
      D32.dslash(hs_qphix1[target_cb],
                 hs_source[source_cb],
                 gauge.u_packed[target_cb],
                 gauge.invclov_packed[target_cb],
                 isign,
                 target_cb);
      hs_qphix1.unpack();

#else
      const int cbsize_in_blocks = rb[0].numSiteTable() / S;

      // Work directly with the QDP-JIT pointers...
      Spinor *chi_targ = (Spinor *)(chi.getFjit()) + target_cb * cbsize_in_blocks;
      Spinor *psi_src = (Spinor *)(psi.getFjit()) + source_cb * cbsize_in_blocks;

      QDPIO::cout << "Direct Interface" << std::endl;
      Spinor *clov_chi_targ =
          (Spinor *)(clov_chi.getFjit()) + target_cb * cbsize_in_blocks;
      Spinor *psi_in_src = (Spinor *)(psi.getFjit()) + source_cb * cbsize_in_blocks;

      D32.dslash(clov_chi_targ,
                 psi_in_src,
                 u_packed[target_cb],
                 invclov_packed[target_cb],
                 isign,
                 target_cb);
#endif
      // Apply QDP Dslash
      dslash(hs_qdp1.qdp(), gauge.u_aniso, hs_source.qdp(), isign, target_cb);
      gauge.invclov_qdp.apply(hs_qdp2.qdp(), hs_qdp1.qdp(), isign, target_cb);

      {
        hs_qdp1.pack();
        hs_qdp2.pack();

        double norm;

        norm2Spinor(norm, hs_qdp1[target_cb], geom, 1);
        masterPrintf("Norm-sq of hs_qdp1[target_cb]: %g\n", norm);

        norm2Spinor(norm, hs_qdp2[target_cb], geom, 1);
        masterPrintf("Norm-sq of hs_qdp2[target_cb]: %g\n", norm);
      }

      expect_near(hs_qdp2.qdp(),
                  hs_qphix1.qdp(),
                  1e-6,
                  geom,
                  target_cb,
                  "Wilson clover Dslash, QPhiX vs. QDP++");

    } // cb
  } // isign

#endif

#if 1
  // Go through the test cases -- apply SSE dslash versus, QDP Dslash
  // Test ax - bDslash y
  QDPIO::cout << "Testing dslashAchiMinusBDPsi" << endl;

  for (int isign = 1; isign >= -1; isign -= 2) {
    for (int cb = 0; cb < 2; cb++) {
      int source_cb = 1 - cb;
      int target_cb = cb;

      hs_qphix1.zero();

      double beta = 0.25; // Always 0.25

      // Apply Optimized Dslash
      D32.dslashAChiMinusBDPsi(hs_qphix1[target_cb],
                               hs_source[source_cb],
                               hs_source[target_cb],
                               gauge.u_packed[target_cb],
                               gauge.clov_packed[target_cb],
                               beta,
                               isign,
                               target_cb);
      hs_qphix1.unpack();

      // Apply QDP Dslash
      dslash(hs_qdp1.qdp(), gauge.u_aniso, hs_source.qdp(), isign, target_cb);
      QdpSpinor res = zero;
      gauge.clov_qdp.apply(res, hs_source.qdp(), isign, target_cb);
      res[rb[target_cb]] -= beta * hs_qdp1.qdp();

      // Check the difference per number in chi vector
      expect_near(res,
                  hs_qphix1.qdp(),
                  1e-6,
                  geom,
                  target_cb,
                  "A chi - b d psi, QPhiX vs. QDP++");
    }
  }

#endif

#if 1
  // Test only Dslash operator.
  // For clover this will be: A^{-1}_(1-cb,1-cb) D_(1-cb, cb)  psi_cb
  QDPIO::cout << "Testing Dslash With antiperiodic BCs \n" << endl;

  RandomGauge<FT, V, S, compress, QdpGauge, QdpSpinor> gauge_antip(geom, -1.0, 0.1);

  // Create Antiperiodic Dslash
  ClovDslash<FT, V, S, compress> D32_ap(&geom,
                                        gauge_antip.t_boundary,
                                        gauge_antip.aniso_fac_s,
                                        gauge_antip.aniso_fac_t);

  for (int isign = 1; isign >= -1; isign -= 2) {
    for (int cb = 0; cb < 2; cb++) {
      int source_cb = 1 - cb;
      int target_cb = cb;

      hs_qphix1.zero();

      // Apply Optimized Dslash
      D32_ap.dslash(hs_qphix1[target_cb],
                    hs_source[source_cb],
                    gauge_antip.u_packed[target_cb],
                    gauge_antip.invclov_packed[target_cb],
                    isign,
                    target_cb);
      hs_qphix1.unpack();

      // Apply QDP Dslash
      dslash(hs_qdp1.qdp(), gauge_antip.u_aniso, hs_source.qdp(), isign, target_cb);
      gauge_antip.invclov_qdp.apply(hs_qdp2.qdp(), hs_qdp1.qdp(), isign, target_cb);

      expect_near(hs_qdp2.qdp(),
                  hs_qphix1.qdp(),
                  1e-6,
                  geom,
                  target_cb,
                  "Dslash with antiperiodic T, QPhiX vs. QDP++");
    } // cb
  } // isign

#endif

#if 1
  // Go through the test cases -- apply SSE dslash versus, QDP Dslash
  // Test ax - bDslash y
  QDPIO::cout << "Testing dslashAchiMinusBDPsi" << endl;

  for (int isign = 1; isign >= -1; isign -= 2) {
    for (int cb = 0; cb < 2; cb++) {
      int source_cb = 1 - cb;
      int target_cb = cb;

      double beta = 0.25; // Always 0.25

      // Apply Optimized Dslash
      hs_qphix1.zero();
      D32_ap.dslashAChiMinusBDPsi(hs_qphix1[target_cb],
                                  hs_source[source_cb],
                                  hs_source[target_cb],
                                  gauge_antip.u_packed[target_cb],
                                  gauge_antip.clov_packed[target_cb],
                                  beta,
                                  isign,
                                  target_cb);
      hs_qphix1.unpack();

      // Apply QDP Dslash
      dslash(hs_qdp1.qdp(), gauge_antip.u_aniso, hs_source.qdp(), isign, target_cb);
      QdpSpinor res = zero;
      gauge_antip.clov_qdp.apply(res, hs_source.qdp(), isign, target_cb);
      res[rb[target_cb]] -= beta * hs_qdp1.qdp();

      expect_near(res,
                  hs_qphix1.qdp(),
                  1e-6,
                  geom,
                  target_cb,
                  "Antiperiodic A chi - b d psi, QPhiX vs. QDP++");
    }
  }

#endif

  Real betaFactor = Real(0.25);
#if 1
  for (int target_cb = 0; target_cb < 2; ++target_cb) {
    int source_cb = 1 - target_cb;
    QDPIO::cout << "Testing Even Odd Operator" << endl;

    EvenOddCloverOperator<FT, V, S, compress> M(
        gauge_antip.u_packed,
        gauge_antip.clov_packed[target_cb],
        gauge_antip.invclov_packed[source_cb],
        &geom,
        gauge_antip.t_boundary,
        gauge_antip.aniso_fac_s,
        gauge_antip.aniso_fac_t);

    for (int isign = 1; isign >= -1; isign -= 2) {
      hs_qphix1.zero();
      M(hs_qphix1[target_cb], hs_source[target_cb], isign, target_cb);
      hs_qphix1.unpack();

      // Apply QDP Dslash
      QdpSpinor ltmp = zero;
      dslash(hs_qdp1.qdp(), gauge_antip.u_aniso, hs_source.qdp(), isign, source_cb);
      gauge_antip.invclov_qdp.apply(hs_qdp2.qdp(), hs_qdp1.qdp(), isign, source_cb);
      dslash(ltmp, gauge_antip.u_aniso, hs_qdp2.qdp(), isign, target_cb);

      gauge_antip.clov_qdp.apply(hs_qdp1.qdp(), hs_source.qdp(), isign, target_cb);

      hs_qdp1.qdp()[rb[target_cb]] -= betaFactor * ltmp;

      expect_near(hs_qdp1.qdp(),
                  hs_qphix1.qdp(),
                  1e-6,
                  geom,
                  target_cb,
                  "Even-odd operator, QPhiX vs. QDP++");
    } // isign loop
  } // cb loop

#endif

#if 1
  {
    for (int cb = 0; cb < 2; ++cb) {
      masterPrintf("Testing Clover CG on cb = %i\n", cb);
      int other_cb = 1 - cb;
      EvenOddCloverOperator<FT, V, S, compress> M(
          gauge_antip.u_packed,
          gauge_antip.clov_packed[cb],
          gauge_antip.invclov_packed[other_cb],
          &geom,
          gauge_antip.t_boundary,
          gauge_antip.aniso_fac_s,
          gauge_antip.aniso_fac_t);

      double rsd_target = rsdTarget<FT>::value;
      int max_iters = 500;
      int niters;
      double rsd_final;
      unsigned long site_flops;
      unsigned long mv_apps;

      InvCG<FT, V, S, compress> solver(M, max_iters);

      hs_qphix1.zero();
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
      QdpSpinor ltmp = zero;
      dslash(hs_qdp1.qdp(), gauge_antip.u_aniso, hs_qphix1.qdp(), 1, other_cb);
      gauge_antip.invclov_qdp.apply(hs_qdp2.qdp(), hs_qdp1.qdp(), 1, other_cb);
      dslash(ltmp, gauge_antip.u_aniso, hs_qdp2.qdp(), 1, cb);

      gauge_antip.clov_qdp.apply(hs_qdp1.qdp(), hs_qphix1.qdp(), 1, cb);
      hs_qdp1.qdp()[rb[cb]] -= betaFactor * ltmp;

      // chi3 = M^\dagger chi2
      QdpSpinor chi3 = zero;
      dslash(chi3, gauge_antip.u_aniso, hs_qdp1.qdp(), -1, other_cb);
      gauge_antip.invclov_qdp.apply(hs_qdp2.qdp(), chi3, -1, other_cb);
      dslash(ltmp, gauge_antip.u_aniso, hs_qdp2.qdp(), -1, cb);

      gauge_antip.clov_qdp.apply(chi3, hs_qdp1.qdp(), -1, cb);
      chi3[rb[cb]] -= betaFactor * ltmp;

      //  dslash(chi3,u,chi2, (-1), 1);
      // dslash(ltmp,u,chi3, (-1), 0);
      // chi3[rb[0]] = massFactor*chi2 - betaFactor*ltmp;

      expect_near(hs_source.qdp(), chi3, 1e-9, geom, cb, "CG");

      int Nxh = Nx / 2;
      unsigned long num_cb_sites = Nxh * Ny * Nz * Nt;

      unsigned long total_flops =
          (site_flops + (1320 + 504 + 1320 + 504 + 48) * mv_apps) * num_cb_sites;

      masterPrintf("GFLOPS=%e\n", 1.0e-9 * (double)(total_flops) / (end - start));
    }
  }
#endif

#if 1
  {
    for (int cb = 0; cb < 2; ++cb) {
      int other_cb = 1 - cb;
      EvenOddCloverOperator<FT, V, S, compress> M(
          gauge_antip.u_packed,
          gauge_antip.clov_packed[cb],
          gauge_antip.invclov_packed[other_cb],
          &geom,
          gauge_antip.t_boundary,
          gauge_antip.aniso_fac_s,
          gauge_antip.aniso_fac_t);

      double rsd_target = rsdTarget<FT>::value;
      int max_iters = 500;
      int niters;
      double rsd_final;
      unsigned long site_flops;
      unsigned long mv_apps;

      InvBiCGStab<FT, V, S, compress> solver(M, max_iters);
      const int isign = 1;
      hs_qphix1.zero();
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
      QdpSpinor ltmp = zero;
      dslash(hs_qdp1.qdp(), gauge_antip.u_aniso, hs_qphix1.qdp(), 1, other_cb);
      gauge_antip.invclov_qdp.apply(hs_qdp2.qdp(), hs_qdp1.qdp(), 1, other_cb);
      dslash(ltmp, gauge_antip.u_aniso, hs_qdp2.qdp(), 1, cb);

      gauge_antip.clov_qdp.apply(hs_qdp1.qdp(), hs_qphix1.qdp(), 1, cb);
      hs_qdp1.qdp()[rb[cb]] -= betaFactor * ltmp;

      QdpSpinor diff = hs_qdp1.qdp() - hs_source.qdp();
      QDPIO::cout << " cb = " << cb << " True norm is: "
                  << sqrt(norm2(diff, rb[cb]) / norm2(hs_source.qdp(), rb[cb]))
                  << endl;
      expect_near(hs_source.qdp(), hs_qdp1.qdp(), 1e-6, geom, cb, "BiCGStab");

      int Nxh = Nx / 2;
      unsigned long num_cb_sites = Nxh * Ny * Nz * Nt;

      unsigned long total_flops =
          (site_flops + (1320 + 504 + 1320 + 504 + 48) * mv_apps) * num_cb_sites;
      masterPrintf("GFLOPS=%e\n", 1.0e-9 * (double)(total_flops) / (end - start));
    } // cb loop
  } // scope
#endif
}
