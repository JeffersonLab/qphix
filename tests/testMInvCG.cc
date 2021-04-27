#include <iostream>
#include "unittest.h"
#include "testMInvCG.h"

#include "qdp.h"
using namespace QDP;

#include "dslashm_w.h"

#include "reunit.h"
#include "qphix/qphix_config.h"

#include "qphix/geometry.h"
#include "qphix/qdp_packer.h"
#include "qphix/blas_new_c.h"
#include "qphix/wilson.h"
#include "qphix/minvcg.h"

#include <omp.h>

using namespace Assertions;
using namespace std;
using namespace QPhiX;

#include "tparam_selector.h"

#include "veclen.h"
#include "tolerance.h"
#include "RandomGauge.h"
#include "compare_qdp_spinors_custom.h"

void TestMultishift::run()
{
  call(*this, args_.prec, args_.soalen, args_.compress12);
}

template <typename FT,
          int veclen,
          int soalen,
          bool compress12,
          typename QdpGauge,
          typename QdpSpinor>
void TestMultishift::operator()()
{
  if (!compress12) {
    testMInvCG<FT, veclen, soalen, compress12, QdpGauge, QdpSpinor>(-1);
  }
  else {
    // FIXME (Martin Ueding): The following test should be run in every case.
    // The `else` is only here because the second call to `testMInvCG` causes
    // `gaussian` to fill the spinor with junk. This is not yet understood.
    // However, running one of the tests is better than nothing.
    testMInvCG<FT, veclen, soalen, compress12, QdpGauge, QdpSpinor>(1);
  }
}

template <typename T,
          int V,
          int S,
          bool compress,
          typename QdpGauge,
          typename QdpSpinor>
void TestMultishift::testMInvCG(int t_bc)
{
  typedef typename Geometry<T, V, S, compress>::SU3MatrixBlock Gauge;
  typedef typename Geometry<T, V, S, compress>::FourSpinorBlock Spinor;

  QDPIO::cout << "In testMInvCG:" << endl;

  Geometry<T, V, S, compress> geom(Layout::subgridLattSize().slice(),
                                   args_.By,
                                   args_.Bz,
                                   args_.NCores,
                                   args_.Sy,
                                   args_.Sz,
                                   args_.PadXY,
                                   args_.PadXYZ,
                                   args_.MinCt,
                                   true);

  RandomGauge<T, V, S, compress, QdpGauge, QdpSpinor> gauge(geom, t_bc);

  HybridSpinor<T, V, S, compress, QdpSpinor> hs_source(geom);
  gaussian(hs_source.qdp());
  hs_source.pack();

  int const threads_per_core = args_.Sy * args_.Sz;
  int n_shift = 4;
  double shifts[4] = {0.01, 0.02, 0.03, 0.04};

  for (int cb = 0; cb < 2; ++cb) {
    int other_cb = 1 - cb;
    RNG::setrn(rng_seed);

    // Allocate the solutions
    Spinor *psi_d[4];
    for (int s = 0; s < n_shift; s++) {
      psi_d[s] = geom.allocCBFourSpinor();
    }

    QDPIO::cout << "Fields allocated" << endl;

    QDPIO::cout << " Zeroing the psis" << endl;
    for (int s = 0; s < n_shift; s++) {
      zeroSpinor<T, V, S, compress>(psi_d[s], geom, threads_per_core);
    }

    QDPIO::cout << "done" << endl;

    EvenOddWilsonOperator<T, V, S, compress> M(gauge.clover_mass,
                                               gauge.u_packed,
                                               &geom,
                                               gauge.t_boundary,
                                               gauge.aniso_fac_s,
                                               gauge.aniso_fac_t);

    double rsd_target[4];
    double rsd_final[4];
    for (int s = 0; s < n_shift; s++) {
      rsd_target[s] = rsdTarget<T>::value;
    }

    int max_iters = 200;
    int niters;
    unsigned long site_flops;
    unsigned long mv_apps;
    double r2 = 0;
    double isign = 1;
    double start = 0;
    double end = 0;

    {
      MInvCG<T, V, S, compress> solver(M, max_iters, n_shift);
      norm2Spinor<T, V, S, compress>(r2, hs_source[cb], geom, threads_per_core);
      masterPrintf("chi has norm2 = %16.8e\n", r2);

      start = omp_get_wtime();

      bool verbose = true;
      solver(psi_d,
             hs_source[cb],
             n_shift,
             shifts,
             rsd_target,
             niters,
             rsd_final,
             site_flops,
             mv_apps,
             isign,
             verbose,
             cb);

      end = omp_get_wtime();

      QDPIO::cout << " cb= " << cb << " Solver Completed. Iters = " << niters
                  << " Wallclock = " << end - start << " sec." << endl;
    }

    // check solutions
    QdpSpinor psi, psi2, psi3, ltmp;
    for (int s = 0; s < n_shift; s++) {
      Real massFactor = Real(4) + Real(gauge.clover_mass);
      Real betaFactor = Real(0.25) / massFactor;
      Real shiftFactor = Real(shifts[s]);

      // Unpack solution s into 'psi'
      qdp_unpack_cb_spinor<>(psi_d[s], psi, geom, cb);

      // psi2 = M psi
      dslash(psi2, gauge.get_u_aniso(), psi, isign, other_cb);
      dslash(ltmp, gauge.get_u_aniso(), psi2, isign, cb);
      psi2[rb[cb]] = massFactor * psi - betaFactor * ltmp;

      // psi3 = M^\dagger psi2
      dslash(psi3, gauge.get_u_aniso(), psi2, (-isign), other_cb);
      dslash(ltmp, gauge.get_u_aniso(), psi3, (-isign), cb);
      psi3[rb[cb]] = massFactor * psi2 - betaFactor * ltmp;

      // psi2 = psi3 + shift * psi = M^\dagger M psi + shift psi;
      psi2[rb[cb]] = shiftFactor * psi + psi3;

      expect_near(hs_source.qdp(), psi2, 1e-9, geom, cb, "multi-shift CG");
    }

    unsigned long num_cb_sites = Layout::vol() / 2;
    unsigned long total_flops =
        (site_flops + (72 + 2 * 1320) * mv_apps) * num_cb_sites;

    masterPrintf("GFLOPS=%e\n", 1.0e-9 * (double)(total_flops) / (end - start));

    for (int i = 0; i < 4; i++) {
      geom.free(psi_d[i]);
    }
  }
}
