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

int Nx, Ny, Nz, Nt, Nxh;
bool verbose = true;

void TestMultishift::run() {
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
  Nx = args_.nrow_in[0];
  Ny = args_.nrow_in[1];
  Nz = args_.nrow_in[2];
  Nt = args_.nrow_in[3];

  QDPIO::cout << "Inititalizing QDP++ gauge field" << endl;
  // Make a random gauge field
  multi1d<QdpGauge> u(4);
  QdpGauge g;
  QdpGauge uf;
  for (int mu = 0; mu < 4; mu++) {
    uf = 1; // Unit gauge

    Real factor = Real(0.09);
    gaussian(g);
    u[mu] = uf + factor * g;
    reunit(u[mu]);
  }

  testMInvCG<FT, veclen, soalen, compress12, QdpGauge, QdpSpinor>(u, 1);
  if (!compress12) {
    testMInvCG<FT, veclen, soalen, compress12, QdpGauge, QdpSpinor>(u, -1);
  }
}

template <typename T, int V, int S, bool compress, typename U, typename Phi>
void TestMultishift::testMInvCG(const multi1d<U> &u, int t_bc)
{
  for (int cb = 0; cb < 2; ++cb) {
    int other_cb = 1 - cb;
    RNG::setrn(rng_seed);
    typedef typename Geometry<T, V, S, compress>::SU3MatrixBlock Gauge;
    typedef typename Geometry<T, V, S, compress>::FourSpinorBlock Spinor;

    int threads_per_core = args_.Sy * args_.Sz;

    QDPIO::cout << "In testMInvCG:" << endl;

    Phi chi;
    QDPIO::cout << "Filling chi with gaussian noise" << endl;
    gaussian(chi);

    QDPIO::cout << "Norm2 || psi || = " << norm2(chi, rb[cb]) << endl;
    QDPIO::cout << "Done" << endl;
    double aniso_fac_s = (double)(0.35);
    double aniso_fac_t = (double)(1.4);
    double t_boundary = (double)(t_bc);

    Geometry<T, V, S, compress> geom(Layout::subgridLattSize().slice(),
                                     args_.By,
                                     args_.Bz,
                                     args_.NCores,
                                     args_.Sy,
                                     args_.Sz,
                                     args_.PadXY,
                                     args_.PadXYZ,
                                     args_.MinCt, true);

    // NEED TO MOVE ALL THIS INTO DSLASH AT SOME POINT
    // -- Allocate the gauge field
    Gauge *packed_gauge_cb0 = (Gauge *)geom.allocCBGauge();
    Gauge *packed_gauge_cb1 = (Gauge *)geom.allocCBGauge();
    Gauge *u_packed[2] = {packed_gauge_cb0, packed_gauge_cb1};

    // Allocate the right hand side
    Spinor *chi_d = (Spinor *)geom.allocCBFourSpinor();
    if (chi_d == 0x0) {
      QDPIO::cout << "Error allocating chi" << endl;
      QDP_abort(1);
    }

    int n_shift = 4;
    double shifts[4] = {0.01, 0.02, 0.03, 0.04};

    // Allocate the solutions
    Spinor *psi_d[4];
    for (int s = 0; s < n_shift; s++) {
      psi_d[s] = (Spinor *)geom.allocCBFourSpinor();
      if (psi_d[s] == 0x0) {
        QDPIO::cout << "Error allocating psi[" << s << "]" << endl;
        QDP_abort(1);
      }
    }

    QDPIO::cout << "Fields allocated" << endl;

    // Pack the gauge field
    QDPIO::cout << "Packing gauge field...";
    qdp_pack_gauge<>(u, packed_gauge_cb0, packed_gauge_cb1, geom);
    QDPIO::cout << "done" << endl;

    QDPIO::cout << " Packing chi";
    qdp_pack_cb_spinor<>(chi, chi_d, geom, cb);

    QDPIO::cout << " Zeroing the psis";
    for (int s = 0; s < n_shift; s++) {
      zeroSpinor<T, V, S, compress>(psi_d[s], geom, threads_per_core);
    }
    QDPIO::cout << "done" << endl;
    QDPIO::cout << "T BCs = " << t_boundary << endl;

    QDPIO::cout << "Applying anisotropy to test gauge field" << endl;
    multi1d<U> u_test(Nd);
    for (int mu = 0; mu < Nd; mu++) {
      Real factor = Real(aniso_fac_s);
      if (mu == Nd - 1) {
        factor = Real(aniso_fac_t);
      }
      u_test[mu] = factor * u[mu];
    }
    // Apply BCs on u-test for QDP++ test (Dslash gets unmodified field)
    u_test[3] *= where(Layout::latticeCoordinate(3) == (Layout::lattSize()[3] - 1),
                       Real(t_boundary),
                       Real(1));

    double Mass = 0.01;
    EvenOddWilsonOperator<T, V, S, compress> M(
        Mass, u_packed, &geom, t_boundary, aniso_fac_s, aniso_fac_t);
    Phi ltmp = zero;
    Real massFactor = Real(4) + Real(Mass);
    Real betaFactor = Real(0.25) / massFactor;
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
      norm2Spinor<T, V, S, compress>(r2, chi_d, geom, threads_per_core);
      masterPrintf("chi has norm2 = %16.8e\n", r2);

      start = omp_get_wtime();

      solver(psi_d,
             chi_d,
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
    Phi psi, psi2, psi3;
    for (int s = 0; s < n_shift; s++) {
      Real shiftFactor = Real(shifts[s]);
      // Unpack solution s into 'psi'
      qdp_unpack_cb_spinor<>(psi_d[s], psi, geom, cb);

      // psi2 = M psi
      dslash(psi2, u_test, psi, isign, other_cb);
      dslash(ltmp, u_test, psi2, isign, cb);
      psi2[rb[cb]] = massFactor * psi - betaFactor * ltmp;

      // psi3 = M^\dagger psi2
      dslash(psi3, u_test, psi2, (-isign), other_cb);
      dslash(ltmp, u_test, psi3, (-isign), cb);
      psi3[rb[cb]] = massFactor * psi2 - betaFactor * ltmp;

      // psi2 = psi3 + shift * psi = M^\dagger M psi + shift psi;
      psi2[rb[cb]] = shiftFactor * psi + psi3;

      Phi diff = chi - psi2;
      Double true_norm = sqrt(norm2(diff, rb[cb]) / norm2(chi, rb[cb]));
      QDPIO::cout << "cb = " << cb << " True norm for shift[" << s
                  << "]=" << shifts[s] << " is: " << true_norm << endl;
      //    assertion( toBool( true_norm < (rsd_target +
      //    tolerance<T>::small) ) );
    }

    unsigned long num_cb_sites = Layout::vol() / 2;
    unsigned long total_flops =
        (site_flops + (72 + 2 * 1320) * mv_apps) * num_cb_sites;

    masterPrintf("GFLOPS=%e\n", 1.0e-9 * (double)(total_flops) / (end - start));

#if 1

    geom.free(packed_gauge_cb0);
    geom.free(packed_gauge_cb1);
    geom.free(chi_d);
    for (int i = 0; i < 4; i++) {
      geom.free(psi_d[i]);
    }

#endif
  }
}
