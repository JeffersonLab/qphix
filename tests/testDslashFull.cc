#include <iostream>
using namespace std;

#include "unittest.h"
#include "testDslashFull.h"
#include "qdp.h"
using namespace QDP;

#include "dslashm_w.h"
#include "reunit.h"

#include "qphix/geometry.h"
#include "qphix/qdp_packer.h"
#include "qphix/blas_new_c.h"
#include "qphix/wilson.h"
#include "qphix/invcg.h"
#include "qphix/invbicgstab.h"
#include "qphix/inv_richardson_multiprec.h"
#include "./invbicgstab_test.h"

#include <omp.h>

using namespace Assertions;
using namespace std;
using namespace QPhiX;

#include "veclen.h"
#include "tolerance.h"
#include "tparam_selector.h"
#include "compare_qdp_spinors_custom.h"
#include "RandomGauge.h"

int Nx, Ny, Nz, Nt, Nxh;
bool verbose = true;

void TestDslash::run(void)
{
  RNG::savern(rng_seed);

  call(*this, args_.prec, args_.soalen, args_.compress12);
}

template <typename FT,
          int veclen,
          int soalen,
          bool compress12,
          typename QdpGauge,
          typename QdpSpinor>
void TestDslash::operator()()
{
  RNG::savern(rng_seed);

  QDPIO::cout << "Inititalizing QDP++ gauge field" << endl;
  // Make a random gauge field
  multi1d<QdpGauge> u(4);
  QdpGauge g;
  QdpGauge uf;
  for (int mu = 0; mu < 4; mu++) {
#if 1
    uf = 1; // Unit gauge

    Real factor = Real(0.09);
    gaussian(g);
    u[mu] = uf + factor * g;
    reunit(u[mu]);
#else
    u[mu] = 1;
#endif
  }

  for (int const t_bc : {1, -1}) {
    testDslash<FT, veclen, soalen, compress12, QdpGauge, QdpSpinor>(u, t_bc);
    testDslashAChiMBDPsi<FT, veclen, soalen, compress12, QdpGauge, QdpSpinor>(u,
                                                                              t_bc);
    testM<FT, veclen, soalen, compress12, QdpGauge, QdpSpinor>(u, t_bc);
    testCG<FT, veclen, soalen, compress12, QdpGauge, QdpSpinor>(u, t_bc);
    testBiCGStab<FT, veclen, soalen, compress12, QdpGauge, QdpSpinor>(u, t_bc);
    
    // unlike the other tests, we run these with fixed soalen
#if defined(QPHIX_MIC_SOURCE) || defined(QPHIX_AVX512_SOURCE)
    if (Nx % 32 == 0) {
      testRichardson<double, VECLEN_DP, 8, true, half, VECLEN_HP, 16, QdpGauge, QdpSpinor>(
          u, t_bc);
    } else {
      if (Nx % 16 == 0) {
        testRichardson<double, VECLEN_DP, 8, true, half, VECLEN_HP, 8, QdpGauge, QdpSpinor>(
            u, t_bc);
      } else {
        masterPrintf("I havent set up that mixed precision solver combination\n");
      }
    }
#elif defined(QPHIX_AVX_SOURCE) || defined(QPHIX_AVX2_SOURCE)
    // AVX: Double SOALEN = 4
    if (Nx % 16 == 0) {
      testRichardson<double, VECLEN_DP, 4, true, float, VECLEN_SP, 8, QdpGauge, QdpSpinor>(
          u, t_bc);
    } else {
      if (Nx % 8 == 0) {
        testRichardson<double, VECLEN_DP, 4, true, float, VECLEN_SP, 4, QdpGauge, QdpSpinor>(
            u, t_bc);
      } else {
        masterPrintf("I havent set up that mixed precision solver combination\n");
      }
    }
#endif
  }
}

template <typename T, int V, int S, bool compress, typename U, typename Phi>
void TestDslash::testDslash(const multi1d<U> &u, int t_bc)
{
  QDPIO::cout << "RNG seeed = " << rng_seed << std::endl;
  RNG::setrn(rng_seed);

  typedef typename Geometry<T, V, S, compress>::SU3MatrixBlock Gauge;
  typedef typename Geometry<T, V, S, compress>::FourSpinorBlock Spinor;

  QDPIO::cout << "In testDslash" << endl;
  double aniso_fac_s = ((double)0.35);
  double aniso_fac_t = ((double)1.4);
  double t_boundary = ((double)t_bc);
   QDPIO::cout << "Testing Dslash...  T BCs = " << t_boundary << endl;

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

  RandomGauge<T, V, S, compress, U, Phi> gauge(geom, t_bc);

  Dslash<T, V, S, compress> D32(
      &geom, gauge.t_boundary, gauge.aniso_fac_s, gauge.aniso_fac_t);

  HybridSpinor<T, V, S, compress, Phi> hs_source(geom), hs_qphix1(geom),
      hs_qphix2(geom), hs_qdp1(geom), hs_qdp2(geom);
  gaussian(hs_source.qdp());
  hs_source.pack();

  int isign = 1;
  for (int isign = 1; isign >= -1; isign -= 2) {
    for (int cb = 0; cb < 2; cb++) {
      int source_cb = 1 - cb;
      int target_cb = cb;

      D32.dslash(hs_qphix1[target_cb],
                 hs_source[source_cb],
                 gauge.u_packed[target_cb],
                 isign,
                 target_cb);
      hs_qphix1.unpack();

      dslash(hs_qdp1.qdp(), gauge.u_aniso, hs_source.qdp(), isign, target_cb);

      expect_near(hs_qdp1.qdp(),
                  hs_qphix1.qdp(),
                  tolerance<T>::value,
                  geom,
                  target_cb,
                  "Wilson::Dslash");
    } // cb
  } // isign
}

template <typename T, int V, int S, bool compress, typename U, typename Phi>
void TestDslash::testDslashAChiMBDPsi(const multi1d<U> &u, int t_bc)
{
  RNG::setrn(rng_seed);
  typedef typename Geometry<T, V, S, compress>::SU3MatrixBlock Gauge;
  typedef typename Geometry<T, V, S, compress>::FourSpinorBlock Spinor;

  QDPIO::cout << "In testAChiMBDPsi " << endl;
  double aniso_fac_s = (double)(0.35);
  double aniso_fac_t = (double)(1.4);
  double t_boundary = (double)(t_bc);

  Phi psi, chi, chi2;
  QDPIO::cout << "Filling psi with random noise" << endl;
  gaussian(psi);

  Geometry<T, V, S, compress> geom(Layout::subgridLattSize().slice(),
                                   args_.By,
                                   args_.Bz,
                                   args_.NCores,
                                   args_.Sy,
                                   args_.Sz,
                                   args_.PadXY,
                                   args_.PadXYZ,
                                   args_.MinCt);

  // NEED TO MOVE ALL THIS INTO DSLASH AT SOME POINT
  Dslash<T, V, S, compress> D32(&geom, t_boundary, aniso_fac_s, aniso_fac_t);

  Gauge *packed_gauge_cb0 = (Gauge *)geom.allocCBGauge();
  Gauge *packed_gauge_cb1 = (Gauge *)geom.allocCBGauge();

  // Over allocate, so that an unsigned load doesnt cause segfault accesing
  // either end...
  Spinor *psi_even = (Spinor *)geom.allocCBFourSpinor();
  Spinor *psi_odd = (Spinor *)geom.allocCBFourSpinor();
  Spinor *chi_even = (Spinor *)geom.allocCBFourSpinor();
  Spinor *chi_odd = (Spinor *)geom.allocCBFourSpinor();

  QDPIO::cout << "Fields allocated" << endl;

  // Pack the gauge field
  QDPIO::cout << "Packing gauge field...";
  qdp_pack_gauge<>(u, packed_gauge_cb0, packed_gauge_cb1, geom);

  QDPIO::cout << "done" << endl;
  Gauge *u_packed[2];
  u_packed[0] = packed_gauge_cb0;
  u_packed[1] = packed_gauge_cb1;

  Spinor *psi_s[2] = {psi_even, psi_odd};
  Spinor *chi_s[2] = {chi_even, chi_odd};

  QDPIO::cout << " Packing fermions...";
  qdp_pack_spinor<>(psi, psi_even, psi_odd, geom);

  QDPIO::cout << "done" << endl;

  QDPIO::cout << "Testing DslashAChiMBDPsi \n" << endl;
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

  for (int isign = 1; isign >= -1; isign -= 2) {
    for (int cb = 0; cb < 2; cb++) {
      int source_cb = 1 - cb;
      int target_cb = cb;
      chi = zero;
      qdp_pack_spinor<>(chi, chi_even, chi_odd, geom);

      double alpha = (double)(4.01); // Nd + M, M=0.01
      double beta = (double)(0.5); // Operator is (Nd+M) - (1/2)Dslash

      // Apply Optimized Dslash
      D32.dslashAChiMinusBDPsi(chi_s[target_cb],
                               psi_s[source_cb],
                               psi_s[target_cb],
                               u_packed[target_cb],
                               alpha,
                               beta,
                               isign,
                               target_cb);

      qdp_unpack_spinor<>(chi_s[0], chi_s[1], chi, geom);

      // Apply QDP Dslash
      chi2 = zero;
      dslash(chi2, u_test, psi, isign, target_cb);
      Phi res;
      res[rb[source_cb]] = zero;
      res[rb[target_cb]] = alpha * psi - beta * chi2;

      // Check the difference per number in chi vector
      Phi diff = res - chi;

      Double diff_norm = sqrt(norm2(diff, rb[target_cb])) /
                         (Real(4 * 3 * 2 * Layout::vol()) / Real(2));

      QDPIO::cout << "\t cb = " << target_cb << "  isign = " << isign
                  << "  diff_norm = " << diff_norm << endl;
      // Assert things are OK...
      if (toBool(diff_norm >= tolerance<T>::small)) {
        for (int i = 0; i < rb[target_cb].siteTable().size(); i++) {
          for (int s = 0; s < Ns; s++) {
            for (int c = 0; c < Nc; c++) {
              QDPIO::cout << "site=" << i << " spin=" << s << " color=" << c
                          << " Diff = "
                          << diff.elem(rb[target_cb].start() + i).elem(s).elem(c)
                          << endl;
            }
          }
        }
      }
      assertion(toBool(diff_norm < tolerance<T>::small));
    }
  }

  geom.free(packed_gauge_cb0);
  geom.free(packed_gauge_cb1);
  geom.free(psi_even);
  geom.free(psi_odd);
  geom.free(chi_even);
  geom.free(chi_odd);
}

template <typename T, int V, int S, bool compress, typename U, typename Phi>
void TestDslash::testM(const multi1d<U> &u, int t_bc)
{
  RNG::setrn(rng_seed);
  typedef typename Geometry<T, V, S, compress>::SU3MatrixBlock Gauge;
  typedef typename Geometry<T, V, S, compress>::FourSpinorBlock Spinor;

  QDPIO::cout << "in testM:" << endl;
  double aniso_fac_s = (double)(0.35);
  double aniso_fac_t = (double)(1.4);
  double t_boundary = (double)(t_bc);

  Phi psi, chi, chi2;
  QDPIO::cout << "Filling source spinor with gaussian noise" << endl;
  gaussian(psi);

  Geometry<T, V, S, compress> geom(Layout::subgridLattSize().slice(),
                                   args_.By,
                                   args_.Bz,
                                   args_.NCores,
                                   args_.Sy,
                                   args_.Sz,
                                   args_.PadXY,
                                   args_.PadXYZ,
                                   args_.MinCt);

  Gauge *packed_gauge_cb0 = (Gauge *)geom.allocCBGauge();
  Gauge *packed_gauge_cb1 = (Gauge *)geom.allocCBGauge();

  // Over allocate, so that an unsigned load doesnt cause segfault accesing
  // either end...
  Spinor *psi_even = (Spinor *)geom.allocCBFourSpinor();
  Spinor *psi_odd = (Spinor *)geom.allocCBFourSpinor();
  Spinor *chi_even = (Spinor *)geom.allocCBFourSpinor();
  Spinor *chi_odd = (Spinor *)geom.allocCBFourSpinor();

  QDPIO::cout << "Fields allocated" << endl;

  // Pack the gauge field
  QDPIO::cout << "Packing gauge field...";
  qdp_pack_gauge<>(u, packed_gauge_cb0, packed_gauge_cb1, geom);

  QDPIO::cout << "done" << endl;
  Gauge *u_packed[2];
  u_packed[0] = packed_gauge_cb0;
  u_packed[1] = packed_gauge_cb1;

  Spinor *psi_s[2] = {psi_even, psi_odd};
  Spinor *chi_s[2] = {chi_even, chi_odd};

  QDPIO::cout << " Packing fermions...";
  qdp_pack_spinor<>(psi, psi_even, psi_odd, geom);

  QDPIO::cout << "done" << endl;

  QDPIO::cout << "Testing M \n" << endl;
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

  double Mass = 0.1;
  EvenOddWilsonOperator<T, V, S, compress> M(
      Mass, u_packed, &geom, t_boundary, aniso_fac_s, aniso_fac_t);

  Phi ltmp = zero;
  Real massFactor = Real(4) + Real(Mass);
  Real betaFactor = Real(0.25) / massFactor;
  // Apply optimized
  for (int isign = 1; isign >= -1; isign -= 2) {
    for (int cb = 0; cb < 2; ++cb) {
      int other_cb = 1 - cb;
      chi = zero;
      qdp_pack_spinor<>(chi, chi_even, chi_odd, geom);

      M(chi_s[cb], psi_s[cb], isign, cb);

      qdp_unpack_spinor<>(chi_s[0], chi_s[1], chi, geom);

      // Apply QDP Dslash
      chi2 = zero;
      dslash(chi2, u_test, psi, isign, other_cb);
      dslash(ltmp, u_test, chi2, isign, cb);
      chi2[rb[cb]] = massFactor * psi - betaFactor * ltmp;

      // Check the difference per number in chi vector
      Phi diff = zero;
      diff[rb[cb]] = chi2 - chi;

      Double diff_norm =
          sqrt(norm2(diff, rb[cb])) / (Real(4 * 3 * 2 * Layout::vol()) / Real(2));

      QDPIO::cout << "\t cb= " << cb << " isign = " << isign
                  << "  diff_norm = " << diff_norm << endl;
      // Assert things are OK...
      if (toBool(diff_norm >= tolerance<T>::small)) {
        for (int i = 0; i < rb[cb].siteTable().size(); i++) {
          for (int s = 0; s < Ns; s++) {
            for (int c = 0; c < Nc; c++) {
              Double re =
                  Double(diff.elem(rb[cb].start() + i).elem(s).elem(c).real());
              Double im =
                  Double(diff.elem(rb[cb].start() + i).elem(s).elem(c).imag());
              if (toBool(fabs(re) > tolerance<T>::small) ||
                  toBool(fabs(im) > tolerance<T>::small)) {
                QDPIO::cout << "site=" << i << " spin=" << s << " color=" << c
                            << " Diff = "
                            << diff.elem(rb[cb].start() + i).elem(s).elem(c) << endl;
              }
            }
          }
        }
      }
      assertion(toBool(diff_norm < tolerance<T>::small));

    } // cb loop
  } // isign loop

  geom.free(packed_gauge_cb0);
  geom.free(packed_gauge_cb1);
  geom.free(psi_even);
  geom.free(psi_odd);
  geom.free(chi_even);
  geom.free(chi_odd);
}

template <typename T, int V, int S, bool compress, typename U, typename Phi>
void TestDslash::testCG(const multi1d<U> &u, int t_bc)
{
  for (int cb = 0; cb < 2; ++cb) {
    int other_cb = 1 - cb;
    RNG::setrn(rng_seed);
    typedef typename Geometry<T, V, S, compress>::SU3MatrixBlock Gauge;
    typedef typename Geometry<T, V, S, compress>::FourSpinorBlock Spinor;

    QDPIO::cout << "In testCG:" << endl;

    Phi psi, chi, chi2, chi3;
    QDPIO::cout << "Filling psi with gaussian noise" << endl;

    gaussian(psi);
    QDPIO::cout << "Norm2 || psi || = " << norm2(psi, rb[cb]) << endl;
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
                                     args_.MinCt);

    // NEED TO MOVE ALL THIS INTO DSLASH AT SOME POINT
    Gauge *packed_gauge_cb0 = (Gauge *)geom.allocCBGauge();
    Gauge *packed_gauge_cb1 = (Gauge *)geom.allocCBGauge();

    // Over allocate, so that an unsigned load doesnt cause segfault
    // accesing either end...
    Spinor *psi_even = (Spinor *)geom.allocCBFourSpinor();
    Spinor *psi_odd = (Spinor *)geom.allocCBFourSpinor();
    Spinor *chi_even = (Spinor *)geom.allocCBFourSpinor();
    Spinor *chi_odd = (Spinor *)geom.allocCBFourSpinor();

    QDPIO::cout << "Fields allocated" << endl;

    // Pack the gauge field
    QDPIO::cout << "Packing gauge field...";
    //  qdp_pack_gauge< T,V,S,compress, U >(u,
    //  packed_gauge_cb0,packed_gauge_cb1, geom);
    qdp_pack_gauge<>(u, packed_gauge_cb0, packed_gauge_cb1, geom);

    QDPIO::cout << "done" << endl;
    Gauge *u_packed[2];
    u_packed[0] = packed_gauge_cb0;
    u_packed[1] = packed_gauge_cb1;

    Spinor *psi_s[2] = {psi_even, psi_odd};
    Spinor *chi_s[2] = {chi_even, chi_odd};

    QDPIO::cout << " Packing fermions...";

    qdp_pack_spinor<>(psi, psi_even, psi_odd, geom);

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

    double Mass = 0.1;
    EvenOddWilsonOperator<T, V, S, compress> M(
        Mass, u_packed, &geom, t_boundary, aniso_fac_s, aniso_fac_t);
    Phi ltmp = zero;
    Real massFactor = Real(4) + Real(Mass);
    Real betaFactor = Real(0.25) / massFactor;

    {
      chi = zero;
      qdp_pack_cb_spinor<>(chi, chi_even, geom, cb);

      double rsd_target = rsdTarget<T>::value;
      int max_iters = 200;
      int niters;
      double rsd_final;
      unsigned long site_flops;
      unsigned long mv_apps;

      InvCG<T, V, S, compress> solver(M, max_iters);
      double r2;
      norm2Spinor<T, V, S, compress>(r2, psi_even, geom, 1);
      masterPrintf("psi has norm2 = %16.8e\n", r2);
      int isign = 1;
      double start = omp_get_wtime();
      solver(chi_s[cb],
             psi_s[cb],
             rsd_target,
             niters,
             rsd_final,
             site_flops,
             mv_apps,
             isign,
             verbose,
             cb);
      double end = omp_get_wtime();

      qdp_unpack_cb_spinor<>(chi_s[cb], chi, geom, cb);

      // Multiply back
      // chi2 = M chi
      int other_cb = 1 - cb;
      dslash(chi2, u_test, chi, isign, other_cb);
      dslash(ltmp, u_test, chi2, isign, cb);
      chi2[rb[cb]] = massFactor * chi - betaFactor * ltmp;

      // chi3 = M^\dagger chi2
      dslash(chi3, u_test, chi2, (-isign), other_cb);
      dslash(ltmp, u_test, chi3, (-isign), cb);
      chi3[rb[cb]] = massFactor * chi2 - betaFactor * ltmp;
      
      Phi diff = chi3 - psi;
      Double true_norm = sqrt(norm2(diff, rb[cb]) / norm2(psi, rb[cb]));
      QDPIO::cout << "Wilson CG Solve isign=" << isign
                  << " True norm is: " << true_norm << endl;

      expect_near(psi, chi3, 1e-9, geom, cb, "Wilson CG");

      unsigned long num_cb_sites = Layout::vol() / 2;
      unsigned long total_flops =
          (site_flops + (72 + 2 * 1320) * mv_apps) * num_cb_sites;

      masterPrintf("GFLOPS=%e\n", 1.0e-9 * (double)(total_flops) / (end - start));
    }

    geom.free(packed_gauge_cb0);
    geom.free(packed_gauge_cb1);
    geom.free(psi_even);
    geom.free(psi_odd);
    geom.free(chi_even);
    geom.free(chi_odd);
  } // cb loop
}

template <typename T, int V, int S, bool compress, typename U, typename Phi>
void TestDslash::testBiCGStab(const multi1d<U> &u, int t_bc)
{
  for (int cb = 0; cb < 2; ++cb) {
    int other_cb = 1 - cb;
    RNG::setrn(rng_seed);
    typedef typename Geometry<T, V, S, compress>::SU3MatrixBlock Gauge;
    typedef typename Geometry<T, V, S, compress>::FourSpinorBlock Spinor;

    QDPIO::cout << "In testBiCGStab" << endl;

    //  double aniso_fac_s = (double)(0.35);
    // double aniso_fac_t = (double)(1.4);
    double aniso_fac_s = (double)(1);
    double aniso_fac_t = (double)(1);

    double t_boundary = (double)(t_bc);

    Phi psi, chi, chi2;
    gaussian(psi);

    Geometry<T, V, S, compress> geom(Layout::subgridLattSize().slice(),
                                     args_.By,
                                     args_.Bz,
                                     args_.NCores,
                                     args_.Sy,
                                     args_.Sz,
                                     args_.PadXY,
                                     args_.PadXYZ,
                                     args_.MinCt);

    // NEED TO MOVE ALL THIS INTO DSLASH AT SOME POINT

    Gauge *packed_gauge_cb0 = (Gauge *)geom.allocCBGauge();
    Gauge *packed_gauge_cb1 = (Gauge *)geom.allocCBGauge();

    // Over allocate, so that an unsigned load doesnt cause segfault
    // accesing either end...
    Spinor *psi_even = (Spinor *)geom.allocCBFourSpinor();
    Spinor *psi_odd = (Spinor *)geom.allocCBFourSpinor();
    Spinor *chi_even = (Spinor *)geom.allocCBFourSpinor();
    Spinor *chi_odd = (Spinor *)geom.allocCBFourSpinor();

    QDPIO::cout << "Fields allocated" << endl;

    // Pack the gauge field
    QDPIO::cout << "Packing gauge field...";
    qdp_pack_gauge<>(u, packed_gauge_cb0, packed_gauge_cb1, geom);

    QDPIO::cout << "done" << endl;
    Gauge *u_packed[2];
    u_packed[0] = packed_gauge_cb0;
    u_packed[1] = packed_gauge_cb1;

    Spinor *psi_s[2] = {psi_even, psi_odd};
    Spinor *chi_s[2] = {chi_even, chi_odd};

    QDPIO::cout << " Packing fermions...";
    qdp_pack_spinor<>(psi, psi_even, psi_odd, geom);

    QDPIO::cout << "done" << endl;
    QDPIO::cout << "Applying anisotropy to test gauge field" << endl;
    multi1d<U> u_test(Nd);
    for (int mu = 0; mu < Nd; mu++) {
      Real factor = Real(aniso_fac_s);
      if (mu == Nd - 1) {
        factor = Real(aniso_fac_t);
      }
      u_test[mu] = factor * u[mu];
    }
    u_test[3] *= where(Layout::latticeCoordinate(3) == (Layout::lattSize()[3] - 1),
                       Real(t_boundary),
                       Real(1));

    double Mass = 0.1;
    EvenOddWilsonOperator<T, V, S, compress> M(
        Mass, u_packed, &geom, t_boundary, aniso_fac_s, aniso_fac_t);
    Phi ltmp = zero;
    Real massFactor = Real(4) + Real(Mass);
    Real betaFactor = Real(0.25) / massFactor;

    {
      float rsd_target = rsdTarget<T>::value;
      int max_iters = 100;
      int niters;
      double rsd_final;
      unsigned long site_flops;
      unsigned long mv_apps;

      InvBiCGStab<T, V, S, compress> solver(M, max_iters);

      for (int isign = 1; isign >= -1; isign -= 2) {
        chi = zero;
        qdp_pack_cb_spinor<>(chi, chi_even, geom, cb);
        masterPrintf("Entering solver\n");

        double start = omp_get_wtime();
        solver(chi_s[cb],
               psi_s[cb],
               rsd_target,
               niters,
               rsd_final,
               site_flops,
               mv_apps,
               isign,
               verbose,
               cb);
        double end = omp_get_wtime();

        qdp_unpack_cb_spinor<>(chi_s[cb], chi, geom, cb);

        // Multiply back
        // chi2 = M chi
        dslash(chi2, u_test, chi, isign, other_cb);
        dslash(ltmp, u_test, chi2, isign, cb);
        chi2[rb[cb]] = massFactor * chi - betaFactor * ltmp;

        Phi diff = chi2 - psi;
        Double true_norm = sqrt(norm2(diff, rb[cb]) / norm2(psi, rb[cb]));
        QDPIO::cout << "BiCGStab Solve isign=" << isign
                    << " True norm is: " << true_norm << endl;
        expect_near(psi, chi2, 1e-9, geom, cb, "BiCGstab: ");
        assertion(toBool(true_norm < (rsd_target + tolerance<T>::small)));
        int Nxh = Nx / 2;
        unsigned long num_cb_sites = Layout::vol() / 2;
        unsigned long total_flops =
            (site_flops + (72 + 2 * 1320) * mv_apps) * num_cb_sites;
        masterPrintf("BICGSTAB Solve isign=%d, args_.iters=%d MV apps=%lu "
                     "Site flops=%lu\n",
                     isign,
                     niters,
                     mv_apps,
                     site_flops);
        masterPrintf("BICGSTAB Solve isign=%d, time=%e (sec) GFLOPS=%e\n",
                     isign,
                     (end - start),
                     1.0e-9 * (double)(total_flops) / (end - start));
      }
    }

    geom.free(packed_gauge_cb0);
    geom.free(packed_gauge_cb1);
    geom.free(psi_even);
    geom.free(psi_odd);
    geom.free(chi_even);
    geom.free(chi_odd);
  } // cb loop
}

template <typename T1,
          int VEC1,
          int SOA1,
          bool compress,
          typename T2,
          int VEC2,
          int SOA2,
          typename U,
          typename Phi>
void TestDslash::testRichardson(const multi1d<U> &u, int t_bc)
{

  for (int cb = 0; cb < 2; ++cb) {
    int other_cb = 1 - cb;

    RNG::setrn(rng_seed);
    typedef typename Geometry<T1, VEC1, SOA1, compress>::SU3MatrixBlock GaugeOuter;
    typedef typename Geometry<T1, VEC1, SOA1, compress>::FourSpinorBlock SpinorOuter;

    typedef typename Geometry<T2, VEC2, SOA2, compress>::SU3MatrixBlock GaugeInner;
    typedef typename Geometry<T2, VEC2, SOA2, compress>::FourSpinorBlock SpinorInner;

    QDPIO::cout << "In testRichardson" << endl;

    double aniso_fac_s = (double)(1);
    double aniso_fac_t = (double)(1);
    double t_boundary = (double)(t_bc);

    Phi psi, chi, chi2;
    gaussian(psi);

    Geometry<T1, VEC1, SOA1, compress> geom_outer(Layout::subgridLattSize().slice(),
                                                  args_.By,
                                                  args_.Bz,
                                                  args_.NCores,
                                                  args_.Sy,
                                                  args_.Sz,
                                                  args_.PadXY,
                                                  args_.PadXYZ,
                                                  args_.MinCt);

    Geometry<T2, VEC2, SOA2, compress> geom_inner(Layout::subgridLattSize().slice(),
                                                  args_.By,
                                                  args_.Bz,
                                                  args_.NCores,
                                                  args_.Sy,
                                                  args_.Sz,
                                                  args_.PadXY,
                                                  args_.PadXYZ,
                                                  args_.MinCt);

    // NEED TO MOVE ALL THIS INTO DSLASH AT SOME POINT

    GaugeOuter *packed_gauge_cb0 = (GaugeOuter *)geom_outer.allocCBGauge();
    GaugeOuter *packed_gauge_cb1 = (GaugeOuter *)geom_outer.allocCBGauge();

    QDPIO::cout << "Packing gauge field...";
    //  qdp_pack_gauge< T,V,S,compress, U >(u,
    //  packed_gauge_cb0,packed_gauge_cb1, geom);
    qdp_pack_gauge<>(u, packed_gauge_cb0, packed_gauge_cb1, geom_outer);

    GaugeOuter *u_packed[2];
    u_packed[0] = packed_gauge_cb0;
    u_packed[1] = packed_gauge_cb1;
    QDPIO::cout << "done" << endl;

    GaugeInner *packed_gauge_cb0_inner = (GaugeInner *)geom_inner.allocCBGauge();
    GaugeInner *packed_gauge_cb1_inner = (GaugeInner *)geom_inner.allocCBGauge();

    QDPIO::cout << "Packing inner gauge field...";
    qdp_pack_gauge<>(u, packed_gauge_cb0_inner, packed_gauge_cb1_inner, geom_inner);
    GaugeInner *u_packed_inner[2];
    u_packed_inner[0] = packed_gauge_cb0_inner;
    u_packed_inner[1] = packed_gauge_cb1_inner;
    QDPIO::cout << "done" << endl;

    // Over allocate, so that an unsigned load doesnt cause segfault
    // accesing either end...
    SpinorOuter *psi_even = (SpinorOuter *)geom_outer.allocCBFourSpinor();
    SpinorOuter *psi_odd = (SpinorOuter *)geom_outer.allocCBFourSpinor();
    SpinorOuter *chi_even = (SpinorOuter *)geom_outer.allocCBFourSpinor();
    SpinorOuter *chi_odd = (SpinorOuter *)geom_outer.allocCBFourSpinor();

    QDPIO::cout << "Fields allocated" << endl;

    SpinorOuter *psi_s[2] = {psi_even, psi_odd};
    SpinorOuter *chi_s[2] = {chi_even, chi_odd};

    QDPIO::cout << " Packing fermions...";
    qdp_pack_spinor<>(psi, psi_even, psi_odd, geom_outer);

    QDPIO::cout << "done" << endl;
    QDPIO::cout << "Applying anisotropy to test gauge field" << endl;
    multi1d<U> u_test(Nd);
    for (int mu = 0; mu < Nd; mu++) {
      Real factor = Real(aniso_fac_s);
      if (mu == Nd - 1) {
        factor = Real(aniso_fac_t);
      }
      u_test[mu] = factor * u[mu];
    }
    u_test[3] *= where(Layout::latticeCoordinate(3) == (Layout::lattSize()[3] - 1),
                       Real(t_boundary),
                       Real(1));

    double Mass = 0.001;
    EvenOddWilsonOperator<T1, VEC1, SOA1, compress> M_outer(
        Mass, u_packed, &geom_outer, t_boundary, aniso_fac_s, aniso_fac_t);

    EvenOddWilsonOperator<T2, VEC2, SOA2, compress> M_inner(
        Mass, u_packed_inner, &geom_inner, t_boundary, aniso_fac_s, aniso_fac_t);

    Phi ltmp = zero;
    Real massFactor = Real(4) + Real(Mass);
    Real betaFactor = Real(0.25) / massFactor;

    {
      float rsd_target = rsdTarget<T1>::value;
      int max_iters = 200;
      int niters;
      double rsd_final;
      unsigned long site_flops;
      unsigned long mv_apps;

      InvBiCGStab<T2, VEC2, SOA2, compress> inner_solver(M_inner, max_iters);
      InvRichardsonMultiPrec<T1, VEC1, SOA1, compress, T2, VEC2, SOA2, compress>
          outer_solver(M_outer, inner_solver, 0.01, max_iters);

      for (int isign = 1; isign >= -1; isign -= 2) {
        // BiCGStab Inner SOlver

        chi = zero;
        qdp_pack_cb_spinor<>(chi, chi_even, geom_outer, 0);
        masterPrintf("Entering solver\n");

        double start = omp_get_wtime();
        outer_solver(chi_s[cb],
                     psi_s[cb],
                     rsd_target,
                     niters,
                     rsd_final,
                     site_flops,
                     mv_apps,
                     isign,
                     verbose,
                     cb);
        double end = omp_get_wtime();

        qdp_unpack_cb_spinor<>(chi_s[cb], chi, geom_outer, cb);

        // Multiply back
        // chi2 = M chi
        dslash(chi2, u_test, chi, isign, other_cb);
        dslash(ltmp, u_test, chi2, isign, cb);
        chi2[rb[cb]] = massFactor * chi - betaFactor * ltmp;

        Phi diff = chi2 - psi;
        Double true_norm = sqrt(norm2(diff, rb[cb]) / norm2(psi, rb[cb]));
        QDPIO::cout << "RICHARDSON Solve isign=" << isign
                    << " True norm is: " << true_norm << endl;
        expect_near(psi, chi2, 1e-9, geom_outer, cb, "BiCGstab: ");

        assertion(toBool(true_norm < (rsd_target + tolerance<T1>::small)));
        unsigned long num_cb_sites = Layout::vol() / 2;
        unsigned long total_flops =
            (site_flops + (72 + 2 * 1320) * mv_apps) * num_cb_sites;

        masterPrintf("RICHARDSON Solve isign=%d time=%e (sec) GFLOPS=%e\n",
                     isign,
                     (end - start),
                     1.0e-9 * (double)(total_flops) / (end - start));
      }
    }

    geom_outer.free(packed_gauge_cb0);
    geom_outer.free(packed_gauge_cb1);

    geom_inner.free(packed_gauge_cb0_inner);
    geom_inner.free(packed_gauge_cb1_inner);

    geom_outer.free(psi_even);
    geom_outer.free(psi_odd);
    geom_outer.free(chi_even);
    geom_outer.free(chi_odd);
  } // cb loop.
}
