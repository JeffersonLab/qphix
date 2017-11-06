// TODO Make an `A chi - b D psi` function for QDP spinors.

#include "unittest.h"
#include "testTWMDslashFull.h"
#include "qdp.h"
using namespace QDP;

#include "dslashm_w.h"
#include "reunit.h"

#include "qphix/geometry.h"
#include "qphix/qdp_packer.h"
#include "qphix/blas_new_c.h"
#include "qphix/twisted_mass.h"
#include "qphix/invcg.h"
#include "qphix/invbicgstab.h"
#include "qphix/inv_richardson_multiprec.h"
#if 0
#include "./invbicgstab_test.h"
#endif

#include <omp.h>
#include <complex>
#if 0
#include "qphix/memmap.h"
#endif

using namespace Assertions;
using namespace std;
using namespace QPhiX;

#include "RandomGauge.h"
#include "veclen.h"
#include "tolerance.h"
#include "tparam_selector.h"
#include "compare_qdp_spinors_custom.h"

int Nx, Ny, Nz, Nt, Nxh;
bool verbose = true;

void TestTMDslash::run(void)
{
  RNG::savern(rng_seed);

  const multi1d<int> &lattSize = Layout::subgridLattSize();
  Nx = lattSize[0];
  Ny = lattSize[1];
  Nz = lattSize[2];
  Nt = lattSize[3];

  call(*this, args_.prec, args_.soalen, args_.compress12);

  /*

  typedef multi1d<UF> MUF;
  typedef multi1d<UD> MUD;

  // Make a random gauge field
  QDPIO::cout << "Inititalizing QDP++ gauge field" << endl;
  multi1d<LatticeColorMatrix> u(4);
  LatticeColorMatrix g;
  LatticeColorMatrix uf;
  for (int mu = 0; mu < 4; mu++) {
    uf = 1; // Unit gauge
    Real factor = Real(0.09);
    gaussian(g);
    u[mu] = uf + factor * g;
    reunit(u[mu]);
  }

  // Save build time
  if (precision == FLOAT_PREC) {

    multi1d<LatticeColorMatrixF> u_in(4);
    for (int mu = 0; mu < Nd; mu++) {
      u_in[mu] = u[mu];
    }

#if defined(QPHIX_SCALAR_SOURCE)
    if (soalen == 1) {
      testTWMBiCGStabWrapper<float, VECLEN_SP, 1, MUF, PhiF>(u_in);
    }
#endif

    if (soalen == 4) {
#if defined(QPHIX_AVX_SOURCE) || defined(QPHIX_AVX2_SOURCE) ||                      \
    defined(QPHIX_MIC_SOURCE) || defined(QPHIX_AVX512_SOURCE) ||                    \
    defined(QPHIX_SSE_SOURCE)
      testTWMBiCGStabWrapper<float, VECLEN_SP, 4, MUF, PhiF>(u_in);
#endif
    }

    if (soalen == 8) {
#if defined(QPHIX_AVX_SOURCE) || defined(QPHIX_AVX2_SOURCE) ||                      \
    defined(QPHIX_MIC_SOURCE) || defined(QPHIX_AVX512_SOURCE)
      testTWMBiCGStabWrapper<float, VECLEN_SP, 8, MUF, PhiF>(u_in);
#endif
    }

    if (soalen == 16) {
#if defined(QPHIX_MIC_SOURCE) || defined(QPHIX_AVX512_SOURCE)
      testTWMBiCGStabWrapper<float, VECLEN_SP, 16, MUF, PhiF>(u_in);
#else
      masterPrintf("SOALEN = 16 not available");
      return;
#endif
    }

  } // FLOAT_PREC

  if (precision == HALF_PREC) {
#if defined(QPHIX_MIC_SOURCE) || defined(QPHIX_AVX2_SOURCE) ||                      \
    defined(QPHIX_AVX512_SOURCE)
    multi1d<LatticeColorMatrixF> u_in(4);
    for (int mu = 0; mu < Nd; mu++) {
      u_in[mu] = u[mu];
    }
    if (soalen == 4) {
      testTWMBiCGStabWrapper<half, VECLEN_HP, 4, MUF, PhiF>(u_in);
    }

    if (soalen == 8) {
      testTWMBiCGStabWrapper<half, VECLEN_HP, 8, MUF, PhiF>(u_in);
    }

    if (soalen == 16) {
      testTWMBiCGStabWrapper<half, VECLEN_HP, 16, MUF, PhiF>(u_in);
    }
#else
    QDPIO::cout << " Half Prec is only supported on MIC and AVX2 Targets just now "
                << endl;
#endif
  } // HALF_PREC

  if (precision == DOUBLE_PREC) {
    MUD u_in(4);
    for (int mu = 0; mu < Nd; mu++) {
      u_in[mu] = u[mu];
    }

    if (soalen == 2) {
#if defined(QPHIX_AVX_SOURCE) || defined(QPHIX_AVX2_SOURCE)
      testTWMBiCGStabWrapper<double, VECLEN_DP, 2, MUD, PhiD>(u_in);
#endif
    }

    if (soalen == 4) {
#if defined(QPHIX_AVX_SOURCE) || defined(QPHIX_AVX2_SOURCE) ||                      \
    defined(QPHIX_MIC_SOURCE) || defined(QPHIX_AVX512_SOURCE)
      testTWMBiCGStabWrapper<double, VECLEN_DP, 4, MUD, PhiD>(u_in);
#endif
    }

    if (soalen == 8) {
#if defined(QPHIX_MIC_SOURCE) || defined(QPHIX_AVX512_SOURCE)
      testTWMBiCGStabWrapper<double, VECLEN_DP, 8, MUD, PhiD>(u_in);
#endif
    }
  } // DOUBLE_PREC
  */
}

template <typename T, int V, int S, bool compress, typename U, typename Phi>
void TestTMDslash::testTWMDslash(int t_bc)
{
  RNG::setrn(rng_seed);

  typedef typename Geometry<T, V, S, compress>::SU3MatrixBlock Gauge;
  typedef typename Geometry<T, V, S, compress>::FourSpinorBlock Spinor;

  const double Mass = 0.2;
  const double TwistedMass = 0.1;

  const double alpha = 4.0 + Mass;
  const double Mu = TwistedMass / alpha;
  const double MuInv = alpha / (alpha * alpha + TwistedMass * TwistedMass);

  QDPIO::cout << endl
              << "TESTING TM Dlash (time boundary condition = " << t_bc << ")"
              << endl
              << "================" << endl
              << endl;

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

  TMDslash<T, V, S, compress> TMD(&geom,
                                  gauge.t_boundary,
                                  gauge.aniso_fac_s,
                                  gauge.aniso_fac_t,
                                  Mass,
                                  TwistedMass);

  HybridSpinor<T, V, S, compress, Phi> hs_source(geom), hs_qphix1(geom),
      hs_qphix2(geom), hs_qdp1(geom), hs_qdp2(geom);
  gaussian(hs_source.qdp());
  hs_source.pack();

  for (int isign = 1; isign >= -1; isign -= 2) {
    for (int cb = 0; cb <= 1; ++cb) {

      int source_cb = 1 - cb;
      int target_cb = cb;

      // 1. Apply QPhiX Dslash
      hs_qphix1.zero();
      TMD.dslash(hs_qphix1[target_cb],
                 hs_source[source_cb],
                 gauge.u_packed[target_cb],
                 isign,
                 target_cb);
      hs_qphix1.unpack();

      // 2. Apply QDP++ Dslash (Plain Wilson Dslash + Twisted Mass)
      hs_qdp1.zero();
      qdp_dslash(hs_qdp1.qdp(),
                 hs_source.qdp(),
                 gauge.u_aniso,
                 Mu,
                 MuInv,
                 isign,
                 target_cb);

      expect_near(hs_qdp1.qdp(), hs_qphix1.qdp(), 1e-6, geom, target_cb, "Dslash");
    } // cb
  } // isign
}

template <typename T, int V, int S, bool compress, typename U, typename Phi>
void TestTMDslash::testTWMDslashAChiMBDPsi(int t_bc)
{
  RNG::setrn(rng_seed);
  typedef typename Geometry<T, V, S, compress>::SU3MatrixBlock Gauge;
  typedef typename Geometry<T, V, S, compress>::FourSpinorBlock Spinor;

  const double aniso_fac_s = 1.0;
  const double aniso_fac_t = 1.0;
  const double t_boundary = (double)t_bc;

  const double Mass = 0.1;
  const double TwistedMass = 0.1;

  const double alpha = 4.0 + Mass;
  const double beta = 0.25;
  const double Mu = TwistedMass / alpha;
  const double MuInv = alpha / (alpha * alpha + TwistedMass * TwistedMass);

  QDPIO::cout << endl
              << "TESTING TM A chi - b D psi (time boundary condition = "
              << t_boundary << ")" << endl
              << "==========================" << endl
              << endl;

  Geometry<T, V, S, compress> geom(Layout::subgridLattSize().slice(),
                                   args_.By,
                                   args_.Bz,
                                   args_.NCores,
                                   args_.Sy,
                                   args_.Sz,
                                   args_.PadXY,
                                   args_.PadXYZ,
                                   args_.MinCt);

  RandomGauge<T, V, S, compress, U, Phi> gauge(geom, t_bc);

  TMDslash<T, V, S, compress> TMD(
      &geom, t_boundary, aniso_fac_s, aniso_fac_t, Mass, TwistedMass);

  HybridSpinor<T, V, S, compress, Phi> hs_source1(geom), hs_source2(geom),
      hs_qphix1(geom), hs_qphix2(geom), hs_qdp1(geom), hs_qdp2(geom);
  gaussian(hs_source1.qdp());
  gaussian(hs_source2.qdp());
  hs_source1.pack();
  hs_source2.pack();

  for (int isign = 1; isign >= -1; isign -= 2) {
    for (int target_cb = 0; target_cb <= 1; ++target_cb) {
      masterPrintf("Target CB: %i\n", target_cb);

      int const source_cb = 1 - target_cb;

      // 1. Apply QPhiX Dslash
      TMD.dslashAChiMinusBDPsi(hs_qphix1[target_cb],
                               hs_source2[source_cb],
                               hs_source1[target_cb],
                               gauge.u_packed[target_cb],
                               alpha,
                               beta,
                               isign,
                               target_cb);
      hs_qphix1.unpack();

      // 2. Apply QDP++ Dslash
      qdp_achimbdpsi(hs_qdp1.qdp(),
                     hs_source1.qdp(),
                     hs_source2.qdp(),
                     gauge.u_aniso,
                     Mu,
                     MuInv,
                     alpha,
                     beta,
                     isign,
                     target_cb);

      expect_near(hs_qdp1.qdp(), hs_qphix1.qdp(), 1e-6, geom, target_cb, "A chi - b D psi");
    } // cb
  } // isign
}

template <typename T, int V, int S, bool compress, typename U, typename Phi>
void TestTMDslash::testTWMM(int t_bc)
{
  RNG::setrn(rng_seed);

  typedef typename Geometry<T, V, S, compress>::SU3MatrixBlock Gauge;
  typedef typename Geometry<T, V, S, compress>::FourSpinorBlock Spinor;

  const double Mass = 0.1;
  const double TwistedMass = 0.1;

  const double alpha = 4.0 + Mass;
  const double beta = 0.25;
  const double Mu = TwistedMass / alpha;
  const double MuInv = alpha / (alpha * alpha + TwistedMass * TwistedMass);

  Geometry<T, V, S, compress> geom(Layout::subgridLattSize().slice(),
                                   args_.By,
                                   args_.Bz,
                                   args_.NCores,
                                   args_.Sy,
                                   args_.Sz,
                                   args_.PadXY,
                                   args_.PadXYZ,
                                   args_.MinCt);

  RandomGauge<T, V, S, compress, U, Phi> gauge(geom, t_bc);

  QDPIO::cout << endl
              << "TESTING TM Fermion Matrix (time boundary condition = "
              << gauge.t_boundary << ")" << endl
              << "=========================" << endl
              << endl;

  HybridSpinor<T, V, S, compress, Phi> hs_source(geom), hs_qphix1(geom),
      hs_qphix2(geom), hs_qdp1(geom), hs_qdp2(geom);
  gaussian(hs_source.qdp());
  hs_source.pack();

  EvenOddTMWilsonOperator<T, V, S, compress> M(Mass,
                                               TwistedMass,
                                               gauge.u_packed,
                                               &geom,
                                               gauge.t_boundary,
                                               gauge.aniso_fac_s,
                                               gauge.aniso_fac_t);

  for (int isign = 1; isign >= -1; isign -= 2) {
    for (int target_cb = 0; target_cb <= 1; ++target_cb) {
      int source_cb = 1 - target_cb;
      QDPIO::cout << "Target CB = " << target_cb << ", isign = " << isign << endl;

      // QPhiX version
      hs_qphix1.zero();
      M(hs_qphix1[target_cb], hs_source[target_cb], isign, target_cb);
      hs_qphix1.unpack();

      // QDP++ version:
      qdp_apply_operator(hs_qdp1.qdp(),
                         hs_source.qdp(),
                         gauge.u_aniso,
                         Mu,
                         MuInv,
                         alpha,
                         beta,
                         isign,
                         target_cb);

      expect_near(hs_qdp1.qdp(), hs_qphix1.qdp(), 1e-6, geom, target_cb, "TM Fermion Matrix");
    }
  }
}

// add twisted mass term.
template <typename T, int V, int S, bool compress, typename U, typename Phi>
void TestTMDslash::testTWMCG(int t_bc)
{
  RNG::setrn(rng_seed);

  typedef typename Geometry<T, V, S, compress>::SU3MatrixBlock Gauge;
  typedef typename Geometry<T, V, S, compress>::FourSpinorBlock Spinor;

  const double Mass = 0.1;
  const double TwistedMass = 0.1;

  const double alpha = 4.0 + Mass;
  const double beta = 0.25;
  const double Mu = TwistedMass / alpha;
  const double MuInv = alpha / (alpha * alpha + TwistedMass * TwistedMass);

  Geometry<T, V, S, compress> geom(Layout::subgridLattSize().slice(),
                                   args_.By,
                                   args_.Bz,
                                   args_.NCores,
                                   args_.Sy,
                                   args_.Sz,
                                   args_.PadXY,
                                   args_.PadXYZ,
                                   args_.MinCt);

  RandomGauge<T, V, S, compress, U, Phi> gauge(geom, t_bc);
  HybridSpinor<T, V, S, compress, Phi> hs_source(geom), hs_qphix1(geom),
      hs_qphix2(geom), hs_qdp1(geom), hs_qdp2(geom);
  gaussian(hs_source.qdp());
  hs_source.pack();

  QDPIO::cout << endl
              << "TESTING TM CG Inversion (time boundary condition = "
              << gauge.t_boundary << ")" << endl
              << "=======================" << endl
              << endl;

  QDPIO::cout << "Constructing TM Wilson Fermion Matrix... " << endl;
  EvenOddTMWilsonOperator<T, V, S, compress> M(Mass,
                                               TwistedMass,
                                               gauge.u_packed,
                                               &geom,
                                               gauge.t_boundary,
                                               gauge.aniso_fac_s,
                                               gauge.aniso_fac_t);
  QDPIO::cout << "done." << endl << endl;

  double r2 = 0.0;
  norm2Spinor<T, V, S, compress>(r2, hs_source[0], geom, 1);
  QDPIO::cout << "||psi||^2 = " << r2 << endl << endl;

  int const source_target_cb = 0;
  int const other_cb = 1 - source_target_cb;

  // 0. Setup Solver
  // ===============
  double rsd_target = rsdTarget<T>::value;
  int max_iters = 200;
  int niters = 0;
  double rsd_final = -1.0;
  unsigned long site_flops = 0;
  unsigned long mv_apps = 0;
  int isign = 1;
  InvCG<T, V, S, compress> solver(M, max_iters);

  // 1. QPhiX CG Solve
  // =================
  double start = omp_get_wtime();
  solver(hs_qphix1[source_target_cb],
         hs_source[source_target_cb],
         rsd_target,
         niters,
         rsd_final,
         site_flops,
         mv_apps,
         isign,
         verbose,
         source_target_cb);
  double end = omp_get_wtime();
  hs_qphix1.unpack();

  // 2. Multiply back with QDP + TM term
  // ===================================
  // chi2 = M chi
  qdp_apply_operator(hs_qdp1.qdp(),
                     hs_qphix1.qdp(),
                     gauge.u_aniso,
                     Mu,
                     MuInv,
                     alpha,
                     beta,
                     1,
                     source_target_cb);

  // chi3 = M^\dagger chi2
  qdp_apply_operator(hs_qdp2.qdp(),
                     hs_qdp1.qdp(),
                     gauge.u_aniso,
                     Mu,
                     MuInv,
                     alpha,
                     beta,
                     -1,
                     source_target_cb);

  unsigned long num_cb_sites = Layout::vol() / 2;
  unsigned long total_flops =
      (site_flops + (72 + 2 * 1320) * mv_apps) * num_cb_sites;
  masterPrintf("GFLOPS=%e\n", 1.0e-9 * (double)(total_flops) / (end - start));

  expect_near(hs_source.qdp(), hs_qdp2.qdp(), 1e-8, geom, source_target_cb, "CG");
}

template <typename T, int V, int S, bool compress, typename U, typename Phi>
void TestTMDslash::testTWMBiCGStab(int t_bc)
{
  RNG::setrn(rng_seed);

  typedef typename Geometry<T, V, S, compress>::SU3MatrixBlock Gauge;
  typedef typename Geometry<T, V, S, compress>::FourSpinorBlock Spinor;

  const double aniso_fac_s = 1.0;
  const double aniso_fac_t = 1.0;
  const double t_boundary = (double)(t_bc);

  const double Mass = 0.1;
  const double TwistedMass = 0.1;

  const double alpha = 4.0 + Mass;
  const double beta = 0.25;
  const double Mu = TwistedMass / alpha;
  const double MuInv = alpha / (alpha * alpha + TwistedMass * TwistedMass);

  QDPIO::cout << endl
              << "TESTING TM BiCGStab Inversion (time boundary condition = "
              << t_boundary << ")" << endl
              << "=============================" << endl
              << endl;

  Geometry<T, V, S, compress> geom(Layout::subgridLattSize().slice(),
                                   args_.By,
                                   args_.Bz,
                                   args_.NCores,
                                   args_.Sy,
                                   args_.Sz,
                                   args_.PadXY,
                                   args_.PadXYZ,
                                   args_.MinCt);

  RandomGauge<T, V, S, compress, U, Phi> gauge(geom, t_bc);
  HybridSpinor<T, V, S, compress, Phi> hs_source(geom), hs_qphix1(geom),
      hs_qphix2(geom), hs_qdp1(geom), hs_qdp2(geom);
  gaussian(hs_source.qdp());
  hs_source.pack();

  QDPIO::cout << "Constructing TM Wilson Fermion Matrix... " << endl;
  EvenOddTMWilsonOperator<T, V, S, compress> M(Mass,
                                               TwistedMass,
                                               gauge.u_packed,
                                               &geom,
                                               gauge.t_boundary,
                                               gauge.aniso_fac_s,
                                               gauge.aniso_fac_t);
  QDPIO::cout << "done." << endl;

  // 0. Setup Solver
  // ===============
  double rsd_target = rsdTarget<T>::value;
  int max_iters = 200;
  int niters = 0;
  double rsd_final = -1.0;
  unsigned long site_flops = 0;
  unsigned long mv_apps = 0;
  InvBiCGStab<T, V, S, compress> solver(M, max_iters);

  int const cb = 0;

  for (int isign = 1; isign >= -1; isign -= 2) {
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
    qdp_apply_operator(hs_qdp1.qdp(),
                       hs_qphix1.qdp(),
                       gauge.u_aniso,
                       Mu,
                       MuInv,
                       alpha,
                       beta,
                       isign,
                       cb);

    expect_near(hs_source.qdp(), hs_qdp1.qdp(), 1e-8, geom, cb, "TM Wilson BiCGStab");

    unsigned long num_cb_sites = Layout::vol() / 2;
    unsigned long total_flops =
        (site_flops + (72 + 2 * 1320) * mv_apps) * num_cb_sites;
    masterPrintf("BICGSTAB Solve isign = %d, GFLOPS = %e\n",
                 isign,
                 1.0e-9 * (double)(total_flops) / (end - start));
    QDPIO::cout << endl;

  } // isign
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
void TestTMDslash::testTWMRichardson(const U &u, int t_bc)
{
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

  // Over allocate, so that an unsigned load doesnt cause segfault accesing
  // either end...
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
  U u_test(Nd);
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

  const double Mass = 0.1;
  const double TwistedMass = 0.1;

  const double alpha = 4.0 + Mass;
  const double beta = 0.25;
  const double Mu = TwistedMass / alpha;
  const double MuInv = alpha / (alpha * alpha + TwistedMass * TwistedMass);

  EvenOddTMWilsonOperator<T1, VEC1, SOA1, compress> M_outer(Mass,
                                                            TwistedMass,
                                                            u_packed,
                                                            &geom_outer,
                                                            t_boundary,
                                                            aniso_fac_s,
                                                            aniso_fac_t);
  EvenOddTMWilsonOperator<T2, VEC2, SOA2, compress> M_inner(Mass,
                                                            TwistedMass,
                                                            u_packed_inner,
                                                            &geom_inner,
                                                            t_boundary,
                                                            aniso_fac_s,
                                                            aniso_fac_t);

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
      //      qdp_pack_spinor<T,V,S, compress,Phi >(chi, chi_even,
      //      chi_odd, geom);
      qdp_pack_cb_spinor<>(chi, chi_even, geom_outer, 0);
      masterPrintf("Entering solver\n");

      double start = omp_get_wtime();
      outer_solver(chi_s[0],
                   psi_s[0],
                   rsd_target,
                   niters,
                   rsd_final,
                   site_flops,
                   mv_apps,
                   isign,
                   verbose);
      double end = omp_get_wtime();

      //      qdp_unpack_spinor<T,V,S, compress, Phi >(chi_s[0], chi_s[1],
      //      chi, geom);
      qdp_unpack_cb_spinor<>(chi_s[0], chi, geom_outer, 0);

      // Multiply back
      // chi2 = M chi
      dslash(chi2, u_test, chi, isign, 1);
      dslash(ltmp, u_test, chi2, isign, 0);
      chi2[rb[0]] = massFactor * chi - betaFactor * ltmp;

      expect_near(psi, chi2, 1e-8, geom_outer, 0, "TM Wilson Richardson");

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
}

template <typename QDPSpinor>
void TestTMDslash::applyTwist(
    QDPSpinor &psi, double Mu, double Alpha, int isign, int target_cb)
{
  int Nxh = Nx / 2;

  for (int t = 0; t < Nt; t++) {
    for (int z = 0; z < Nz; z++) {
      for (int y = 0; y < Ny; y++) {
        for (int x = 0; x < Nxh; x++) {

          // These are unpadded QDP++ indices
          int index = x + Nxh * (y + Ny * (z + Nz * t));

          for (int spin = 0; spin < Ns; ++spin) {
            for (int color = 0; color < Nc; ++color) {

              double signed_mu = (spin < 2) ? 1.0 : -1.0;
              signed_mu *= isign * Mu;

              REAL x_re = psi.elem(rb[target_cb].start() + index)
                              .elem(spin)
                              .elem(color)
                              .real();
              REAL x_im = psi.elem(rb[target_cb].start() + index)
                              .elem(spin)
                              .elem(color)
                              .imag();

              REAL y_re = x_re - signed_mu * x_im;
              REAL y_im = x_im + signed_mu * x_re;

              psi.elem(rb[target_cb].start() + index).elem(spin).elem(color).real() =
                  Alpha * y_re;
              psi.elem(rb[target_cb].start() + index).elem(spin).elem(color).imag() =
                  Alpha * y_im;
            }
          }

        } // x
      } // y
    } // z
  } // t
}

template <typename QDPSpinor>
void TestTMDslash::applyInvTwist(
    QDPSpinor &psi, double Mu, double MuInv, int isign, int target_cb)
{
  int Nxh = Nx / 2;

  for (int t = 0; t < Nt; t++) {
    for (int z = 0; z < Nz; z++) {
      for (int y = 0; y < Ny; y++) {
        for (int x = 0; x < Nxh; x++) {

          // These are unpadded QDP++ indices
          int index = x + Nxh * (y + Ny * (z + Nz * t));

          for (int spin = 0; spin < Ns; ++spin) {
            for (int color = 0; color < Nc; ++color) {

              double signed_mu = (spin < 2) ? 1.0 : -1.0;
              signed_mu *= isign * Mu;

              REAL x_re = psi.elem(rb[target_cb].start() + index)
                              .elem(spin)
                              .elem(color)
                              .real();
              REAL x_im = psi.elem(rb[target_cb].start() + index)
                              .elem(spin)
                              .elem(color)
                              .imag();

              REAL y_re = x_re + signed_mu * x_im;
              REAL y_im = x_im - signed_mu * x_re;

              psi.elem(rb[target_cb].start() + index).elem(spin).elem(color).real() =
                  MuInv * y_re;
              psi.elem(rb[target_cb].start() + index).elem(spin).elem(color).imag() =
                  MuInv * y_im;
            }
          }

        } // x
      } // y
    } // z
  } // t
}

template <typename QdpGauge, typename QdpSpinor>
void TestTMDslash::qdp_dslash(QdpSpinor &out,
                                   QdpSpinor const &in,
                                   QDP::multi1d<QdpGauge> const &u_aniso,
                                   double const Mu,
                                   double const MuInv,
                                   int const isign,
                                   int const target_cb)
{
  dslash(out, u_aniso, in, isign, target_cb);
  applyInvTwist<>(out, Mu, MuInv, isign, target_cb);
}

template <typename QdpGauge, typename QdpSpinor>
void TestTMDslash::qdp_achimbdpsi(QdpSpinor &out,
                                  QdpSpinor const &chi,
                                  QdpSpinor const &psi,
                                  QDP::multi1d<QdpGauge> const &u_aniso,
                                  double const Mu,
                                  double const MuInv,
                                  double const alpha,
                                  double const beta,
                                  int const isign,
                                  int const target_cb)
{
  int const other_cb = 1 - target_cb;

  QdpSpinor D_psi = zero;
  dslash(D_psi, u_aniso, psi, isign, target_cb);

  QdpSpinor A_chi = chi;
  applyTwist<>(A_chi, Mu, alpha, isign, target_cb);

  out[rb[target_cb]] = A_chi - beta * D_psi;
}

template <typename QdpGauge, typename QdpSpinor>
void TestTMDslash::qdp_apply_operator(QdpSpinor &out,
                                      QdpSpinor const &in,
                                      QDP::multi1d<QdpGauge> const &u_aniso,
                                      double const Mu,
                                      double const MuInv,
                                      double const alpha,
                                      double const beta,
                                      int const isign,
                                      int const target_cb)
{
  int const other_cb = 1 - target_cb;

  QdpSpinor A_inv_D_psi = zero;
  qdp_dslash(A_inv_D_psi, in, u_aniso, Mu, MuInv, isign, other_cb);

  qdp_achimbdpsi(
      out, in, A_inv_D_psi, u_aniso, Mu, MuInv, alpha, beta, isign, target_cb);
}
