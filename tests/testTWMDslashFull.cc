#include "unittest.h"
#include "testTWMDslashFull.h"
#undef SEEK_SET
#undef SEEK_CUR
#undef SEEK_END
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

// What we consider to be small enough...
int Nx, Ny, Nz, Nt, Nxh;
bool verbose = true;

void testTWMDslashFull::run(void)
{
  RNG::savern(rng_seed);

  typedef LatticeColorMatrixF UF;
  typedef LatticeColorMatrixD UD;
  typedef LatticeDiracFermionF PhiF;
  typedef LatticeDiracFermionD PhiD;

  typedef multi1d<UF> MUF;
  typedef multi1d<UD> MUD;

  // Diagnostic information:
  const multi1d<int> &lattSize = Layout::subgridLattSize();
  Nx = lattSize[0];
  Ny = lattSize[1];
  Nz = lattSize[2];
  Nt = lattSize[3];

  QDPIO::cout << "Lattice Size: ";
  for (int mu = 0; mu < lattSize.size(); mu++) {
    QDPIO::cout << " " << lattSize[mu];
  }
  QDPIO::cout << endl;

  QDPIO::cout << "Block Sizes: By = " << By << ", Bz = " << Bz << endl;
  QDPIO::cout << "Number of Cores = " << NCores << endl;
  QDPIO::cout << "SMT Grid: Sy = " << Sy << ", Sz = " << Sz << endl;
  QDPIO::cout << "Pad Factors: PadXY = " << PadXY << ", PadXYZ = " << PadXYZ << endl;
  QDPIO::cout << "MinCt = " << MinCt << endl;
  QDPIO::cout << "Threads_per_core = " << N_simt << endl;
  QDPIO::cout << "Inititalizing QDP++ gauge field" << endl;

  // Make a random gauge field
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

    QDPIO::cout << "SINGLE PRECISION TESTING:" << endl;
    multi1d<LatticeColorMatrixF> u_in(4);
    for (int mu = 0; mu < Nd; mu++) {
      u_in[mu] = u[mu];
    }

    if (soalen == 1) {
      QDPIO::cout << "VECLEN = " << VECLEN_SP << ", SOALEN = 4" << endl;
      testTWMDslashWrapper<float, VECLEN_SP, 1, UF, PhiF>();
      testTWMDslashAChiMBDPsiWrapper<float, VECLEN_SP, 1, UF, PhiF>();
      testTWMMWrapper<float, VECLEN_SP, 1, UF, PhiF>();
      testTWMCGWrapper<float, VECLEN_SP, 1, MUF, PhiF>(u_in);
      testTWMBiCGStabWrapper<float, VECLEN_SP, 1, MUF, PhiF>(u_in);
    }

    if (soalen == 4) {
      QDPIO::cout << "VECLEN = " << VECLEN_SP << ", SOALEN = 4" << endl;
      testTWMDslashWrapper<float, VECLEN_SP, 4, UF, PhiF>();
      testTWMDslashAChiMBDPsiWrapper<float, VECLEN_SP, 4, UF, PhiF>();
      testTWMMWrapper<float, VECLEN_SP, 4, UF, PhiF>();
      testTWMCGWrapper<float, VECLEN_SP, 4, MUF, PhiF>(u_in);
      testTWMBiCGStabWrapper<float, VECLEN_SP, 4, MUF, PhiF>(u_in);
    }

    if (soalen == 8) {
      QDPIO::cout << "VECLEN = " << VECLEN_SP << ", SOALEN = 8" << endl;
      testTWMDslashWrapper<float, VECLEN_SP, 8, UF, PhiF>();
      testTWMDslashAChiMBDPsiWrapper<float, VECLEN_SP, 8, UF, PhiF>();
      testTWMMWrapper<float, VECLEN_SP, 8, UF, PhiF>();
      testTWMCGWrapper<float, VECLEN_SP, 8, MUF, PhiF>(u_in);
      testTWMBiCGStabWrapper<float, VECLEN_SP, 8, MUF, PhiF>(u_in);
    }

    if (soalen == 16) {
#if defined(QPHIX_MIC_SOURCE)
      QDPIO::cout << "VECLEN = " << VECLEN_SP << ", SOALEN = 16" << endl;
      testTWMDslashWrapper<float, VECLEN_SP, 16, UF, PhiF>();
      testTWMDslashAChiMBDPsiWrapper<float, VECLEN_SP, 16, UF, PhiF>();
      testTWMMWrapper<float, VECLEN_SP, 16, UF, PhiF>();
      testTWMCGWrapper<float, VECLEN_SP, 16, MUF, PhiF>(u_in);
      testTWMBiCGStabWrapper<float, VECLEN_SP, 16, MUF, PhiF>(u_in);
#else
      masterPrintf("SOALEN = 16 not available");
      return;
#endif
    }

  } // FLOAT_PREC

  if (precision == HALF_PREC) {
#if defined(QPHIX_MIC_SOURCE) || defined(QPHIX_AVX2_SOURCE)
    QDPIO::cout << "HALF PRECISION TESTING:" << endl;
    multi1d<LatticeColorMatrixF> u_in(4);
    for (int mu = 0; mu < Nd; mu++) {
      u_in[mu] = u[mu];
    }
    if (soalen == 4) {
      QDPIO::cout << "VECLEN = " << VECLEN_HP << ", SOALEN=4 " << endl;
      testTWMDslashWrapper<half, VECLEN_HP, 4, UF, PhiF>();
      testTWMDslashAChiMBDPsiWrapper<half, VECLEN_HP, 4, UF, PhiF>();
      testTWMCGWrapper<half, VECLEN_HP, 4, MUF, PhiF>(u_in);
      testTWMBiCGStabWrapper<half, VECLEN_HP, 4, MUF, PhiF>(u_in);
    }

    if (soalen == 8) {
      QDPIO::cout << "VECLEN = " << VECLEN_HP << ", SOALEN=8 " << endl;
      testTWMDslashWrapper<half, VECLEN_HP, 8, UF, PhiF>();
      testTWMDslashAChiMBDPsiWrapper<half, VECLEN_HP, 8, UF, PhiF>();
      testTWMCGWrapper<half, VECLEN_HP, 8, MUF, PhiF>(u_in);
      testTWMBiCGStabWrapper<half, VECLEN_HP, 8, MUF, PhiF>(u_in);
    }

    if (soalen == 16) {
#if defined(QPHIX_MIC_SOURCE)
      QDPIO::cout << "VECLEN = " << VECLEN_HP << ", SOALEN=16 " << endl;
      testTWMDslashWrapper<half, VECLEN_HP, 16, UF, PhiF>();
      testTWMDslashAChiMBDPsiWrapper<half, VECLEN_HP, 16, UF, PhiF>();
      testTWMCGWrapper<half, VECLEN_HP, 16, MUF, PhiF>(u_in);
      testTWMBiCGStabWrapper<half, VECLEN_HP, 16, MUF, PhiF>(u_in);
#endif
    }
#else
    QDPIO::cout << " Half Prec is only supported on MIC and AVX2 Targets just now "
                << endl;
#endif
  } // HALF_PREC

  if (precision == DOUBLE_PREC) {
    QDPIO::cout << "DOUBLE PRECISION TESTING:" << endl;
    MUD u_in(4);
    for (int mu = 0; mu < Nd; mu++) {
      u_in[mu] = u[mu];
    }

    if (soalen == 2) {
#if defined(QPHIX_AVX_SOURCE) || defined(QPHIX_AVX2_SOURCE)
      QDPIO::cout << "VECLEN = " << VECLEN_DP << ", SOALEN = 2 " << endl;
      testTWMDslashWrapper<double, VECLEN_DP, 2, UD, PhiD>();
      testTWMDslashAChiMBDPsiWrapper<double, VECLEN_DP, 2, UD, PhiD>();
      testTWMMWrapper<double, VECLEN_DP, 2, UD, PhiD>();
      testTWMCGWrapper<double, VECLEN_DP, 2, MUD, PhiD>(u_in);
      testTWMBiCGStabWrapper<double, VECLEN_DP, 2, MUD, PhiD>(u_in);
#endif
    }

    if (soalen == 4) {
      QDPIO::cout << "VECLEN = " << VECLEN_DP << ", SOALEN = 4 " << endl;
      testTWMDslashWrapper<double, VECLEN_DP, 4, UD, PhiD>();
      testTWMDslashAChiMBDPsiWrapper<double, VECLEN_DP, 4, UD, PhiD>();
      testTWMMWrapper<double, VECLEN_DP, 4, UD, PhiD>();
      testTWMCGWrapper<double, VECLEN_DP, 4, MUD, PhiD>(u_in);
      testTWMBiCGStabWrapper<double, VECLEN_DP, 4, MUD, PhiD>(u_in);
    }

    if (soalen == 8) {
#if defined(QPHIX_MIC_SOURCE)
      QDPIO::cout << "VECLEN = " << VECLEN_DP << ", SOALEN = 8 " << endl;
      testTWMDslashWrapper<double, VECLEN_DP, 8, UD, PhiD>();
      testTWMDslashAChiMBDPsiWrapper<double, VECLEN_DP, 8, UD, PhiD>();
      testTWMMWrapper<double, VECLEN_DP, 8, UD, PhiD>();
      testTWMCGWrapper<double, VECLEN_DP, 8, MUD, PhiD>(u_in);
      testTWMBiCGStabWrapper<double, VECLEN_DP, 8, MUD, PhiD>(u_in);
#endif
    }
  } // DOUBLE_PREC
}

template <typename T, int V, int S, bool compress, typename U, typename Phi>
void testTWMDslashFull::testTWMDslash(int t_bc)
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
                                   By,
                                   Bz,
                                   NCores,
                                   Sy,
                                   Sz,
                                   PadXY,
                                   PadXYZ,
                                   MinCt);

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

  for (auto isign : {1, -1}) {
    for (auto cb : {0, 1}) {

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
      dslash(hs_qdp1.qdp(), gauge.get_u_aniso(), hs_source.qdp(), isign, target_cb);
      applyInvTwist<>(hs_qdp1.qdp(), Mu, MuInv, isign, target_cb);

      expect_near(hs_qphix1, hs_qdp1, 1e-6, geom, target_cb, "Dslash");
    } // cb
  } // isign
}

template <typename T, int V, int S, bool compress, typename U, typename Phi>
void testTWMDslashFull::testTWMDslashAChiMBDPsi(int t_bc)
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
                                   By,
                                   Bz,
                                   NCores,
                                   Sy,
                                   Sz,
                                   PadXY,
                                   PadXYZ,
                                   MinCt);

  RandomGauge<T, V, S, compress, U, Phi> gauge(geom, t_bc);

  TMDslash<T, V, S, compress> TMD(
      &geom, t_boundary, aniso_fac_s, aniso_fac_t, Mass, TwistedMass);

  HybridSpinor<T, V, S, compress, Phi> hs_source1(geom), hs_source2(geom),
      hs_qphix1(geom), hs_qphix2(geom), hs_qdp1(geom), hs_qdp2(geom);
  gaussian(hs_source1.qdp());
  gaussian(hs_source2.qdp());
  hs_source1.pack();
  hs_source2.pack();

  for (auto isign : {1, -1}) {
    for (auto target_cb : {1, 0}) {
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
      dslash(hs_qdp1.qdp(), gauge.u_aniso, hs_source2.qdp(), isign, target_cb);
      Phi source_copy = hs_source1.qdp();
      applyTwist<>(source_copy, Mu, alpha, isign, target_cb);

      hs_qdp2.zero();
      hs_qdp2.qdp()[rb[target_cb]] = source_copy - beta * hs_qdp1.qdp();

      expect_near(hs_qphix1, hs_qdp2, 1e-6, geom, target_cb, "A chi - b D psi");
    } // cb
  } // isign
}

template <typename T, int V, int S, bool compress, typename U, typename Phi>
void testTWMDslashFull::testTWMM(int t_bc)
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
                                   By,
                                   Bz,
                                   NCores,
                                   Sy,
                                   Sz,
                                   PadXY,
                                   PadXYZ,
                                   MinCt);

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

  for (auto isign : {1, -1}) {
    for (int target_cb : {1, 0}) {
      int source_cb = 1 - target_cb;
      QDPIO::cout << "Target CB = " << target_cb << endl;

      // QPhiX version
      hs_qphix1.zero();
      M(hs_qphix1[target_cb], hs_source[target_cb], isign);
      hs_qphix1.unpack();

      // QDP++ version:
      // (a) hs_qdp2 := D A^{-1} D hs_source
      dslash(hs_qdp1.qdp(), gauge.u_aniso, hs_source.qdp(), isign, source_cb);
      applyInvTwist<>(hs_qdp1.qdp(), Mu, MuInv, isign, source_cb);
      dslash(hs_qdp2.qdp(), gauge.u_aniso, hs_qdp1.qdp(), isign, target_cb);

      // (b) tmp := A hs_source
      Phi tmp = hs_source.qdp();
      expect_near(
          hs_source.qdp(), tmp, 1e-6, geom, target_cb, "QDP Spinor copy-ctor");
      applyTwist<>(tmp, Mu, alpha, isign, target_cb);

      // (c) result := A hs_source - beta D A^{-1} D hs_source
      Phi qdp_result = zero;
      qdp_result[rb[target_cb]] = tmp - beta * hs_qdp2.qdp();

      expect_near(
          hs_qphix1.qdp(), qdp_result, 1e-6, geom, target_cb, "TM Fermion Matrix");
    }
  }
}

// add twisted mass term.
template <typename T, int V, int S, bool compress, typename U, typename Phi>
void testTWMDslashFull::testTWMCG(const U &u, int t_bc)
{
  RNG::setrn(rng_seed);

  typedef typename Geometry<T, V, S, compress>::SU3MatrixBlock Gauge;
  typedef typename Geometry<T, V, S, compress>::FourSpinorBlock Spinor;

  const double aniso_fac_s = 0.35;
  const double aniso_fac_t = 1.4;
  const double t_boundary = (double)(t_bc);

  const double Mass = 0.1;
  const double TwistedMass = 0.1;

  const double alpha = 4.0 + Mass;
  const double beta = 0.25;
  const double Mu = TwistedMass / alpha;
  const double MuInv = alpha / (alpha * alpha + TwistedMass * TwistedMass);

  QDPIO::cout << endl
              << "TESTING TM CG Inversion (time boundary condition = " << t_boundary
              << ")" << endl
              << "=======================" << endl
              << endl;

  Geometry<T, V, S, compress> geom(Layout::subgridLattSize().slice(),
                                   By,
                                   Bz,
                                   NCores,
                                   Sy,
                                   Sz,
                                   PadXY,
                                   PadXYZ,
                                   MinCt);

  QDPIO::cout << "Allocating Fields... ";
  Gauge *packed_gauge_cb0 = (Gauge *)geom.allocCBGauge();
  Gauge *packed_gauge_cb1 = (Gauge *)geom.allocCBGauge();
  Spinor *psi_even = (Spinor *)geom.allocCBFourSpinor();
  Spinor *psi_odd = (Spinor *)geom.allocCBFourSpinor();
  Spinor *chi_even = (Spinor *)geom.allocCBFourSpinor();
  Spinor *chi_odd = (Spinor *)geom.allocCBFourSpinor();
  Gauge *u_packed[2];
  u_packed[0] = packed_gauge_cb0;
  u_packed[1] = packed_gauge_cb1;
  Spinor *psi_s[2] = {psi_even, psi_odd};
  Spinor *chi_s[2] = {chi_even, chi_odd};
  QDPIO::cout << "done." << endl;

  QDPIO::cout << "Filling psi with gaussian noise... ";
  Phi psi;
  gaussian(psi);
  QDPIO::cout << "done." << endl;

  QDPIO::cout << "Packing gauge field... ";
  qdp_pack_gauge<>(u, packed_gauge_cb0, packed_gauge_cb1, geom);
  QDPIO::cout << "done." << endl;

  QDPIO::cout << "Packing fermions... ";
  qdp_pack_spinor<>(psi, psi_even, psi_odd, geom);
  QDPIO::cout << "done." << endl;

  QDPIO::cout << "Applying anisotropy to test gauge field... ";
  U u_test(Nd);
  for (int mu = 0; mu < Nd; mu++) {
    Real factor = Real(aniso_fac_s);
    if (mu == Nd - 1) {
      factor = Real(aniso_fac_t);
    }
    u_test[mu] = factor * u[mu];
  }
  QDPIO::cout << "done." << endl;

  // Apply BCs on u-test for QDP++ test (Dslash gets unmodified field)
  QDPIO::cout << "Applying BCs to test gauge field... ";
  u_test[3] *= where(Layout::latticeCoordinate(3) == (Layout::lattSize()[3] - 1),
                     Real(t_boundary),
                     Real(1));
  QDPIO::cout << "done." << endl;

  QDPIO::cout << "Constructing TM Wilson Fermion Matrix... " << endl;
  EvenOddTMWilsonOperator<T, V, S, compress> M(
      Mass, TwistedMass, u_packed, &geom, t_boundary, aniso_fac_s, aniso_fac_t);
  QDPIO::cout << "done." << endl << endl;

  double r2 = 0.0;
  norm2Spinor<T, V, S, compress>(r2, psi_even, geom, 1);
  QDPIO::cout << "||psi||^2 = " << r2 << endl << endl;

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
  Phi chi = zero;
  qdp_pack_spinor<>(chi, chi_even, chi_odd, geom);
  double start = omp_get_wtime();
  solver(chi_s[0],
         psi_s[0],
         rsd_target,
         niters,
         rsd_final,
         site_flops,
         mv_apps,
         isign,
         verbose);
  double end = omp_get_wtime();
  qdp_unpack_spinor<>(chi_s[0], chi_s[1], chi, geom);

  // 2. Multiply back with QDP + TM term
  // ===================================
  // chi2 = M chi
  Phi chi2 = zero;
  Phi ltmp = zero;
  dslash(chi2, u_test, chi, 1, 1);
  applyInvTwist<>(chi2, Mu, MuInv, 1, 1);
  dslash(ltmp, u_test, chi2, 1, 0);
  applyTwist<>(chi, Mu, alpha, 1, 0);
  chi2[rb[0]] = chi - beta * ltmp;

  // chi3 = M^\dagger chi2
  Phi chi3 = zero;
  dslash(chi3, u_test, chi2, (-1), 1);
  applyInvTwist<>(chi3, Mu, MuInv, (-1), 1);
  dslash(ltmp, u_test, chi3, (-1), 0);
  applyTwist<>(chi2, Mu, alpha, (-1), 0);
  chi3[rb[0]] = chi2 - beta * ltmp;

  // 3. Assert difference & GFLOP/s
  // ===================================
  Phi diff = chi3 - psi;
  Double true_norm = sqrt(norm2(diff, rb[0]) / norm2(psi, rb[0]));
  QDPIO::cout << "True norm of residuum is: " << true_norm << endl;

  unsigned long num_cb_sites = Layout::vol() / 2;
  unsigned long total_flops =
      (site_flops + (72 + 2 * 1320) * mv_apps) * num_cb_sites;
  masterPrintf("GFLOPS=%e\n", 1.0e-9 * (double)(total_flops) / (end - start));
  assertion(toBool(true_norm < (rsd_target + tolerance<T>::small)));

  geom.free(packed_gauge_cb0);
  geom.free(packed_gauge_cb1);
  geom.free(psi_even);
  geom.free(psi_odd);
  geom.free(chi_even);
  geom.free(chi_odd);
}

template <typename T, int V, int S, bool compress, typename U, typename Phi>
void testTWMDslashFull::testTWMBiCGStab(const U &u, int t_bc)
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
                                   By,
                                   Bz,
                                   NCores,
                                   Sy,
                                   Sz,
                                   PadXY,
                                   PadXYZ,
                                   MinCt);

  Gauge *packed_gauge_cb0 = (Gauge *)geom.allocCBGauge();
  Gauge *packed_gauge_cb1 = (Gauge *)geom.allocCBGauge();
  Spinor *psi_even = (Spinor *)geom.allocCBFourSpinor();
  Spinor *psi_odd = (Spinor *)geom.allocCBFourSpinor();
  Spinor *chi_even = (Spinor *)geom.allocCBFourSpinor();
  Spinor *chi_odd = (Spinor *)geom.allocCBFourSpinor();
  Gauge *u_packed[2];
  u_packed[0] = packed_gauge_cb0;
  u_packed[1] = packed_gauge_cb1;
  Spinor *psi_s[2] = {psi_even, psi_odd};
  Spinor *chi_s[2] = {chi_even, chi_odd};
  QDPIO::cout << "Fields allocated" << endl;

  QDPIO::cout << "Filling psi with gaussian noise... ";
  Phi psi;
  gaussian(psi);
  QDPIO::cout << "done." << endl;

  QDPIO::cout << "Packing gauge field... ";
  qdp_pack_gauge<>(u, packed_gauge_cb0, packed_gauge_cb1, geom);
  QDPIO::cout << "done." << endl;

  QDPIO::cout << "Packing fermions... ";
  qdp_pack_spinor<>(psi, psi_even, psi_odd, geom);
  QDPIO::cout << "done." << endl;

  QDPIO::cout << "Applying anisotropy to test gauge field... ";
  U u_test(Nd);
  for (int mu = 0; mu < Nd; mu++) {
    Real factor = Real(aniso_fac_s);
    if (mu == Nd - 1) {
      factor = Real(aniso_fac_t);
    }
    u_test[mu] = factor * u[mu];
  }
  QDPIO::cout << "done." << endl;

  QDPIO::cout << "Applying BCs to test gauge field... ";
  u_test[3] *= where(Layout::latticeCoordinate(3) == (Layout::lattSize()[3] - 1),
                     Real(t_boundary),
                     Real(1));
  QDPIO::cout << "done." << endl;

  QDPIO::cout << "Constructing TM Wilson Fermion Matrix... " << endl;
  EvenOddTMWilsonOperator<T, V, S, compress> M(
      Mass, TwistedMass, u_packed, &geom, t_boundary, aniso_fac_s, aniso_fac_t);
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

  for (auto isign : {1, -1}) {

    Phi chi = zero;
    qdp_pack_cb_spinor<>(chi, chi_even, geom, 0);
    double start = omp_get_wtime();
    solver(chi_s[0],
           psi_s[0],
           rsd_target,
           niters,
           rsd_final,
           site_flops,
           mv_apps,
           isign,
           verbose);
    double end = omp_get_wtime();
    qdp_unpack_cb_spinor<>(chi_s[0], chi, geom, 0);

    // Multiply back
    // chi2 = M chi
    Phi chi2 = zero;
    Phi ltmp = zero;
    dslash(chi2, u_test, chi, isign, 1);
    applyInvTwist<>(chi2, Mu, MuInv, isign, 1);
    dslash(ltmp, u_test, chi2, isign, 0);
    applyTwist<>(chi, Mu, alpha, isign, 0);
    chi2[rb[0]] = chi - beta * ltmp;

    Phi diff = chi2 - psi;
    Double true_norm = sqrt(norm2(diff, rb[0]) / norm2(psi, rb[0]));
    QDPIO::cout << "BiCGStab Solve isign = " << isign
                << ": True norm of residual is " << true_norm << endl;

    unsigned long num_cb_sites = Layout::vol() / 2;
    unsigned long total_flops =
        (site_flops + (72 + 2 * 1320) * mv_apps) * num_cb_sites;
    masterPrintf("BICGSTAB Solve isign = %d, GFLOPS = %e\n",
                 isign,
                 1.0e-9 * (double)(total_flops) / (end - start));
    QDPIO::cout << endl;
    assertion(toBool(true_norm < (rsd_target + tolerance<T>::small)));

  } // isign

  geom.free(packed_gauge_cb0);
  geom.free(packed_gauge_cb1);
  geom.free(psi_even);
  geom.free(psi_odd);
  geom.free(chi_even);
  geom.free(chi_odd);
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
void testTWMDslashFull::testTWMRichardson(const U &u, int t_bc)
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
                                                By,
                                                Bz,
                                                NCores,
                                                Sy,
                                                Sz,
                                                PadXY,
                                                PadXYZ,
                                                MinCt);

  Geometry<T2, VEC2, SOA2, compress> geom_inner(Layout::subgridLattSize().slice(),
                                                By,
                                                Bz,
                                                NCores,
                                                Sy,
                                                Sz,
                                                PadXY,
                                                PadXYZ,
                                                MinCt);

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

      Phi diff = chi2 - psi;
      Double true_norm = sqrt(norm2(diff, rb[0]) / norm2(psi, rb[0]));
      QDPIO::cout << "RICHARDSON Solve isign=" << isign
                  << " True norm is: " << true_norm << endl;
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
}

template <typename QDPSpinor>
void testTWMDslashFull::applyTwist(
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
void testTWMDslashFull::applyInvTwist(
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
