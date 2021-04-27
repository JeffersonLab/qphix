
#include "unittest.h"
#include "testNDTMvsTMCloverFull.h"

#include "qphix/ndtm_clov_dslash_def.h"
#include "qphix/ndtm_clov_dslash_body.h"

#include "qdp.h"
using namespace QDP;

#ifndef DSLASH_M_W_H
#include "dslashm_w.h"
#endif

#ifndef REUNIT_H
#include "reunit.h"
#endif

#include "qphix/qdp_packer.h"
#include "./clover_term.h"

#include "qphix/twisted_mass_clover.h"
#include "qphix/ndtm_operator_clover.h"
#include "qphix/invcg.h"
#include "qphix/invbicgstab.h"

#include <omp.h>

using namespace Assertions;
using namespace std;
using namespace QPhiX;

#include "soalen.h"
#include "tolerance_type_traits.h"
#include "assertion_helper.h"
#include "packed_spinor.h"
#include "packed_clover.h"

/**
  Add a twisted mass to an existing clover term.

  Essentially we are given \f$ c_\mathrm{SW} \sigma_{\mu\nu} F_{\mu\nu} \f$ and
  want to add \f$ \mathrm i \mu \gamma_5 \tau_3 \f$ where the \f$ \tau_3 \f$ is
  just boiled down to \f$ \pm \f$ because the current implementation is not
  aware of multiple flavors and just takes on the two flavors via the sign
  change.

  @tparam FT Floating point type
  @tparam V Length of SIMD vectors in terms of the floating type
  @tparam soalen Length of the SoA
  @tparam compress Whether to use gauge compression

  @param[in] twisted_mass_mu Absolute value of twisted mass to add, \f$ \mu \f$
  @param[in,out] clov_packed Clover term with the following structure:
  - Checkerboard, length 2.
  - Plus/Minus (normal or hermitian conjugate), length 2.
  - One combined spacetime index.
  - Then structure like described in \ref QPhiX::Geometry::FullCloverBlock
  */
template <typename FT, int V, int soalen, bool compress>
void add_twisted_mass_to_full_clover(
    const double twisted_mass_mu,
    const Geometry<FT, V, soalen, compress> &geometry,
    typename Geometry<FT, V, soalen, compress>::FullCloverBlock *clov_packed[2][2])
{
  typedef typename Geometry<FT, V, soalen, compress>::FullCloverBlock FullClover;

  const auto Nt = geometry.Nt();
  const auto Ny = geometry.Ny();
  const auto Nz = geometry.Nz();
  const auto Pxy = geometry.getPadXY();
  const auto Pxyz = geometry.getPadXYZ();
  const auto nvecs = geometry.nVecs();
  const auto nyg = geometry.nGY();

  for (int cb : {0, 1}) {
    for (int pm : {0, 1}) {
#pragma omp parallel for collapse(4)
      for (int64_t t = 0; t < Nt; t++) {
        for (int64_t z = 0; z < Nz; z++) {
          for (int64_t y = 0; y < Ny; y++) {
            for (int64_t s = 0; s < nvecs; s++) {
              int block = (t * Pxyz + z * Pxy) / nyg + (y / nyg) * nvecs + s;
              auto &clov_block = clov_packed[cb][pm][block];
              auto &clov_block1 = clov_block.block1;
              auto &clov_block2 = clov_block.block2;

              for (int c = 0; c < 3; c++) {
                for (int spin = 0; spin < 2; spin++) {
                  const auto cs = 3 * spin + c;
                  for (int64_t x = 0; x < soalen; x++) {
                    int xx = (y % nyg) * soalen + x;
                    clov_block1[cs][cs][IM][xx] += -twisted_mass_mu;
                    clov_block2[cs][cs][IM][xx] += twisted_mass_mu;
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
void testNDTMvsTMCloverFull::runTest(void)
{
  typedef typename Geometry<FT, V, S, compress>::FourSpinorBlock Spinor;
  typedef typename Geometry<FT, V, S, compress>::SU3MatrixBlock Gauge;
  typedef typename Geometry<FT, V, S, compress>::FullCloverBlock FullClover;

  bool verbose = false;
  QDPIO::cout << endl << "ENTERING NDTM vs TM CLOVER TEST" << endl;
  QDPIO::cout << "(testing NDTM vs TM clover with no mass splitting)" << endl
              << endl;

  // Diagnostic information:
  const multi1d<int> &lattSize = Layout::subgridLattSize();
  int Nx = lattSize[0];
  int Ny = lattSize[1];
  int Nz = lattSize[2];
  int Nt = lattSize[3];

  QDPIO::cout << "VECLEN = " << V << ", SOALEN = " << S << endl;
  QDPIO::cout << "Lattice Size: ";
  for (int mu = 0; mu < lattSize.size(); mu++) {
    QDPIO::cout << " " << lattSize[mu];
  }
  QDPIO::cout << endl;

  QDPIO::cout << "Block Sizes: By = " << By << ", Bz = " << Bz << endl;
  QDPIO::cout << "N Cores: " << NCores << endl;
  QDPIO::cout << "SMT Grid: Sy = " << Sy << ", Sz = " << Sz << endl;
  QDPIO::cout << "Pad Factors: PadXY = " << PadXY << ", PadXYZ = " << PadXYZ << endl;
  QDPIO::cout << "MinCt = " << MinCt << endl;
  QDPIO::cout << "Threads_per_core = " << N_simt << endl;

  QDPIO::cout << "Inititalizing QDP++ gauge field." << endl;

  // Make a random gauge field
  // Start up the gauge field somehow
  // We choose: u = reunit(  1 + factor*gaussian(u) );
  // Adjust factor to get reasonable inversion times in invert test.
  // Bug gauge field is always unitary
  multi1d<U> u(4);
  U g;
  U uf;
  for (int mu = 0; mu < 4; mu++) {
    uf = 1; // Unit gauge
    Real factor = Real(0.08);
    gaussian(g);
    u[mu] = uf + factor * g;
    reunit(u[mu]);
  }

  // Set anisotropy parameters -- pick some random numbers
  double xi_0_f = 0.3;
  double nu_f = 1.4;
  // double xi_0_f = 1;
  // double nu_f = 1;

  // This is what makeFermCoeffs does under the hood.
  // Get the direct spatial and temporal anisotropy factors
  double aniso_fac_s = (double)nu_f / xi_0_f;
  double aniso_fac_t = (double)(1);

  QDPIO::cout << "Setting Clover Term Parameters" << endl;

  CloverFermActParams clparam;
  AnisoParam_t aniso;

  // Aniso prarams
  aniso.anisoP = true;
  aniso.xi_0 = xi_0_f;
  aniso.nu = nu_f;
  aniso.t_dir = 3;

  // Set up the Clover params
  clparam.anisoParam = aniso;

  // Some mass
  clparam.Mass = Real(0.1);

  // Some random clover coeffs
  // This should be wilson...
  clparam.clovCoeffR = Real(1.2);
  clparam.clovCoeffT = Real(0.9);

  // Set up the 'periodic BC dslash'
  QDPIO::cout << "Dslash will run with " << omp_get_max_threads() << " threads"
              << endl;

  double t_boundary = (double)(1);
  QDPIO::cout << "Instantiating ClovDslash<FT," << V << "," << S << ">"
              << " with t_boundary = " << t_boundary << std::endl;
  Geometry<FT, V, S, compress> geom(Layout::subgridLattSize().slice(),
                                    By,
                                    Bz,
                                    NCores,
                                    Sy,
                                    Sz,
                                    PadXY,
                                    PadXYZ,
                                    MinCt);
  TMClovDslash<FT, V, S, compress> tmdslash(
      &geom, t_boundary, aniso_fac_s, aniso_fac_t);
  NDTMClovDslash<FT, V, S, compress> ndtmdslash(
      &geom, t_boundary, aniso_fac_s, aniso_fac_t);

  // Make a random source
  QDPIO::cout << "Initializing QDP++ input spinor..." << endl;
  Phi chi, clov_chi;
  Phi psi, chi2, clov_chi2;
  QDPIO::cout << "Filling psi with gaussian noise..." << endl;
  gaussian(psi);

  // Allocating fields
  QDPIO::cout << "Allocating packged gauge fields..." << endl;
  Gauge *packed_gauge_cb0 = (Gauge *)geom.allocCBGauge();
  Gauge *packed_gauge_cb1 = (Gauge *)geom.allocCBGauge();

  QDPIO::cout << "Allocating packed spinor fields..." << endl;
  Spinor *psi_even = (Spinor *)geom.allocCBFourSpinor();
  Spinor *psi_odd = (Spinor *)geom.allocCBFourSpinor();
  Spinor *chi_even = (Spinor *)geom.allocCBFourSpinor();
  Spinor *chi_odd = (Spinor *)geom.allocCBFourSpinor();
  Spinor *chi2_even = (Spinor *)geom.allocCBFourSpinor();
  Spinor *chi2_odd = (Spinor *)geom.allocCBFourSpinor();

  QDPIO::cout << "Allocate Packed (Inverse) Clover Terms..." << endl;
  FullClover *A_cb0_plus = (FullClover *)geom.allocCBFullClov();
  FullClover *A_cb0_minus = (FullClover *)geom.allocCBFullClov();
  FullClover *A_cb1_plus = (FullClover *)geom.allocCBFullClov();
  FullClover *A_cb1_minus = (FullClover *)geom.allocCBFullClov();
  FullClover *A_inv_cb0_plus = (FullClover *)geom.allocCBFullClov();
  FullClover *A_inv_cb0_minus = (FullClover *)geom.allocCBFullClov();
  FullClover *A_inv_cb1_plus = (FullClover *)geom.allocCBFullClov();
  FullClover *A_inv_cb1_minus = (FullClover *)geom.allocCBFullClov();

  FullClover *clov_packed[2][2];
  FullClover *invclov_packed[2][2];

  clov_packed[0][0] = A_cb0_plus;
  clov_packed[0][1] = A_cb0_minus;
  clov_packed[1][0] = A_cb1_plus;
  clov_packed[1][1] = A_cb1_minus;
  invclov_packed[0][0] = A_inv_cb0_plus;
  invclov_packed[0][1] = A_inv_cb0_minus;
  invclov_packed[1][0] = A_inv_cb1_plus;
  invclov_packed[1][1] = A_inv_cb1_minus;

  QDPIO::cout << "Fields allocated." << endl;

  // Packing fields to QPhiX layout
  // Pack the gauge field
  QDPIO::cout << "Packing gauge field...";
  qdp_pack_gauge<>(u, packed_gauge_cb0, packed_gauge_cb1, geom);
  Gauge *u_packed[2];
  u_packed[0] = packed_gauge_cb0;
  u_packed[1] = packed_gauge_cb1;
  QDPIO::cout << " Done." << std::endl;

  // Pack the spinor field
  QDPIO::cout << "Packing fermions..." << std::flush;
  Spinor *psi_s[2] = {psi_even, psi_odd};
  Spinor *chi_s[2] = {chi_even, chi_odd};
  Spinor *chi2_s[2] = {chi2_even, chi2_odd};
  qdp_pack_spinor<>(psi, psi_even, psi_odd, geom);
  QDPIO::cout << " Done." << endl;

  ////////////////////////////////////////////////////////////

  // Create a QDP++ clover term from the gauge field. Anisotropy is already
  // taken care of by that clover term class.
  CloverTermT<Phi, U> clov_qdp;
  clov_qdp.create(u, clparam);

  // Let QDP++ invert that.
  CloverTermT<Phi, U> invclov_qdp(clov_qdp);
  for (int cb : {0, 1}) {
    invclov_qdp.choles(cb);
  }

  // Pack the clover term onto the checkerboard, also copy it to the plus and
  // minus case.
  for (int cb : {0, 1}) {
    for (int pm : {0, 1}) {
      qdp_pack_full_clover<>(clov_qdp, clov_packed[cb][pm], geom, cb);
    }
  }

  // The clover term now is in a data structure that supports twisted mass.
  // The QDP++ data structure assumes that it is hermitian and therefore only
  // stores half the data.
  add_twisted_mass_to_full_clover(0.3845, geom, clov_packed);

  // Also copy the inverse clover term onto the checkerboard.
  for (int cb : {0, 1}) {
    for (int pm : {0, 1}) {
      qdp_pack_full_clover<>(invclov_qdp, invclov_packed[cb][pm], geom, cb);
    }
  }

  ////////////////////////////////////////////////////////////

  QDPIO::cout << "Multiplying aniso coeffs into gauge field for testing" << endl;
  multi1d<U> u_test(4);
  for (int mu = 0; mu < Nd; mu++) {
    Real factor;
    if (mu == 3) {
      factor = Real(aniso_fac_t);
    } else {
      factor = Real(aniso_fac_s);
    }
    u_test[mu] = factor * u[mu];
  }

#if 1 // This tests dslash
  {

    typedef PackedSpinor<FT, V, S, compress, Phi> MySpinor;

    MySpinor spinor_psi_up(geom);
    MySpinor spinor_psi_down(geom);
    MySpinor spinor_result_base_up(geom);
    MySpinor spinor_result_base_down(geom);
    MySpinor spinor_result_slow_up(geom);
    MySpinor spinor_result_slow_down(geom);

    typedef PackedClover<FT, V, S, compress> MyClover;

    MyClover zero_clover(geom);
    zero_packed_clover(zero_clover);

    // 1. Test only Dslash operator.
    // For clover this will be: A^{-1}_(1-cb,1-cb) D_(1-cb, cb)  psi_cb
    QDPIO::cout << std::endl << "TESTING TWISTED-MASS CLOVER DSLASH:" << std::endl;

    const int cbsize_in_blocks = rb[0].numSiteTable() / S;

    for (auto isign : {1, -1}) {
      for (auto target_cb : {0, 1}) {
        int source_cb = 1 - target_cb;

        spinor_result_base_up.all = zero;
        spinor_result_base_up.pack();
        spinor_result_base_down.all = zero;
        spinor_result_base_down.pack();
        spinor_result_slow_up.all = zero;
        spinor_result_slow_up.pack();
        spinor_result_slow_down.all = zero;
        spinor_result_slow_down.pack();

        // (a) Apply the base version.
        tmdslash.dslash(spinor_result_base_up[target_cb],
                        spinor_psi_up[source_cb],
                        u_packed[target_cb],
                        (const FullClover **)invclov_packed[target_cb],
                        isign,
                        target_cb);
        tmdslash.dslash(spinor_result_base_down[target_cb],
                        spinor_psi_down[source_cb],
                        u_packed[target_cb],
                        (const FullClover **)invclov_packed[target_cb],
                        isign,
                        target_cb);

        // (b) Apply Martin's slow version.
        Spinor *two_flav_res[2] = {spinor_result_slow_up[target_cb],
                                   spinor_result_slow_down[target_cb]};
        const Spinor *two_flav_psi[2] = {spinor_psi_up[source_cb],
                                         spinor_psi_down[source_cb]};
        const FullClover *two_flav_invclov[2][2][2] = {
            {{invclov_packed[target_cb][0], invclov_packed[target_cb][1]},
             {zero_clover.invclov_packed[target_cb][0],
              zero_clover.invclov_packed[target_cb][1]}},
            {{zero_clover.invclov_packed[target_cb][0],
              zero_clover.invclov_packed[target_cb][1]},
             {invclov_packed[target_cb][0], invclov_packed[target_cb][1]}}};

        tmdslash.two_flav_dslash(two_flav_res,
                                 two_flav_psi,
                                 u_packed[target_cb],
                                 two_flav_invclov,
                                 isign,
                                 target_cb);

        // (c) Compare.
        spinor_result_base_up.unpack();
        spinor_result_base_down.unpack();
        spinor_result_slow_up.unpack();
        spinor_result_slow_down.unpack();

        print_spinor_differences(spinor_result_base_up.all,
                                 spinor_result_slow_up.all,
                                 geom,
                                 target_cb,
                                 isign,
                                 50);
        print_spinor_differences(spinor_result_base_down.all,
                                 spinor_result_slow_down.all,
                                 geom,
                                 target_cb,
                                 isign,
                                 50);
      } // cb
    } // isign
  }
#endif

#if 1 // This tests achimbdpsi
  {
    // 2. Test ax - bDslash y
    QDPIO::cout << std::endl << "TESTING A CHI - b D PSI" << std::endl;

    typedef PackedSpinor<FT, V, S, compress, Phi> MySpinor;

    MySpinor spinor_psi_up(geom);
    MySpinor spinor_psi_down(geom);
    MySpinor spinor_chi_up(geom);
    MySpinor spinor_chi_down(geom);
    MySpinor spinor_result_base_up(geom);
    MySpinor spinor_result_base_down(geom);
    MySpinor spinor_result_slow_up(geom);
    MySpinor spinor_result_slow_down(geom);

    for (auto isign : {1, -1}) {
      for (auto target_cb : {0, 1}) {
        int source_cb = 1 - target_cb;
        double beta = 0.25; // Always 0.25

        spinor_result_base_up.all = zero;
        spinor_result_base_up.pack();
        spinor_result_base_down.all = zero;
        spinor_result_base_down.pack();
        spinor_result_slow_up.all = zero;
        spinor_result_slow_up.pack();
        spinor_result_slow_down.all = zero;
        spinor_result_slow_down.pack();

        // (a) Apply the base version.
        tmdslash.dslashAChiMinusBDPsi(spinor_result_base_up[target_cb],
                                      spinor_chi_up[source_cb],
                                      spinor_psi_up[source_cb],
                                      u_packed[target_cb],
                                      (const FullClover **)clov_packed[target_cb],
                                      beta,
                                      isign,
                                      target_cb);

        tmdslash.dslashAChiMinusBDPsi(spinor_result_base_down[target_cb],
                                      spinor_chi_down[source_cb],
                                      spinor_psi_down[source_cb],
                                      u_packed[target_cb],
                                      (const FullClover **)clov_packed[target_cb],
                                      beta,
                                      isign,
                                      target_cb);

        // (b) Apply Martin's slow version.
        Spinor *two_flav_res[2] = {spinor_result_slow_up[target_cb],
                                   spinor_result_slow_down[target_cb]};
        const Spinor *two_flav_chi[2] = {spinor_chi_up[source_cb],
                                         spinor_chi_down[source_cb]};
        const Spinor *two_flav_psi[2] = {spinor_psi_up[source_cb],
                                         spinor_psi_down[source_cb]};
        const FullClover *two_flav_clov[2][2] = {
            {clov_packed[target_cb][0], clov_packed[target_cb][1]},
            {clov_packed[target_cb][0], clov_packed[target_cb][1]}};

        const auto epsilon = 0.0;
        tmdslash.two_flav_achimbdpsi(two_flav_res,
                                     two_flav_chi,
                                     two_flav_psi,
                                     u_packed[target_cb],
                                     two_flav_clov,
                                     beta,
                                     epsilon,
                                     isign,
                                     target_cb);

        // (c) Compare.
        spinor_result_base_up.unpack();
        spinor_result_base_down.unpack();
        spinor_result_slow_up.unpack();
        spinor_result_slow_down.unpack();

        print_spinor_differences(spinor_result_base_up.all,
                                 spinor_result_slow_up.all,
                                 geom,
                                 target_cb,
                                 isign,
                                 50);
        print_spinor_differences(spinor_result_base_down.all,
                                 spinor_result_slow_down.all,
                                 geom,
                                 target_cb,
                                 isign,
                                 50);

      } // cb
    } // isign
  }
#endif

#if 0 // This tests dslash with anti-periodic BCs
    {
        // Test only Dslash operator.
        // For clover this will be: A^{-1}_(1-cb,1-cb) D_(1-cb, cb)  psi_cb
        QDPIO::cout << std::endl
                    << std::endl
                    << "TESTING DSLASH WITH ANTI-PERIODIC BOUNDARY CONDITIONS"
                    << std::endl;

        // Create Antiperiodic Dslash
        t_boundary = (double)(-1);
        TMClovDslash<FT, V, S, compress> tmdslash_antip(
            &geom, t_boundary, aniso_fac_s, aniso_fac_t);
        NDTMClovDslash<FT, V, S, compress> ndtmdslash_antip(
            &geom, t_boundary, aniso_fac_s, aniso_fac_t);

        // Step 1: Convert u_test into one with antiperiodic BCs.
        // NB: This alone does not need a re-pack for the gauge field.
        QDPIO::cout << "Applying BCs to original gauge field" << std::endl;
        int mu_t = Nd - 1;

        // Let us make u_test to be u with antiperiodic_BCs
        for (int mu = 0; mu < Nd; mu++) {
            u_test[mu] = u[mu];
        }

        u_test[mu_t] *= where(Layout::latticeCoordinate(mu_t) ==
                                  (Layout::lattSize()[mu_t] - 1),
                              Real(t_boundary),
                              Real(1));

        QDPIO::cout << "Creating QPD++ Clover term using Gauge field with "
                       "antiperiodic BCs"
                    << std::endl;
        CloverTermT<Phi, U> clov_qdp_ap;
        clov_qdp_ap.create(u_test, clparam);
        QDPIO::cout << "Inverting QPD++ Clover Term..." << std::endl;
        CloverTermT<Phi, U> invclov_qdp_ap(clov_qdp_ap);
        for (int cb = 0; cb < 2; cb++) {
            invclov_qdp_ap.choles(cb);
        }
        QDPIO::cout << "Done." << std::endl;

        // Now we need to repack this.
        QDPIO::cout << "Re-Packing QPhiX Twisted-Mass Clover Term..."
                    << std::endl;
        for (int cb = 0; cb < 2; cb++) {
            for (int pm = 0; pm < 2; pm++) {
                qdp_pack_full_clover<>(
                    clov_qdp_ap, clov_packed[cb][pm], geom, cb);
            }
        }
        QDPIO::cout << "...and QPhiX Clover Term Inverse..." << std::endl;
        for (int cb = 0; cb < 2; cb++) {
            for (int pm = 0; pm < 2; pm++) {
                qdp_pack_full_clover<>(
                    invclov_qdp_ap, invclov_packed[cb][pm], geom, cb);
            }
        }
        QDPIO::cout << "Done." << std::endl;

        QDPIO::cout << "Folding aniso factors into gauge field for testing"
                    << endl;
        for (int mu = 0; mu < Nd; mu++) {
            Real factor;
            if (mu == 3) {
                factor = Real(aniso_fac_t);
            } else {
                factor = Real(aniso_fac_s);
            }
            u_test[mu] *= factor;
        }

        // NB: Gauge field doesn't need to be repacked. It is the original 'u'
        // without the aniso factors or boundaries
        // As these are now handled in the Dslash.

        QDPIO::cout << std::endl
                    << "TESTING DSLASH w/ ANTI-PERIODIC BCs:" << std::endl;

        for (auto isign : {1, -1}) {
            for (auto cb : {0, 1}) {
                int source_cb = 1 - cb;
                int target_cb = cb;

                // (a) Apply TM Dslash, Pack & Unpack
                clov_chi = zero;
                qdp_pack_spinor<>(chi, chi_even, chi_odd, geom);
                tmdslash_antip.dslash(chi_s[target_cb],
                                psi_s[source_cb],
                                u_packed[target_cb],
                                (const FullClover **)invclov_packed[target_cb],
                                isign,
                                target_cb);
                qdp_unpack_spinor<>(chi_even, chi_odd, clov_chi, geom);

                // (b) Apply NDTM Dslash, Pack & Unpack
                clov_chi2 = zero;
                qdp_pack_spinor<>(chi2, chi_even, chi_odd, geom);
                ndtmdslash_antip.dslash(
                    chi_s[target_cb],
                    psi_s[source_cb],
                    u_packed[target_cb],
                    (const FullClover **)invclov_packed[target_cb],
                    isign,
                    target_cb);
                qdp_unpack_spinor<>(chi_even, chi_odd, clov_chi2, geom);

                print_spinor_differences(
                    clov_chi, clov_chi2, geom, target_cb, isign);
            } // cb
        } // isign
    }
#endif

#if 0 // This tests achimbdpsi with anti-periodic BCs
    {
        QDPIO::cout << std::endl
                    << "TESTING A CHI - b D PSI w/ ANTI-PERIODIC BCs"
                    << std::endl;

        for (auto isign : {1, -1}) {
            for (auto cb : {0, 1}) {

                int source_cb = 1 - cb;
                int target_cb = cb;
                double beta = (double)(0.25); // Always 0.25

                // (a) Apply Optimized Dslash, Pack & Unpack
                chi = zero;
                qdp_pack_spinor<>(chi, chi_even, chi_odd, geom);
                tmdslash.dslashAChiMinusBDPsi(
                    chi_s[target_cb],
                    psi_s[source_cb],
                    psi_s[target_cb],
                    u_packed[target_cb],
                    (const FullClover **)clov_packed[target_cb],
                    beta,
                    isign,
                    target_cb);
                qdp_unpack_spinor<>(chi_s[0], chi_s[1], chi, geom);

                // (a) Apply Optimized Dslash, Pack & Unpack
                chi2 = zero;
                qdp_pack_spinor<>(chi2, chi_even, chi_odd, geom);
                ndtmdslash.dslashAChiMinusBDPsi(
                    chi_s[target_cb],
                    psi_s[source_cb],
                    psi_s[target_cb],
                    u_packed[target_cb],
                    (const FullClover **)clov_packed[target_cb],
                    beta,
                    isign,
                    target_cb);
                qdp_unpack_spinor<>(chi_s[0], chi_s[1], chi2, geom);

                print_spinor_differences(chi, chi2, geom, target_cb, isign);
            } // cb
        } // isign
    }
#endif

  {
#if 0 // This tests the full even/odd operators
        QDPIO::cout << std::endl
                    << std::endl
                    << "TESTING TWISTED MASS EVEN/ODD (DEGENERATE) DSLASH "
                       "CLOVER OPERATOR"
                    << std::endl;
        t_boundary = (double)(-1);
        EvenOddTMCloverOperator<FT, V, S, compress> tm_op(u_packed,
                                                          clov_packed[1],
                                                          invclov_packed[0],
                                                          &geom,
                                                          t_boundary,
                                                          aniso_fac_s,
                                                          aniso_fac_t);
        EvenOddNDTMCloverOperator<FT, V, S, compress> ndtm_op(u_packed,
                                                              clov_packed[1],
                                                              invclov_packed[0],
                                                              &geom,
                                                              t_boundary,
                                                              aniso_fac_s,
                                                              aniso_fac_t);
        Phi ltmp = zero;
        Real betaFactor = Real(0.25);

        for (auto isign : {1, -1}) {
            // (a) Apply TM Operator
            chi = zero;
            qdp_pack_spinor<>(chi, chi_even, chi_odd, geom);
            tm_op(chi_s[1], psi_s[1], isign);
            qdp_unpack_spinor<>(chi_s[0], chi_s[1], chi, geom);

            // (b) Apply NDTM Operator
            chi2 = zero;
            qdp_pack_spinor<>(chi2, chi_even, chi_odd, geom);
            ndtm_op(chi_s[1], psi_s[1], isign);
            qdp_unpack_spinor<>(chi_s[0], chi_s[1], chi2, geom);

            print_spinor_differences(chi, chi2, geom, -1, isign);
        } // isign
#endif

#if 0 // This tests the CG
        {
            QDPIO::cout << std::endl
                        << std::endl
                        << "TESTING CG SOLVER W/ TWISTED MASS CLOVER OPERATOR"
                        << std::endl;

            chi = zero;
            qdp_pack_cb_spinor<>(chi, chi_s[1], geom, 1);

            double rsd_target = rsdTarget<FT>::value;
            int max_iters = 500;
            int niters;
            double rsd_final;
            unsigned long site_flops;
            unsigned long mv_apps;

            InvCG<FT, V, S, compress> solver(tm_op, max_iters);
            solver.tune();

            double start = omp_get_wtime();
            solver(chi_s[1],
                   psi_s[1],
                   rsd_target,
                   niters,
                   rsd_final,
                   site_flops,
                   mv_apps,
                   1,
                   verbose);
            double end = omp_get_wtime();

            qdp_unpack_cb_spinor<>(chi_s[1], chi, geom, 1);

            CloverTermT<Phi, U> clov_qdp_ap;
            clov_qdp_ap.create(u_test, clparam);
            CloverTermT<Phi, U> invclov_qdp_ap(clov_qdp_ap);

            // Multiply back
            // chi2 = M chi
            dslash(chi2, u_test, chi, 1, 0);
            invclov_qdp_ap.apply(clov_chi2, chi2, 1, 0);
            dslash(ltmp, u_test, clov_chi2, 1, 1);

            clov_qdp_ap.apply(chi2, chi, 1, 1);
            chi2[rb[1]] -= betaFactor * ltmp;

            // chi3 = M^\dagger chi2
            Phi chi3 = zero;
            dslash(chi3, u_test, chi2, -1, 0);
            invclov_qdp_ap.apply(clov_chi2, chi3, -1, 0);
            dslash(ltmp, u_test, clov_chi2, -1, 1);

            clov_qdp_ap.apply(chi3, chi2, -1, 1);
            chi3[rb[1]] -= betaFactor * ltmp;

            //  dslash(chi3,u,chi2, (-1), 1);
            // dslash(ltmp,u,chi3, (-1), 0);
            // chi3[rb[0]] = massFactor*chi2 - betaFactor*ltmp;

            Phi diff = chi3 - psi;
            QDPIO::cout << "True norm is: "
                        << sqrt(norm2(diff, rb[1]) / norm2(psi, rb[1])) << endl;

            int Nxh = Nx / 2;
            unsigned long num_cb_sites = Nxh * Ny * Nz * Nt;

            unsigned long total_flops =
                (site_flops + (1320 + 504 + 1320 + 504 + 48) * mv_apps) *
                num_cb_sites;

            masterPrintf("GFLOPS=%e\n",
                         1.0e-9 * (double)(total_flops) / (end - start));
        }
#endif

// This tests the BiCG
#if 0
        {
            QDPIO::cout
                << std::endl
                << std::endl
                << "TESTING BICGSTAB SOLVER W/ TWISTED MASS CLOVER OPERATOR"
                << std::endl;

            chi = zero;
            qdp_pack_spinor<>(chi, chi_even, chi_odd, geom);

            double rsd_target = rsdTarget<FT>::value;
            int max_iters = 500;
            int niters;
            double rsd_final;
            unsigned long site_flops;
            unsigned long mv_apps;

            InvBiCGStab<FT, V, S, compress> solver(tm_op, max_iters);
            solver.tune();
            const int isign = 1;
            double start = omp_get_wtime();
            solver(chi_s[1],
                   psi_s[1],
                   rsd_target,
                   niters,
                   rsd_final,
                   site_flops,
                   mv_apps,
                   isign,
                   verbose);
            double end = omp_get_wtime();

            qdp_unpack_cb_spinor<>(chi_s[1], chi, geom, 1);

            CloverTermT<Phi, U> clov_qdp_ap;
            clov_qdp_ap.create(u_test, clparam);
            CloverTermT<Phi, U> invclov_qdp_ap(clov_qdp_ap);

            // Multiply back
            // chi2 = M chi
            dslash(chi2, u_test, chi, 1, 0);
            invclov_qdp_ap.apply(clov_chi2, chi2, 1, 0);
            dslash(ltmp, u_test, clov_chi2, 1, 1);

            clov_qdp_ap.apply(chi2, chi, 1, 1);
            chi2[rb[1]] -= betaFactor * ltmp;

            Phi diff = chi2 - psi;
            QDPIO::cout << "True norm is: "
                        << sqrt(norm2(diff, rb[1]) / norm2(psi, rb[1])) << endl;

            int Nxh = Nx / 2;
            unsigned long num_cb_sites = Nxh * Ny * Nz * Nt;

            unsigned long total_flops =
                (site_flops + (1320 + 504 + 1320 + 504 + 48) * mv_apps) *
                num_cb_sites;
            masterPrintf("GFLOPS=%e\n",
                         1.0e-9 * (double)(total_flops) / (end - start));
        }
#endif
  }

  QDPIO::cout << std::endl;

  geom.free(packed_gauge_cb0);
  geom.free(packed_gauge_cb1);
  geom.free(psi_even);
  geom.free(psi_odd);
  geom.free(chi_even);
  geom.free(chi_odd);
  geom.free(A_cb0_plus);
  geom.free(A_cb0_minus);
  geom.free(A_cb1_plus);
  geom.free(A_cb1_minus);
  geom.free(A_inv_cb0_plus);
  geom.free(A_inv_cb0_minus);
  geom.free(A_inv_cb1_plus);
  geom.free(A_inv_cb1_minus);
}

void testNDTMvsTMCloverFull::run(void)
{
  typedef LatticeColorMatrixF UF;
  typedef LatticeDiracFermionF PhiF;

  typedef LatticeColorMatrixD UD;
  typedef LatticeDiracFermionD PhiD;

  QDPIO::cout << std::endl;

#if defined(QPHIX_SCALAR_SOURCE)
  if (precision == FLOAT_PREC) {
    QDPIO::cout << "SINGLE PRECISION TESTING " << endl;
    if (compress12) {
      runTest<float, 1, 1, true, UF, PhiF>();
    } else {
      runTest<float, 1, 1, false, UF, PhiF>();
    }
  }
  if (precision == DOUBLE_PREC) {
    QDPIO::cout << "DOUBLE PRECISION TESTING " << endl;
    if (compress12) {
      runTest<double, 1, 1, true, UF, PhiF>();
    } else {
      runTest<double, 1, 1, false, UF, PhiF>();
    }
  }
#else

  if (precision == FLOAT_PREC) {
#if defined(QPHIX_AVX_SOURCE) || defined(QPHIX_AVX2_SOURCE) ||                      \
    defined(QPHIX_MIC_SOURCE) || defined(QPHIX_AVX512_SOURCE) ||                    \
    defined(QPHIX_SSE_SOURCE)
    QDPIO::cout << "SINGLE PRECISION TESTING " << endl;
    if (compress12) {
      runTest<float, VECLEN_SP, 4, true, UF, PhiF>();
    } else {
      runTest<float, VECLEN_SP, 4, false, UF, PhiF>();
    }

#if !defined(QPHIX_SSE_SOURCE)
    if (compress12) {
      runTest<float, VECLEN_SP, 8, true, UF, PhiF>();
    } else {
      runTest<float, VECLEN_SP, 8, false, UF, PhiF>();
    }
#endif

#if defined(QPHIX_MIC_SOURCE) || defined(QPHIX_AVX512_SOURCE)
    if (compress12) {
      runTest<float, VECLEN_SP, 16, true, UF, PhiF>();
    } else {
      runTest<float, VECLEN_SP, 16, false, UF, PhiF>();
    }
#endif
#endif // QPHIX_AVX_SOURCE|| QPHIX_AVX2_SOURCE|| QPHIX_MIC_SOURCE |
    // QPHIX_AVX512_SOURCE
  }

  if (precision == HALF_PREC) {
#if defined(QPHIX_MIC_SOURCE) || defined(QPHIX_AVX512_SOURCE) ||                    \
    defined(QPHIX_AVX2_SOURCE)
    QDPIO::cout << "SINGLE PRECISION TESTING " << endl;
    if (compress12) {
      runTest<half, VECLEN_HP, 4, true, UF, PhiF>();
    } else {
      runTest<half, VECLEN_HP, 4, false, UF, PhiF>();
    }

    if (compress12) {
      runTest<half, VECLEN_HP, 8, true, UF, PhiF>();
    } else {
      runTest<half, VECLEN_HP, 8, false, UF, PhiF>();
    }
#if defined(QPHIX_MIC_SOURCE) || defined(QPHIX_AVX512_SOURCE)
    if (compress12) {
      runTest<half, VECLEN_HP, 16, true, UF, PhiF>();
    } else {
      runTest<half, VECLEN_HP, 16, false, UF, PhiF>();
    }
#endif
#else
    QDPIO::cout << "Half precision tests are not available in this build. "
                   "Currently only in MIC builds"
                << endl;
#endif
  } // HALF_PREC

  if (precision == DOUBLE_PREC) {
    QDPIO::cout << "DOUBLE PRECISION TESTING" << endl;

#if defined(QPHIX_AVX_SOURCE) || defined(QPHIX_AVX2_SOURCE)
    // Only AVX can do DP 2
    if (compress12) {
      runTest<double, VECLEN_DP, 2, true, UD, PhiD>();
    } else {
      runTest<double, VECLEN_DP, 2, false, UD, PhiD>();
    }
#endif

    if (compress12) {
      runTest<double, VECLEN_DP, 4, true, UD, PhiD>();
    } else {
      runTest<double, VECLEN_DP, 4, false, UD, PhiD>();
    }

#if defined(QPHIX_MIC_SOURCE) || defined(QPHIX_AVX512_SOURCE)
    // Only MIC can do DP 8
    if (compress12) {
      runTest<double, VECLEN_DP, 8, true, UD, PhiD>();
    } else {
      runTest<double, VECLEN_DP, 8, false, UD, PhiD>();
    }
#endif
  } // DOUBLE PREC

#endif // SCALAR SOURCE
}
