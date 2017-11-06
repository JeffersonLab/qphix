#include <qphix/qphix_config_internal.h>
#include "../tests/prec.h"
#include "../tests/tparam_selector.h"

using namespace QPhiX;

#include "../tests/veclen.h"
#include "../tests/tolerance.h"
#include "../tests/clover_fermact_params_w.h"
#include "../tests/clover_term.h"
#include "../tests/compare_qdp_spinors_gtest.h"

#include <qphix/clover.h>
#include <qphix/clover_product.h>
#include <qphix/geometry.h>
#include <qphix/qdp_packer.h>
#include <qphix/qphix_cli_args.h>
#include <qphix/wilson.h>

#include <qdp.h>

#include <iostream>

void reunit(QDP::LatticeColorMatrixF &a);
void reunit(QDP::LatticeColorMatrixD &a);

extern int nrow_in[];
extern QPhiX::QPhiXCLIArgs some_user_args;

template <typename G>
class CloverProductTest : public ::testing::Test
{
 public:
  typedef QDP::LatticeColorMatrixD U;
  typedef QDP::LatticeDiracFermionD Phi;
  typedef G Geom;
  typedef typename Geom::FourSpinorBlock Spinor;
  typedef typename Geom::SU3MatrixBlock Gauge;
  typedef typename Geom::CloverBlock Clover;
  typedef typename Geom::FullCloverBlock FullClover;

  typedef typename Geom::FT FT;
  static int constexpr veclen = Geom::veclen;
  static int constexpr soalen = Geom::soalen;
  static bool constexpr compress12 = Geom::compress12;

  double xi_0_f;
  double nu_f;
  double aniso_fac_s;
  double aniso_fac_t;

  CloverProductTest()
      : geom(Layout::subgridLattSize().slice(),
             some_user_args.getBy(),
             some_user_args.getBz(),
             some_user_args.getNCores(),
             some_user_args.getSy(),
             some_user_args.getSz(),
             some_user_args.getPxy(),
             some_user_args.getPxyz(),
             some_user_args.getMinCt()),
        xi_0_f(0.3), nu_f(1.4), aniso_fac_s(nu_f / xi_0_f),
        D32(&geom, 1, aniso_fac_s, aniso_fac_t),
        cD32(&geom, 1, aniso_fac_s, aniso_fac_t), packed_gauge_cb0(geom),
        packed_gauge_cb1(geom), psi_even(geom), psi_odd(geom), chi_even(geom),
        chi_odd(geom), tmp(geom), A_cb0(geom), A_cb1(geom), A_inv_cb0(geom),
        A_inv_cb1(geom), full_A_cb0(geom), full_A_cb1(geom), full_A_inv_cb0(geom),
        full_A_inv_cb1(geom)
  {
    // Make a random gauge field
    // Start up the gauge field somehow
    // We choose: u = reunit(  1 + factor*gaussian(u) );
    // Adjust factor to get reasonable inversion times in invert test.
    // Bug gauge field is always unitary
    QDP::multi1d<U> u(4);
    {
      U g;
      U uf;
      for (int mu = 0; mu < 4; mu++) {
#if 1
        uf = 1; // Unit gauge

        QDP::Real factor = QDP::Real(0.08);
        QDP::gaussian(g);
        u[mu] = uf + factor * g;
        reunit(u[mu]);
#else
        // Only Unit Gauge: testing....
        u[mu] = 1;
#endif
      }
    }

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
    clparam.Mass = QDP::Real(0.1);

    // Some random clover coeffs
    // This should be wilson...
    clparam.clovCoeffR = QDP::Real(1.2);
    clparam.clovCoeffT = QDP::Real(0.9);

    // Make a random source
    gaussian(psi);

    invclov_packed[0] = A_inv_cb0.get();
    invclov_packed[1] = A_inv_cb1.get();
    clov_packed[0] = A_cb0.get();
    clov_packed[1] = A_cb1.get();

    full_invclov_packed[0] = full_A_inv_cb0.get();
    full_invclov_packed[1] = full_A_inv_cb1.get();
    full_clov_packed[0] = full_A_cb0.get();
    full_clov_packed[1] = full_A_cb1.get();

    // Pack the gauge field
    QPhiX::qdp_pack_gauge<>(u, packed_gauge_cb0.get(), packed_gauge_cb1.get(), geom);
    u_packed[0] = packed_gauge_cb0.get();
    u_packed[1] = packed_gauge_cb1.get();

    psi_s[0] = psi_even.get();
    psi_s[1] = psi_odd.get();
    chi_s[0] = chi_even.get();
    chi_s[1] = chi_odd.get();

    QPhiX::qdp_pack_spinor<>(psi, psi_even.get(), psi_odd.get(), geom);

    // Clover term deals with anisotropy internally -- so use original u
    // field.
    clov_qdp.create(u, clparam);
    invclov_qdp = clov_qdp;

    for (int cb = 0; cb < 2; cb++) {
      invclov_qdp.choles(cb);
    }

    for (int cb = 0; cb < 2; cb++) {
      QPhiX::qdp_pack_clover<>(invclov_qdp, invclov_packed[cb], geom, cb);
      QPhiX::qdp_pack_full_clover<>(invclov_qdp, full_invclov_packed[cb], geom, cb);
    }

    for (int cb = 0; cb < 2; cb++) {
      QPhiX::qdp_pack_clover<>(clov_qdp, clov_packed[cb], geom, cb);
      QPhiX::qdp_pack_full_clover<>(clov_qdp, full_clov_packed[cb], geom, cb);
    }

    multi1d<U> u_test(4);
    for (int mu = 0; mu < Nd; mu++) {
      QDP::Real factor;
      if (mu == 3) {
        factor = QDP::Real(aniso_fac_t);
      } else {
        factor = QDP::Real(aniso_fac_s);
      }
      u_test[mu] = factor * u[mu];
    }

    const int cbsize_in_blocks = rb[0].numSiteTable() / Geom::soalen;
  }

  Geom geom;

  QPhiX::ClovDslash<typename Geom::FT, Geom::veclen, Geom::soalen, Geom::compress12>
      D32;

  QPhiX::Dslash<typename Geom::FT, Geom::veclen, Geom::soalen, Geom::compress12>
      cD32;

  Phi chi, clov_chi;
  Phi psi, chi2, clov_chi2;

  QPhiX::GaugeHandle<FT, veclen, soalen, compress12> packed_gauge_cb0;
  QPhiX::GaugeHandle<FT, veclen, soalen, compress12> packed_gauge_cb1;

  QPhiX::FourSpinorHandle<FT, veclen, soalen, compress12> psi_even;
  QPhiX::FourSpinorHandle<FT, veclen, soalen, compress12> psi_odd;
  QPhiX::FourSpinorHandle<FT, veclen, soalen, compress12> chi_even;
  QPhiX::FourSpinorHandle<FT, veclen, soalen, compress12> chi_odd;
  QPhiX::FourSpinorHandle<FT, veclen, soalen, compress12> tmp;

  QPhiX::CloverHandle<FT, veclen, soalen, compress12> A_cb0;
  QPhiX::CloverHandle<FT, veclen, soalen, compress12> A_cb1;
  QPhiX::CloverHandle<FT, veclen, soalen, compress12> A_inv_cb0;
  QPhiX::CloverHandle<FT, veclen, soalen, compress12> A_inv_cb1;

  QPhiX::FullCloverHandle<FT, veclen, soalen, compress12> full_A_cb0;
  QPhiX::FullCloverHandle<FT, veclen, soalen, compress12> full_A_cb1;
  QPhiX::FullCloverHandle<FT, veclen, soalen, compress12> full_A_inv_cb0;
  QPhiX::FullCloverHandle<FT, veclen, soalen, compress12> full_A_inv_cb1;

  Clover *invclov_packed[2];
  Clover *clov_packed[2];
  FullClover *full_invclov_packed[2];
  FullClover *full_clov_packed[2];

  Gauge *u_packed[2];

  Spinor *psi_s[2];
  Spinor *chi_s[2];

  QPhiX::CloverTermT<Phi, U> clov_qdp;
  QPhiX::CloverTermT<Phi, U> invclov_qdp;
};

typedef ::testing::Types<QPhiX::Geometry<double, VECLEN_DP, VECLEN_DP, true>>
    MyTypes;

TYPED_TEST_CASE(CloverProductTest, MyTypes);

TYPED_TEST(CloverProductTest, CombinedDslashVersusTwoParts)
{
  for (int isign = 1; isign >= -1; isign -= 2) {
    for (int cb = 0; cb < 2; cb++) {
      int source_cb = 1 - cb;
      int target_cb = cb;

      // Apply the Wilson Clover Dslash.
      this->chi = zero;
      QPhiX::qdp_pack_spinor<>(
          this->chi, this->chi_even.get(), this->chi_odd.get(), this->geom);
      this->D32.dslash(this->chi_s[target_cb],
                       this->psi_s[source_cb],
                       this->u_packed[target_cb],
                       this->invclov_packed[target_cb],
                       isign,
                       target_cb);
      QPhiX::qdp_unpack_spinor<>(
          this->chi_even.get(), this->chi_odd.get(), this->clov_chi, this->geom);

      // Apply the Wilson Dslash and then apply the new clover inverse
      // multiplication..
      this->chi = zero;
      QPhiX::qdp_pack_spinor<>(
          this->chi, this->chi_even.get(), this->chi_odd.get(), this->geom);
      this->cD32.dslash(this->tmp.get(),
                        this->psi_s[source_cb],
                        this->u_packed[target_cb],
                        isign,
                        target_cb);
      clover_product(this->chi_s[target_cb],
                     this->tmp.get(),
                     this->invclov_packed[target_cb],
                     this->geom);
      QPhiX::qdp_unpack_spinor<>(
          this->chi_even.get(), this->chi_odd.get(), this->clov_chi2, this->geom);

      expect_near(this->clov_chi,
                  this->clov_chi2,
                  toDouble(tolerance<typename TestFixture::FT>::small),
                  this->geom,
                  target_cb);
    }
  }
}

TYPED_TEST(CloverProductTest, CombinedDslashVersusTwoPartsFull)
{
  for (int isign = 1; isign >= -1; isign -= 2) {
    for (int cb = 0; cb < 2; cb++) {
      int source_cb = 1 - cb;
      int target_cb = cb;

      // Apply the Wilson Clover Dslash.
      this->chi = zero;
      QPhiX::qdp_pack_spinor<>(
          this->chi, this->chi_even.get(), this->chi_odd.get(), this->geom);
      this->D32.dslash(this->chi_s[target_cb],
                       this->psi_s[source_cb],
                       this->u_packed[target_cb],
                       this->invclov_packed[target_cb],
                       isign,
                       target_cb);
      QPhiX::qdp_unpack_spinor<>(
          this->chi_even.get(), this->chi_odd.get(), this->clov_chi, this->geom);

      // Apply the Wilson Dslash and then apply the new clover inverse
      // multiplication..
      this->chi = zero;
      QPhiX::qdp_pack_spinor<>(
          this->chi, this->chi_even.get(), this->chi_odd.get(), this->geom);
      this->cD32.dslash(this->tmp.get(),
                        this->psi_s[source_cb],
                        this->u_packed[target_cb],
                        isign,
                        target_cb);
      clover_product(this->chi_s[target_cb],
                     this->tmp.get(),
                     this->full_invclov_packed[target_cb],
                     this->geom);
      QPhiX::qdp_unpack_spinor<>(
          this->chi_even.get(), this->chi_odd.get(), this->clov_chi2, this->geom);

      expect_near(this->clov_chi,
                  this->clov_chi2,
                  toDouble(tolerance<typename TestFixture::FT>::small),
                  this->geom,
                  target_cb);
    }
  }
}

TYPED_TEST(CloverProductTest, AinvA)
{
  for (int isign = 1; isign >= -1; isign -= 2) {
    for (int cb = 0; cb < 2; cb++) {
      int source_cb = 1 - cb;
      int target_cb = cb;

      this->chi = zero;
      QPhiX::qdp_pack_spinor<>(
          this->chi, this->chi_even.get(), this->chi_odd.get(), this->geom);
      clover_product(this->tmp.get(),
                     this->psi_s[target_cb],
                     this->clov_packed[target_cb],
                     this->geom);
      clover_product(this->chi_s[target_cb],
                     this->tmp.get(),
                     this->invclov_packed[target_cb],
                     this->geom);
      QPhiX::qdp_unpack_spinor<>(
          this->chi_even.get(), this->chi_odd.get(), this->clov_chi2, this->geom);

      expect_near(this->psi,
                  this->clov_chi2,
                  toDouble(tolerance<typename TestFixture::FT>::small),
                  this->geom,
                  target_cb);
    }
  }
}

TYPED_TEST(CloverProductTest, AinvAFull)
{
  for (int isign = 1; isign >= -1; isign -= 2) {
    for (int cb = 0; cb < 2; cb++) {
      int source_cb = 1 - cb;
      int target_cb = cb;

      this->chi = zero;
      QPhiX::qdp_pack_spinor<>(
          this->chi, this->chi_even.get(), this->chi_odd.get(), this->geom);
      clover_product(this->tmp.get(),
                     this->psi_s[target_cb],
                     this->full_clov_packed[target_cb],
                     this->geom);
      clover_product(this->chi_s[target_cb],
                     this->tmp.get(),
                     this->full_invclov_packed[target_cb],
                     this->geom);
      QPhiX::qdp_unpack_spinor<>(
          this->chi_even.get(), this->chi_odd.get(), this->clov_chi2, this->geom);

      expect_near(this->psi,
                  this->clov_chi2,
                  toDouble(tolerance<typename TestFixture::FT>::small),
                  this->geom,
                  target_cb);
    }
  }
}
