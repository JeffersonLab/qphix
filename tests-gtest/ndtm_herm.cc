#include <qphix/qphix_config_internal.h>
#include "../tests/prec.h"
#include "../tests/tparam_selector.h"

using namespace QPhiX;

#include "../tests/veclen.h"
#include "../tests/tolerance.h"
#include "../tests/compare_qdp_spinors_gtest.h"

#include <qphix/blas_new_c.h>
#include <qphix/ndtm_reuse_operator.h>
#include <qphix/geometry.h>
#include <qphix/qdp_packer.h>
#include <qphix/qphix_cli_args.h>
#include <qphix/wilson.h>

#include <qdp.h>


#include <iostream>
#include <vector>

void reunit(QDP::LatticeColorMatrixF &a);
void reunit(QDP::LatticeColorMatrixD &a);

extern int nrow_in[];
extern QPhiX::QPhiXCLIArgs some_user_args;

template <typename G>
class NDTMReuseOperatorHermTest : public ::testing::Test
{
 public:
  typedef QDP::LatticeColorMatrixD QDPGauge;
  typedef QDP::LatticeDiracFermionD QDPSpinor;
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

  NDTMReuseOperatorHermTest()
      : geom(Layout::subgridLattSize().slice(),
             some_user_args.getBy(),
             some_user_args.getBz(),
             some_user_args.getNCores(),
             some_user_args.getSy(),
             some_user_args.getSz(),
             some_user_args.getPxy(),
             some_user_args.getPxyz(),
             some_user_args.getMinCt()),
        xi_0_f(1.0), nu_f(1.0), aniso_fac_s(nu_f / xi_0_f), aniso_fac_t(1.0),
        packed_gauge_cb0(geom), packed_gauge_cb1(geom),
        u_packed{packed_gauge_cb0.get(), packed_gauge_cb1.get()}, NDTMEO(nullptr),
        wilson(nullptr)
  {
    NDTMEO = new EvenOddNDTMWilsonReuseOperator<typename Geom::FT,
                                                Geom::veclen,
                                                Geom::soalen,
                                                Geom::compress12>(
        -1, 0.1, 0.05, u_packed, &geom, 1, aniso_fac_s, aniso_fac_t);

    wilson = new EvenOddWilsonOperator<typename Geom::FT,
                                       Geom::veclen,
                                       Geom::soalen,
                                       Geom::compress12>(
        -1, u_packed, &geom, 1, aniso_fac_s, aniso_fac_t);

    constexpr int cb_even = 0;
    constexpr int cb_odd = 1;
    // Make a random gauge field
    // Start up the gauge field somehow
    // We choose: u = reunit(  1 + factor*gaussian(u) );
    // Adjust factor to get reasonable inversion times in invert test.
    // Bug gauge field is always unitary
    Seed rng_seed = 123;
    RNG::setrn(rng_seed);

    QDP::multi1d<QDPGauge> u(4);
    {
      QDPGauge g;
      QDPGauge uf;
      for (int mu = 0; mu < 4; mu++) {
        uf = 1; // Unit gauge

        QDP::Real factor = QDP::Real(0.08);
        QDP::gaussian(g);
#if 1
        u[mu] = uf + factor * g;
        reunit(u[mu]);
#else
        u[mu] = uf;
#endif
      }
    }

    // Make a random source
    gaussian(chi_qdp_f0);
    gaussian(chi_qdp_f1);

    gaussian(chi2_qdp_f0);
    gaussian(chi2_qdp_f1);

    // Pack the gauge field
    QPhiX::qdp_pack_gauge<>(u, packed_gauge_cb0.get(), packed_gauge_cb1.get(), geom);

    // flavour index on the inside to allow passing to two-flavour functions
    for (int cb : {cb_even, cb_odd}) {
      qhandles.push_back(makeFourSpinorHandle(geom));
      MwChi[cb] = qhandles.back().get();
      qhandles.push_back(makeFourSpinorHandle(geom));
      MwChi2[cb] = qhandles.back().get();
      qhandles.push_back(makeFourSpinorHandle(geom));
      MwDagChi2[cb] = qhandles.back().get();
      qhandles.push_back(makeFourSpinorHandle(geom));
      MwDagMwChi2[cb] = qhandles.back().get();
      for (int f : {0, 1}) {
        qhandles.push_back(makeFourSpinorHandle(geom));
        chi[cb][f] = qhandles.back().get();
        qhandles.push_back(makeFourSpinorHandle(geom));
        chi2[cb][f] = qhandles.back().get();

        qhandles.push_back(makeFourSpinorHandle(geom));
        Mchi[cb][f] = qhandles.back().get();
        qhandles.push_back(makeFourSpinorHandle(geom));
        Mchi2[cb][f] = qhandles.back().get();
        qhandles.push_back(makeFourSpinorHandle(geom));
        MdagMchi2[cb][f] = qhandles.back().get();
        qhandles.push_back(makeFourSpinorHandle(geom));
        MdagChi2[cb][f] = qhandles.back().get();
      }
    }

    // pack the gaussian sources
    QPhiX::qdp_pack_spinor<>(chi_qdp_f0, chi[cb_even][0], chi[cb_odd][0], geom);
    QPhiX::qdp_pack_spinor<>(chi_qdp_f1, chi[cb_even][1], chi[cb_odd][1], geom);
    QPhiX::qdp_pack_spinor<>(chi2_qdp_f0, chi2[cb_even][0], chi2[cb_odd][0], geom);
    QPhiX::qdp_pack_spinor<>(chi2_qdp_f1, chi2[cb_even][1], chi2[cb_odd][1], geom);

    const int cbsize_in_blocks = rb[0].numSiteTable() / Geom::soalen;
  }

  ~NDTMReuseOperatorHermTest()
  {
    delete NDTMEO;
    delete wilson;
  }

  Geom geom;

  QPhiX::EvenOddNDTMWilsonReuseOperator<typename Geom::FT,
                                        Geom::veclen,
                                        Geom::soalen,
                                        Geom::compress12> *NDTMEO;
  QPhiX::EvenOddWilsonOperator<typename Geom::FT,
                               Geom::veclen,
                               Geom::soalen,
                               Geom::compress12> *wilson;

  QDPSpinor chi_qdp_f0;
  QDPSpinor chi_qdp_f1;
  QDPSpinor chi2_qdp_f0;
  QDPSpinor chi2_qdp_f1;
  QDPSpinor psi_qdp_f0;
  QDPSpinor psi_qdp_f1;

  QPhiX::GaugeHandle<FT, veclen, soalen, compress12> packed_gauge_cb0;
  QPhiX::GaugeHandle<FT, veclen, soalen, compress12> packed_gauge_cb1;

  std::vector<QPhiX::FourSpinorHandle<FT, veclen, soalen, compress12>> qhandles;

  Gauge *u_packed[2];

  Spinor *chi[2][2];
  Spinor *chi2[2][2];
  Spinor *Mchi[2][2];
  Spinor *Mchi2[2][2];
  Spinor *MdagMchi2[2][2];
  Spinor *MdagChi2[2][2];

  Spinor *MwChi[2];
  Spinor *MwChi2[2];
  Spinor *MwDagChi2[2];
  Spinor *MwDagMwChi2[2];
};

typedef ::testing::Types<QPhiX::Geometry<double, VECLEN_DP, VECLEN_DP, true>> MyTypes;

TYPED_TEST_CASE(NDTMReuseOperatorHermTest, MyTypes);

TYPED_TEST(NDTMReuseOperatorHermTest, GaussianSource)
{
  constexpr int cb_even = 0;
  constexpr int cb_odd = 1;
  for (int isign = 1; isign >= -1; isign -= 2) {
    for (int cb = 0; cb < 2; cb++) {
      this->psi_qdp_f0 = zero;
      this->psi_qdp_f1 = zero;
      // zero output spinors on both checkerboards
      QPhiX::qdp_pack_spinor<>(
          this->psi_qdp_f0, this->MwChi[cb_even], this->MwChi[cb_odd], this->geom);
      QPhiX::qdp_pack_spinor<>(
          this->psi_qdp_f0, this->MwChi2[cb_even], this->MwChi2[cb_odd], this->geom);
      QPhiX::qdp_pack_spinor<>(this->psi_qdp_f0,
                               this->MwDagChi2[cb_even],
                               this->MwDagChi2[cb_odd],
                               this->geom);
      QPhiX::qdp_pack_spinor<>(this->psi_qdp_f0,
                               this->MwDagMwChi2[cb_even],
                               this->MwDagMwChi2[cb_odd],
                               this->geom);

      QPhiX::qdp_pack_spinor<>(
          this->psi_qdp_f0, this->Mchi[cb_even][0], this->Mchi[cb_odd][0], this->geom);
      QPhiX::qdp_pack_spinor<>(
          this->psi_qdp_f1, this->Mchi[cb_even][1], this->Mchi[cb_odd][1], this->geom);

      QPhiX::qdp_pack_spinor<>(
          this->psi_qdp_f0, this->Mchi2[cb_even][0], this->Mchi2[cb_odd][0], this->geom);
      QPhiX::qdp_pack_spinor<>(
          this->psi_qdp_f1, this->Mchi2[cb_even][1], this->Mchi2[cb_odd][1], this->geom);

      QPhiX::qdp_pack_spinor<>(this->psi_qdp_f0,
                               this->MdagMchi2[cb_even][0],
                               this->MdagMchi2[cb_odd][0],
                               this->geom);
      QPhiX::qdp_pack_spinor<>(this->psi_qdp_f1,
                               this->MdagMchi2[cb_even][1],
                               this->MdagMchi2[cb_odd][1],
                               this->geom);

      QPhiX::qdp_pack_spinor<>(this->psi_qdp_f0,
                               this->MdagChi2[cb_even][0],
                               this->MdagChi2[cb_odd][0],
                               this->geom);
      QPhiX::qdp_pack_spinor<>(this->psi_qdp_f1,
                               this->MdagChi2[cb_even][1],
                               this->MdagChi2[cb_odd][1],
                               this->geom);

      // apply the Wilson EO operator to check hermiticity of this
      (this->wilson)->operator()(this->MwChi[cb], this->chi[cb][0], isign, cb);
      (this->wilson)->operator()(this->MwChi2[cb], this->chi2[cb][0], isign, cb);
      (this->wilson)->operator()(this->MwDagChi2[cb], this->chi2[cb][0], -isign, cb);
      (this->wilson)->operator()(this->MwDagMwChi2[cb], this->MwChi2[cb], -isign, cb);

      // apply two-flavour operator
      (this->NDTMEO)->operator()(this->Mchi[cb], this->chi[cb], isign, cb);
      (this->NDTMEO)->operator()(this->Mchi2[cb], this->chi2[cb], isign, cb);
      // apply daggered operator
      (this->NDTMEO)->operator()(this->MdagChi2[cb], this->chi2[cb], -isign, cb);
      (this->NDTMEO)->operator()(this->MdagMchi2[cb], this->Mchi2[cb], -isign, cb);

      double inner_MwChi_Chi2[2];
      double inner_Chi_MwDagChi2[2];

      double inner_MwChi_MwChi2[2];
      double inner_Chi_MwDagMwChi2[2];

      QPhiX::innerProduct<double, VECLEN_DP, VECLEN_DP, true>(inner_MwChi_Chi2,
                                                                 this->MwChi[cb],
                                                                 this->chi2[cb][0],
                                                                 this->geom,
                                                                 this->geom.getNSIMT());

      QPhiX::innerProduct<double, VECLEN_DP, VECLEN_DP, true>(inner_Chi_MwDagChi2,
                                                                 this->chi[cb][0],
                                                                 this->MwDagChi2[cb],
                                                                 this->geom,
                                                                 this->geom.getNSIMT());

      QPhiX::innerProduct<double, VECLEN_DP, VECLEN_DP, true>(inner_MwChi_MwChi2,
                                                                 this->MwChi[cb],
                                                                 this->MwChi2[cb],
                                                                 this->geom,
                                                                 this->geom.getNSIMT());

      QPhiX::innerProduct<double, VECLEN_DP, VECLEN_DP, true>(inner_Chi_MwDagMwChi2,
                                                                 this->chi[cb][0],
                                                                 this->MwDagMwChi2[cb],
                                                                 this->geom,
                                                                 this->geom.getNSIMT());

      QPhiX::masterPrintf(
          "%20s: %lf %lf\n %20s: %lf %lf \n %20s %lf %lf \n %20s %lf %lf \n\n",
          "inner_MwChi_Chi2",
          inner_MwChi_Chi2[0],
          inner_MwChi_Chi2[1],
          "inner_Chi_MwDagChi2",
          inner_Chi_MwDagChi2[0],
          inner_Chi_MwDagChi2[1],
          "inner_MwChi_MwChi2",
          inner_MwChi_MwChi2[0],
          inner_MwChi_MwChi2[1],
          "inner_Chi_MwDagMwChi2",
          inner_Chi_MwDagMwChi2[0],
          inner_Chi_MwDagMwChi2[1]);

      double inner_Mchi_Chi2[2];
      double inner_Chi_MdagChi2[2];

      double inner_Mchi_Mchi2[2];
      double inner_Chi_MdagMchi2[2];

      QPhiX::innerProduct<double, VECLEN_DP, VECLEN_DP, true, 2>(
          inner_Mchi_Chi2,
          this->Mchi[cb],
          this->chi2[cb],
          this->geom,
          this->geom.getNSIMT());

      QPhiX::innerProduct<double, VECLEN_DP, VECLEN_DP, true, 2>(
          inner_Chi_MdagChi2,
          this->chi[cb],
          this->MdagChi2[cb],
          this->geom,
          this->geom.getNSIMT());

      QPhiX::innerProduct<double, VECLEN_DP, VECLEN_DP, true, 2>(
          inner_Mchi_Mchi2,
          this->Mchi[cb],
          this->Mchi2[cb],
          this->geom,
          this->geom.getNSIMT());

      QPhiX::innerProduct<double, VECLEN_DP, VECLEN_DP, true, 2>(
          inner_Chi_MdagMchi2,
          this->chi[cb],
          this->MdagMchi2[cb],
          this->geom,
          this->geom.getNSIMT());

      QPhiX::masterPrintf(
          "%20s: %lf %lf\n %20s: %lf %lf \n %20s %lf %lf \n %20s %lf %lf \n\n",
          "inner_Mchi_Chi2",
          inner_Mchi_Chi2[0],
          inner_Mchi_Chi2[1],
          "inner_Chi_MdagChi2",
          inner_Chi_MdagChi2[0],
          inner_Chi_MdagChi2[1],
          "inner_Mchi_Mchi2",
          inner_Mchi_Mchi2[0],
          inner_Mchi_Mchi2[1],
          "inner_Chi_MdagMchi2",
          inner_Chi_MdagMchi2[0],
          inner_Chi_MdagMchi2[1]);
    }
  }
}
