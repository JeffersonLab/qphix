#include <qphix/qphix_config_internal.h>
#include "prec.h"
#include "tparam_selector.h"

using namespace QPhiX;

#include "./soalen.h"
//#include "../tests/clover_fermact_params_w.h"
//#include "../tests/clover_term.h"

#include <qphix/ndtm_reuse_operator.h>
#include <qphix/geometry.h>
#include <qphix/qdp_packer.h>
#include <qphix/qphix_cli_args.h>
#include <qphix/wilson.h>
#include <qphix/blas_new_c.h>

#include <qdp.h>
#include <gtest/gtest.h>

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
        xi_0_f(1.0), nu_f(1.0), aniso_fac_s(nu_f / xi_0_f),
        aniso_fac_t(1.0),
        packed_gauge_cb0(geom), packed_gauge_cb1(geom),
        u_packed{ packed_gauge_cb0.get(), packed_gauge_cb1.get() },
        NDTMEO(nullptr)
  {
    NDTMEO = new EvenOddNDTMWilsonReuseOperator<typename Geom::FT, Geom::veclen, Geom::soalen, Geom::compress12>
                   (-0.37, 0.1, 0.01, u_packed, &geom, 1, aniso_fac_s, aniso_fac_t);

    constexpr int cb_even = 0;
    constexpr int cb_odd = 1;
    // Make a random gauge field
    // Start up the gauge field somehow
    // We choose: u = reunit(  1 + factor*gaussian(u) );
    // Adjust factor to get reasonable inversion times in invert test.
    // Bug gauge field is always unitary
    QDP::multi1d<QDPGauge> u(4);
    {
      QDPGauge g;
      QDPGauge uf;
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

    // Make a random source
    gaussian(chi_qdp_f0);
    gaussian(chi_qdp_f1);

    // Pack the gauge field
    QPhiX::qdp_pack_gauge<>(u, packed_gauge_cb0.get(), packed_gauge_cb1.get(), geom);

    // flavour index on the inside to allow passing to two-flavour functions
    for( int cb : {cb_even, cb_odd} ){
      for( int f : {0, 1} ){
        qhandles.push_back( makeFourSpinorHandle(geom) );
        chi[cb][f] = qhandles.back().get();
        qhandles.push_back( makeFourSpinorHandle(geom) );
        Mchi[cb][f] = qhandles.back().get();
        qhandles.push_back( makeFourSpinorHandle(geom) );
        MdagMchi[cb][f] = qhandles.back().get();

        QPhiX::masterPrintf("%p %p %p\n", (void*)chi[cb][f], (void*)Mchi[cb][f], (void*)MdagMchi[cb][f]);
      }
    }
    QPhiX::masterPrintf("%p %p\n", (void*)u_packed[0], (void*)u_packed[1]);

    // pack the gaussian sources
    QPhiX::qdp_pack_spinor<>(chi_qdp_f0, chi[cb_even][0], chi[cb_odd][0], geom);
    QPhiX::qdp_pack_spinor<>(chi_qdp_f1, chi[cb_even][1], chi[cb_odd][1], geom);

    const int cbsize_in_blocks = rb[0].numSiteTable() / Geom::soalen;
  }

  ~NDTMReuseOperatorHermTest(){
    delete NDTMEO;
  }

  Geom geom;

  QPhiX::EvenOddNDTMWilsonReuseOperator<typename Geom::FT, Geom::veclen, Geom::soalen, Geom::compress12> *NDTMEO;

  QDPSpinor chi_qdp_f0;
  QDPSpinor chi_qdp_f1;
  QDPSpinor psi_qdp_f0;
  QDPSpinor psi_qdp_f1;

  QPhiX::GaugeHandle<FT, veclen, soalen, compress12> packed_gauge_cb0;
  QPhiX::GaugeHandle<FT, veclen, soalen, compress12> packed_gauge_cb1;

  std::vector< QPhiX::FourSpinorHandle< FT, veclen, soalen, compress12> > qhandles;

  Gauge *u_packed[2];

  Spinor *Mchi[2][2];
  Spinor *chi[2][2];
  Spinor *MdagMchi[2][2];

};

typedef ::testing::Types<QPhiX::Geometry<double, VECLEN_DP, QPHIX_SOALEN, true>>
    MyTypes;

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
          this->psi_qdp_f0, this->Mchi[cb_even][0], this->Mchi[cb_odd][0], this->geom);
      QPhiX::qdp_pack_spinor<>(
          this->psi_qdp_f1, this->Mchi[cb_even][1], this->Mchi[cb_odd][1], this->geom);
      QPhiX::qdp_pack_spinor<>(
          this->psi_qdp_f0, this->MdagMchi[cb_even][0], this->MdagMchi[cb_odd][0], this->geom);
      QPhiX::qdp_pack_spinor<>(
          this->psi_qdp_f1, this->MdagMchi[cb_even][1], this->MdagMchi[cb_odd][1], this->geom);

      // apply two-flavour operator
      (this->NDTMEO)->operator()(this->Mchi[cb],
                   this->chi[cb],
                   isign,
                   cb);
      // apply daggered operator
      (this->NDTMEO)->operator()(this->MdagMchi[cb],
                   this->chi[cb],
                   -isign,
                   cb);
      
      double inner_Chi_MdagMchi[2];
      double inner_Mchi_Mchi[2];
      QPhiX::innerProduct<double, VECLEN_DP, QPHIX_SOALEN, true, 2>(
          inner_Mchi_Mchi, this->Mchi[cb], this->Mchi[cb], this->geom, this->geom.getNSIMT() );
      QPhiX::innerProduct<double, VECLEN_DP, QPHIX_SOALEN, true, 2>(
          inner_Chi_MdagMchi, this->chi[cb], this->MdagMchi[cb], this->geom, this->geom.getNSIMT() );

      QPhiX::masterPrintf("%20s: %lf %lf\n %20s: %lf %lf \n\n",
                          "inner_Mchi_Mchi", inner_Mchi_Mchi[0], inner_Mchi_Mchi[1],
                          "inner_Chi_MdagMchi", inner_Chi_MdagMchi[0], inner_Chi_MdagMchi[1]);  
    }
  }
}

