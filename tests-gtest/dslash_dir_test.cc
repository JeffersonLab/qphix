/*
 * dslash_dir_test.cc
 *
 *  Created on: Oct 12, 2017
 *      Author: bjoo
 */
#include "qdpxx_test_env.h"
#include <qphix/qphix_config.h>
#include <qphix/qphix_cli_args.h>
#include <qphix/geometry.h>
#include <qphix/clover.h>
#include <qphix/wilson.h>
#include <qphix/invbicgstab.h>
#include <qphix/invmr.h>
#include <qphix/qdp_packer.h>
#include <qphix/full_spinor.h>
#include <qphix/unprec_solver_wrapper.h>

#include "../tests/veclen.h"
#include "../tests/dslashm_w.h"
#include "../tests/reunit.h"
#include "../tests/clover_fermact_params_w.h"
#include "../tests/clover_term.h"


using namespace QPhiXTesting;
using namespace QDP;
using namespace QPhiX;

using Geom = QPhiX::Geometry<double, VECLEN_DP, QPHIX_SOALEN, false>;
using ClovOp = QPhiX::EvenOddCloverOperator<double, VECLEN_DP, QPHIX_SOALEN, false>;
using WilsOp = QPhiX::EvenOddWilsonOperator<double, VECLEN_DP, QPHIX_SOALEN, false>;
using DslashOp = QPhiX::Dslash<double, VECLEN_DP,QPHIX_SOALEN,false>;
using QPhiXCBSpinor = QPhiX::FourSpinorHandle<double,VECLEN_DP,QPHIX_SOALEN,false>;
using QPhiXCBGauge = QPhiX::GaugeHandle<double,VECLEN_DP,QPHIX_SOALEN, false>;
using QPhiXCBClover = QPhiX::CloverHandle<double, VECLEN_DP,QPHIX_SOALEN, false>;
using QPhiXFullSpinor = QPhiX::FullSpinor<double,VECLEN_DP,QPHIX_SOALEN,false>;

#if 0
using BiCGStab = QPhiX::InvBiCGStab<double,VECLEN_DP,QPHIX_SOALEN,false>;
using MRBase = QPhiX::InvMRBase<double,VECLEN_DP,QPHIX_SOALEN,false>;
using MR = QPhiX::InvMR<double,VECLEN_DP,QPHIX_SOALEN,false>;
using MRSmoother = QPhiX::InvMR<double,VECLEN_DP,QPHIX_SOALEN,false>;
using EOPrecOp = QPhiX::EvenOddLinearOperator<double,VECLEN_DP,QPHIX_SOALEN,false>;
using QPhiXUnprecSolver = QPhiX::UnprecSolverWrapper<double,VECLEN_DP,
                                                    QPHIX_SOALEN,false,EOPrecOp>;
#endif

// Apply a single direction of Dslash
void DslashDirQDPXX(LatticeFermion& out, const multi1d<LatticeColorMatrix>& u, const LatticeFermion& in, int dir)
{
  switch(dir) {
  case 0: // Dir 0, Forward
    out = spinReconstructDir0Minus(u[0] * shift(spinProjectDir0Minus(in), FORWARD, 0));
    break;
  case 1: // Dir 0, Backward
    out = spinReconstructDir0Plus(shift(adj(u[0]) * spinProjectDir0Plus(in), BACKWARD, 0));
    break;
  case 2: // Dir 1, Forward
    out = spinReconstructDir1Minus(u[1] * shift(spinProjectDir1Minus(in), FORWARD, 1));
    break;
  case 3: // Dir 1, Backward
    out = spinReconstructDir1Plus(shift(adj(u[1]) * spinProjectDir1Plus(in), BACKWARD, 1));
    break;
  case 4: // Dir 2, Forward
    out = spinReconstructDir2Minus(u[2] * shift(spinProjectDir2Minus(in), FORWARD, 2));
    break;
  case 5: // Dir 2, Backward
    out = spinReconstructDir2Plus(shift(adj(u[2]) * spinProjectDir2Plus(in), BACKWARD, 2));
    break;
  case 6: // Dir 3, Forward,
    out = spinReconstructDir3Minus(u[3] * shift(spinProjectDir3Minus(in), FORWARD, 3));
    break;
  case 7: // Dir 3, Backward
    out = spinReconstructDir3Plus(shift(adj(u[3]) * spinProjectDir3Plus(in), BACKWARD, 3));
    break;
  default:
    QDPIO::cout << "Unknown direction. You oughtnt call this\n" << std::endl;
    QDP_abort(1);
  }
}

TEST(QPhiXMGTests, DslashDirTest)
{
  // This sets up the environment and the lattice
   const QDPXXTestEnv& testEnv = getQDPXXTestEnv();

   // Get QPhiX Command Line args
   const QPhiX::QPhiXCLIArgs& CLI = testEnv.getCLIArgs();

   // Set up the params
   int tBc = -1;
   double tBoundary= static_cast<double>(tBoundary);
   double mQ = 0.01;
   double xi0 = 1.3;
   double nu = 0.78;
   double anisoFacS=nu/xi0;
   double anisoFacT=1.0;
   int oddCb = 1;
   int evenCb = 1-oddCb;

   // Set geometry
   Geom geom(Layout::subgridLattSize().slice(),
       CLI.getBy(),
       CLI.getBz(),
       CLI.getNCores(),
       CLI.getSy(),
       CLI.getSz(),
       CLI.getPxy(),
       CLI.getPxyz(),
       CLI.getMinCt(),
       true);

   // Create lattice
   multi1d<LatticeColorMatrix> u(Nd);
   multi1d<LatticeColorMatrix> uAniso(Nd);
   for(int mu=0 ; mu < Nd; ++mu ) {
     LatticeColorMatrix g;
     gaussian(g);
     u[mu] = 1 + 0.1*g;
     reunit(u[mu]);
   }

   // Create lattice with ansitropy and BC applied
   // This will be used by QDP++
   // Anisotropy:
   for (int mu = 0; mu < Nd; mu++) {
     Real factor;
     if (mu == 3) {
       factor = Real(anisoFacT);
     } else {
       factor = Real(anisoFacS);
     }
     uAniso[mu] = factor * u[mu];
   }

   // Antiperiodic BCs
   int const muT = QDP::Nd - 1;
   uAniso[muT] *= QDP::where(QDP::Layout::latticeCoordinate(muT) ==
       (QDP::Layout::lattSize()[muT] - 1),
       QDP::Real(tBoundary),
       QDP::Real(1));

   // Pack U for QPhiX
   QPhiXCBGauge uCb0(geom);
   QPhiXCBGauge uCb1(geom);
   qdp_pack_gauge<>(u, uCb0.get(),uCb1.get(),geom);
   Geom::SU3MatrixBlock *qphixGauge[2] = { uCb0.get(),uCb1.get() };

   DslashOp D(&geom, tBoundary, anisoFacS, anisoFacT);

   // QDP++ Test vector:
   LatticeFermion source;
   gaussian(source);

   LatticeFermion result=zero;
   LatticeFermion q_result=zero;
   LatticeFermion diff = zero;
   Double normDiff=Double(-1);
   Double perSiteDiff=Double(-1);

   QPhiXCBSpinor qphixSource(geom);
   QPhiXCBSpinor qphixResult(geom);

   for(int dir=0; dir < 8; ++dir) {
     DslashDirQDPXX(result, uAniso, source, dir);

     for(int target_cb=0; target_cb < 2; ++target_cb) {

       // We will shift from 1-taret_cb onto target_cb
       QPhiX::qdp_pack_cb_spinor<>(source, qphixSource.get(),geom, 1-target_cb);
       D.dslashDir(qphixResult.get(),qphixSource.get(),qphixGauge[target_cb],target_cb,dir);
       QPhiX::qdp_unpack_cb_spinor(qphixResult.get(),q_result,geom,target_cb);
       diff[rb[target_cb]] = result - q_result;
       normDiff = sqrt(norm2(diff,rb[target_cb]));
       Double normCB = Double(Layout::vol()/2);
       perSiteDiff = normDiff/normCB;
       QDPIO::cout << "Dir = " << dir <<" cb = " << target_cb << " || diff || = " << normDiff
             << "|| diff || / site = " << perSiteDiff <<"\n";
       ASSERT_LT( toDouble(perSiteDiff), 5.0e-17);
     }
   }
}

TEST(UnprecOpTest, TestUnprecCloverOpsDslashDir)
{
  // This sets up the environment and the lattice
   const QDPXXTestEnv& testEnv = getQDPXXTestEnv();

   // Get QPhiX Command Line args
   const QPhiX::QPhiXCLIArgs& CLI = testEnv.getCLIArgs();

   // Set up the params
   int tBc = -1;
   double tBoundary= static_cast<double>(tBoundary);
   double mQ = 0.01;
   double cSwS = 1.23;
   double cSwT = 0.95;
   double xi0 = 1.3;
   double nu = 0.78;
   double anisoFacS=nu/xi0;
   double anisoFacT=1.0;
   int oddCb = 1;
   int evenCb = 1-oddCb;


   // Set geometry
   Geom geom(Layout::subgridLattSize().slice(),
       CLI.getBy(),
       CLI.getBz(),
       CLI.getNCores(),
       CLI.getSy(),
       CLI.getSz(),
       CLI.getPxy(),
       CLI.getPxyz(),
       CLI.getMinCt(),
       true);

   // Create lattice
   multi1d<LatticeColorMatrix> u(Nd);
   multi1d<LatticeColorMatrix> uAniso(Nd);
   for(int mu=0 ; mu < Nd; ++mu ) {
     LatticeColorMatrix g;
     gaussian(g);
     u[mu] = 1 + 0.1*g;
     reunit(u[mu]);

   }

   // Create lattice with ansitropy and BC applied
   for (int mu = 0; mu < Nd; mu++) {
     Real factor;
     if (mu == 3) {
       factor = Real(anisoFacT);
     } else {
       factor = Real(anisoFacS);
     }
     uAniso[mu] = factor * u[mu];
   }

   int const muT = QDP::Nd - 1;
   uAniso[muT] *= QDP::where(QDP::Layout::latticeCoordinate(muT) ==
       (QDP::Layout::lattSize()[muT] - 1),
       QDP::Real(tBoundary),
       QDP::Real(1));

   QPhiXCBGauge uCb0(geom);
   QPhiXCBGauge uCb1(geom);
   qdp_pack_gauge<>(u, uCb0.get(),uCb1.get(),geom);
   Geom::SU3MatrixBlock *qphixGauge[2] = { uCb0.get(),uCb1.get() };


   // Need to make a clover term and inverse to pack.
   CloverFermActParams param;
   param.Mass = Real(mQ);
   param.clovCoeffR = Real(cSwS);
   param.clovCoeffT = Real(cSwT);
   param.u0 = Real(1);
   param.anisoParam.anisoP = true;
   param.anisoParam.t_dir = 3;
   param.anisoParam.xi_0 = Real(xi0);
   param.anisoParam.nu = Real(nu);

   CloverTermT<LatticeFermion, LatticeColorMatrix>  clov;
   clov.create(u, param);
   CloverTermT<LatticeFermion, LatticeColorMatrix>  invclov;
   invclov.create(u, param, clov);
   invclov.choles(evenCb); // Even Even Inv


   QPhiXCBClover A_0(geom);
   QPhiXCBClover A_1(geom);
   QPhiXCBClover AInvEe(geom);
   qdp_pack_clover<>(clov, A_0.get(),geom,0);
   qdp_pack_clover<>(clov, A_1.get(),geom,1);
   Geom::CloverBlock *qphixClov[2] = { A_0.get(), A_1.get() };

   qdp_pack_clover<>(invclov, AInvEe.get(), geom, evenCb);

   // Create QPhiXClovOp
   ClovOp QPhiXEOClov(qphixGauge,
                      qphixClov,
                      AInvEe.get(),
                      &geom,
                      tBoundary,
                      anisoFacS,
                      anisoFacT);


   LatticeFermion source;
   gaussian(source);

   LatticeFermion result=zero;
   LatticeFermion q_result=zero;
   LatticeFermion diff = zero;
   Double normDiff=Double(-1);
   Double perSiteDiff=Double(-1);

   QPhiXFullSpinor qphixSource(geom);
   QPhiXFullSpinor qphixResult(geom);

   for(int dir=0; dir < 8; ++dir) {
     DslashDirQDPXX(result, uAniso, source, dir);

     for(int cb=0; cb < 2; ++cb ) {
       QPhiX::qdp_pack_cb_spinor<>(source, qphixSource.getCBData(cb),geom, cb);
     }
     QPhiXEOClov.DslashDir(qphixResult,qphixSource,dir);

     for(int cb=0; cb < 2; ++cb ) {
       QPhiX::qdp_unpack_cb_spinor(qphixResult.getCBData(cb),q_result,geom,cb);
     }

     for(int cb=0; cb < 2; ++cb) {
         diff[rb[cb]] = result - q_result;
         normDiff = sqrt(norm2(diff,rb[cb]));
         Double normCB = Double(Layout::vol()/2);
         perSiteDiff = normDiff/normCB;
         QDPIO::cout << "Dir = " << dir <<" cb = " << cb << " || diff || = " << normDiff
               << " || diff || / site = " << perSiteDiff <<"\n";
         ASSERT_LT( toDouble(perSiteDiff), 5.0e-17);
     }
     diff = result - q_result;
     normDiff = sqrt(norm2(diff));
     perSiteDiff = normDiff/Double(Layout::vol());
     QDPIO::cout << "Dir = " << dir << " || diff || = " << normDiff
         << " || diff || / site = " << perSiteDiff << "\n";
     ASSERT_LT( toDouble(perSiteDiff), 5.0e-17);

   }

}



int main(int argc, char *argv[])
{
  return QDPXXTestMain(&argc, argv);
}
