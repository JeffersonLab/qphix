/*
 * clover_mult_test.cc
 *
 *  Created on: May 16, 2018
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

#include "../tests/compare_qdp_spinors_gtest.h"

using namespace QPhiXTesting;
using namespace QDP;
using namespace QPhiX;

using Geom = QPhiX::Geometry<double, VECLEN_DP, QPHIX_SOALEN, false>;
using ClovOp = QPhiX::EvenOddCloverOperator<double, VECLEN_DP, QPHIX_SOALEN, false>;
using WilsOp = QPhiX::EvenOddWilsonOperator<double, VECLEN_DP, QPHIX_SOALEN, false>;
using DslashOp = QPhiX::Dslash<double, VECLEN_DP,QPHIX_SOALEN,false>;
using ClovDslashOp = QPhiX::ClovDslash<double,VECLEN_DP,QPHIX_SOALEN,false>;

using QPhiXCBSpinor = QPhiX::FourSpinorHandle<double,VECLEN_DP,QPHIX_SOALEN,false>;
using QPhiXCBGauge = QPhiX::GaugeHandle<double,VECLEN_DP,QPHIX_SOALEN, false>;
using QPhiXCBClover = QPhiX::CloverHandle<double, VECLEN_DP,QPHIX_SOALEN, false>;
using QPhiXFullSpinor = QPhiX::FullSpinor<double,VECLEN_DP,QPHIX_SOALEN,false>;

TEST(CloverMultTest, TestGeneratedCloverTerm)
{
  // This sets up the environment and the lattice
   const QDPXXTestEnv& testEnv = getQDPXXTestEnv();

   // Get QPhiX Command Line args
   const QPhiX::QPhiXCLIArgs& CLI = testEnv.getCLIArgs();

   // Set up the params
   int tBc = -1;
   double tBoundary= static_cast<double>(tBc);
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

   QPhiXCBGauge uCb0(geom);
   QPhiXCBGauge uCb1(geom);
   qdp_pack_gauge<>(u, uCb0.get(),uCb1.get(),geom);
   Geom::SU3MatrixBlock *qphixGauge[2] = { uCb0.get(),uCb1.get() };


   int const muT = QDP::Nd - 1;

   // Apply boundary for clover term.

   u[muT] *= QDP::where(QDP::Layout::latticeCoordinate(muT) ==
       (QDP::Layout::lattSize()[muT] - 1),
       QDP::Real(tBoundary),
       QDP::Real(1));



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

   ClovDslashOp D_clov(&geom,
		   	   	   	   	   	   	   	 tBoundary,
		   	   	   	   	   	   	    anisoFacS,
									anisoFacT);



   LatticeFermion source;
   gaussian(source);

   LatticeFermion q_result=zero;
   LatticeFermion result=zero;
   LatticeFermion diff;
   gaussian(diff);


   Double normDiff=Double(-1);
   Double perSiteDiff=Double(-1);

   QPhiXFullSpinor qphixSource(geom);
   QPhiXFullSpinor qphixResult(geom);

   for(int cb=0; cb < 2; ++cb) {

	   	   QDPIO::cout << "Applying QPhiX Operator\n" ;

	        QPhiX::qdp_pack_cb_spinor<>(source, qphixSource.getCBData(cb),geom, cb);

	        D_clov.clovMult(qphixResult.getCBData(cb), qphixSource.getCBData(cb),
	        			qphixClov[cb]);

	        QPhiX::qdp_unpack_cb_spinor(qphixResult.getCBData(cb),q_result,geom,cb);

	        QDPIO::cout << "Applying QDP++ operator\n";
	        clov.apply(result, source, 1, cb);

	        QDPIO::cout << "Norm QPhix Result [" << cb << "]=" << sqrt(norm2(q_result,rb[cb])) << std::endl;
	        QDPIO::cout << "Norm QDP++ result [ " << cb << "]=" << sqrt(norm2(result,rb[cb])) << std::endl;

	        diff[rb[cb]] = result - q_result;
	         normDiff = sqrt(norm2(diff,rb[cb]));
	         Double normCB = Double(Layout::vol()/2);
	         perSiteDiff = normDiff/normCB;
	         QDPIO::cout << " cb = " << cb << " || diff || = " << normDiff
	               << " || diff || / site = " << perSiteDiff <<"\n";

	         expect_near<>(result, q_result, 1.0e-8, geom, cb, "Foo");
//	ASSERT_LT( toDouble(perSiteDiff), 5.0e-17);

   }
}

#if 1
TEST(CloverMultTest, TestGeneratedCloverTermEOOp)
{
  // This sets up the environment and the lattice
   const QDPXXTestEnv& testEnv = getQDPXXTestEnv();

   // Get QPhiX Command Line args
   const QPhiX::QPhiXCLIArgs& CLI = testEnv.getCLIArgs();

   // Set up the params
   int tBc = -1;
   double tBoundary= static_cast<double>(tBc);
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

     for(int cb=0; cb < 2; ++cb ) {
       QPhiX::qdp_pack_cb_spinor<>(source, qphixSource.getCBData(cb),geom, cb);

       QPhiXEOClov.M_diag(qphixResult,qphixSource,1,cb);

       QPhiX::qdp_unpack_cb_spinor(qphixResult.getCBData(cb),q_result,geom,cb);

       QDPIO::cout << "Applying QDP++ operator\n";
       clov.apply(result, source, 1, cb);

       diff[rb[cb]] = result - q_result;
       normDiff = sqrt(norm2(diff,rb[cb]));
       Double normCB = Double(Layout::vol()/2);
       perSiteDiff = normDiff/normCB;
       QDPIO::cout << " cb = " << cb << " || diff || = " << normDiff
    		   << " || diff || / site = " << perSiteDiff <<"\n";

         expect_near(result,q_result, 1.0e-8, geom, cb, "M_diag_check");
     }
     diff = result - q_result;
     normDiff = sqrt(norm2(diff));
     perSiteDiff = normDiff/Double(Layout::vol());
     QDPIO::cout << "|| diff || = " << normDiff
         << " || diff || / site = " << perSiteDiff << "\n";




     QPhiXFullSpinor qphixInvApply(geom);
     {
    	 int cb=0;
    	 QPhiXEOClov.M_diag_inv(qphixInvApply.getCBData(cb),qphixResult.getCBData(cb),1);
    	 QPhiX::qdp_unpack_cb_spinor(qphixInvApply.getCBData(cb),q_result,geom,cb);
    	 expect_near<>(source, q_result, 1.0e-8, geom, cb, "M_diag_Inv_check");

     }


}
#endif


TEST(CloverMultTest, TestGeneratedCloverTermSpeed)
{
  // This sets up the environment and the lattice
   const QDPXXTestEnv& testEnv = getQDPXXTestEnv();

   // Get QPhiX Command Line args
   const QPhiX::QPhiXCLIArgs& CLI = testEnv.getCLIArgs();

   // Set up the params
   int tBc = -1;
   double tBoundary= static_cast<double>(tBc);
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

     for(int cb=0; cb < 2; ++cb ) {
       QPhiX::qdp_pack_cb_spinor<>(source, qphixSource.getCBData(cb),geom, cb);
       int timing_iters=1000;
       StopWatch swatch;
       swatch.reset();
       swatch.start();
       for(int i=0; i < timing_iters; ++i ) {
    	   QPhiXEOClov.M_diag(qphixResult,qphixSource,1,cb);
       }
       swatch.stop();
       double optimized_time = swatch.getTimeInSeconds();

       swatch.reset();
       swatch.start();
       for(int i=0; i < timing_iters; ++i ) {
     	   clover_product(qphixResult.getCBData(cb),qphixSource.getCBData(cb),qphixClov[cb],geom);
        }
        swatch.stop();
        double clover_product_time = swatch.getTimeInSeconds();

        QDPIO::cout << "cb = " << cb << " : Original clover_product time=" << clover_product_time << " sec for " << timing_iters << " iters."<< std::endl;
        QDPIO::cout << "cb = " << cb << " : Codegen based cloverMult time=" << optimized_time << " sec for " << timing_iters << " iters. " << std::endl;
        QDPIO::cout << "cb = " << cb << " : Speedup from codegen = " << clover_product_time/optimized_time << " x\n";



     }
     diff = result - q_result;
     normDiff = sqrt(norm2(diff));
     perSiteDiff = normDiff/Double(Layout::vol());
     QDPIO::cout << "|| diff || = " << normDiff
         << " || diff || / site = " << perSiteDiff << "\n";

}

int main(int argc, char *argv[])
{
  return QDPXXTestMain(&argc, argv);
}
