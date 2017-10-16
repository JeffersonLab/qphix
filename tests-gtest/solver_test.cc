/*
 * bicgstab_solver_test.cc
 *
 *  Created on: Jul 10, 2017
 *      Author: traveluser
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
using BiCGStab = QPhiX::InvBiCGStab<double,VECLEN_DP,QPHIX_SOALEN,false>;
using MRBase = QPhiX::InvMRBase<double,VECLEN_DP,QPHIX_SOALEN,false>;
using MR = QPhiX::InvMR<double,VECLEN_DP,QPHIX_SOALEN,false>;
using MRSmoother = QPhiX::InvMR<double,VECLEN_DP,QPHIX_SOALEN,false>;
using QPhiXCBSpinor = QPhiX::FourSpinorHandle<double,VECLEN_DP,QPHIX_SOALEN,false>;
using QPhiXCBGauge = QPhiX::GaugeHandle<double,VECLEN_DP,QPHIX_SOALEN, false>;
using QPhiXCBClover = QPhiX::CloverHandle<double, VECLEN_DP,QPHIX_SOALEN, false>;
using EOPrecOp = QPhiX::EvenOddLinearOperator<double,VECLEN_DP,QPHIX_SOALEN,false>;

using QPhiXFullSpinor = QPhiX::FullSpinor<double,VECLEN_DP,QPHIX_SOALEN,false>;

// Weird thing: If I add the LinOp type to be specifically ClovOp or WilsOp
// the virtual functions from EOPrecOp are not found, and I'd need to implemented
// in the actual ClovOp or WilsOp
// So basically casting those back to the base class is the way forward.
//
using QPhiXUnprecSolver = QPhiX::UnprecSolverWrapper<double,VECLEN_DP,
                                                    QPHIX_SOALEN,false,EOPrecOp>;

TEST(BiCGStabSolverTest, TestBiCGStabSolver)
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
	QPhiX::ResiduumType residType = QPhiX::RELATIVE;

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
	invclov.choles(0); // Even Even Inv

	QPhiXCBClover aOo(geom);
	QPhiXCBClover aInvEe(geom);
	qdp_pack_clover<>(clov, aOo.get(),geom,oddCb);
	qdp_pack_clover<>(invclov, aInvEe.get(), geom, evenCb);



	LatticeFermion source;
	gaussian(source);

	LatticeFermion solution=zero;

	QPhiXCBSpinor qphixSource(geom);
	QPhiXCBSpinor qphixSolution(geom);

	// Now pack fields.

	QPhiX::qdp_pack_cb_spinor<>(source, qphixSource.get(),geom, oddCb);
	QPhiX::qdp_pack_cb_spinor<>(solution, qphixSolution.get(),geom,oddCb);

	// Create QPhiXClovOp
	ClovOp QPhiXEOClov(qphixGauge,
			aOo.get(),
			aInvEe.get(),
			&geom,
			tBoundary,
			anisoFacS,
			anisoFacT);


	// make a BiCGStab Solver
	BiCGStab solver(QPhiXEOClov, 5000);



	QPhiXCBSpinor::ValueType* soln[1] = { qphixSolution.get() };
	QPhiXCBSpinor::ValueType* rhs[1] = { qphixSource.get() };

	int nIters;
	double rsdSqFinal;
	unsigned long siteFlops;
	unsigned long mvApps;

	solver(soln,rhs, 1.0e-7, nIters, rsdSqFinal,siteFlops, mvApps,1,true, oddCb,residType);

	// Solution[odd] is the same between the transformed and untransformed system
	QPhiX::qdp_unpack_cb_spinor<>(qphixSolution.get(),solution,geom,oddCb);

    LatticeFermion tmp1 = zero;
    LatticeFermion tmp2 = zero;
    LatticeFermion handMSoln=zero;
    Real betaFactor=Real(0.25); // 2x0.5 from the Dslash-es

    dslash(tmp1, uAniso, solution, 1, evenCb);
    invclov.apply(tmp2, tmp1, 1, evenCb);
    QPhiX::dslash(tmp1, uAniso, tmp2, 1, oddCb);

    clov.apply(handMSoln, solution, 1, oddCb);
    handMSoln[rb[oddCb]] -= betaFactor * tmp1;


    LatticeFermionD diff;
    diff[rb[oddCb]] = source - handMSoln;
    Double norm_diff = sqrt(norm2(diff,rb[oddCb]));

    if( residType == RELATIVE ) {
    	Double norm_src = sqrt(norm2(source,rb[oddCb]));
    	QDPIO::cout << "|| r || = " << norm_diff << " || r || / || b || = " << norm_diff/norm_src << "\n";
    	double rel_resid = toDouble(norm_diff/norm_src);
    	ASSERT_LT( rel_resid, 1.0e-7);

    }
    else {
    	QDPIO::cout << "|| r || = " << norm_diff << "\n";
    }
}


TEST(InvMRTest, TestInvMR)
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
	QPhiX::ResiduumType residType = QPhiX::RELATIVE;

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
	invclov.choles(0); // Even Even Inv

	QPhiXCBClover aOo(geom);
	QPhiXCBClover aInvEe(geom);
	qdp_pack_clover<>(clov, aOo.get(),geom,oddCb);
	qdp_pack_clover<>(invclov, aInvEe.get(), geom, evenCb);



	LatticeFermion source;
	gaussian(source);

	LatticeFermion solution=zero;

	QPhiXCBSpinor qphixSource(geom);
	QPhiXCBSpinor qphixSolution(geom);

	// Now pack fields.

	QPhiX::qdp_pack_cb_spinor<>(source, qphixSource.get(),geom, oddCb);
	QPhiX::qdp_pack_cb_spinor<>(solution, qphixSolution.get(),geom,oddCb);

	// Create QPhiXClovOp
	ClovOp QPhiXEOClov(qphixGauge,
			aOo.get(),
			aInvEe.get(),
			&geom,
			tBoundary,
			anisoFacS,
			anisoFacT);


	// make a BiCGStab Solver
	double omega=1.1;
	MR mr_solver(QPhiXEOClov, 5000, omega);


	QPhiXCBSpinor::ValueType* soln[1] = { qphixSolution.get() };
	QPhiXCBSpinor::ValueType* rhs[1] = { qphixSource.get() };

	int nIters;
	double rsdSqFinal;
	unsigned long siteFlops;
	unsigned long mvApps;

	mr_solver(soln, rhs, 1.0e-7, nIters, rsdSqFinal, siteFlops, mvApps, 1, true, oddCb,residType);

	// Solution[odd] is the same between the transformed and untransformed system
	QPhiX::qdp_unpack_cb_spinor<>(qphixSolution.get(),solution,geom,oddCb);

    LatticeFermion tmp1 = zero;
    LatticeFermion tmp2 = zero;
    LatticeFermion handMSoln=zero;
    Real betaFactor=Real(0.25); // 2x0.5 from the Dslash-es

    dslash(tmp1, uAniso, solution, 1, evenCb);
    invclov.apply(tmp2, tmp1, 1, evenCb);
    QPhiX::dslash(tmp1, uAniso, tmp2, 1, oddCb);

    clov.apply(handMSoln, solution, 1, oddCb);
    handMSoln[rb[oddCb]] -= betaFactor * tmp1;


    LatticeFermionD diff;
    diff[rb[oddCb]] = source - handMSoln;
    Double norm_diff = sqrt(norm2(diff,rb[oddCb]));

    if( residType == RELATIVE ) {
    	Double norm_src = sqrt(norm2(source,rb[oddCb]));
    	QDPIO::cout << "|| r || = " << norm_diff << " || r || / || b || = " << norm_diff/norm_src << "\n";
    	double rel_resid = toDouble(norm_diff/norm_src);
    	ASSERT_LT( rel_resid, 1.0e-7);

    }
    else {
    	QDPIO::cout << "|| r || = " << norm_diff << "\n";
    }
}

TEST(InvMRSmootherTest, TestInvMRSmoother)
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
	QPhiX::ResiduumType residType = QPhiX::RELATIVE;

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
	invclov.choles(0); // Even Even Inv

	QPhiXCBClover aOo(geom);
	QPhiXCBClover aInvEe(geom);
	qdp_pack_clover<>(clov, aOo.get(),geom,oddCb);
	qdp_pack_clover<>(invclov, aInvEe.get(), geom, evenCb);



	LatticeFermion source;
	gaussian(source);

	LatticeFermion solution=zero;

	QPhiXCBSpinor qphixSource(geom);
	QPhiXCBSpinor qphixSolution(geom);

	// Now pack fields.

	QPhiX::qdp_pack_cb_spinor<>(source, qphixSource.get(),geom, oddCb);
	QPhiX::qdp_pack_cb_spinor<>(solution, qphixSolution.get(),geom,oddCb);

	// Create QPhiXClovOp
	ClovOp QPhiXEOClov(qphixGauge,
			aOo.get(),
			aInvEe.get(),
			&geom,
			tBoundary,
			anisoFacS,
			anisoFacT);


	// make a BiCGStab Solver
	double omega=1.1;
	MRSmoother  mr_solver(QPhiXEOClov, 5, omega);


	QPhiXCBSpinor::ValueType* soln[1] = { qphixSolution.get() };
	QPhiXCBSpinor::ValueType* rhs[1] = { qphixSource.get() };

	int niters;
	double rsdSqFinal;
	unsigned long siteFlops;
	unsigned long mvApps;

	mr_solver(soln, rhs, 1.0e-7, niters, rsdSqFinal, siteFlops, mvApps, 1, true, oddCb,residType);

	// Solution[odd] is the same between the transformed and untransformed system
	QPhiX::qdp_unpack_cb_spinor<>(qphixSolution.get(),solution,geom,oddCb);

    LatticeFermion tmp1 = zero;
    LatticeFermion tmp2 = zero;
    LatticeFermion handMSoln=zero;
    Real betaFactor=Real(0.25); // 2x0.5 from the Dslash-es

    dslash(tmp1, uAniso, solution, 1, evenCb);
    invclov.apply(tmp2, tmp1, 1, evenCb);
    QPhiX::dslash(tmp1, uAniso, tmp2, 1, oddCb);

    clov.apply(handMSoln, solution, 1, oddCb);
    handMSoln[rb[oddCb]] -= betaFactor * tmp1;


    LatticeFermionD diff;
    diff[rb[oddCb]] = source - handMSoln;
    Double norm_diff = sqrt(norm2(diff,rb[oddCb]));

    if( residType == RELATIVE ) {
    	Double norm_src = sqrt(norm2(source,rb[oddCb]));
    	QDPIO::cout << "|| r || = " << norm_diff << " || r || / || b || = " << norm_diff/norm_src << "\n";
    	double rel_resid = toDouble(norm_diff/norm_src);
    }
    else {
    	QDPIO::cout << "|| r || = " << norm_diff << "\n";
    }
    ASSERT_EQ(niters,5);

}

TEST(UnprecBiCGStabSolverTest, TestUnprecBiCGStabSolver)
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
  QPhiX::ResiduumType residType = QPhiX::RELATIVE;

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

  QPhiXCBClover aOo(geom);
  QPhiXCBClover aInvEe(geom);
  qdp_pack_clover<>(clov, aOo.get(),geom,oddCb);
  qdp_pack_clover<>(invclov, aInvEe.get(), geom, evenCb);



  LatticeFermion source;
  gaussian(source);

  LatticeFermion solution=zero;

  QPhiXFullSpinor qphixSource(geom);
  QPhiXFullSpinor qphixSolution(geom);

  // Now pack fields.
  QPhiX::qdp_pack_cb_spinor<>(source, qphixSource.getCBData(evenCb),geom, evenCb);
  QPhiX::qdp_pack_cb_spinor<>(source, qphixSource.getCBData(oddCb),geom, oddCb);

  QPhiX::qdp_pack_cb_spinor<>(solution, qphixSolution.getCBData(evenCb),geom,evenCb);
  QPhiX::qdp_pack_cb_spinor<>(solution, qphixSolution.getCBData(oddCb),geom,oddCb);

  // Create QPhiXClovOp
  ClovOp QPhiXEOClov(qphixGauge,
      aOo.get(),
      aInvEe.get(),
      &geom,
      tBoundary,
      anisoFacS,
      anisoFacT);


  // make a BiCGStab Solver
  BiCGStab solver(QPhiXEOClov, 5000);
  QPhiXUnprecSolver unprec_wrapper(solver,QPhiXEOClov);




  //QPhiXFullSpinor* soln[1] = { &qphixSolution };
  //QPhiXFullSpinor* rhs[1] = { &qphixSource };
  QPhiXFullSpinor* soln =  &qphixSolution;
  QPhiXFullSpinor* rhs = &qphixSource;
  int nIters;
  double rsdSqFinal;
  unsigned long siteFlops;
  unsigned long mvApps;

  unprec_wrapper(soln,rhs, 1.0e-7, nIters, rsdSqFinal,siteFlops, mvApps,1,true, oddCb,residType);

  // Solution[odd] is the same between the transformed and untransformed system
  QPhiX::qdp_unpack_cb_spinor<>(qphixSolution.getCBData(evenCb),solution,geom,evenCb);
  QPhiX::qdp_unpack_cb_spinor<>(qphixSolution.getCBData(oddCb),solution,geom,oddCb);

  QDPIO::cout << "Checking Odd Solution" << std::endl;
  LatticeFermion tmp1 = zero;
  LatticeFermion tmp2 = zero;
  LatticeFermion handMSoln=zero;
  Real betaFactor=Real(0.25); // 2x0.5 from the Dslash-es

  dslash(tmp1, uAniso, solution, 1, evenCb);
  invclov.apply(tmp2, tmp1, 1, evenCb);
  dslash(tmp1, uAniso, tmp2, 1, oddCb);

  clov.apply(handMSoln, solution, 1, oddCb);
  handMSoln[rb[oddCb]] -= betaFactor * tmp1;

  // construct transformed source:
  LatticeFermion transf_source = source;

  {
    LatticeFermion t1,t2;
    invclov.apply(t1, source, 1, evenCb);
    dslash( t2, uAniso, t1,1,oddCb);
    transf_source[ rb[oddCb] ] -= Real(-0.5)*t2;
   }


  LatticeFermionD diff;
  diff[rb[oddCb]] = transf_source - handMSoln;
  Double norm_diff = sqrt(norm2(diff,rb[oddCb]));

  Double norm_src=sqrt(norm2(transf_source,rb[oddCb]));
  if( residType == RELATIVE ) {
    QDPIO::cout << "|| b || = " << norm_src << " || r || = " << norm_diff << " || r || / || b || = " << norm_diff/norm_src << "\n";
    double rel_resid = toDouble(norm_diff/norm_src);
    ASSERT_LT( rel_resid, 1.0e-7);

  }
  else {
    QDPIO::cout << "|| b || = " << norm_src << " || r || = " << norm_diff << "\n";
    ASSERT_LT( toDouble(norm_diff), 1.0e-7);
  }

  QDPIO::cout << "Checking Full Solution" << std::endl;
  for(int target_cb=0; target_cb < 2; ++target_cb) {
    clov.apply(handMSoln,solution,1,target_cb);
    dslash(tmp1,uAniso,solution, 1, target_cb);
    handMSoln[rb[target_cb]] -= Real(0.5)*tmp1;
  }

  for(int target_cb=0; target_cb < 2; ++target_cb) {


    diff[rb[target_cb]] = source - handMSoln;
    norm_diff = sqrt(norm2(diff,rb[target_cb]));
    norm_src = sqrt(norm2(source,rb[target_cb]));
    if ( residType == RELATIVE ) {
      QDPIO::cout << "CB="<<target_cb<<": || b || = "<<norm_src
                  << "|| r || = " << norm_diff<< " || r || / || b || = " << norm_diff/norm_src
                  << std::endl;
      double rel_resid = toDouble( norm_diff/norm_src );
      ASSERT_LT( rel_resid, 1.0e-7);
    }
    else {
      QDPIO::cout << "CB="<<target_cb<<": || norm_src || = "<<norm_src
                        << " || r || = " << norm_diff << std::endl;
      ASSERT_LT( toDouble(norm_diff), 1.0e-7);
    }

  }

  // On whole lattice
  diff = source - handMSoln;
  norm_diff = sqrt(norm2(diff));
  norm_src = sqrt(norm2(source));

  if( residType == RELATIVE ) {

    double rel_resid = toDouble(norm_diff/norm_src);
    QDPIO::cout << "Full: || b || = " << norm_src
                << " || r || = " << norm_diff
                << " || r || / || b || = " << rel_resid << "\n";
  }
  else {
    QDPIO::cout << "Full: || b || = " << norm_src
              <<  " || r || = " << norm_diff <<"\n";
  }
}


TEST(UnprecMRSolverTest, TestUnprecBiCGStabSolverWilson)
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
  QPhiX::ResiduumType residType = QPhiX::RELATIVE;

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

  LatticeFermion source;
  gaussian(source);

  LatticeFermion solution=zero;

  QPhiXFullSpinor qphixSource(geom);
  QPhiXFullSpinor qphixSolution(geom);

  // Now pack fields.
  QPhiX::qdp_pack_cb_spinor<>(source, qphixSource.getCBData(evenCb),geom, evenCb);
  QPhiX::qdp_pack_cb_spinor<>(source, qphixSource.getCBData(oddCb),geom, oddCb);

  QPhiX::qdp_pack_cb_spinor<>(solution, qphixSolution.getCBData(evenCb),geom,evenCb);
  QPhiX::qdp_pack_cb_spinor<>(solution, qphixSolution.getCBData(oddCb),geom,oddCb);

  // Create QPhiXClovOp
  WilsOp QPhiXEOWils(mQ,qphixGauge,
            &geom,
            tBoundary,
            anisoFacS,
            anisoFacT);


  // make a BiCGStab Solver
  MR solver(QPhiXEOWils, 5000,1.1);
  QPhiXUnprecSolver unprec_wrapper(solver,QPhiXEOWils);

  QPhiXFullSpinor* soln =  &qphixSolution;
  QPhiXFullSpinor* rhs = &qphixSource;
  int nIters;
  double rsdSqFinal;
  unsigned long siteFlops;
  unsigned long mvApps;

  unprec_wrapper(soln,rhs, 1.0e-7, nIters, rsdSqFinal,siteFlops, mvApps,1,true, oddCb,residType);

  // Solution[odd] is the same between the transformed and untransformed system
  QPhiX::qdp_unpack_cb_spinor<>(qphixSolution.getCBData(evenCb),solution,geom,evenCb);
  QPhiX::qdp_unpack_cb_spinor<>(qphixSolution.getCBData(oddCb),solution,geom,oddCb);

  QDPIO::cout << "Checking Odd Solution" << std::endl;
  LatticeFermion tmp1 = zero;
  LatticeFermion tmp2 = zero;
  LatticeFermion handMSoln=zero;
  Real alphaFactor=Real(Nd)+mQ;
  Real betaFactor=Real(0.25)/alphaFactor; // 2x0.5 from the Dslash-es

  dslash(tmp1, uAniso, solution, 1, evenCb);
  dslash(tmp2, uAniso, tmp1, 1, oddCb);
  handMSoln[rb[oddCb]] = alphaFactor*solution-betaFactor*tmp2;

  // construct transformed source:
  LatticeFermion transf_source = source;

  {
    LatticeFermion t1;
    dslash( t1, uAniso, source,1,oddCb);
    transf_source[ rb[oddCb] ] -= (Real(-0.5)/alphaFactor)*t1;
   }


  LatticeFermionD diff;
  diff[rb[oddCb]] = transf_source - handMSoln;
  Double norm_diff = sqrt(norm2(diff,rb[oddCb]));

  Double norm_src=sqrt(norm2(transf_source,rb[oddCb]));
  if( residType == RELATIVE ) {
    QDPIO::cout << "|| b || = " << norm_src << " || r || = " << norm_diff << " || r || / || b || = " << norm_diff/norm_src << "\n";
    double rel_resid = toDouble(norm_diff/norm_src);
    ASSERT_LT( rel_resid, 1.0e-7);

  }
  else {
    QDPIO::cout << "|| b || = " << norm_src << " || r || = " << norm_diff << "\n";
    ASSERT_LT( toDouble(norm_diff), 1.0e-7);
  }

  QDPIO::cout << "Checking Full Solution" << std::endl;
  for(int target_cb=0; target_cb < 2; ++target_cb) {
    dslash(tmp1,uAniso,solution, 1, target_cb);
    handMSoln[rb[target_cb]] = alphaFactor*solution - Real(0.5)*tmp1;
  }

  for(int target_cb=0; target_cb < 2; ++target_cb) {


    diff[rb[target_cb]] = source - handMSoln;
    norm_diff = sqrt(norm2(diff,rb[target_cb]));
    norm_src = sqrt(norm2(source,rb[target_cb]));
    if ( residType == RELATIVE ) {
      QDPIO::cout << "CB="<<target_cb<<": || b || = "<<norm_src
                  << "|| r || = " << norm_diff<< " || r || / || b || = " << norm_diff/norm_src
                  << std::endl;
      double rel_resid = toDouble( norm_diff/norm_src );
      ASSERT_LT( rel_resid, 5.0e-7);
    }
    else {
      QDPIO::cout << "CB="<<target_cb<<": || norm_src || = "<<norm_src
                        << "|| r || = " << norm_diff << std::endl;
      ASSERT_LT( toDouble(norm_diff), 5.0e-7);
    }

  }

  // On whole lattice
  diff = source - handMSoln;
  norm_diff = sqrt(norm2(diff));
  norm_src = sqrt(norm2(source));

  if( residType == RELATIVE ) {

    double rel_resid = toDouble(norm_diff/norm_src);
    QDPIO::cout << "Full: || b || = " << norm_src
                << " || r || = " << norm_diff
                << " || r || / || b || = " << rel_resid << "\n";
  }
  else {
    QDPIO::cout << "Full: || b || = " << norm_src
              <<  " || r || = " << norm_diff <<"\n";
  }
}

TEST(UnprecOpTest, TestUnprecCloverOps)
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

   LatticeFermion source;
   gaussian(source);


   QPhiXFullSpinor qphixSource(geom);
   QPhiXFullSpinor qphixResult(geom);

   // Now pack fields.
   QPhiX::qdp_pack_cb_spinor<>(source, qphixSource.getCBData(evenCb),geom, evenCb);
   QPhiX::qdp_pack_cb_spinor<>(source, qphixSource.getCBData(oddCb),geom, oddCb);


   // Create QPhiXClovOp
   ClovOp QPhiXEOClov(qphixGauge,
       qphixClov,
       AInvEe.get(),
       &geom,
       tBoundary,
       anisoFacS,
       anisoFacT);

   LatticeFermion result;
   LatticeFermion qdp_result;
   LatticeFermion diff=zero;
   Double norm_diff = Double(0);
   Double vol = Double(Layout::vol());
   Double vol_cb = Double(Layout::vol()/2);

   for(int cb =0; cb < 2; ++cb) {
     QPhiXEOClov.M_diag(qphixResult,qphixSource,1,cb);
     QPhiX::qdp_unpack_cb_spinor<>(qphixResult.getCBData(cb),result,geom,cb);
     clov.apply(qdp_result,source,1,cb);

     diff[rb[cb]]=qdp_result - result;
     norm_diff = sqrt(norm2(diff,rb[cb]));
     QDPIO::cout << "cb=" << cb << " || diff || = " << norm_diff << " || diff || / site = " << norm_diff/vol_cb
              << std::endl;
     ASSERT_LT( toDouble(norm_diff), 5.0e-14);
   }
     // Apply the diagonal operator
   QPhiXEOClov.M_diag(qphixResult,qphixSource,1);
   QPhiX::qdp_unpack_cb_spinor<>(qphixResult.getCBData(evenCb),result,geom,evenCb);
   QPhiX::qdp_unpack_cb_spinor<>(qphixResult.getCBData(oddCb),result,geom,oddCb);

   for(int cb=0; cb < 2; ++cb) {
     clov.apply(qdp_result,source,1,cb);
   }

   QDPIO::cout << "Checking Clover Op Diagonal Part" << std::endl;

   for(int cb=0; cb < 2; ++cb) {
     norm_diff = sqrt(norm2(diff,rb[cb]));
     QDPIO::cout << "cb=" << cb << " || diff || = " << norm_diff << " || diff || / site = " << norm_diff/vol_cb
         << std::endl;
     ASSERT_LT( toDouble(norm_diff), 5.0e-14);
   }
   norm_diff = sqrt(norm2(diff));
   QDPIO::cout << "Full: || diff || = " << norm_diff << " || diff || / site = " << norm_diff/vol
       << std::endl;
   ASSERT_LT( toDouble(norm_diff), 7.0e-14);


  QDPIO::cout << "Checking UnprecOp Part" << std::endl;
  for(int isign=-1; isign < +2; isign+=2) {
   QPhiXEOClov.M_unprec(qphixResult,qphixSource,isign);
   QPhiX::qdp_unpack_cb_spinor<>(qphixResult.getCBData(evenCb),result,geom,evenCb);
   QPhiX::qdp_unpack_cb_spinor<>(qphixResult.getCBData(oddCb),result,geom,oddCb);

     for(int cb=0; cb < 2; ++cb) {
       clov.apply(qdp_result,source,isign,cb);
       LatticeFermion tmp_dsl;
       dslash(tmp_dsl,uAniso,source,isign,cb);
       qdp_result[rb[cb]] -= Real(0.5)*tmp_dsl;
     }
     diff = qdp_result - result;
     for(int cb=0; cb < 2; ++cb) {
       norm_diff = sqrt(norm2(diff,rb[cb]));
       QDPIO::cout << "cb=" << cb << " isign = " << isign <<" || diff || = " << norm_diff << " || diff || / site = " << norm_diff/vol_cb
           << std::endl;
       ASSERT_LT( toDouble(norm_diff), 5.0e-14);
     }
     norm_diff = sqrt(norm2(diff));
     QDPIO::cout << "Full: isign = " << isign << " || diff || = " << norm_diff << " || diff || / site = " << norm_diff/vol
         << std::endl;
     ASSERT_LT( toDouble(norm_diff), 7.0e-14);

  }
}


int main(int argc, char *argv[])
{
	return QDPXXTestMain(&argc, argv);
}


