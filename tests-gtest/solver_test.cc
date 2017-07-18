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
#include <qphix/invbicgstab.h>
#include <qphix/qdp_packer.h>


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
using BiCGStab = QPhiX::InvBiCGStab<double,VECLEN_DP,QPHIX_SOALEN,false>;
using QPhiXCBSpinor = QPhiX::FourSpinorHandle<double,VECLEN_DP,QPHIX_SOALEN,false>;
using QPhiXCBGauge = QPhiX::GaugeHandle<double,VECLEN_DP,QPHIX_SOALEN, false>;
using QPhiXCBClover = QPhiX::CloverHandle<double, VECLEN_DP,QPHIX_SOALEN, false>;


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

int main(int argc, char *argv[])
{
	return QDPXXTestMain(&argc, argv);
}


