/*
 * full_blas_test.cc
 *
 *  Created on: Oct 16, 2017
 *      Author: bjoo
 */
#include "qdpxx_test_env.h"
#include <qphix/qphix_config.h>
#include <qphix/qphix_cli_args.h>
#include <qphix/geometry.h>
#include <qphix/qdp_packer.h>
#include <qphix/full_spinor.h>
#include <qphix/blas_full_spinor.h>
#include <qphix/unprec_solver_wrapper.h>

#include "../tests/veclen.h"


using namespace QPhiXTesting;
using namespace QDP;
using namespace QPhiX;

using Geom = QPhiX::Geometry<double, VECLEN_DP, QPHIX_SOALEN, false>;
using QPhiXCBSpinor = QPhiX::FourSpinorHandle<double,VECLEN_DP,QPHIX_SOALEN,false>;
using QPhiXFullSpinor = QPhiX::FullSpinor<double,VECLEN_DP,QPHIX_SOALEN,false>;

template<typename FT, int V, int S, bool compress>
void compareSpinors(const LatticeFermion& reference, const QPhiXFullSpinor& qphix_test, QPhiX::Geometry<FT,V,S,compress>& geom, double tol)
{
  LatticeFermion test; qdp_unpack_spinor<>(qphix_test, test, geom);
  LatticeFermion diff = reference - test;
  Double norm_diff=sqrt(norm2(diff));
  Double vol = QDP::Layout::vol();
  Double norm_diff_per_site = norm_diff / vol;
  Double vol_cb = vol/2;
  for(int cb=0; cb < 2; ++cb) {
    Double norm_diff_cb=sqrt(norm2(diff,rb[cb]));

    QDPIO::cout << "CB "<< cb <<": || r || = " << norm_diff_cb
        << " || r || / site = " << norm_diff_cb / vol_cb << "\n";
  }
  QDPIO::cout << "Full: || r || = " << norm_diff << " || r || / site = " << norm_diff_per_site << "\n";
  ASSERT_LT( toDouble(norm_diff), tol);
}

TEST(QPhiXFullBlasTest, FullBlasTestCopySpinor)
{
  // This sets up the environment and the lattice
   const QDPXXTestEnv& testEnv = getQDPXXTestEnv();

   // Get QPhiX Command Line args
   const QPhiX::QPhiXCLIArgs& CLI = testEnv.getCLIArgs();

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
   int n_blas_simt = geom.getNSIMT();

   LatticeFermion  source; gaussian(source);
   QPhiXFullSpinor qphix_source(geom);
   QPhiXFullSpinor qphix_copy(geom);
   qdp_pack_spinor<>(source,qphix_source,geom);

   copySpinor(qphix_copy,qphix_source,geom,n_blas_simt);
   compareSpinors<>(source,qphix_copy,geom, 5.0e-14 );
}

TEST(QPhiXFullBlasTest, FullBlasTestZeroSpinor)
{
  // This sets up the environment and the lattice
   const QDPXXTestEnv& testEnv = getQDPXXTestEnv();

   // Get QPhiX Command Line args
   const QPhiX::QPhiXCLIArgs& CLI = testEnv.getCLIArgs();

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
   int n_blas_simt = geom.getNSIMT();

   LatticeFermion  source=zero;
   QPhiXFullSpinor qphix_source(geom);

   zeroSpinor(qphix_source,geom,n_blas_simt);
   compareSpinors<>(source,qphix_source,geom, 5.0e-14 );
}

TEST(QPhiXFullBlasTest, FullBlasTestAXPYSpinor)
{
  // This sets up the environment and the lattice
   const QDPXXTestEnv& testEnv = getQDPXXTestEnv();

   // Get QPhiX Command Line args
   const QPhiX::QPhiXCLIArgs& CLI = testEnv.getCLIArgs();

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
   int n_blas_simt = geom.getNSIMT();

   LatticeFermion  x; gaussian(x);
   LatticeFermion  y; gaussian(y);

   QPhiXFullSpinor qphix_x(geom);
   QPhiXFullSpinor qphix_y(geom);
   qdp_pack_spinor<>(x,qphix_x,geom);
   qdp_pack_spinor<>(y,qphix_y,geom);
   Real a(0.2);

   // DO QDP++
   y += a*x;


   // DO QPhiX
   double alpha = 0.2;
   axpySpinor(alpha, qphix_x, qphix_y, geom,n_blas_simt);
   compareSpinors<>(y , qphix_y,geom, 5.0e-14 );

}

TEST(QPhiXFullBlasTest, FullBlasTestCAXPYSpinor)
{
  // This sets up the environment and the lattice
   const QDPXXTestEnv& testEnv = getQDPXXTestEnv();

   // Get QPhiX Command Line args
   const QPhiX::QPhiXCLIArgs& CLI = testEnv.getCLIArgs();

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
   int n_blas_simt = geom.getNSIMT();

   LatticeFermion  x; gaussian(x);
   LatticeFermion  y; gaussian(y);

   QPhiXFullSpinor qphix_x(geom);
   QPhiXFullSpinor qphix_y(geom);
   qdp_pack_spinor<>(x,qphix_x,geom);
   qdp_pack_spinor<>(y,qphix_y,geom);
   Complex a(Real(0.2));
   a +=  timesI(Real(1.5));

   // DO QDP++
   y += a*x;


   // DO QPhiX
   double alpha[2] = { 0.2, 1.5 };
   caxpySpinor(alpha, qphix_x, qphix_y, geom,n_blas_simt);

   compareSpinors<>(y , qphix_y,geom, 5.0e-14 );

}
TEST(QPhiXFullBlasTest, FullBlasTestXmYNorm2Spinor)
{
  // This sets up the environment and the lattice
   const QDPXXTestEnv& testEnv = getQDPXXTestEnv();

   // Get QPhiX Command Line args
   const QPhiX::QPhiXCLIArgs& CLI = testEnv.getCLIArgs();

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
   int n_blas_simt = geom.getNSIMT();

   LatticeFermion  x; gaussian(x);
   LatticeFermion  y; gaussian(y);

   QPhiXFullSpinor qphix_x(geom);
   QPhiXFullSpinor qphix_y(geom);
   qdp_pack_spinor<>(x,qphix_x,geom);
   qdp_pack_spinor<>(y,qphix_y,geom);

   // DO QDP++
   y -= x;
   Double norm2_y = norm2(y);

   double norm2_qphix_y;

   xmy2Norm2Spinor( qphix_x, qphix_y, norm2_qphix_y, geom,n_blas_simt);

   compareSpinors<>(y , qphix_y,geom, 5.0e-14 );
   double norm_diff = abs( toDouble(norm2_y - norm2_qphix_y));
   ASSERT_LT(norm_diff, 1.0e-10);
}

TEST(QPhiXFullBlasTest, FullBlasInnerProduct)
{
  // This sets up the environment and the lattice
   const QDPXXTestEnv& testEnv = getQDPXXTestEnv();

   // Get QPhiX Command Line args
   const QPhiX::QPhiXCLIArgs& CLI = testEnv.getCLIArgs();

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
   int n_blas_simt = geom.getNSIMT();

   LatticeFermion  x; gaussian(x);
   LatticeFermion  y; gaussian(y);

   QPhiXFullSpinor qphix_x(geom);
   QPhiXFullSpinor qphix_y(geom);
   qdp_pack_spinor<>(x,qphix_x,geom);
   qdp_pack_spinor<>(y,qphix_y,geom);

   DComplex iprod = innerProduct(x,y);
   double iprod_qphix[2];
   innerProductSpinor<>(iprod_qphix,
                        qphix_x,qphix_y,
                        geom, n_blas_simt);

   double diff_real = abs( toDouble(real(iprod))-iprod_qphix[0]);
   double diff_imag = abs( toDouble(imag(iprod))-iprod_qphix[1]);
   ASSERT_LT(diff_real, 5.0e-13);
   ASSERT_LT(diff_imag, 5.0e-13);
}

TEST(QPhiXFullBlasTest, FullBlasNorm)
{
  // This sets up the environment and the lattice
   const QDPXXTestEnv& testEnv = getQDPXXTestEnv();

   // Get QPhiX Command Line args
   const QPhiX::QPhiXCLIArgs& CLI = testEnv.getCLIArgs();

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
   int n_blas_simt = geom.getNSIMT();

   LatticeFermion  x; gaussian(x);

   QPhiXFullSpinor qphix_x(geom);
   qdp_pack_spinor<>(x,qphix_x,geom);

   Double n = norm2(x);

   double n_qphix;
   norm2Spinor<>(n_qphix,
       qphix_x,
       geom, n_blas_simt);

   double diff = abs( toDouble(n)-n_qphix);

   ASSERT_LT(diff, 1.0e-10);
}


TEST(QPhiXFullBlasTest, FullAx)
{
  // This sets up the environment and the lattice
   const QDPXXTestEnv& testEnv = getQDPXXTestEnv();

   // Get QPhiX Command Line args
   const QPhiX::QPhiXCLIArgs& CLI = testEnv.getCLIArgs();

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
   int n_blas_simt = geom.getNSIMT();

   LatticeFermion  x; gaussian(x);

   QPhiXFullSpinor qphix_x(geom);


   qdp_pack_spinor<>(x,qphix_x,geom);
   double alpha = 6.2;
   Real qalpha(alpha);


   axSpinor<>(alpha,qphix_x,geom,n_blas_simt);
   x *= qalpha;
   compareSpinors<>(x , qphix_x,geom, 5.0e-14 );

}

TEST(QPhiXFullBlasTest, FullYPeqX)
{
  // This sets up the environment and the lattice
   const QDPXXTestEnv& testEnv = getQDPXXTestEnv();

   // Get QPhiX Command Line args
   const QPhiX::QPhiXCLIArgs& CLI = testEnv.getCLIArgs();

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
   int n_blas_simt = geom.getNSIMT();

   LatticeFermion  x; gaussian(x);
   LatticeFermion  y; gaussian(y);

   QPhiXFullSpinor qphix_x(geom);
   QPhiXFullSpinor qphix_y(geom);


   qdp_pack_spinor<>(x,qphix_x,geom);

   qdp_pack_spinor<>(y,qphix_y,geom);


   ypeqxSpinor<>(qphix_x,qphix_y,geom,n_blas_simt);
   y += x;

   compareSpinors<>(x , qphix_x,geom, 5.0e-14 );
   compareSpinors<>(y , qphix_y,geom, 5.0e-14 );

}

TEST(QPhiXFullBlasTest, FullYMeqX)
{
  // This sets up the environment and the lattice
   const QDPXXTestEnv& testEnv = getQDPXXTestEnv();

   // Get QPhiX Command Line args
   const QPhiX::QPhiXCLIArgs& CLI = testEnv.getCLIArgs();

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
   int n_blas_simt = geom.getNSIMT();

   LatticeFermion  x; gaussian(x);
   LatticeFermion  y; gaussian(y);

   QPhiXFullSpinor qphix_x(geom);
   QPhiXFullSpinor qphix_y(geom);


   qdp_pack_spinor<>(x,qphix_x,geom);

   qdp_pack_spinor<>(y,qphix_y,geom);


   ymeqxSpinor<>(qphix_x,qphix_y,geom,n_blas_simt);
   y -= x;

   compareSpinors<>(x , qphix_x,geom, 5.0e-14 );
   compareSpinors<>(y , qphix_y,geom, 5.0e-14 );

}
int main(int argc, char *argv[])
{
  return QDPXXTestMain(&argc, argv);
}
