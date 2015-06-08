#include <iostream>
#include "unittest.h"
#include "testMInvCG.h"

#undef SEEK_SET
#undef SEEK_CUR
#undef SEEK_END
#include "qdp.h"
using namespace QDP;

#ifndef DSLASH_M_W_H
#include "dslashm_w.h"
#endif

#ifndef REUNIT_H
#include "reunit.h"
#endif
#include "qphix/qphix_config.h"

#include "qphix/geometry.h"
#include "qphix/qdp_packer.h"
#include "qphix/blas_new_c.h"
#include "qphix/wilson.h"
#include "qphix/minvcg.h"

#include <omp.h>

using namespace Assertions;
using namespace std;
using namespace QPhiX;

#ifndef QPHIX_SOALEN
#define QPHIX_SOALEN 4
#endif

#if defined(QPHIX_MIC_SOURCE)

#define VECLEN_SP 16 
#define VECLEN_HP 16
#define VECLEN_DP 8
#include <immintrin.h>

#elif defined(QPHIX_AVX_SOURCE) 

#define VECLEN_SP 8
#define VECLEN_DP 4

#elif defined(QPHIX_SCALAR_SOURCE)
#warning SCALAR_SOURCE
#define VECLEN_SP 1
#define VECLEN_DP 1

#elif defined(QPHIX_QPX_SOURCE)
#warning QPX_SOURCE
#define VECLEN_SP 4
#define VECLEN_DP 4

#endif

  // What we consider to be small enough...
int Nx, Ny, Nz, Nt, Nxh;
bool verbose = false;

template<typename F> 
struct tolerance { 
  static const Double small; // Always fail
};

template<>
const Double tolerance<half>::small = Double(5.0e-3);

template<>
const Double tolerance<float>::small = Double(1.0e-6);


template<>
const Double tolerance<double>::small = Double(1.0e-13);

template<typename T>
struct rsdTarget { 
  static const double value;
};

template<>
const double rsdTarget<half>::value = (double)(1.0e-4);

template<>
const double rsdTarget<float>::value = (double)(1.0e-7);


template<>
const double rsdTarget<double>::value = (double)(1.0e-12);

template<typename T, int V, int S, bool compress, typename U, typename Phi>
void 
MInvCGTester::testMInvCG(const U& u, int t_bc)
{
  RNG::setrn(rng_seed);
  typedef typename Geometry<T,V,S,compress>::SU3MatrixBlock   Gauge;
  typedef typename Geometry<T,V,S,compress>::FourSpinorBlock  Spinor;

  int threads_per_core = Sy*Sz;

  QDPIO::cout << "In testMInvCG:" << endl;

  Phi chi;
  QDPIO::cout << "Filling chi with gaussian noise" << endl;
  gaussian(chi);
  
  QDPIO::cout << "Norm2 || psi || = " << norm2(chi,rb[0]) << endl;
  QDPIO::cout << "Done" << endl;
  double aniso_fac_s = (double)(0.35);
  double aniso_fac_t = (double)(1.4);
  double t_boundary = (double)(t_bc);
  

  Geometry<T,V,S,compress> geom(Layout::subgridLattSize().slice(), By, Bz, NCores, Sy, Sz, PadXY, PadXYZ, MinCt);
  
  // NEED TO MOVE ALL THIS INTO DSLASH AT SOME POINT
  // -- Allocate the gauge field 
  Gauge* packed_gauge_cb0 = (Gauge*)geom.allocCBGauge();
  Gauge* packed_gauge_cb1 = (Gauge*)geom.allocCBGauge();
  Gauge* u_packed[2] = { packed_gauge_cb0, packed_gauge_cb1};

  // Allocate the right hand side
  Spinor *chi_d = (Spinor *)geom.allocCBFourSpinor();
  if( chi_d == 0x0 ) { 
    QDPIO::cout << "Error allocating chi" << endl;
    QDP_abort(1);
  }
  
  int n_shift=4;
  double shifts[4] = {0.01,0.02,0.03,0.04};

  // Allocate the solutions
  Spinor *psi_d[4];
  for(int s=0; s < n_shift; s++) { 
    psi_d[s] = (Spinor*)geom.allocCBFourSpinor();
    if ( psi_d[s] == 0x0 ) { 
      QDPIO::cout << "Error allocating psi[" << s <<"]" << endl;
      QDP_abort(1);
    }
  }

  QDPIO::cout << "Fields allocated" << endl;
 
  // Pack the gauge field
  QDPIO::cout << "Packing gauge field..." ;
  qdp_pack_gauge<>(u, packed_gauge_cb0,packed_gauge_cb1, geom);
  QDPIO::cout << "done" << endl;

  QDPIO::cout << " Packing chi" ;	
  qdp_pack_cb_spinor<>(chi, chi_d, geom,0);

  QDPIO::cout << " Zeroing the psis"; 
  for(int s=0; s < n_shift; s++) { 
    zeroSpinor<T,V,S,compress>(psi_d[s],geom, threads_per_core);
  }
  QDPIO::cout << "done" << endl; 
  QDPIO::cout << "T BCs = " << t_boundary << endl;

  QDPIO::cout << "Applying anisotropy to test gauge field" << endl;
  U u_test(Nd);
  for(int mu=0; mu < Nd; mu++) { 
    Real factor = Real(aniso_fac_s);
    if( mu == Nd -1) { 
      factor = Real(aniso_fac_t);
    }
    u_test[mu] = factor*u[mu];
  }
  // Apply BCs on u-test for QDP++ test (Dslash gets unmodified field)
  u_test[3] *= where(Layout::latticeCoordinate(3) == (Layout::lattSize()[3]-1),
  			Real(t_boundary), Real(1));

  double Mass=0.01;
  EvenOddWilsonOperator<T,V,S,compress> M(Mass, u_packed, &geom, t_boundary, aniso_fac_s, aniso_fac_t);
  Phi ltmp=zero;
  Real massFactor=Real(4) + Real(Mass);
  Real betaFactor=Real(0.25)/massFactor;
  double rsd_target[4];
  double rsd_final[4];
  
  for(int s=0; s < n_shift; s++) { 
    rsd_target[s]=rsdTarget<T>::value;
  }
  int max_iters=200;
  int niters;
  unsigned long site_flops;
  unsigned long mv_apps;
  
  MInvCG<T,V,S,compress> solver(M, max_iters, n_shift);
  //  solver.tune();
  double r2;
  norm2Spinor<T,V,S,compress>(r2,chi_d,geom,threads_per_core);
  masterPrintf("chi has norm2 = %16.8e\n", r2);

  int isign=1;
  double start = omp_get_wtime();
  
  solver(psi_d, chi_d, n_shift,shifts, rsd_target, niters, rsd_final, site_flops, mv_apps, isign, verbose);

  
  double end = omp_get_wtime();

  QDPIO::cout << "Solver Completed. Iters = " << niters << " Wallclock = " << end -start << " sec." << endl;

  // check solutions
  Phi psi,psi2,psi3;
  for(int s=0; s < n_shift; s++) {
    Real shiftFactor = Real(shifts[s]);
    // Unpack solution s into 'psi'
    qdp_unpack_cb_spinor<>(psi_d[s], psi, geom, 0);
 

    // psi2 = M psi
    dslash(psi2,u_test,psi, isign, 1);
    dslash(ltmp,u_test,psi2, isign, 0);
    psi2[rb[0]] = massFactor*psi - betaFactor*ltmp;
    
    // psi3 = M^\dagger psi2
    dslash(psi3,u_test,psi2, (-isign), 1);
    dslash(ltmp,u_test,psi3, (-isign), 0);
    psi3[rb[0]] = massFactor*psi2 - betaFactor*ltmp;

    // psi2 = psi3 + shift * psi = M^\dagger M psi + shift psi;
    psi2[rb[0]] = shiftFactor*psi + psi3;
    
    Phi diff = chi - psi2;
    Double true_norm = sqrt(norm2(diff, rb[0])/norm2(chi,rb[0]));
    QDPIO::cout << "True norm for shift["<<s<<"]="<< shifts[s] << " is: " << true_norm << endl;
    //    assertion( toBool( true_norm < (rsd_target + tolerance<T>::small) ) );
  }

  unsigned long num_cb_sites=Layout::vol()/2;
  unsigned long total_flops = (site_flops + (72+2*1320)*mv_apps)*num_cb_sites;
  
  masterPrintf("GFLOPS=%e\n", 1.0e-9*(double)(total_flops)/(end -start));
  
  geom.free(packed_gauge_cb0);
  geom.free(packed_gauge_cb1);
  geom.free(chi_d);
  for(int i=0; i < 4; i++){ 
    geom.free(psi_d[i]);
  }


}

void
MInvCGTester::run(void) 
{
  RNG::savern(rng_seed);

  typedef multi1d<LatticeColorMatrixF> UF;
  typedef multi1d<LatticeColorMatrixD> UD;
  typedef LatticeDiracFermionF PhiF;
  typedef LatticeDiracFermionD PhiD;

  // Diagnostic information:
  const multi1d<int>& lattSize = Layout::subgridLattSize();
  Nx = lattSize[0];
  Ny = lattSize[1];
  Nz = lattSize[2];
  Nt = lattSize[3];

  QDPIO::cout << "Lattice Size: ";
  for(int mu=0; mu < lattSize.size(); mu++){ 
    QDPIO::cout << " " << lattSize[mu];
  }
  QDPIO::cout << endl;

  QDPIO::cout << "Block Sizes: By=" << By << " Bz=" << Bz << endl;
  QDPIO::cout << "N Cores" << NCores << endl;
  QDPIO::cout << "SMT Grid: Sy=" << Sy << " Sz=" << Sz << endl;
  QDPIO::cout << "Pad Factors: PadXY=" << PadXY << " PadXYZ=" << PadXYZ << endl;
  QDPIO::cout << "MinCt=" << MinCt << endl;
  QDPIO::cout << "Threads_per_core = " << N_simt << endl;

 


  QDPIO::cout << "Inititalizing QDP++ gauge field"  << endl;
  // Make a random gauge field 
  multi1d<LatticeColorMatrix> u(4);
  LatticeColorMatrix g;
  LatticeColorMatrix uf;
  for(int mu=0; mu < 4; mu++) { 
    uf = 1;   // Unit gauge

    Real factor=Real(0.09);
    gaussian(g);
    u[mu] = uf + factor*g;
    reunit(u[mu]);
  }



#if 1 // Save build time
  if( precision == FLOAT_PREC ) {
    
    QDPIO::cout << "SINGLE PRECISION TESTING:" << endl;
    multi1d<LatticeColorMatrixF> u_in(4);
    for(int mu=0; mu < Nd; mu++) {
      u_in[mu] = u[mu];
    }

    {

      if( soalen == 1 ) { 
#if defined(QPHIX_SCALAR_SOURCE)
	QDPIO::cout << "VECLEN = " << VECLEN_SP << " SOALEN=1 " << endl;
	testMInvCGWrapper<float,VECLEN_SP,1,UF,PhiF>(u_in);
#endif
      }

#if 0
      if( soalen == 4 ) { 
#if defined (QPHIX_AVX_SOURCE) || defined(QPHIX_MIC_SOURCE)
	QDPIO::cout << "VECLEN = " << VECLEN_SP << " SOALEN=4 " << endl;
	testMInvCGWrapper<float,VECLEN_SP,4,UF,PhiF>(u_in);
#endif
      }
#endif

      if( soalen == 8 ) {
	QDPIO::cout << "In the SOALEN =8 branch" << endl;

#if defined (QPHIX_AVX_SOURCE) || defined(QPHIX_MIC_SOURCE)
	QDPIO::cout << "VECLEN = " << VECLEN_SP << " SOALEN=8"  << endl;
	testMInvCGWrapper<float,VECLEN_SP,8,UF,PhiF>(u_in);
#endif
      }


      if ( soalen == 16 ) { 
#if defined(QPHIX_MIC_SOURCE)
	QDPIO::cout << "VECLEN = " << VECLEN_SP << " SOALEN=16 " << endl;
	testMInvCGWrapper<float,VECLEN_SP,16,UF,PhiF>(u_in);
#else 
	masterPrintf("SOALEN=16 not available");
	return;
#endif
      }


    }
  }

#endif // FLOAT PREC

#if 1
  if (precision == HALF_PREC ) { 
#if defined(QPHIX_MIC_SOURCE)
    QDPIO::cout << "HALF PRECISION TESTING:" << endl;
    multi1d<LatticeColorMatrixF> u_in(4);
    for(int mu=0; mu < Nd; mu++) {
      u_in[mu] = u[mu];
    }
    {
      if( soalen == 4 ) { 
	QDPIO::cout << "VECLEN = " << VECLEN_HP << " SOALEN=4 " << endl;
	testMInvCGWrapper<half,VECLEN_HP,4,UF,PhiF>(u_in);
      }

      if (soalen == 8 ) {
	QDPIO::cout << "VECLEN = " << VECLEN_HP << " SOALEN=8 " << endl;
	testMInvCGWrapper<half,VECLEN_HP,8,UF,PhiF>(u_in);
      }

      if( soalen == 16 ) {
	QDPIO::cout << "VECLEN = " << VECLEN_HP << " SOALEN=16 " << endl;
	testMInvCGWrapper<half,VECLEN_HP,16,UF,PhiF>(u_in);
      }
    }
#else
    QDPIO::cout << " Half Prec is only supported on MIC Targets just now " << endl;
#endif
  }
#endif // If HALF-PREC

#if 1
  if( precision == DOUBLE_PREC ) { 
    QDPIO::cout << "DOUBLE PRECISION TESTING:" << endl;
    UD u_in(4);
    for(int mu=0; mu < Nd; mu++) {
      u_in[mu] = u[mu];
    }
    
    {

      if( soalen == 1) {
#if defined (QPHIX_SCALAR_SOURCE)
	QDPIO::cout << "VECLEN = " << VECLEN_DP << " SOALEN=1 " << endl;
	testMInvCGWrapper<double,VECLEN_DP,1,UD,PhiD>(u_in);
#endif
      }


      if( soalen == 2) {
#if defined (QPHIX_AVX_SOURCE)
	QDPIO::cout << "VECLEN = " << VECLEN_DP << " SOALEN=2 " << endl;
	testMInvCGWrapper<double,VECLEN_DP,2,UD,PhiD>(u_in);
#endif
      }

      if( soalen == 4 ) { 
	testMInvCGWrapper<double,VECLEN_DP,4,UD,PhiD>(u_in);
      }

      if( soalen == 8 ) { 
#if defined(QPHIX_MIC_SOURCE)
	QDPIO::cout << "VECLEN = " << VECLEN_DP << " SOALEN=8 " << endl;
	testMInvCGWrapper<double,VECLEN_DP,8,UD,PhiD>(u_in);
#endif
      }
      
    }
  }
#endif // If DOUBLE-PREC

}



