#ifndef QPHIX_MINVCG2_H
#define QPHIX_MINVCG2_H

///// USER SERVICEABLE OPTIONS -- SHOULD BE MOVED TO AUTOCONF CONTROL
#include "qphix/linearOp.h"
#include "qphix/print_utils.h"
#include "qphix/blas_new_c.h"
#include "qphix/tsc.h"
#include "qphix/abs_solver.h"

namespace  QPhiX
{


  // That will be a second refactor step when everything works
  template<typename FT, int veclen, int soalen, bool compress12>
  class MInvCG : public AbstractMultiSolver<FT,veclen,soalen,compress12> {
  public: 
    typedef typename Geometry<FT,veclen,soalen,compress12>::FourSpinorBlock Spinor;
    MInvCG(EvenOddLinearOperator<FT,veclen,soalen,compress12>& M_,
	   int MaxIters_,
	   int MaxShifts_) : M(M_), geom(M_.getGeometry()), MaxIters(MaxIters_), MaxShifts(MaxShifts_)
    {
      
      // Length of vectors in floats, including PADs
      masterPrintf("Initializing Multi-Shift CG Solver with n_shift = %d: Nvec=%d Ny=%d Nz=%d Nt=%d\n", MaxShifts, geom.nVecs(), geom.Ny(), geom.Nz(), geom.Nt());
      
      // Initialie norm2_threads to all
      norm2_threads=geom.getNSIMT();
      copy_threads=geom.getNSIMT();
      axpy_threads=geom.getNSIMT();
      axpby_threads=geom.getNSIMT();
      axy_threads=geom.getNSIMT();
      axpynorm2_threads=geom.getNSIMT();

      allocateSpace();
    }

    ~MInvCG() {
      freeSpace();
    }
    

    void operator()(Spinor **x, 
		    const Spinor *rhs, 
		    const int n_shift,
		    const double* shifts,
		    const double* RsdTarget,
		    int& n_iters, 
		    double* rsd_sq_final, 
		    unsigned long& site_flops,
		    unsigned long& mv_apps, 
		    int isign,
		    bool verboseP) 
    {
      
      mv_apps = 0;
      site_flops=0;
      if( n_shift > MaxShifts) { 
	masterPrintf("Too many shifts requested: n_shift=%d MaxShifts=%d\n", 
		     n_shift, MaxShifts);
	exit(-1);
      }

#if 0 
      // Find index of the smallest shift
      int smallest=0;
      for(int s=1; s < n_shift; s++) { 
	if ( shifts[s] < shifts[smallest] ) smallest = s;
      }
#endif
      // Get the norm of the rhs
      double chi_norm_sq, chi_norm;
      norm2Spinor(chi_norm_sq, rhs,geom, norm2_threads);
      site_flops += 47;

      double cp = chi_norm_sq;
      for(int s=0; s < n_shift; s++) { 
	rsd_sq[s] = RsdTarget[s]*RsdTarget[s]*cp;
      }
      
      // r[0] = p_0 = p[s] = rhs
      copySpinor(r,rhs, geom, copy_threads);
      copySpinor(p_0,rhs,geom, copy_threads);
      for(int s=0; s < n_shift; s++) { 
	copySpinor(p[s],rhs,geom,copy_threads);
      }

      // MMp = M^\dagger M p_0 
      double d;
      M(mp,p_0,+1);
      M(mmp,mp,-1);  
      mv_apps +=2;

      // d = norm2(Mp) = < M p | M p > =  <p | M^\dagger M p>

      norm2Spinor(d,mp,geom, norm2_threads);
      site_flops += 47;  // 3*Nc*Ns + Nc*Ns-1

      double b = -cp / d;
      double c;
      // r1 = r + b MM p
      // c = norm2(r)
      //      axpy(b,mmp,r, geom, axpy_threads);
      //norm2Spinor(c,r,geom, norm2_threads);
      axpyNorm2(b,mmp,r,c,geom,axpynorm2_threads);
      site_flops += 95; 

      // Initialize zetas and betas
      for(int s=0; s < n_shift; s++) { 
	zeta_prev[s] = (double)1;
	zeta[s] = (double)1/ ( (double)1 - shifts[s]*b );
	beta[s] = b/( (double)1 - shifts[s]*b );
	convsP[s]=false;
      }      
     

      for(int s=0; s < n_shift; s++) { 
	double mbeta = -beta[s];
	axy(mbeta, rhs, x[s], geom,axy_threads); // x[s] = -beta[s]*rhs
      }
      site_flops += 24*n_shift;

      bool convP =false;
      double z0, z1;
      double ztmp;
      double cs;
      double a;
      double as;
      double bp;
      int k;
      for(k=1; k < MaxIters && !convP; ++k) {

	// a[k+1] = || r ||^2 / || r[k-1] ||^2 = c/p
	a = c/cp;   

	// p_0 = r + a p_0  -- aypx op: a = a, x = r, y = p0
	aypx(a,r,p_0, geom, axpy_threads);


	// update shifted p-s
	for(int s=0; s < n_shift; s++) { 
	  if ( ! convsP[s] ) { 
	    // only update unconverged systems
	    as = a*zeta[s]*beta[s]/ (zeta_prev[s]*b);
	    axpby(zeta[s], r, as, p[s], geom, axpby_threads);
	    site_flops+=72;
	  }
	}

	cp = c;
	bp = b;

	M(mp,p_0,+1);
	M(mmp,mp,-1);  
	mv_apps += 2;

	// d = norm2(Mp) = < M p | M p > =  <p | M^\dagger M p>
	norm2Spinor(d,mp,geom, norm2_threads);
	b = -cp/d;	
	// r = r + b MM p
	// c = norm2(r)
	//	axpy(b,mmp,r, geom, axpy_threads);
	// norm2Spinor(c,r,geom, norm2_threads);
	axpyNorm2(b,mmp,r,c,geom, axpynorm2_threads);

	site_flops += 190;    // 2 norms + 2 axpy

	convP = true;

	for(int s=0; s < n_shift; s++) { 
	  if( !convsP[s] ) { 
	    z0 = zeta[s];
	    z1 = zeta_prev[s];
	    zeta_prev[s] = z0;

	    zeta[s] = z0*z1*bp;
	    zeta[s] /=(  b*a*(z1-z0) + z1*bp*( (double)1 - shifts[s]*b ) );
	    beta[s] =  b*zeta[s] / z0;
	    
	    double mbeta = -beta[s];
	    axpy( mbeta,p[s],x[s], geom, axpy_threads);
	    site_flops += 48;

	    double css = c*zeta[s]*zeta[s];

	    if( verboseP ) {
	      masterPrintf("iter=%d shift=%d css=%16.8e target=%16.8e\n",
			   k, s, css, rsd_sq[s]);
	    }

	    convsP[s] = ( css < rsd_sq[s] ) ;
	    if( convsP[s] ) { 
	      rsd_sq_final[s] = css;
	    }
	  }
	  convP &= convsP[s];
	}
	n_iters = k; // Note iteration in case converged.
      }
    


      // freeSpace();
      return;
    }

    

    void tune()
    {
    }

    Geometry<FT,veclen,soalen,compress12>& getGeometry() {
      return geom;
    }


    private:

    EvenOddLinearOperator<FT, veclen,soalen,compress12>& M;
    Geometry<FT, veclen,soalen,compress12>& geom;
    int MaxIters;
    int MaxShifts;

    Spinor *mp;
    Spinor *mmp;
    Spinor *p_0;
    Spinor **p;
    Spinor *r;
    double *rsd_sq;
    double *beta;
    double *zeta_prev;
    double *zeta;
    bool *convsP;

    int norm2_threads;
    int copy_threads;
    int axpy_threads;
    int axpby_threads;
    int axy_threads;
    int axpynorm2_threads;

    void allocateSpace(void) {
      mp = (Spinor*)geom.allocCBFourSpinor();
      if( mp == 0x0 ) { 
	masterPrintf("MInvCG Failed to allocate mp\n");
	exit(-1);
      }

      mmp = (Spinor*)geom.allocCBFourSpinor();
      if( mmp == 0x0 ) { 
	masterPrintf("MInvCG Failed to allocate mmp\n");
	exit(-1);
      }

      p =(Spinor**)ALIGNED_MALLOC(sizeof(Spinor*)*MaxShifts,QPHIX_LLC_CACHE_ALIGN);
      if( p == 0x0 ) { 
	masterPrintf("MInvCG Failed to allocate p-array\n");
	exit(-1);
      }

      for(int s=0; s < MaxShifts; s++) { 
	p[s] = (Spinor*)geom.allocCBFourSpinor();
	if( p[s] == 0x0 ) { 
	  masterPrintf("MInvCG Failed to allocate p[%d]\n",s);
	  exit(-1);
	}
      
      }
      r = (Spinor *)geom.allocCBFourSpinor();
      if( r == 0x0 ) { 
	masterPrintf("MInvCG Failed to allocate r\n");
	exit(-1);
      }
      p_0 = (Spinor *)geom.allocCBFourSpinor();
      if( p_0 == 0x0 ) { 
	masterPrintf("MInvCG Failed to allocate p0\n");
	exit(-1);
      }

      rsd_sq = new double[MaxShifts];
      if ( rsd_sq == 0x0 ) { 
	masterPrintf("Unable to allocate rsd_sq array\n");
	exit(-1);
      }

      beta = (double *)ALIGNED_MALLOC(MaxShifts*sizeof(double), QPHIX_LLC_CACHE_ALIGN);
      if( beta == 0x0 ) { 
	masterPrintf("Unable to allocate beta array\n");
	exit(-1);
      }

      zeta_prev = (double *)ALIGNED_MALLOC(MaxShifts*sizeof(double), QPHIX_LLC_CACHE_ALIGN);
      if( zeta_prev == 0x0 ) { 
	masterPrintf("Unable to allocate zeta_prev\n");
	exit(-1);
      }
      

      zeta = (double *)ALIGNED_MALLOC(MaxShifts*sizeof(double), QPHIX_LLC_CACHE_ALIGN);
      if( zeta == 0x0 ) { 
	masterPrintf("Unable to allocate zeta\n");
	exit(-1);
      }

      convsP = (bool *)ALIGNED_MALLOC(MaxShifts*sizeof(bool), QPHIX_LLC_CACHE_ALIGN);
      if( zeta == 0x0 ) { 
	masterPrintf("Unable to allocate convsP\n");
	exit(-1);
      }
      
    }

    void freeSpace(void) {
      geom.free(mp);
      geom.free(mmp);
      geom.free(r);
      geom.free(p_0);
      for(int s=0; s < MaxShifts; s++) { 
	geom.free(p[s]);
      }
      ALIGNED_FREE(p);
      delete [] rsd_sq;
      ALIGNED_FREE(beta);
      ALIGNED_FREE(zeta_prev);
      ALIGNED_FREE(zeta);
      ALIGNED_FREE(convsP);
    }

  };



};



#endif
