#ifndef QPHIX_MINVCG2_H
#define QPHIX_MINVCG2_H

///// USER SERVICEABLE OPTIONS -- SHOULD BE MOVED TO AUTOCONF CONTROL
#include "qphix/linearOp.h"
#include "qphix/print_utils.h"
#include "qphix/blas_new_c.h"
#include "qphix/tsc.h"
#include "qphix/abs_solver.h"

#include <vector>

#undef DEBUG_MINVCG

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

      masterPrintf("MinvCG::MinvCG: Allocating Space");
      allocateSpace();
      masterPrintf("MinvCG::~MinvCG: Done\n");
    }

    ~MInvCG() {
      masterPrintf("MInvCG::~MInvCG: Freeing Space\n");
      freeSpace();
      masterPrintf("MinvCG::~MInvCG: Done\n");

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
		    bool verboseP,
		    int cb=1)  const
    {
      
      mv_apps = 0;
      site_flops=0;
      if( n_shift > MaxShifts) { 
	masterPrintf("Too many shifts requested: n_shift=%d MaxShifts=%d\n", 
		     n_shift, MaxShifts);
	exit(-1);
      }

#ifdef DEBUG_MINVCG
      /* Sanity Check */
      {
	double n2_tmp = 0;
	norm2Spinor(n2_tmp, rhs, geom, norm2_threads);
	masterPrintf("rhs has norm=%16.8e\n", n2_tmp);
	for(int s=0; s < n_shift; ++s) { 
	  norm2Spinor(n2_tmp, x[s], geom, norm2_threads);
	  masterPrintf("x[%d] has norm=%16.8e\n",s, n2_tmp);
        }
	for(int s=0; s < n_shift; ++s) { 
	  masterPrintf("shift[%d]= %16.8e   RsdTarget[%d]=%16.8e\n", s, shifts[s], s, RsdTarget[s]);
	}
      }
#endif

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

#ifdef DEBUG_MINVCG
      {
	masterPrintf("Initial ChiSqNorm= %16.8e\n", chi_norm_sq);
      }
#endif

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

#ifdef DEBUG_MINVCG 
      {
	double n2_tmp;
	norm2Spinor(n2_tmp, r, geom, norm2_threads);
	masterPrintf("r has norm=%16.8e\n", n2_tmp);

	norm2Spinor(n2_tmp, p_0, geom, norm2_threads);
	masterPrintf("p_0, has norm=%16.8e\n", n2_tmp);

	for(int s=0; s < n_shift; ++s) { 
	  norm2Spinor(n2_tmp, p[s], geom, norm2_threads);
	  masterPrintf("p[%d] has norm=%16.8e\n", s, n2_tmp);
	}
      }
#endif

      // MMp = M^\dagger M p_0 
      double d;
      zeroSpinor(mp, geom, norm2_threads);
      zeroSpinor(mmp,geom, norm2_threads);

#ifdef DEBUG_MINVCG
      {
	double n2_tmp;                                                                                                                    
        norm2Spinor(n2_tmp,mp, geom, norm2_threads);                                                                                      
        masterPrintf("mp has norm=%16.8e\n", n2_tmp);   

        norm2Spinor(n2_tmp,mmp, geom, norm2_threads);                                                                                      
        masterPrintf("mmp has norm=%16.8e\n", n2_tmp);   
      }
#endif

      M(mp,p_0,+1,cb);


#ifdef DEBUG_MINVCG
      {
	double n2_tmp;                                                                                                                    
        norm2Spinor(n2_tmp,mp, geom, norm2_threads);                                                                                      
        masterPrintf("mp has norm=%16.8e\n", n2_tmp);   
      }
#endif


      M(mmp,mp,-1,cb);  

#ifdef DEBUG_MINVCG
      {
	double n2_tmp;                                                                                                                    
        norm2Spinor(n2_tmp,mmp, geom, norm2_threads);                                                                                      
        masterPrintf("mmp has norm=%16.8e\n", n2_tmp);   
      }
      
#endif
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

#ifdef DEBUG_MINVCG 
      for(int s=0; s < n_shift;++s) { 
	masterPrintf("zeta_prev[%d]=%16.8e zeta[%d]=%16.8e beta[%d]=%16.8e convsP[%d]=%d\n",
		     s,zeta_prev[s],s,zeta[s],s,beta[s],s,(int)(convsP[s]));
      }
#endif

      for(int s=0; s < n_shift; s++) { 
	double mbeta = -beta[s];
	axy(mbeta, rhs, x[s], geom,axy_threads); // x[s] = -beta[s]*rhs
      }
      site_flops += 24*n_shift;

#ifdef DEBUG_MINVCG
      for(int s=0; s < n_shift; ++s) {
	double n2_tmp = 0;
	norm2Spinor(n2_tmp, x[s], geom,norm2_threads);
	masterPrintf("x[%d] has norm = ", s, n2_tmp);
      }
#endif

      bool convP =false;
      double z0, z1;
      double ztmp;
      double cs;
      double a;
      double as;
      double bp;
      int k=0;
      for(k=1; k < MaxIters && !convP; ++k) {

	// a[k+1] = || r ||^2 / || r[k-1] ||^2 = c/p
	a = c/cp;   

#ifdef DEBUG_MINVCG
	masterPrintf("a=%16.8e\n", a);
#endif

	// p_0 = r + a p_0  -- aypx op: a = a, x = r, y = p0
	aypx(a,r,p_0, geom, axpy_threads);

#ifdef DEBUG_MINVCG
	{
	  double n2_tmp=0;
	  norm2Spinor(n2_tmp, p_0, geom, norm2_threads);
	  masterPrintf("p_0 has norm=", n2_tmp);
	}
#endif

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

	M(mp,p_0,+1,cb);
	M(mmp,mp,-1,cb);  
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

#ifdef DEBUG_MINVCG
	    masterPrintf("z0=%16.8e, z1=%16.8e, bp=%16.8e\n", z0,z1,bp);
#endif

	    zeta[s] = z0*z1*bp;


#ifdef DEBUG_MINVCG
	    masterPrintf("a=%16.8e b=%16.8e, bp=%16.8e, shifts[s]=%16.8e\n", 
			 a,b,bp,shifts[s]);
#endif

	    zeta[s] /=(  b*a*(z1-z0) + z1*bp*( (double)1 - shifts[s]*b ) );

#ifdef DEBUG_MINVCG
	    masterPrintf("zeta[s]=%16.8e\n",zeta[s]);
#endif
	    beta[s] =  b*zeta[s] / z0;
	    
	    double mbeta = -beta[s];
	    axpy( mbeta,p[s],x[s], geom, axpy_threads);
	    site_flops += 48;


#ifdef DEBUG_MINVCG 
	    masterPrintf("c=%16.8e zeta[%d]=%16.8e\n", c, s, zeta[s]);
#endif

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
    


      return;
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
    mutable std::vector<double> rsd_sq;
    mutable std::vector<double> beta;
    mutable std::vector<double> zeta_prev;
    mutable std::vector<double> zeta;
    mutable std::vector<bool> convsP;

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

      p =new Spinor*[MaxShifts];
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

      rsd_sq.resize(MaxShifts);	
      beta.resize(MaxShifts);
      zeta_prev.resize(MaxShifts);
      zeta.resize(MaxShifts);
      convsP.resize(MaxShifts);
      
    }

    void freeSpace(void) {
      if( mp != 0x0 ) geom.free(mp);
      if ( mmp != 0x0 ) geom.free(mmp);
      if ( r != 0x0 ) geom.free(r);
      if ( p_0 != 0x0 ) geom.free(p_0);
      if ( p != 0x0 ) { 
	for(int s=0; s < MaxShifts; s++) { 
	  if ( p[s] != 0x0 ) geom.free(p[s]);
        }
       
        delete [] p;
      }
    }

  };



};



#endif
