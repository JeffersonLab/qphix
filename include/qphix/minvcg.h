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
	   int n_shift_) : M(M_), geom(M_.getGeometry()), MaxIters(MaxIters_), n_shift(n_shift_)
    {
      
      // Length of vectors in floats, including PADs
      masterPrintf("Initializing CG Solver: Nvec=%d Ny=%d Nz=%d Nt=%d\n", geom.nVecs(), geom.Ny(), geom.Nz(), geom.Nt());
      
    }

    ~MInvCG() {

    }
    

    void operator()(Spinor *x[], 
		    const Spinor *rhs, 
		    const double shifts[],
		    const double RsdTarget[],
		    int& n_iters, 
		    double rsd_sq_final[], 
		    unsigned long& site_flops,
		    unsigned long& mv_apps, 
		    int isign,
		    bool verboseP) 
    {
      int threads_per_core=2;
     
      return;
    }

    

    void tune()
    {
    }

    Geometry<FT,veclen,soalen,compress12>& getGeometry() {
      return geom;
    }

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

      p =(Spinor**)ALIGNED_MALLOC(sizeof(Spinor*)*n_shift,QPHIX_LLC_CACHE_ALIGN);
      if( p == 0x0 ) { 
	masterPrintf("MInvCG Failed to allocate p-array\n");
	exit(-1);
      }

      for(int s=0; s < n_shift; s++) { 
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

    }

    void freeSpace(void) {
      geom.free(mp);
      geom.free(mmp);
      geom.free(r);
      for(int s=0; s < n_shift; s++) { 
	geom.free(p[s]);
      }
      ALIGNED_FREE(p);
    }

    private:

    EvenOddLinearOperator<FT, veclen,soalen,compress12>& M;
    Geometry<FT, veclen,soalen,compress12>& geom;
    int MaxIters;
    int n_shift;

    Spinor *mp;
    Spinor *mmp;
    Spinor **p;
    Spinor *r;



  };



};



#endif
