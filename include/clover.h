#ifndef CLOVER_H
#define CLOVER_H

#include "linearOp.h"
#include "cpp_clovdslash_scalar.h"

namespace CPlusPlusWilsonDslash { 

  template<typename FT, int veclen, int soalen, bool compress12>
    class EvenOddCloverOperator : public EvenOddLinearOperator<FT,veclen,soalen,compress12>{
  public: 
    typedef typename Geometry<FT,veclen,soalen,compress12>::CloverBlock CloverBlock;
    typedef typename Geometry<FT,veclen,soalen,compress12>::FourSpinorBlock FourSpinorBlock;
    typedef typename Geometry<FT,veclen,soalen,compress12>::SU3MatrixBlock SU3MatrixBlock;
    
    // Constructor
    // No anisotropy, all boundaries periodic for now.
    EvenOddCloverOperator(const int *latt_size, 
			  SU3MatrixBlock* u_[2],
			  CloverBlock* clov_,
			  CloverBlock* invclov_,
			  int By, 
			  int Bz, 
			  int NCores,
			  int Sy,
			  int Sz,
			  int PadXY, 
			  int PadXYZ, 
			  int MinCt,
			  double t_boundary,
			  double aniso_coeff_s,
			  double aniso_coeff_t): D(new ClovDslash<FT, veclen,soalen,compress12>(latt_size,By,Bz,NCores,Sy,Sz,PadXY,PadXYZ,MinCt,t_boundary, aniso_coeff_s, aniso_coeff_t)), 
      geom(D->getGeometry()) 
	{

	  u[0] = u_[0];
	  u[1] = u_[1];
	  clov = clov_;
	  invclov = invclov_;
	  tmp = (FourSpinorBlock *)geom.allocCBFourSpinor();
	}
    
  EvenOddCloverOperator(const int *latt_size, 
			int By, 
			int Bz, 
			int NCores,
			int Sy,
			int Sz,
			int PadXY, 
			int PadXYZ, 
			int MinCt,
			bool compress12_,
			double t_boundary,
			double aniso_coeff_s,
			double aniso_coeff_t): D(new ClovDslash<FT, veclen,soalen>(latt_size,By,Bz,NCores,Sy,Sz,PadXY,PadXYZ,MinCt,compress12_, t_boundary, aniso_coeff_s, aniso_coeff_t)), geom(D->getGeometry())
      {
	tmp = (FourSpinorBlock *)geom.allocCBFourSpinor();
      }
    
    void setFields(SU3MatrixBlock* u_[2], CloverBlock* clov_, CloverBlock* invclov_)
    {
      u[0] = u_[0];
      u[1] = u_[1];
      clov = clov_;
      invclov = invclov_;
    }

    ~EvenOddCloverOperator() {
      geom.free(tmp);
      delete D;
    }


    void operator()(FourSpinorBlock *res, const FourSpinorBlock* in, int isign) {
      double beta=(double)0.25;

      D->dslash(tmp, in, u[0], invclov,  isign, 0);
      D->dslashAChiMinusBDPsi(res, tmp, in, u[1], clov, beta, isign, 1);
    }


    Geometry<FT,veclen, soalen, compress12>& getGeometry() { return D->getGeometry(); }
  private:
    double Mass;
    ClovDslash<FT, veclen,soalen,compress12>* D;
    Geometry<FT,veclen,soalen,compress12>& geom;
    SU3MatrixBlock *u[2];
    FourSpinorBlock *tmp;
    CloverBlock* clov;
    CloverBlock* invclov;

						   }; // Class
}; // Namespace

    

#endif
