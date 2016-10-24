#ifndef QPHIX_TM_CLOVER_H
#define QPHIX_TM_CLOVER_H

#include "qphix/linearOp.h"
#include "qphix/tm_clov_dslash_def.h"

namespace QPhiX {

  template <typename FT, int veclen, int soalen, bool compress12>
    class EvenOddTMCloverOperator : public EvenOddLinearOperator<FT,veclen,soalen,compress12> {

      public:

        typedef typename Geometry<FT,veclen,soalen,compress12>::FullCloverBlock FullCloverBlock;
        typedef typename Geometry<FT,veclen,soalen,compress12>::FourSpinorBlock FourSpinorBlock;
        typedef typename Geometry<FT,veclen,soalen,compress12>::SU3MatrixBlock  SU3MatrixBlock;

        // Constructor
        // No anisotropy, all boundaries periodic for now.
        EvenOddTMCloverOperator( SU3MatrixBlock* u_[2],
            FullCloverBlock* clov_,
            FullCloverBlock* invclov_,
            Geometry<FT,veclen,soalen,compress12>* geom_,
            double t_boundary,
            double aniso_coeff_s,
            double aniso_coeff_t )
          : D( new TMClovDslash<FT,veclen,soalen,compress12>(geom_, t_boundary, aniso_coeff_s, aniso_coeff_t) )
        {
          Geometry<FT,veclen,soalen,compress12>& geom = D->getGeometry();
          u[0] = u_[0];
          u[1] = u_[1];
          clov = clov_;
          invclov = invclov_;
          tmp = (FourSpinorBlock*) geom.allocCBFourSpinor();
        }

        // Destructor
        ~EvenOddTMCloverOperator()
        {
          Geometry<FT,veclen,soalen,compress12>& geom = D->getGeometry();
          geom.free(tmp);
          delete D;
        }

        EvenOddTMCloverOperator( Geometry<FT,veclen,soalen,compress12>* geom_,
            double t_boundary, double aniso_coeff_s, double aniso_coeff_t )
          : D( new TMClovDslash<FT,veclen,soalen,compress12>(geom_, t_boundary, aniso_coeff_s, aniso_coeff_t) )
        {
          Geometry<FT,veclen,soalen,compress12>& geom = D->getGeometry();
          tmp = (FourSpinorBlock *)geom.allocCBFourSpinor();
        }

        void setFields(SU3MatrixBlock* u_[2], FullCloverBlock* clov_, FullCloverBlock* invclov_)
        {
          u[0] = u_[0];
          u[1] = u_[1];
          clov = clov_;
          invclov = invclov_;
        }

        void operator()(FourSpinorBlock *res, const FourSpinorBlock* in, int isign)
        {
          double beta=(double)0.25;

          D->dslash(tmp, in, u[0], invclov,  isign, 0);
          D->dslashAChiMinusBDPsi(res, tmp, in, u[1], clov, beta, isign, 1);
        }

        Geometry<FT,veclen,soalen,compress12>& getGeometry() { return D->getGeometry(); }


      private:

        double Mass;
        TMClovDslash<FT,veclen,soalen,compress12>* D;
        mutable SU3MatrixBlock *u[2]; // Mutable because of setFields
        mutable FullCloverBlock* clov;    // Mutable because of setFields
        mutable FullCloverBlock* invclov; // Mutable because of setFields
        FourSpinorBlock *tmp;

    }; // Class
}; // Namespace

#endif // Include guard
