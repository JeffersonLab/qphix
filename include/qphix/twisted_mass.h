#ifndef QPHIX_TWISTED_MASS_H
#define QPHIX_TWISTED_MASS_H

#include "qphix/linearOp.h"
#include "qphix/tm_dslash_def.h"

namespace QPhiX
{

  template<typename FT, int veclen, int soalen, bool compress12>
    class EvenOddTMWilsonOperator : public EvenOddLinearOperator<FT, veclen, soalen, compress12> {

      public:

        typedef typename Geometry<FT,veclen,soalen,compress12>::FourSpinorBlock FourSpinorBlock;
        typedef typename Geometry<FT,veclen,soalen,compress12>::SU3MatrixBlock  SU3MatrixBlock;


        EvenOddTMWilsonOperator(const double Mass_,
            const double TwistedMass_,
            SU3MatrixBlock* u_[2],
            Geometry<FT,veclen,soalen,compress12>* geom_,
            double t_boundary_,
            double aniso_fac_s_,
            double aniso_fac_t_)
          :
          Mass(Mass_),
          TwistedMass(TwistedMass_),
          mass_factor_alpha(4.0+Mass),
          mass_factor_beta(0.25),
          D(
              new TMDslash<FT, veclen, soalen, compress12>(
                geom_,
                t_boundary_,
                aniso_fac_s_,
                aniso_fac_t_,
                Mass,
                TwistedMass
              )
          )
        {
          u[0] = u_[0];
          u[1] = u_[1];
          Geometry<FT, veclen, soalen, compress12>& geom = D->getGeometry();
          tmp = (FourSpinorBlock *)geom.allocCBFourSpinor();
        }

        ~EvenOddTMWilsonOperator()
        {
          Geometry<FT, veclen, soalen, compress12>& geom = D->getGeometry();
          geom.free(tmp);
          delete D;
        }

        Geometry<FT,veclen, soalen, compress12>& getGeometry() { return D->getGeometry(); }

        // TODO: Add detailed description of the operator!
        // This is the even-odd preconditioned (odd-odd) fermion matrix
        void operator()(
            FourSpinorBlock *res,      // result spinor field
            const FourSpinorBlock* in, // input spinor field
            int isign                  // non-conjugate = 1, hermitian conjugate = -1
        )
        {
          D->dslash(tmp, in, u[1], isign, 1);
          D->dslashAChiMinusBDPsi(res, tmp, in, u[0], mass_factor_alpha, mass_factor_beta, isign, 0);
        }


      private:

        double Mass;
        double TwistedMass;

        double mass_factor_alpha;
        double mass_factor_beta;

        const SU3MatrixBlock *u[2];
        FourSpinorBlock *tmp;
        TMDslash<FT, veclen, soalen, compress12>* D;

    }; // Class
}; // Namespace

#endif
