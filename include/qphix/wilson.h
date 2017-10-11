#pragma once

#include "qphix/linearOp.h"
#include "qphix/dslash_def.h"
#include "qphix/blas.h"

namespace QPhiX
{

template <typename FT, int veclen, int soalen, bool compress12>
class EvenOddWilsonOperator
    : public EvenOddLinearOperator<FT, veclen, soalen, compress12>
{
 public:
  typedef typename Geometry<FT, veclen, soalen, compress12>::FourSpinorBlock
      FourSpinorBlock;
  typedef typename Geometry<FT, veclen, soalen, compress12>::SU3MatrixBlock
      SU3MatrixBlock;

  EvenOddWilsonOperator(const double Mass_,
                        SU3MatrixBlock *u_[2],
                        Geometry<FT, veclen, soalen, compress12> *geom_,
                        double t_boundary_,
                        double aniso_fac_s_,
                        double aniso_fac_t_,
                        bool use_tbc_[4] = nullptr,
                        double tbc_phases_[4][2] = nullptr)
      : Mass(Mass_),
        D(new Dslash<FT, veclen, soalen, compress12>(
            geom_, t_boundary_, aniso_fac_s_, aniso_fac_t_, use_tbc_, tbc_phases_))
  {

    Geometry<FT, veclen, soalen, compress12> &geom = D->getGeometry();
    u[0] = u_[0];
    u[1] = u_[1];
    tmp = (FourSpinorBlock *)geom.allocCBFourSpinor();

    mass_factor_alpha = (double)4 + Mass;
    mass_factor_beta = (double)(0.25) / mass_factor_alpha;
  }

  ~EvenOddWilsonOperator()
  {
    Geometry<FT, veclen, soalen, compress12> &geom = D->getGeometry();
    geom.free(tmp);
    delete D;
  }

  inline void operator()(FourSpinorBlock *res,
                         const FourSpinorBlock *in,
                         int isign,
                         int target_cb = 1) const override
  {

    int other_cb = 1 - target_cb;
    D->dslash(tmp, in, u[other_cb], isign, other_cb);
    D->dslashAChiMinusBDPsi(res,
                            tmp,
                            in,
                            u[target_cb],
                            mass_factor_alpha,
                            mass_factor_beta,
                            isign,
                            target_cb);
  }

  // Offdiag is just the dslash
  inline void M_offdiag(FourSpinorBlock *res,
                          FourSpinorBlock const *in,
                          int isign,
                          int target_cb) const override {
    Geometry<FT, veclen, soalen, compress12> &geom = D->getGeometry();
    double beta = -0.5;
    D->dslash(tmp,in, u[target_cb],isign,target_cb);
    axy(beta,tmp,res,geom,geom.getNSIMT());
  }

   // EE-inv is just the identity
   inline void M_diag_inv(FourSpinorBlock *res,
                         FourSpinorBlock const *in,
                         int isign) const override {
     double beta = static_cast<double>(1)/mass_factor_alpha;
     Geometry<FT, veclen, soalen, compress12> &geom = D->getGeometry();
     axy(beta,in,res, geom, geom.getNSIMT());
   }

  Geometry<FT, veclen, soalen, compress12> &getGeometry() override
  {
    return D->getGeometry();
  }

 private:
  double Mass;
  Dslash<FT, veclen, soalen, compress12> *D;
  const SU3MatrixBlock *u[2];
  FourSpinorBlock *tmp;
  double mass_factor_alpha;
  double mass_factor_beta;
}; // Class
}; // Namespace
