#pragma once

#include "qphix/linearOp.h"
#include "qphix/dslash_def.h"
#include "qphix/clover_dslash_def.h"
#include "qphix/clover_product.h"
#include "qphix/blas_new_c.h"
#include "qphix/full_spinor.h"
#include <cassert>

namespace QPhiX
{

template <typename FT, int veclen, int soalen, bool compress12>
class EvenOddCloverOperator
    : public EvenOddLinearOperator<FT, veclen, soalen, compress12>
{
public:
  typedef typename Geometry<FT, veclen, soalen, compress12>::CloverBlock CloverBlock;
  typedef typename Geometry<FT, veclen, soalen, compress12>::FourSpinorBlock
      FourSpinorBlock;
  typedef typename Geometry<FT, veclen, soalen, compress12>::SU3MatrixBlock
      SU3MatrixBlock;


  // Constructor
  // No anisotropy, all boundaries periodic for now.
  // User passes only a clov_oo, and invclov_ee
  // Sufficient for use of props...

  EvenOddCloverOperator(SU3MatrixBlock *u_[2],
      CloverBlock *clov_,
      CloverBlock *invclov_,
      Geometry<FT, veclen, soalen, compress12> *geom_,
      double t_boundary,
      double aniso_coeff_s,
      double aniso_coeff_t,
      bool use_tbc_[4] = nullptr,
      double tbc_phases_[4][2] = nullptr,
      double const prec_mass_rho = 0.0)
  : D(new ClovDslash<FT, veclen, soalen, compress12>(geom_,
      t_boundary,
      aniso_coeff_s,
      aniso_coeff_t,
      use_tbc_,
      tbc_phases_,
      prec_mass_rho)),
      D_dslash(new Dslash<FT, veclen, soalen, compress12>(geom_,
          t_boundary,
          aniso_coeff_s,
          aniso_coeff_t,
          use_tbc_,
          tbc_phases_))
  {
    Geometry<FT, veclen, soalen, compress12> &geom = D->getGeometry();
    u[0] = u_[0];
    u[1] = u_[1];
    clov[0] = nullptr;
    clov[1] = clov_;
    invclov = invclov_;
    tmp = (FourSpinorBlock *)geom.allocCBFourSpinor();
  }

  // Constructor
  // No anisotropy, all boundaries periodic for now.
  // BOTH CLOVER PArities are passed
  EvenOddCloverOperator(SU3MatrixBlock *u_[2],
      CloverBlock *clov_[2],
      CloverBlock *invclov_,
      Geometry<FT, veclen, soalen, compress12> *geom_,
      double t_boundary,
      double aniso_coeff_s,
      double aniso_coeff_t,
      bool use_tbc_[4] = nullptr,
      double tbc_phases_[4][2] = nullptr,
      double const prec_mass_rho = 0.0)
  : D(new ClovDslash<FT, veclen, soalen, compress12>(geom_,
      t_boundary,
      aniso_coeff_s,
      aniso_coeff_t,
      use_tbc_,
      tbc_phases_,
      prec_mass_rho)),
      D_dslash(new Dslash<FT, veclen, soalen, compress12>(geom_,
          t_boundary,
          aniso_coeff_s,
          aniso_coeff_t,
          use_tbc_,
          tbc_phases_))
  {
    Geometry<FT, veclen, soalen, compress12> &geom = D->getGeometry();
    u[0] = u_[0];
    u[1] = u_[1];
    clov[0] = clov_[0];
    clov[1] = clov_[1];
    invclov = invclov_;
    tmp = (FourSpinorBlock *)geom.allocCBFourSpinor();
  }

  // FIelds are not set, use set fields later to set them.
  EvenOddCloverOperator(Geometry<FT, veclen, soalen, compress12> *geom_,
      double t_boundary,
      double aniso_coeff_s,
      double aniso_coeff_t,
      bool use_tbc_[4] = nullptr,
      double tbc_phases_[4][2] = nullptr,
      double const prec_mass_rho = 0.0)
  : D(new ClovDslash<FT, veclen, soalen, compress12>(geom_,
      t_boundary,
      aniso_coeff_s,
      aniso_coeff_t,
      use_tbc_,
      tbc_phases_,
      prec_mass_rho)),
      D_dslash(new Dslash<FT, veclen, soalen, compress12>(geom_,
          t_boundary,
          aniso_coeff_s,
          aniso_coeff_t,
          use_tbc_,
          tbc_phases_)),
          u({nullptr,nullptr}), clov({nullptr,nullptr}),invclov(nullptr)
  {
    Geometry<FT, veclen, soalen, compress12> &geom = D->getGeometry();
    tmp = (FourSpinorBlock *)geom.allocCBFourSpinor();
  }

  // Single Parity of clover is passed
  void setFields(SU3MatrixBlock *u_[2], CloverBlock *clov_, CloverBlock *invclov_)
  {
    u[0] = u_[0];
    u[1] = u_[1];
    clov[0] = nullptr;
    clov[1] = clov_;
    invclov = invclov_;
  }

  // Both Parities of clover are passed
  void setFields(SU3MatrixBlock *u_[2], CloverBlock *clov_[2], CloverBlock *invclov_)
   {
     u[0] = u_[0];
     u[1] = u_[1];

     clov[0] = clov_[0];
     clov[1] = clov_[1];
     invclov = invclov_;
   }

  ~EvenOddCloverOperator()
  {
    Geometry<FT, veclen, soalen, compress12> &geom = D->getGeometry();
    geom.free(tmp);
    delete D;
    delete D_dslash;
  }

  inline void operator()(FourSpinorBlock *res,
      const FourSpinorBlock *in,
      int isign,
      int target_cb = 1) const override
      {
    double beta = 0.25;
    int other_cb = 1 - target_cb;
    D->dslash(tmp, in, u[other_cb], invclov, isign, other_cb);
    D->dslashAChiMinusBDPsi(
        res, tmp, in, u[target_cb], clov[1], beta, isign, target_cb);
      }

  // Offdiag is just the dslash
  inline void M_offdiag(FourSpinorBlock *res,
      FourSpinorBlock const *in,
      int isign,
      int target_cb) const override {
    Geometry<FT, veclen, soalen, compress12> &geom = D->getGeometry();
    double beta = -0.5;
    D_dslash->dslash(tmp,in, u[target_cb],isign,target_cb);
    axy(beta,tmp,res,geom,geom.getNSIMT());
  }

  // EE-inv is just the identity
  inline void M_diag_inv(FourSpinorBlock *res,
      FourSpinorBlock const *in,
      int isign) const override {
    Geometry<FT, veclen, soalen, compress12> &geom = D->getGeometry();
    clover_product(res,in,invclov,geom);
  }

  using SpinorFull = FullSpinor<FT,veclen,soalen,compress12>;

  // Extra operations: Clover Term
  inline void M_diag(SpinorFull& res,
                     const SpinorFull& in,
                     int isign) const

  {
    Geometry<FT, veclen, soalen, compress12> &geom = D->getGeometry();
    assert(clov[0]!=nullptr);
    assert(clov[1]!=nullptr);
    for(int cb=0; cb < 2;++cb)
      clover_product( res.getCBData(cb), in.getCBData(cb), clov[cb], geom);
  }
  // Extra operations: Clover Term
   inline void M_diag(SpinorFull& res,
                      const SpinorFull& in,
                      int isign,int cb) const

   {
     Geometry<FT, veclen, soalen, compress12> &geom = D->getGeometry();
     assert(clov[cb]!=nullptr);
     clover_product( res.getCBData(cb), in.getCBData(cb), clov[cb], geom);
   }

  inline void M_unprec(SpinorFull& res,
                       const SpinorFull& in,
                       int isign) const
  {
    assert(clov[0] != nullptr);
    assert(clov[1] != nullptr);

    double beta = 0.5;
    for(int cb=0; cb < 2; ++cb) {
    D->dslashAChiMinusBDPsi(
            res.getCBData(cb),
            in.getCBData(1-cb),
            in.getCBData(cb),
            u[cb],
            clov[cb],
            beta,
            isign,
            cb);
  } //cb
  }

  inline void DslashDir(SpinorFull& res,
                        const SpinorFull& in,
                        int dir)
     {
       for(int cb=0; cb < 2; ++cb) {
         D_dslash->dslashDir(
               res.getCBData(cb),
               in.getCBData(1-cb),
               u[cb],
               cb,
               dir);
     } //cb
  }

  Geometry<FT, veclen, soalen, compress12> &getGeometry() override
  {
    return D->getGeometry();
  }

private:
  //double Mass;
  ClovDslash<FT, veclen, soalen, compress12> *D;
  Dslash<FT, veclen, soalen, compress12> *D_dslash;
  SU3MatrixBlock *u[2];
  CloverBlock *clov[2];
  CloverBlock *invclov;
  FourSpinorBlock *tmp;
}; // Class
}; // Namespace
