#pragma once

#include <qphix/geometry.h>

namespace QPhiX
{

template <typename FT, int veclen, int soalen, bool compress>
class EvenOddLinearOperator
{
 public:
  virtual ~EvenOddLinearOperator(){};

  typedef typename Geometry<FT, veclen, soalen, compress>::FourSpinorBlock
      FourSpinorBlock;

  static constexpr int num_flav = 1;

  // This is the full Schur Operator also known as M_oo
  virtual void operator()(FourSpinorBlock *res,
                          FourSpinorBlock const *in,
                          int isign,
                          int target_cb = 1) const = 0;


  virtual void M_offdiag(FourSpinorBlock *res,
                         FourSpinorBlock const *in,
                         int isign,
                         int target_cb) const = 0;

  // M_ee_inv is always Hermitian so no need for isign?
  // for wilson it is the identity, for clover it is hermitian, for TWM it is gamma_5?
  virtual void M_diag_inv(FourSpinorBlock *res,
                        FourSpinorBlock const *in,
                        int isign) const = 0;

#ifdef __INTEL_COMPILER
  virtual void operator()(FourSpinorBlock *const res[1],
                          FourSpinorBlock *const in[1],
                          int isign,
                          int target_cb = 1) const
  {
    (*this)(res[0], in[0], isign, target_cb);
  }

  virtual void M_offdiag(FourSpinorBlock *const res[1],
                         FourSpinorBlock *const in[1],
                         int isign,
                         int target_cb) const
  {
    (*this).M_offdiag(res[0], in[0], isign, target_cb);
  }

   // M_ee_inv is always Hermitian so no need for isign?
   // for wilson it is the identity, for clover it is hermitian, for TWM it is gamma_5?
   virtual void M_diag_inv(FourSpinorBlock *const res[1],
                           FourSpinorBlock *const in[1],
                           int isign) const
   {
     (*this).M_diag_inv(res[0],in[0],isign);
   }

#endif

  virtual void operator()(FourSpinorBlock *const res[1],
                          FourSpinorBlock const *const in[1],
                          int isign,
                          int target_cb = 1) const
  {
    (*this)(res[0], in[0], isign, target_cb);
  };

   virtual void M_offdiag(FourSpinorBlock *const res[1],
                          FourSpinorBlock const *const in[1],
                          int isign,
                          int target_cb) const
  {
    this->M_offdiag(res[0], in[0], isign, target_cb);
  }

   // M_ee_inv is always Hermitian so no need for isign?
   // for wilson it is the identity, for clover it is hermitian, for TWM it is gamma_5?
   virtual void M_diag_inv(FourSpinorBlock *const res[1],
                         FourSpinorBlock const *const in[1],
                         int isign) const
   {
     this->M_diag_inv(res[0],in[0],isign);
   }


  virtual Geometry<FT, veclen, soalen, compress> &getGeometry(void) = 0;
};

template <typename FT, int veclen, int soalen, bool compress>
class TwoFlavEvenOddLinearOperator
{
 public:
  virtual ~TwoFlavEvenOddLinearOperator(){};

  typedef typename Geometry<FT, veclen, soalen, compress>::FourSpinorBlock
      FourSpinorBlock;

  static constexpr int num_flav = 2;

#ifdef __INTEL_COMPILER
  virtual void
  operator()(FourSpinorBlock * const res[2], FourSpinorBlock *const in[2], int isign, int target_cb = 1)
  {
    (*this)(res, const_cast<FourSpinorBlock const *const *>(in), isign, target_cb);
  }

  virtual void M_offdiag(FourSpinorBlock *res[2],
                         FourSpinorBlock *const in[2],
                            int isign,
                            int target_cb) const
  {
    this->M_offdiag(res, const_cast<FourSpinorBlock const *const *>(in), isign);
  }

  // M_ee_inv is always Hermitian so no need for isign?
  // for wilson it is the identity, for clover it is hermitian, for TWM it is gamma_5?
  virtual void M_diag_inv(FourSpinorBlock *res[2],
                          FourSpinorBlock *const in[2],
                          int isign) const
  {
    this->M_diag_inv(res,const_cast<FourSpinorBlock const *const *>(in),isign);
  }
#endif

  virtual void operator()(FourSpinorBlock * const res[2],
                          FourSpinorBlock const *const in[2],
                          int isign,
                          int target_cb = 1) = 0;

  virtual void M_offdiag(FourSpinorBlock *res[2],
                           const FourSpinorBlock *const in[2],
                           int isign,
                           int target_cb) const = 0;

    // M_ee_inv is always Hermitian so no need for isign?
    // for wilson it is the identity, for clover it is hermitian, for TWM it is gamma_5?
    virtual void M_diag_inv(FourSpinorBlock *res[2],
                          const FourSpinorBlock *const in[2],
                          int isign) const  = 0;

  virtual Geometry<FT, veclen, soalen, compress> &getGeometry(void) = 0;
};
}

