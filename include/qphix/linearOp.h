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

  virtual void operator()(FourSpinorBlock *res,
                          FourSpinorBlock const *in,
                          int isign,
                          int target_cb = 1) const = 0;

#ifdef __INTEL_COMPILER
  virtual void operator()(FourSpinorBlock *const res[1],
                          FourSpinorBlock *const in[1],
                          int isign,
                          int target_cb = 1) const
  {
    (*this)(res[0], in[0], isign, target_cb);
  };
#endif

  virtual void operator()(FourSpinorBlock *const res[1],
                          FourSpinorBlock const *const in[1],
                          int isign,
                          int target_cb = 1) const
  {
    (*this)(res[0], in[0], isign, target_cb);
  };

  virtual Geometry<FT, veclen, soalen, compress> &getGeometry(void) = 0;
};

template <typename FT, int veclen, int soalen, bool compress>
class TwoFlavEvenOddLinearOperator
{
 public:
  typedef typename Geometry<FT, veclen, soalen, compress>::FourSpinorBlock
      FourSpinorBlock;

  static constexpr int num_flav = 2;

#ifdef __INTEL_COMPILER
  virtual void
  operator()(FourSpinorBlock *res[2], FourSpinorBlock *const in[2], int isign, int target_cb = 1) const
  {
    (*this)(res, const_cast<FourSpinorBlock const *const *>(in), isign, target_cb);
  }
#endif

  virtual void operator()(FourSpinorBlock *res[2],
                          const FourSpinorBlock *const in[2],
                          int isign,
                          int target_cb = 1) const = 0;

  virtual Geometry<FT, veclen, soalen, compress> &getGeometry(void) = 0;
};
}

