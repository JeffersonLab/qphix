#ifndef LINEAR_OPERATOR_H
#define LINEAR_OPERATOR_H

#include "cpp_dslash.h"


namespace CPlusPlusWilsonDslash { 

  template<typename FT, int veclen, int soalen, bool compress> 
  class EvenOddLinearOperator {
  public:
    typedef typename Geometry<FT,veclen,soalen,compress>::FourSpinorBlock FourSpinorBlock;
    virtual void operator()(FourSpinorBlock *res, const FourSpinorBlock* in, int isign) = 0;

    virtual Geometry<FT,veclen,soalen,compress>& getGeometry(void)=0;

  };






}



#endif
