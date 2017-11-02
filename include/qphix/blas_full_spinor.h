/*
 * blas_full_spinor.h
 *
 *  Created on: Oct 16, 2017
 *      Author: bjoo
 */

#ifndef INCLUDE_QPHIX_BLAS_FULL_SPINOR_H_
#define INCLUDE_QPHIX_BLAS_FULL_SPINOR_H_

#include "qphix/full_spinor.h"
#include "qphix/blas_new_c.h"
namespace QPhiX {

template<typename FT, int V, int S, bool compress>
void copySpinor(FullSpinor<FT,V,S,compress>& result,
                const FullSpinor<FT,V,S,compress>& src,
                const Geometry<FT,V,S,compress>& geom,
                int n_blas_simt)
{
  for(int cb=0; cb < 2; ++cb) {
    copySpinor<FT,V,S,compress>(result.getCBData(cb),src.getCBData(cb),geom,n_blas_simt);
  }
}

template<typename FT, int V, int S, bool compress>
void zeroSpinor(FullSpinor<FT,V,S,compress>& result,
                const Geometry<FT,V,S,compress>& geom,
                int n_blas_simt)
{
  for(int cb=0; cb < 2; ++cb) {
    zeroSpinor<FT,V,S,compress>(result.getCBData(cb),geom,n_blas_simt);
  }

}
// Complex AXPY
template<typename FT, int V, int S, bool compress>
void caxpySpinor(double alpha[2],
    const FullSpinor<FT,V,S,compress>& x,
    FullSpinor<FT,V,S,compress>& y,
    const Geometry<FT,V,S,compress>& geom,
    int n_blas_simt)
{
  for(int cb=0; cb < 2; ++cb) {
    caxpySpinor<FT,V,S,compress>(alpha, x.getCBData(cb), y.getCBData(cb), geom, n_blas_simt);
  }
}

// Real AXPY
template<typename FT, int V, int S, bool compress>
void axpySpinor(double alpha,
    const FullSpinor<FT,V,S,compress>& x,
    FullSpinor<FT,V,S,compress>& y,
    const Geometry<FT,V,S,compress>& geom,
    int n_blas_simt)
{
  for(int cb=0; cb < 2; ++cb) {
    axpy<FT,V,S,compress>(alpha, x.getCBData(cb), y.getCBData(cb), geom, n_blas_simt);
  }
}

// XmyNorm2
template<typename FT, int V, int S, bool compress>
void xmy2Norm2Spinor( const FullSpinor<FT,V,S,compress>& x,
                     FullSpinor<FT,V,S,compress>& y,
                     double &r_norm,
                     const Geometry<FT,V,S,compress>& geom,
                     int n_blas_simt)
{
    double res_cb[2] = { 0,0 };
    for(int cb=0; cb < 2; ++cb) {
      xmy2Norm2Spinor<FT,V,S,compress>(x.getCBData(cb),y.getCBData(cb), res_cb[cb],geom,n_blas_simt,
          false);
    }

    r_norm = res_cb[0] + res_cb[1];
    CommsUtils::sumDouble(&r_norm);

}


// InnerProduct
template <typename FT, int V, int S, bool compress>
void innerProductSpinor(double results[2],
                  const FullSpinor<FT,V,S,compress>& x,
                  const FullSpinor<FT,V,S,compress>& y,
                  const Geometry<FT, V, S, compress> &geom,
                  int n_blas_simt)
{
  results[0] = (double) 0;
  results[1] = (double) 0;

  double cb_results_0[2] = {0,0};
  double cb_results_1[2] = {0,0};

  innerProduct<FT,V,S,compress>(cb_results_0,
               x.getCBData(0),
               y.getCBData(0),
               geom,
               n_blas_simt,
               false);

  innerProduct<FT,V,S,compress>(cb_results_1,
                 x.getCBData(1),
                 y.getCBData(1),
                 geom,
                 n_blas_simt,
                 false);

  results[0] = cb_results_0[0]+cb_results_1[0];
  results[1] = cb_results_0[1]+cb_results_1[1];

  CommsUtils::sumDoubleArray(results,2);

}


// Norm2
template<typename FT, int V, int S, bool compress>
void norm2Spinor(double &n2,
                 const FullSpinor<FT,V,S,compress>& x,
                 const Geometry<FT, V, S, compress> &geom,
                 int n_blas_simt)
{
  double norm_cb[2] = {0,0};
  for(int cb=0; cb < 2; ++cb )  {
    norm2Spinor<FT,V,S,compress>(norm_cb[cb], x.getCBData(cb),geom,n_blas_simt,false );
  }
  n2 = norm_cb[0] + norm_cb[1];
  CommsUtils::sumDouble(&n2);
}


}



#endif /* INCLUDE_QPHIX_BLAS_FULL_SPINOR_H_ */
