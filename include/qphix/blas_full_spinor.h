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
		int n_blas_simt,
		const int min_cb=0,
		const int max_cb=2)
{
	for(int cb=min_cb; cb < max_cb; ++cb) {
		copySpinor<FT,V,S,compress>(result.getCBData(cb),src.getCBData(cb),geom,n_blas_simt);
	}
}

template<typename FT, int V, int S, bool compress>
void zeroSpinor(FullSpinor<FT,V,S,compress>& result,
		const Geometry<FT,V,S,compress>& geom,
		int n_blas_simt,
		const int min_cb=0,
		const int max_cb=2)
{
	for(int cb=min_cb; cb < max_cb; ++cb) {
		zeroSpinor<FT,V,S,compress>(result.getCBData(cb),geom,n_blas_simt);
	}

}

template<typename FT, int V, int S, bool compress>
void axSpinor(const double alpha,
		FullSpinor<FT,V,S,compress>& x,
		const Geometry<FT,V,S,compress>& geom,
		int n_blas_simt,
		const int min_cb=0,
		const int max_cb=2)
{
	for(int cb=min_cb; cb < max_cb; ++cb) {
		ax<FT,V,S,compress>(alpha, x.getCBData(cb),geom,n_blas_simt);
	}

}
// Complex AXPY
template<typename FT, int V, int S, bool compress>
void caxpySpinor(double alpha[2],
		const FullSpinor<FT,V,S,compress>& x,
		FullSpinor<FT,V,S,compress>& y,
		const Geometry<FT,V,S,compress>& geom,
		int n_blas_simt,
		const int min_cb=0,
		const int max_cb=2)
{
	for(int cb=min_cb; cb < max_cb; ++cb) {
		caxpySpinor<FT,V,S,compress>(alpha, x.getCBData(cb), y.getCBData(cb), geom, n_blas_simt);
	}
}

// Real AXPY
template<typename FT, int V, int S, bool compress>
void axpySpinor(double alpha,
		const FullSpinor<FT,V,S,compress>& x,
		FullSpinor<FT,V,S,compress>& y,
		const Geometry<FT,V,S,compress>& geom,
		int n_blas_simt,
		const int min_cb=0,
		const int max_cb=2)
{
	for(int cb=min_cb; cb < max_cb; ++cb) {
		axpy<FT,V,S,compress>(alpha, x.getCBData(cb), y.getCBData(cb), geom, n_blas_simt);
	}
}

// Real AXPY
template<typename FT, int V, int S, bool compress>
void ypeqxSpinor(const FullSpinor<FT,V,S,compress>& x,
		FullSpinor<FT,V,S,compress>& y,
		const Geometry<FT,V,S,compress>& geom,
		int n_blas_simt,
		const int min_cb=0,
		const int max_cb=2)
{
	for(int cb=min_cb; cb < max_cb; ++cb) {
		ypeqx<FT,V,S,compress>(x.getCBData(cb), y.getCBData(cb), geom, n_blas_simt);
	}
}

// Real AXPY
template<typename FT, int V, int S, bool compress>
void ymeqxSpinor(const FullSpinor<FT,V,S,compress>& x,
		FullSpinor<FT,V,S,compress>& y,
		const Geometry<FT,V,S,compress>& geom,
		int n_blas_simt,
		const int min_cb=0,
		const int max_cb=2)
{
	for(int cb=min_cb; cb < max_cb; ++cb) {
		ymeqx<FT,V,S,compress>(x.getCBData(cb), y.getCBData(cb), geom, n_blas_simt);
	}
}

// XmyNorm2
template<typename FT, int V, int S, bool compress>
void xmy2Norm2Spinor( const FullSpinor<FT,V,S,compress>& x,
		FullSpinor<FT,V,S,compress>& y,
		double &r_norm,
		const Geometry<FT,V,S,compress>& geom,
		int n_blas_simt,
		const int min_cb=0,
		const int max_cb=2)
{
	r_norm = (double)0;
	for(int cb=min_cb; cb < max_cb; ++cb) {
		double res_cb=(double)0;

		xmy2Norm2Spinor<FT,V,S,compress>(x.getCBData(cb),y.getCBData(cb), res_cb,geom,n_blas_simt,
				false);

		r_norm += res_cb;
	}
	CommsUtils::sumDouble(&r_norm);

}


// InnerProduct
template <typename FT, int V, int S, bool compress>
void innerProductSpinor(double results[2],
		const FullSpinor<FT,V,S,compress>& x,
		const FullSpinor<FT,V,S,compress>& y,
		const Geometry<FT, V, S, compress> &geom,
		int n_blas_simt,
		const int min_cb=0,
		const int max_cb=2)
{
	results[0] = (double) 0;
	results[1] = (double) 0;


	for(int cb=min_cb; cb < max_cb; ++cb) {
		double cb_results[2] = {0,0};

		innerProduct<FT,V,S,compress>(cb_results,
				x.getCBData(cb),
				y.getCBData(cb),
				geom,
				n_blas_simt,
				false);

		results[0] += cb_results[0];
		results[1] += cb_results[1];

	}

	CommsUtils::sumDoubleArray(results,2);

}


// Norm2
template<typename FT, int V, int S, bool compress>
void norm2Spinor(double &n2,
		const FullSpinor<FT,V,S,compress>& x,
		const Geometry<FT, V, S, compress> &geom,
		int n_blas_simt,
		const int min_cb=0,
		const int max_cb=2)
{
	n2 = (double)0;
	for(int cb=min_cb; cb < max_cb; ++cb )  {
		double norm_cb = (double)0;
		norm2Spinor<FT,V,S,compress>(norm_cb, x.getCBData(cb),geom,n_blas_simt,false );
		n2 += norm_cb;
	}
	CommsUtils::sumDouble(&n2);
}


}



#endif /* INCLUDE_QPHIX_BLAS_FULL_SPINOR_H_ */
