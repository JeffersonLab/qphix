#ifndef QPHIX_SSE_UTILS_H
#define QPHIX_SSE_UTILS_H

#include <smmintrin.h>

inline __m128d _mm_int2mask_pd(unsigned int msk) {
	static __m128d allOne = _mm_set1_pd(-1.0);
	__m128d ret = _mm_setzero_pd();
	msk = msk & 0x03;

	switch (msk) {
		case 0: ret = _mm_blend_pd(ret, allOne, 0); break;
		case 1: ret = _mm_blend_pd(ret, allOne, 1); break;
		case 2: ret = _mm_blend_pd(ret, allOne, 2); break;
		case 3: ret = _mm_blend_pd(ret, allOne, 3); break;

	}
	return ret;
}


inline __m128 _mm_int2mask_ps(unsigned int msk) {
	static __m128 allOne = _mm_set1_ps(-1.0);
	__m128 ret = _mm_setzero_ps();
	msk = msk & 0x0F;
	switch (msk) {
                case 0: ret = _mm_blend_ps(ret, allOne, 0); break;
                case 1: ret = _mm_blend_ps(ret, allOne, 1); break;
                case 2: ret = _mm_blend_ps(ret, allOne, 2); break;
                case 3: ret = _mm_blend_ps(ret, allOne, 3); break;
                case 4: ret = _mm_blend_ps(ret, allOne, 4); break;
                case 5: ret = _mm_blend_ps(ret, allOne, 5); break;
                case 6: ret = _mm_blend_ps(ret, allOne, 6); break;
                case 7: ret = _mm_blend_ps(ret, allOne, 7); break;
                case 8: ret = _mm_blend_ps(ret, allOne, 8); break;
                case 9: ret = _mm_blend_ps(ret, allOne, 9); break;
                case 10: ret = _mm_blend_ps(ret, allOne, 10); break;
                case 11: ret = _mm_blend_ps(ret, allOne, 11); break;
                case 12: ret = _mm_blend_ps(ret, allOne, 12); break;
                case 13: ret = _mm_blend_ps(ret, allOne, 13); break;
                case 14: ret = _mm_blend_ps(ret, allOne, 14); break;
                case 15: ret = _mm_blend_ps(ret, allOne, 15); break;
	}
	return ret;
}

#endif
