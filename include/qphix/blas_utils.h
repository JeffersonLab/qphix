#ifndef QPHIX_BLAS_UTILS_H
#define QPHIX_BLAS_UTILS_H

namespace QPhiX {

	namespace BLASUtils { 
    
		//  res = alpha x + y
		//  alpha is complex
		template<typename FT, int S>
		inline void 
			cmadd(FT res[2][S], FT alpha[2], FT x[2][S], FT y[2][S]) 
		{
			//  (a[RE] x[RE] - a[IM] y[IM])  + res[RE]
			//  (a[RE] y[IM] + a[IM] y[RE])  + res[IM]
#pragma omp simd
			for(int s=0; s < S; s++) { 
				res[0][s] = alpha[0]*x[0][s] - alpha[1]*x[1][s]  + y[0][s];
				res[1][s] = alpha[0]*x[1][s] + alpha[1]*x[0][s]  + y[1][s];
			}
		}
    
		// res = -alpha x + y
		// res = y - alpha x
		// alpha is complex
		template<typename FT, int S>
		inline void 
			cnmadd(FT res[2][S], FT alpha[2], FT x[2][S], FT y[2][S]) 
		{
			//  res[RE] -(a[RE] x[RE] - a[IM] y[IM])
			// =res[RE] - a[RE] x[RE] + a[IM] y[IM]
      
			//  res[IM] -(a[RE] y[IM] + a[IM] y[RE])
			// =res[IM] -a[RE]y[IM] - a[IM] y[RE]
#pragma omp simd 
			for(int s=0; s < S; s++) { 
				res[0][s] = y[0][s] - alpha[0]*x[0][s] + alpha[1]*x[1][s];
				res[1][s] = y[1][s] - alpha[0]*x[1][s] - alpha[1]*x[0][s];
			}
		}
    
   
		// Generic stream in
		template<typename FT, int V>
		inline void
		streamInSpinor(volatile FT* restrict dst, const volatile FT* restrict src, int numvec) { 
      
#if defined(__MIC__)
			//Intel MIC
			const int prefdist1 = 12;
			const int prefdist2 = 64;
      
			const char* prefl1base = (const char *)src+prefdist1*64;
			const char* prefl2base = (const char *)src+prefdist2*64;
			
			for(int v=0; v < numvec; v++) { 
				_mm_prefetch(&prefl1base[v*V*sizeof(FT)], _MM_HINT_T0);
				_mm_prefetch(&prefl2base[v*V*sizeof(FT)], _MM_HINT_T1);
#pragma omp simd
				for(int s=0; s < V; s++) {
					dst[v*V+s]=src[v*V+s];
				}
			}
#else
			//Generic
#pragma omp simd collapse(2) aligned(dst,src:V)
			for(int v=0; v < numvec; v++) {
				for(int s=0; s < V; s++) {
					dst[v*V+s]=src[v*V+s];
				}
			}
#endif
		}
    
		// Generic write  out
		template<typename FT, int V>
		inline void
		writeSpinor(volatile FT* restrict dst, const volatile FT* restrict src, int numvec) { 
      
			//#pragma vector aligned(src)
			//#pragma vector aligned(dst)
#pragma omp simd collapse(2) aligned(dst,src:V) //by thorsten
			for(int v=0; v < numvec; v++) { 
				//#pragma simd
				for(int s=0; s < V; s++) {
					dst[v*V+s]=src[v*V+s];
				}
			}
		}
    

		// Generic stream out
		template<typename FT, int V>
		inline void
		streamOutSpinor(volatile FT* restrict dst, const volatile FT* restrict src, int numvec) { 
      
			//#pragma vector aligned(src)
			//#pragma vector aligned(dst)
			//#pragma vector nontemporal(dst)
#pragma omp simd collapse(2) aligned(dst,src:V) //by thorsten
			for(int v=0; v < numvec; v++) { 
				//#pragma simd
				for(int s=0; s < V; s++) {
					dst[v*V+s]=src[v*V+s];
				}
			}
		}
    

		// Stream In to a different type
		template<typename FT, int V>
		inline void
		streamInSpinor(volatile typename ArithType<FT>::Type* restrict dst, const volatile  FT* restrict src, int numvec) { 
      
#if defined(__MIC__)
			//Intel MIC
			const int prefdist1 = 12;
			const int prefdist2 = 64;
      
			const char* prefl1base = (const char *)src+prefdist1*64;
			const char* prefl2base = (const char *)src+prefdist2*64;
			
			for(int v=0; v < numvec; v++) { 
				_mm_prefetch(&prefl1base[v*V*sizeof(FT)], _MM_HINT_T0);
				_mm_prefetch(&prefl2base[v*V*sizeof(FT)], _MM_HINT_T1);
				
#pragma simd
				for(int s=0; s < V; s++) {
					dst[v*V+s]=src[v*V+s];
				}
			}
#else
			//Generic
#pragma omp simd collapse(2) aligned(dst,src:V)
			for(int v=0; v < numvec; v++) { 
				for(int s=0; s < V; s++) {
					dst[v*V+s]=src[v*V+s];
				}
			}
#endif
		}
    
		// Write out to a different type
		template<typename FT, int V>
		inline void
		writeSpinor(volatile FT* restrict dst, const volatile  typename ArithType<FT>::Type* restrict src, int numvec) { 
      
			//#pragma vector aligned(src)
			//#pragma vector aligned(dst)
#pragma omp simd collapse(2) aligned(dst,src:V)
			for(int v=0; v < numvec; v++) { 
				//#pragma simd
				for(int s=0; s < V; s++) {
					dst[v*V+s]=src[v*V+s];
				}
			}
		}
    

		// Stream out to a different type
		template<typename FT, int V>
		inline void
		streamOutSpinor(volatile FT* restrict dst, const volatile typename ArithType<FT>::Type* restrict src, int numvec) { 
      
			//#pragma vector aligned(src)
			//#pragma vector aligned(dst)
			//#pragma vector nontemporal(dst)
#pragma omp simd collapse(2) aligned(dst,src:V)
			for(int v=0; v < numvec; v++) { 
				//#pragma simd
				for(int s=0; s < V; s++) {
					dst[v*V+s]=src[v*V+s];
				}
			}
		}


#if defined(__MIC__)
#include <immintrin.h>

		// Half prec specicialize 
		template<>
		inline void
		streamInSpinor<half,16>(volatile typename ArithType<half>::Type* restrict dst, const volatile half* restrict src, int numvec) { 
     
			const int prefdist1 = 12;
			const int prefdist2 = 64;
      
			const char* prefl1base = (const char *)src+prefdist1*64;
			const char* prefl2base = (const char *)src+prefdist2*64;

			for(int v=0; v < numvec; v++) { 
				_mm_prefetch(&prefl1base[v*16*sizeof(half)], _MM_HINT_T0);
				_mm_prefetch(&prefl2base[v*16*sizeof(half)], _MM_HINT_T1);
	
				__m512 r = _mm512_extload_ps((void *)&src[v*16], _MM_UPCONV_PS_FLOAT16, _MM_BROADCAST32_NONE, _MM_HINT_T0);
				_mm512_store_ps((void *)&dst[v*16], r);

			}
		}
    
   
		template<>
		inline void
		writeSpinor<half,16>(volatile half* restrict dst, const volatile typename ArithType<half>::Type* restrict src, int numvec) { 
			const int prefdist1 = 12;
			const int prefdist2 = 64;
      
			const char* prefl1base = (const char *)src+prefdist1*64;
			const char* prefl2base = (const char *)src+prefdist2*64;

			const char* prefl1baseo = (const char *)dst+prefdist1*64;
			const char* prefl2baseo = (const char *)dst+prefdist2*64;
      
			for(int v=0; v < numvec; v++) { 
				_mm_prefetch(&prefl1baseo[v*16*sizeof(half)], _MM_HINT_T0); // Prefetch for write
				_mm_prefetch(&prefl2baseo[v*16*sizeof(half)], _MM_HINT_T1); // Prefetch for write

				__m512 r = _mm512_load_ps((void *)&src[v*16]);
				_mm512_extstore_ps((void *)&dst[v*16], r, _MM_DOWNCONV_PS_FLOAT16, _MM_HINT_T0);
			}
		}


		template<>
		inline void
		streamOutSpinor<half,16>(volatile half* restrict dst, const volatile typename ArithType<half>::Type* restrict src, int numvec) { 
			const int prefdist1 = 12;
			const int prefdist2 = 64;
      
			const char* prefl1base = (const char *)src+prefdist1*64;
			const char* prefl2base = (const char *)src+prefdist2*64;
      
			for(int v=0; v < numvec; v++) { 
				_mm_prefetch(&prefl1base[v*16*sizeof(half)], _MM_HINT_T0);
				_mm_prefetch(&prefl2base[v*16*sizeof(half)], _MM_HINT_T1);

				__m512 r = _mm512_load_ps((void *)&src[v*16]);
				_mm512_extstore_ps((void *)&dst[v*16], r, _MM_DOWNCONV_PS_FLOAT16, _MM_HINT_NT);
			}
		}

#endif // defined MIC    

    
	}; // Namespace BLAS UTILS
  
};



#endif
