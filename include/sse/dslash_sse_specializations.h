
#include <mmintrin.h>


/**** FIXME: TODO Implement compression ******/
template<>
  inline 
void matMultVec<float,4>( float ub[2][3][2][4], const float b[2][3][2][4], const float  u[3][3][2][4], const int nsites, const bool compress12)
{
__m128 b_S0_C0_RE = _mm_setzero_ps(); 
__m128 b_S0_C0_IM = _mm_setzero_ps(); 
__m128 b_S0_C1_RE = _mm_setzero_ps(); 
__m128 b_S0_C1_IM = _mm_setzero_ps(); 
__m128 b_S0_C2_RE = _mm_setzero_ps(); 
__m128 b_S0_C2_IM = _mm_setzero_ps(); 
__m128 b_S1_C0_RE = _mm_setzero_ps(); 
__m128 b_S1_C0_IM = _mm_setzero_ps(); 
__m128 b_S1_C1_RE = _mm_setzero_ps(); 
__m128 b_S1_C1_IM = _mm_setzero_ps(); 
__m128 b_S1_C2_RE = _mm_setzero_ps(); 
__m128 b_S1_C2_IM = _mm_setzero_ps(); 
__m128 ub_S0_C0_RE = _mm_setzero_ps(); 
__m128 ub_S0_C0_IM = _mm_setzero_ps(); 
__m128 ub_S0_C1_RE = _mm_setzero_ps(); 
__m128 ub_S0_C1_IM = _mm_setzero_ps(); 
__m128 ub_S0_C2_RE = _mm_setzero_ps(); 
__m128 ub_S0_C2_IM = _mm_setzero_ps(); 
__m128 ub_S1_C0_RE = _mm_setzero_ps(); 
__m128 ub_S1_C0_IM = _mm_setzero_ps(); 
__m128 ub_S1_C1_RE = _mm_setzero_ps(); 
__m128 ub_S1_C1_IM = _mm_setzero_ps(); 
__m128 ub_S1_C2_RE = _mm_setzero_ps(); 
__m128 ub_S1_C2_IM = _mm_setzero_ps(); 
__m128 u_00_re = _mm_setzero_ps(); 
__m128 u_00_im = _mm_setzero_ps(); 
__m128 u_10_re = _mm_setzero_ps(); 
__m128 u_10_im = _mm_setzero_ps(); 
__m128 u_20_re = _mm_setzero_ps(); 
__m128 u_20_im = _mm_setzero_ps(); 
__m128 u_01_re = _mm_setzero_ps(); 
__m128 u_01_im = _mm_setzero_ps(); 
__m128 u_11_re = _mm_setzero_ps(); 
__m128 u_11_im = _mm_setzero_ps(); 
__m128 u_21_re = _mm_setzero_ps(); 
__m128 u_21_im = _mm_setzero_ps(); 
__m128 u_02_re = _mm_setzero_ps(); 
__m128 u_02_im = _mm_setzero_ps(); 
__m128 u_12_re = _mm_setzero_ps(); 
__m128 u_12_im = _mm_setzero_ps(); 
__m128 u_22_re = _mm_setzero_ps(); 
__m128 u_22_im = _mm_setzero_ps(); 
b_S0_C0_RE = _mm_load_ps(&(b[0][0][0][0]));

b_S0_C0_IM = _mm_load_ps(&(b[0][0][1][0]));

b_S0_C1_RE = _mm_load_ps(&(b[0][1][0][0]));

b_S0_C1_IM = _mm_load_ps(&(b[0][1][1][0]));

b_S0_C2_RE = _mm_load_ps(&(b[0][2][0][0]));

b_S0_C2_IM = _mm_load_ps(&(b[0][2][1][0]));

b_S1_C0_RE = _mm_load_ps(&(b[1][0][0][0]));

b_S1_C0_IM = _mm_load_ps(&(b[1][0][1][0]));

b_S1_C1_RE = _mm_load_ps(&(b[1][1][0][0]));

b_S1_C1_IM = _mm_load_ps(&(b[1][1][1][0]));

b_S1_C2_RE = _mm_load_ps(&(b[1][2][0][0]));

b_S1_C2_IM = _mm_load_ps(&(b[1][2][1][0]));

u_00_re = _mm_load_ps(u[0][0][0]);

u_00_im = _mm_load_ps(u[0][0][1]);

u_01_re = _mm_load_ps(u[0][1][0]);

u_01_im = _mm_load_ps(u[0][1][1]);

u_02_re = _mm_load_ps(u[0][2][0]);

u_02_im = _mm_load_ps(u[0][2][1]);

u_10_re = _mm_load_ps(u[1][0][0]);

u_10_im = _mm_load_ps(u[1][0][1]);

u_11_re = _mm_load_ps(u[1][1][0]);

u_11_im = _mm_load_ps(u[1][1][1]);

u_12_re = _mm_load_ps(u[1][2][0]);

u_12_im = _mm_load_ps(u[1][2][1]);

 if( ! compress12 ) { 
   u_20_re = _mm_load_ps(u[2][0][0]);
   
   u_20_im = _mm_load_ps(u[2][0][1]);
   
   u_21_re = _mm_load_ps(u[2][1][0]);
   
   u_21_im = _mm_load_ps(u[2][1][1]);
   
   u_22_re = _mm_load_ps(u[2][2][0]);
   
   u_22_im = _mm_load_ps(u[2][2][1]);
 }
 else { 
   u_20_re = _mm_mul_ps( u_01_re , u_12_re );
   u_20_re = _mm_sub_ps( u_20_re,  _mm_mul_ps(  u_01_im, u_12_im) );
   u_20_re = _mm_sub_ps( u_20_re,  _mm_mul_ps(  u_11_re, u_02_re) );
   u_20_re = _mm_add_ps( u_20_re, _mm_mul_ps( u_11_im, u_02_im ) ); 
   u_20_im = _mm_mul_ps( u_11_re , u_02_im );
   u_20_im = _mm_add_ps( u_20_im, _mm_mul_ps( u_11_im, u_02_re ) ); 
   u_20_im = _mm_sub_ps( u_20_im,  _mm_mul_ps(  u_01_re, u_12_im) );
   u_20_im = _mm_sub_ps( u_20_im,  _mm_mul_ps(  u_01_im, u_12_re) );
   u_21_re = _mm_mul_ps( u_10_re , u_02_re );
   u_21_re = _mm_sub_ps( u_21_re,  _mm_mul_ps(  u_10_im, u_02_im) );
   u_21_re = _mm_sub_ps( u_21_re,  _mm_mul_ps(  u_00_re, u_12_re) );
   u_21_re = _mm_add_ps( u_21_re, _mm_mul_ps( u_00_im, u_12_im ) ); 
   u_21_im = _mm_mul_ps( u_00_re , u_12_im );
   u_21_im = _mm_add_ps( u_21_im, _mm_mul_ps( u_00_im, u_12_re ) ); 
   u_21_im = _mm_sub_ps( u_21_im,  _mm_mul_ps(  u_10_re, u_02_im) );
   u_21_im = _mm_sub_ps( u_21_im,  _mm_mul_ps(  u_10_im, u_02_re) );
   u_22_re = _mm_mul_ps( u_00_re , u_11_re );
   u_22_re = _mm_sub_ps( u_22_re,  _mm_mul_ps(  u_00_im, u_11_im) );
   u_22_re = _mm_sub_ps( u_22_re,  _mm_mul_ps(  u_10_re, u_01_re) );
   u_22_re = _mm_add_ps( u_22_re, _mm_mul_ps( u_10_im, u_01_im ) ); 
   u_22_im = _mm_mul_ps( u_10_re , u_01_im );
   u_22_im = _mm_add_ps( u_22_im, _mm_mul_ps( u_10_im, u_01_re ) ); 
   u_22_im = _mm_sub_ps( u_22_im,  _mm_mul_ps(  u_00_re, u_11_im) );
   u_22_im = _mm_sub_ps( u_22_im,  _mm_mul_ps(  u_00_im, u_11_re) );
 }

ub_S0_C0_RE = _mm_mul_ps( u_00_re , b_S0_C0_RE );
ub_S0_C0_RE = _mm_sub_ps( ub_S0_C0_RE,  _mm_mul_ps(  u_00_im, b_S0_C0_IM) );
ub_S0_C0_IM = _mm_mul_ps( u_00_re , b_S0_C0_IM );
ub_S0_C0_IM = _mm_add_ps( ub_S0_C0_IM, _mm_mul_ps( u_00_im, b_S0_C0_RE ) ); 
ub_S0_C0_RE = _mm_add_ps( ub_S0_C0_RE, _mm_mul_ps( u_10_re, b_S0_C1_RE ) ); 
ub_S0_C0_RE = _mm_sub_ps( ub_S0_C0_RE,  _mm_mul_ps(  u_10_im, b_S0_C1_IM) );
ub_S0_C0_IM = _mm_add_ps( ub_S0_C0_IM, _mm_mul_ps( u_10_re, b_S0_C1_IM ) ); 
ub_S0_C0_IM = _mm_add_ps( ub_S0_C0_IM, _mm_mul_ps( u_10_im, b_S0_C1_RE ) ); 
ub_S0_C0_RE = _mm_add_ps( ub_S0_C0_RE, _mm_mul_ps( u_20_re, b_S0_C2_RE ) ); 
ub_S0_C0_RE = _mm_sub_ps( ub_S0_C0_RE,  _mm_mul_ps(  u_20_im, b_S0_C2_IM) );
ub_S0_C0_IM = _mm_add_ps( ub_S0_C0_IM, _mm_mul_ps( u_20_re, b_S0_C2_IM ) ); 
ub_S0_C0_IM = _mm_add_ps( ub_S0_C0_IM, _mm_mul_ps( u_20_im, b_S0_C2_RE ) ); 
ub_S0_C1_RE = _mm_mul_ps( u_01_re , b_S0_C0_RE );
ub_S0_C1_RE = _mm_sub_ps( ub_S0_C1_RE,  _mm_mul_ps(  u_01_im, b_S0_C0_IM) );
ub_S0_C1_IM = _mm_mul_ps( u_01_re , b_S0_C0_IM );
ub_S0_C1_IM = _mm_add_ps( ub_S0_C1_IM, _mm_mul_ps( u_01_im, b_S0_C0_RE ) ); 
ub_S0_C1_RE = _mm_add_ps( ub_S0_C1_RE, _mm_mul_ps( u_11_re, b_S0_C1_RE ) ); 
ub_S0_C1_RE = _mm_sub_ps( ub_S0_C1_RE,  _mm_mul_ps(  u_11_im, b_S0_C1_IM) );
ub_S0_C1_IM = _mm_add_ps( ub_S0_C1_IM, _mm_mul_ps( u_11_re, b_S0_C1_IM ) ); 
ub_S0_C1_IM = _mm_add_ps( ub_S0_C1_IM, _mm_mul_ps( u_11_im, b_S0_C1_RE ) ); 
ub_S0_C1_RE = _mm_add_ps( ub_S0_C1_RE, _mm_mul_ps( u_21_re, b_S0_C2_RE ) ); 
ub_S0_C1_RE = _mm_sub_ps( ub_S0_C1_RE,  _mm_mul_ps(  u_21_im, b_S0_C2_IM) );
ub_S0_C1_IM = _mm_add_ps( ub_S0_C1_IM, _mm_mul_ps( u_21_re, b_S0_C2_IM ) ); 
ub_S0_C1_IM = _mm_add_ps( ub_S0_C1_IM, _mm_mul_ps( u_21_im, b_S0_C2_RE ) ); 
ub_S0_C2_RE = _mm_mul_ps( u_02_re , b_S0_C0_RE );
ub_S0_C2_RE = _mm_sub_ps( ub_S0_C2_RE,  _mm_mul_ps(  u_02_im, b_S0_C0_IM) );
ub_S0_C2_IM = _mm_mul_ps( u_02_re , b_S0_C0_IM );
ub_S0_C2_IM = _mm_add_ps( ub_S0_C2_IM, _mm_mul_ps( u_02_im, b_S0_C0_RE ) ); 
ub_S0_C2_RE = _mm_add_ps( ub_S0_C2_RE, _mm_mul_ps( u_12_re, b_S0_C1_RE ) ); 
ub_S0_C2_RE = _mm_sub_ps( ub_S0_C2_RE,  _mm_mul_ps(  u_12_im, b_S0_C1_IM) );
ub_S0_C2_IM = _mm_add_ps( ub_S0_C2_IM, _mm_mul_ps( u_12_re, b_S0_C1_IM ) ); 
ub_S0_C2_IM = _mm_add_ps( ub_S0_C2_IM, _mm_mul_ps( u_12_im, b_S0_C1_RE ) ); 
ub_S0_C2_RE = _mm_add_ps( ub_S0_C2_RE, _mm_mul_ps( u_22_re, b_S0_C2_RE ) ); 
ub_S0_C2_RE = _mm_sub_ps( ub_S0_C2_RE,  _mm_mul_ps(  u_22_im, b_S0_C2_IM) );
ub_S0_C2_IM = _mm_add_ps( ub_S0_C2_IM, _mm_mul_ps( u_22_re, b_S0_C2_IM ) ); 
ub_S0_C2_IM = _mm_add_ps( ub_S0_C2_IM, _mm_mul_ps( u_22_im, b_S0_C2_RE ) ); 
ub_S1_C0_RE = _mm_mul_ps( u_00_re , b_S1_C0_RE );
ub_S1_C0_RE = _mm_sub_ps( ub_S1_C0_RE,  _mm_mul_ps(  u_00_im, b_S1_C0_IM) );
ub_S1_C0_IM = _mm_mul_ps( u_00_re , b_S1_C0_IM );
ub_S1_C0_IM = _mm_add_ps( ub_S1_C0_IM, _mm_mul_ps( u_00_im, b_S1_C0_RE ) ); 
ub_S1_C0_RE = _mm_add_ps( ub_S1_C0_RE, _mm_mul_ps( u_10_re, b_S1_C1_RE ) ); 
ub_S1_C0_RE = _mm_sub_ps( ub_S1_C0_RE,  _mm_mul_ps(  u_10_im, b_S1_C1_IM) );
ub_S1_C0_IM = _mm_add_ps( ub_S1_C0_IM, _mm_mul_ps( u_10_re, b_S1_C1_IM ) ); 
ub_S1_C0_IM = _mm_add_ps( ub_S1_C0_IM, _mm_mul_ps( u_10_im, b_S1_C1_RE ) ); 
ub_S1_C0_RE = _mm_add_ps( ub_S1_C0_RE, _mm_mul_ps( u_20_re, b_S1_C2_RE ) ); 
ub_S1_C0_RE = _mm_sub_ps( ub_S1_C0_RE,  _mm_mul_ps(  u_20_im, b_S1_C2_IM) );
ub_S1_C0_IM = _mm_add_ps( ub_S1_C0_IM, _mm_mul_ps( u_20_re, b_S1_C2_IM ) ); 
ub_S1_C0_IM = _mm_add_ps( ub_S1_C0_IM, _mm_mul_ps( u_20_im, b_S1_C2_RE ) ); 
ub_S1_C1_RE = _mm_mul_ps( u_01_re , b_S1_C0_RE );
ub_S1_C1_RE = _mm_sub_ps( ub_S1_C1_RE,  _mm_mul_ps(  u_01_im, b_S1_C0_IM) );
ub_S1_C1_IM = _mm_mul_ps( u_01_re , b_S1_C0_IM );
ub_S1_C1_IM = _mm_add_ps( ub_S1_C1_IM, _mm_mul_ps( u_01_im, b_S1_C0_RE ) ); 
ub_S1_C1_RE = _mm_add_ps( ub_S1_C1_RE, _mm_mul_ps( u_11_re, b_S1_C1_RE ) ); 
ub_S1_C1_RE = _mm_sub_ps( ub_S1_C1_RE,  _mm_mul_ps(  u_11_im, b_S1_C1_IM) );
ub_S1_C1_IM = _mm_add_ps( ub_S1_C1_IM, _mm_mul_ps( u_11_re, b_S1_C1_IM ) ); 
ub_S1_C1_IM = _mm_add_ps( ub_S1_C1_IM, _mm_mul_ps( u_11_im, b_S1_C1_RE ) ); 
ub_S1_C1_RE = _mm_add_ps( ub_S1_C1_RE, _mm_mul_ps( u_21_re, b_S1_C2_RE ) ); 
ub_S1_C1_RE = _mm_sub_ps( ub_S1_C1_RE,  _mm_mul_ps(  u_21_im, b_S1_C2_IM) );
ub_S1_C1_IM = _mm_add_ps( ub_S1_C1_IM, _mm_mul_ps( u_21_re, b_S1_C2_IM ) ); 
ub_S1_C1_IM = _mm_add_ps( ub_S1_C1_IM, _mm_mul_ps( u_21_im, b_S1_C2_RE ) ); 
ub_S1_C2_RE = _mm_mul_ps( u_02_re , b_S1_C0_RE );
ub_S1_C2_RE = _mm_sub_ps( ub_S1_C2_RE,  _mm_mul_ps(  u_02_im, b_S1_C0_IM) );
ub_S1_C2_IM = _mm_mul_ps( u_02_re , b_S1_C0_IM );
ub_S1_C2_IM = _mm_add_ps( ub_S1_C2_IM, _mm_mul_ps( u_02_im, b_S1_C0_RE ) ); 
ub_S1_C2_RE = _mm_add_ps( ub_S1_C2_RE, _mm_mul_ps( u_12_re, b_S1_C1_RE ) ); 
ub_S1_C2_RE = _mm_sub_ps( ub_S1_C2_RE,  _mm_mul_ps(  u_12_im, b_S1_C1_IM) );
ub_S1_C2_IM = _mm_add_ps( ub_S1_C2_IM, _mm_mul_ps( u_12_re, b_S1_C1_IM ) ); 
ub_S1_C2_IM = _mm_add_ps( ub_S1_C2_IM, _mm_mul_ps( u_12_im, b_S1_C1_RE ) ); 
ub_S1_C2_RE = _mm_add_ps( ub_S1_C2_RE, _mm_mul_ps( u_22_re, b_S1_C2_RE ) ); 
ub_S1_C2_RE = _mm_sub_ps( ub_S1_C2_RE,  _mm_mul_ps(  u_22_im, b_S1_C2_IM) );
ub_S1_C2_IM = _mm_add_ps( ub_S1_C2_IM, _mm_mul_ps( u_22_re, b_S1_C2_IM ) ); 
ub_S1_C2_IM = _mm_add_ps( ub_S1_C2_IM, _mm_mul_ps( u_22_im, b_S1_C2_RE ) ); 
_mm_stream_ps(&(ub[0][0][0][0]),ub_S0_C0_RE);

_mm_stream_ps(&(ub[0][0][1][0]),ub_S0_C0_IM);

_mm_stream_ps(&(ub[0][1][0][0]),ub_S0_C1_RE);

_mm_stream_ps(&(ub[0][1][1][0]),ub_S0_C1_IM);

_mm_stream_ps(&(ub[0][2][0][0]),ub_S0_C2_RE);

_mm_stream_ps(&(ub[0][2][1][0]),ub_S0_C2_IM);

_mm_stream_ps(&(ub[1][0][0][0]),ub_S1_C0_RE);

_mm_stream_ps(&(ub[1][0][1][0]),ub_S1_C0_IM);

_mm_stream_ps(&(ub[1][1][0][0]),ub_S1_C1_RE);

_mm_stream_ps(&(ub[1][1][1][0]),ub_S1_C1_IM);

_mm_stream_ps(&(ub[1][2][0][0]),ub_S1_C2_RE);

_mm_stream_ps(&(ub[1][2][1][0]),ub_S1_C2_IM);



}

/**** FIXME: TODO Implement compression ******/
template<>
  inline 
void adjMatMultVec( float ub[2][3][2][4], const float b[2][3][2][4], const float  u[3][3][2][4], const int nsites, const bool compress12)
{
  __m128 b_S0_C0_RE = _mm_setzero_ps(); 
__m128 b_S0_C0_IM = _mm_setzero_ps(); 
__m128 b_S0_C1_RE = _mm_setzero_ps(); 
__m128 b_S0_C1_IM = _mm_setzero_ps(); 
__m128 b_S0_C2_RE = _mm_setzero_ps(); 
__m128 b_S0_C2_IM = _mm_setzero_ps(); 
__m128 b_S1_C0_RE = _mm_setzero_ps(); 
__m128 b_S1_C0_IM = _mm_setzero_ps(); 
__m128 b_S1_C1_RE = _mm_setzero_ps(); 
__m128 b_S1_C1_IM = _mm_setzero_ps(); 
__m128 b_S1_C2_RE = _mm_setzero_ps(); 
__m128 b_S1_C2_IM = _mm_setzero_ps(); 
__m128 ub_S0_C0_RE = _mm_setzero_ps(); 
__m128 ub_S0_C0_IM = _mm_setzero_ps(); 
__m128 ub_S0_C1_RE = _mm_setzero_ps(); 
__m128 ub_S0_C1_IM = _mm_setzero_ps(); 
__m128 ub_S0_C2_RE = _mm_setzero_ps(); 
__m128 ub_S0_C2_IM = _mm_setzero_ps(); 
__m128 ub_S1_C0_RE = _mm_setzero_ps(); 
__m128 ub_S1_C0_IM = _mm_setzero_ps(); 
__m128 ub_S1_C1_RE = _mm_setzero_ps(); 
__m128 ub_S1_C1_IM = _mm_setzero_ps(); 
__m128 ub_S1_C2_RE = _mm_setzero_ps(); 
__m128 ub_S1_C2_IM = _mm_setzero_ps(); 
__m128 u_00_re = _mm_setzero_ps(); 
__m128 u_00_im = _mm_setzero_ps(); 
__m128 u_10_re = _mm_setzero_ps(); 
__m128 u_10_im = _mm_setzero_ps(); 
__m128 u_20_re = _mm_setzero_ps(); 
__m128 u_20_im = _mm_setzero_ps(); 
__m128 u_01_re = _mm_setzero_ps(); 
__m128 u_01_im = _mm_setzero_ps(); 
__m128 u_11_re = _mm_setzero_ps(); 
__m128 u_11_im = _mm_setzero_ps(); 
__m128 u_21_re = _mm_setzero_ps(); 
__m128 u_21_im = _mm_setzero_ps(); 
__m128 u_02_re = _mm_setzero_ps(); 
__m128 u_02_im = _mm_setzero_ps(); 
__m128 u_12_re = _mm_setzero_ps(); 
__m128 u_12_im = _mm_setzero_ps(); 
__m128 u_22_re = _mm_setzero_ps(); 
__m128 u_22_im = _mm_setzero_ps(); 
b_S0_C0_RE = _mm_load_ps(&(b[0][0][0][0]));

b_S0_C0_IM = _mm_load_ps(&(b[0][0][1][0]));

b_S0_C1_RE = _mm_load_ps(&(b[0][1][0][0]));

b_S0_C1_IM = _mm_load_ps(&(b[0][1][1][0]));

b_S0_C2_RE = _mm_load_ps(&(b[0][2][0][0]));

b_S0_C2_IM = _mm_load_ps(&(b[0][2][1][0]));

b_S1_C0_RE = _mm_load_ps(&(b[1][0][0][0]));

b_S1_C0_IM = _mm_load_ps(&(b[1][0][1][0]));

b_S1_C1_RE = _mm_load_ps(&(b[1][1][0][0]));

b_S1_C1_IM = _mm_load_ps(&(b[1][1][1][0]));

b_S1_C2_RE = _mm_load_ps(&(b[1][2][0][0]));

b_S1_C2_IM = _mm_load_ps(&(b[1][2][1][0]));

u_00_re = _mm_load_ps(u[0][0][0]);

u_00_im = _mm_load_ps(u[0][0][1]);

u_01_re = _mm_load_ps(u[0][1][0]);

u_01_im = _mm_load_ps(u[0][1][1]);

u_02_re = _mm_load_ps(u[0][2][0]);

u_02_im = _mm_load_ps(u[0][2][1]);

u_10_re = _mm_load_ps(u[1][0][0]);

u_10_im = _mm_load_ps(u[1][0][1]);

u_11_re = _mm_load_ps(u[1][1][0]);

u_11_im = _mm_load_ps(u[1][1][1]);

u_12_re = _mm_load_ps(u[1][2][0]);

u_12_im = _mm_load_ps(u[1][2][1]);

 if( !compress12)  {
   u_20_re = _mm_load_ps(u[2][0][0]);
   
   u_20_im = _mm_load_ps(u[2][0][1]);
   
   u_21_re = _mm_load_ps(u[2][1][0]);
   
   u_21_im = _mm_load_ps(u[2][1][1]);
   
   u_22_re = _mm_load_ps(u[2][2][0]);
   
   u_22_im = _mm_load_ps(u[2][2][1]);
 }
 else { 

u_20_re = _mm_mul_ps( u_01_re , u_12_re );
u_20_re = _mm_sub_ps( u_20_re,  _mm_mul_ps(  u_01_im, u_12_im) );
u_20_re = _mm_sub_ps( u_20_re,  _mm_mul_ps(  u_11_re, u_02_re) );
u_20_re = _mm_add_ps( u_20_re, _mm_mul_ps( u_11_im, u_02_im ) ); 
u_20_im = _mm_mul_ps( u_11_re , u_02_im );
u_20_im = _mm_add_ps( u_20_im, _mm_mul_ps( u_11_im, u_02_re ) ); 
u_20_im = _mm_sub_ps( u_20_im,  _mm_mul_ps(  u_01_re, u_12_im) );
u_20_im = _mm_sub_ps( u_20_im,  _mm_mul_ps(  u_01_im, u_12_re) );
u_21_re = _mm_mul_ps( u_10_re , u_02_re );
u_21_re = _mm_sub_ps( u_21_re,  _mm_mul_ps(  u_10_im, u_02_im) );
u_21_re = _mm_sub_ps( u_21_re,  _mm_mul_ps(  u_00_re, u_12_re) );
u_21_re = _mm_add_ps( u_21_re, _mm_mul_ps( u_00_im, u_12_im ) ); 
u_21_im = _mm_mul_ps( u_00_re , u_12_im );
u_21_im = _mm_add_ps( u_21_im, _mm_mul_ps( u_00_im, u_12_re ) ); 
u_21_im = _mm_sub_ps( u_21_im,  _mm_mul_ps(  u_10_re, u_02_im) );
u_21_im = _mm_sub_ps( u_21_im,  _mm_mul_ps(  u_10_im, u_02_re) );
u_22_re = _mm_mul_ps( u_00_re , u_11_re );
u_22_re = _mm_sub_ps( u_22_re,  _mm_mul_ps(  u_00_im, u_11_im) );
u_22_re = _mm_sub_ps( u_22_re,  _mm_mul_ps(  u_10_re, u_01_re) );
u_22_re = _mm_add_ps( u_22_re, _mm_mul_ps( u_10_im, u_01_im ) ); 
u_22_im = _mm_mul_ps( u_10_re , u_01_im );
u_22_im = _mm_add_ps( u_22_im, _mm_mul_ps( u_10_im, u_01_re ) ); 
u_22_im = _mm_sub_ps( u_22_im,  _mm_mul_ps(  u_00_re, u_11_im) );
u_22_im = _mm_sub_ps( u_22_im,  _mm_mul_ps(  u_00_im, u_11_re) );

 }
ub_S0_C0_RE = _mm_mul_ps( u_00_re , b_S0_C0_RE );
ub_S0_C0_IM = _mm_mul_ps( u_00_re , b_S0_C0_IM );
ub_S0_C0_RE = _mm_add_ps( ub_S0_C0_RE, _mm_mul_ps( u_00_im, b_S0_C0_IM ) ); 
ub_S0_C0_IM = _mm_sub_ps( ub_S0_C0_IM,  _mm_mul_ps(  u_00_im, b_S0_C0_RE) );
ub_S0_C0_RE = _mm_add_ps( ub_S0_C0_RE, _mm_mul_ps( u_01_re, b_S0_C1_RE ) ); 
ub_S0_C0_IM = _mm_add_ps( ub_S0_C0_IM, _mm_mul_ps( u_01_re, b_S0_C1_IM ) ); 
ub_S0_C0_RE = _mm_add_ps( ub_S0_C0_RE, _mm_mul_ps( u_01_im, b_S0_C1_IM ) ); 
ub_S0_C0_IM = _mm_sub_ps( ub_S0_C0_IM,  _mm_mul_ps(  u_01_im, b_S0_C1_RE) );
ub_S0_C0_RE = _mm_add_ps( ub_S0_C0_RE, _mm_mul_ps( u_02_re, b_S0_C2_RE ) ); 
ub_S0_C0_IM = _mm_add_ps( ub_S0_C0_IM, _mm_mul_ps( u_02_re, b_S0_C2_IM ) ); 
ub_S0_C0_RE = _mm_add_ps( ub_S0_C0_RE, _mm_mul_ps( u_02_im, b_S0_C2_IM ) ); 
ub_S0_C0_IM = _mm_sub_ps( ub_S0_C0_IM,  _mm_mul_ps(  u_02_im, b_S0_C2_RE) );
ub_S0_C1_RE = _mm_mul_ps( u_10_re , b_S0_C0_RE );
ub_S0_C1_IM = _mm_mul_ps( u_10_re , b_S0_C0_IM );
ub_S0_C1_RE = _mm_add_ps( ub_S0_C1_RE, _mm_mul_ps( u_10_im, b_S0_C0_IM ) ); 
ub_S0_C1_IM = _mm_sub_ps( ub_S0_C1_IM,  _mm_mul_ps(  u_10_im, b_S0_C0_RE) );
ub_S0_C1_RE = _mm_add_ps( ub_S0_C1_RE, _mm_mul_ps( u_11_re, b_S0_C1_RE ) ); 
ub_S0_C1_IM = _mm_add_ps( ub_S0_C1_IM, _mm_mul_ps( u_11_re, b_S0_C1_IM ) ); 
ub_S0_C1_RE = _mm_add_ps( ub_S0_C1_RE, _mm_mul_ps( u_11_im, b_S0_C1_IM ) ); 
ub_S0_C1_IM = _mm_sub_ps( ub_S0_C1_IM,  _mm_mul_ps(  u_11_im, b_S0_C1_RE) );
ub_S0_C1_RE = _mm_add_ps( ub_S0_C1_RE, _mm_mul_ps( u_12_re, b_S0_C2_RE ) ); 
ub_S0_C1_IM = _mm_add_ps( ub_S0_C1_IM, _mm_mul_ps( u_12_re, b_S0_C2_IM ) ); 
ub_S0_C1_RE = _mm_add_ps( ub_S0_C1_RE, _mm_mul_ps( u_12_im, b_S0_C2_IM ) ); 
ub_S0_C1_IM = _mm_sub_ps( ub_S0_C1_IM,  _mm_mul_ps(  u_12_im, b_S0_C2_RE) );
ub_S0_C2_RE = _mm_mul_ps( u_20_re , b_S0_C0_RE );
ub_S0_C2_IM = _mm_mul_ps( u_20_re , b_S0_C0_IM );
ub_S0_C2_RE = _mm_add_ps( ub_S0_C2_RE, _mm_mul_ps( u_20_im, b_S0_C0_IM ) ); 
ub_S0_C2_IM = _mm_sub_ps( ub_S0_C2_IM,  _mm_mul_ps(  u_20_im, b_S0_C0_RE) );
ub_S0_C2_RE = _mm_add_ps( ub_S0_C2_RE, _mm_mul_ps( u_21_re, b_S0_C1_RE ) ); 
ub_S0_C2_IM = _mm_add_ps( ub_S0_C2_IM, _mm_mul_ps( u_21_re, b_S0_C1_IM ) ); 
ub_S0_C2_RE = _mm_add_ps( ub_S0_C2_RE, _mm_mul_ps( u_21_im, b_S0_C1_IM ) ); 
ub_S0_C2_IM = _mm_sub_ps( ub_S0_C2_IM,  _mm_mul_ps(  u_21_im, b_S0_C1_RE) );
ub_S0_C2_RE = _mm_add_ps( ub_S0_C2_RE, _mm_mul_ps( u_22_re, b_S0_C2_RE ) ); 
ub_S0_C2_IM = _mm_add_ps( ub_S0_C2_IM, _mm_mul_ps( u_22_re, b_S0_C2_IM ) ); 
ub_S0_C2_RE = _mm_add_ps( ub_S0_C2_RE, _mm_mul_ps( u_22_im, b_S0_C2_IM ) ); 
ub_S0_C2_IM = _mm_sub_ps( ub_S0_C2_IM,  _mm_mul_ps(  u_22_im, b_S0_C2_RE) );
ub_S1_C0_RE = _mm_mul_ps( u_00_re , b_S1_C0_RE );
ub_S1_C0_IM = _mm_mul_ps( u_00_re , b_S1_C0_IM );
ub_S1_C0_RE = _mm_add_ps( ub_S1_C0_RE, _mm_mul_ps( u_00_im, b_S1_C0_IM ) ); 
ub_S1_C0_IM = _mm_sub_ps( ub_S1_C0_IM,  _mm_mul_ps(  u_00_im, b_S1_C0_RE) );
ub_S1_C0_RE = _mm_add_ps( ub_S1_C0_RE, _mm_mul_ps( u_01_re, b_S1_C1_RE ) ); 
ub_S1_C0_IM = _mm_add_ps( ub_S1_C0_IM, _mm_mul_ps( u_01_re, b_S1_C1_IM ) ); 
ub_S1_C0_RE = _mm_add_ps( ub_S1_C0_RE, _mm_mul_ps( u_01_im, b_S1_C1_IM ) ); 
ub_S1_C0_IM = _mm_sub_ps( ub_S1_C0_IM,  _mm_mul_ps(  u_01_im, b_S1_C1_RE) );
ub_S1_C0_RE = _mm_add_ps( ub_S1_C0_RE, _mm_mul_ps( u_02_re, b_S1_C2_RE ) ); 
ub_S1_C0_IM = _mm_add_ps( ub_S1_C0_IM, _mm_mul_ps( u_02_re, b_S1_C2_IM ) ); 
ub_S1_C0_RE = _mm_add_ps( ub_S1_C0_RE, _mm_mul_ps( u_02_im, b_S1_C2_IM ) ); 
ub_S1_C0_IM = _mm_sub_ps( ub_S1_C0_IM,  _mm_mul_ps(  u_02_im, b_S1_C2_RE) );
ub_S1_C1_RE = _mm_mul_ps( u_10_re , b_S1_C0_RE );
ub_S1_C1_IM = _mm_mul_ps( u_10_re , b_S1_C0_IM );
ub_S1_C1_RE = _mm_add_ps( ub_S1_C1_RE, _mm_mul_ps( u_10_im, b_S1_C0_IM ) ); 
ub_S1_C1_IM = _mm_sub_ps( ub_S1_C1_IM,  _mm_mul_ps(  u_10_im, b_S1_C0_RE) );
ub_S1_C1_RE = _mm_add_ps( ub_S1_C1_RE, _mm_mul_ps( u_11_re, b_S1_C1_RE ) ); 
ub_S1_C1_IM = _mm_add_ps( ub_S1_C1_IM, _mm_mul_ps( u_11_re, b_S1_C1_IM ) ); 
ub_S1_C1_RE = _mm_add_ps( ub_S1_C1_RE, _mm_mul_ps( u_11_im, b_S1_C1_IM ) ); 
ub_S1_C1_IM = _mm_sub_ps( ub_S1_C1_IM,  _mm_mul_ps(  u_11_im, b_S1_C1_RE) );
ub_S1_C1_RE = _mm_add_ps( ub_S1_C1_RE, _mm_mul_ps( u_12_re, b_S1_C2_RE ) ); 
ub_S1_C1_IM = _mm_add_ps( ub_S1_C1_IM, _mm_mul_ps( u_12_re, b_S1_C2_IM ) ); 
ub_S1_C1_RE = _mm_add_ps( ub_S1_C1_RE, _mm_mul_ps( u_12_im, b_S1_C2_IM ) ); 
ub_S1_C1_IM = _mm_sub_ps( ub_S1_C1_IM,  _mm_mul_ps(  u_12_im, b_S1_C2_RE) );
ub_S1_C2_RE = _mm_mul_ps( u_20_re , b_S1_C0_RE );
ub_S1_C2_IM = _mm_mul_ps( u_20_re , b_S1_C0_IM );
ub_S1_C2_RE = _mm_add_ps( ub_S1_C2_RE, _mm_mul_ps( u_20_im, b_S1_C0_IM ) ); 
ub_S1_C2_IM = _mm_sub_ps( ub_S1_C2_IM,  _mm_mul_ps(  u_20_im, b_S1_C0_RE) );
ub_S1_C2_RE = _mm_add_ps( ub_S1_C2_RE, _mm_mul_ps( u_21_re, b_S1_C1_RE ) ); 
ub_S1_C2_IM = _mm_add_ps( ub_S1_C2_IM, _mm_mul_ps( u_21_re, b_S1_C1_IM ) ); 
ub_S1_C2_RE = _mm_add_ps( ub_S1_C2_RE, _mm_mul_ps( u_21_im, b_S1_C1_IM ) ); 
ub_S1_C2_IM = _mm_sub_ps( ub_S1_C2_IM,  _mm_mul_ps(  u_21_im, b_S1_C1_RE) );
ub_S1_C2_RE = _mm_add_ps( ub_S1_C2_RE, _mm_mul_ps( u_22_re, b_S1_C2_RE ) ); 
ub_S1_C2_IM = _mm_add_ps( ub_S1_C2_IM, _mm_mul_ps( u_22_re, b_S1_C2_IM ) ); 
ub_S1_C2_RE = _mm_add_ps( ub_S1_C2_RE, _mm_mul_ps( u_22_im, b_S1_C2_IM ) ); 
ub_S1_C2_IM = _mm_sub_ps( ub_S1_C2_IM,  _mm_mul_ps(  u_22_im, b_S1_C2_RE) );
_mm_stream_ps(&(ub[0][0][0][0]),ub_S0_C0_RE);

_mm_stream_ps(&(ub[0][0][1][0]),ub_S0_C0_IM);

_mm_stream_ps(&(ub[0][1][0][0]),ub_S0_C1_RE);

_mm_stream_ps(&(ub[0][1][1][0]),ub_S0_C1_IM);

_mm_stream_ps(&(ub[0][2][0][0]),ub_S0_C2_RE);

_mm_stream_ps(&(ub[0][2][1][0]),ub_S0_C2_IM);

_mm_stream_ps(&(ub[1][0][0][0]),ub_S1_C0_RE);

_mm_stream_ps(&(ub[1][0][1][0]),ub_S1_C0_IM);

_mm_stream_ps(&(ub[1][1][0][0]),ub_S1_C1_RE);

_mm_stream_ps(&(ub[1][1][1][0]),ub_S1_C1_IM);

_mm_stream_ps(&(ub[1][2][0][0]),ub_S1_C2_RE);

_mm_stream_ps(&(ub[1][2][1][0]),ub_S1_C2_IM);




}

// NEED PROJ PLUS X FORW SSE  here
