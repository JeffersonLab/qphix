#include <mmintrin.h>
#if 1
template<> 
inline void 
dslash_plus_vec( const int curblock,
		 const int fblock_x_last, const int fblock_x,
		 const int fsite_x_last, const int fsite_x,
		 const int bblock_x_last, const int bblock_x,
		 const int bsite_x_last, const int bsite_x,
		 const int fblock_y, const int bblock_y,
		 const int fblock_z, const int bblock_z,
		 const int fblock_t, const int bblock_t,
		 const float  (*psi)[4][3][2][4], 
		 const float  (*u)[8][3][3][2][4],
		 float (*out)[4][3][2][4],
		 const int nsites,
		 const bool accumulate[6],  
		 const int gauge_prefdist,
		 const int bblock_y_next,
		 const int fblock_y_next,
                 const int bblock_z_next,
                 const int fblock_z_next,
                 const int fblock_t_next,
		 const int spinor_out_prefdist,
		 const bool compress12)
{
  
  if (compress12) {
#include "sse/dslash_plus_body_sse_float_12"
  }
  else {
#include "sse/dslash_plus_body_sse_float_18"
  }
}



template<> 
inline void 
dslash_minus_vec(const int curblock,
		 const int fblock_x_last, const int fblock_x,
		 const int fsite_x_last, const int fsite_x,
		 const int bblock_x_last, const int bblock_x,
		 const int bsite_x_last, const int bsite_x,
		 const int fblock_y, const int bblock_y,
		 const int fblock_z, const int bblock_z,
		 const int fblock_t, const int bblock_t,
		 const float  (*psi)[4][3][2][4], 
		 const float  (*u)[8][3][3][2][4],
		 float (*out)[4][3][2][4],
		 const int nsites,
		 const bool accumulate[6], 
		 const int gauge_prefdist,
		 const int bblock_y_next,
		 const int fblock_y_next,
                 const int bblock_z_next,
                 const int fblock_z_next,
                 const int fblock_t_next,
		 const int spinor_out_prefdist,
		 const bool compress12)
{
  if( compress12 ) { 
#include "sse/dslash_minus_body_sse_float_12"
  }
  else{ 
#include "sse/dslash_minus_body_sse_float_18"
  }
}


#if 0
template<> 
inline void 
dslash_plus_vec_aChiMinusBDPsi( const int curblock,
				const int fblock_x_last, const int fblock_x,
				const int fsite_x_last, const int fsite_x,
				const int bblock_x_last, const int bblock_x,
				const int bsite_x_last, const int bsite_x,
				const int fblock_y, const int bblock_y,
				const int fblock_z, const int bblock_z,
				const int fblock_t, const int bblock_t,
				const float  (*psi)[4][3][2][4], 
				const float  (*chi)[4][3][2][4], 
				const float  (*u)[8][3][3][2][4],
				const float alpha, 
				const float beta,
				float (*out)[4][3][2][4],
				const int nsites,
				const bool accumulate[6],  
				const int gauge_prefdist,
				const int bblock_y_next,
				const int fblock_y_next,
				const int bblock_z_next,
				const int fblock_z_next,
				const int fblock_t_next,
				const int spinor_out_prefdist,
				const bool compress12)
{
  
  if (compress12) {
#include "sse/dslash_achimbdpsi_plus_body_sse_float_12"
  }
  else {
#include "sse/dslash_achimbdpsi_plus_body_sse_float_18"
  }
}



template<> 
inline void 
dslash_minus_vec_aChiMinusBDPsi(const int curblock,
				const int fblock_x_last, const int fblock_x,
				const int fsite_x_last, const int fsite_x,
				const int bblock_x_last, const int bblock_x,
				const int bsite_x_last, const int bsite_x,
				const int fblock_y, const int bblock_y,
				const int fblock_z, const int bblock_z,
				const int fblock_t, const int bblock_t,
				const float  (*psi)[4][3][2][4], 
				const float  (*chi)[4][3][2][4], 
				const float  (*u)[8][3][3][2][4],
				const float alpha, 
				const float beta,
				float (*out)[4][3][2][4],
				const int nsites,
				const bool accumulate[6], 
				const int gauge_prefdist,
				const int bblock_y_next,
				const int fblock_y_next,
				const int bblock_z_next,
				const int fblock_z_next,
				const int fblock_t_next,
				const int spinor_out_prefdist,
				const bool compress12)
{
  if( compress12 ) { 
#include "sse/dslash_achimbdpsi_minus_body_sse_float_12"
  }
  else{ 
#include "sse/dslash_achimbdpsi_minus_body_sse_float_18"
  }
}
#endif

#endif


// THESE ARE HERE SO WE CAN LOOK AT THE ASSEMBLER OUTPUT
//
extern "C" {
__attribute__((noinline))
void 
dslash_plus_vec4_noinline(
			  const int curblock,
			  const int fblock_x_last, const int fblock_x,
			  const int fsite_x_last, const int fsite_x,
			  const int bblock_x_last, const int bblock_x,
			  const int bsite_x_last, const int bsite_x,
			  const int fblock_y, const int bblock_y,
			  const int fblock_z, const int bblock_z,
			  const int fblock_t, const int bblock_t,
			  const float  (*psi)[4][3][2][4], 
			  const float  (*u)[8][3][3][2][4],
			  float (*out)[4][3][2][4],
			  const int nsites,
			  const bool accumulate[6],
			  const int gauge_prefdist,
			  const int bblock_y_next,
			  const int fblock_y_next,
			  const int bblock_z_next,
			  const int fblock_z_next,
			  const int fblock_t_next,
			  const int spinor_out_prefdist, 
			  const bool compress12)
{

  if( compress12 ) { 
#include "sse/dslash_plus_body_sse_float_12"
  }
  else { 
#include "sse/dslash_plus_body_sse_float_18"
  }

}


__attribute__((noinline))
void 
dslash_minus_vec4_noinline(
			   const int curblock, 
			   const int fblock_x_last, const int fblock_x,
			   const int fsite_x_last, const int fsite_x,
			   const int bblock_x_last, const int bblock_x,
			   const int bsite_x_last, const int bsite_x,
			   const int fblock_y, const int bblock_y,
			   const int fblock_z, const int bblock_z,
			   const int fblock_t, const int bblock_t,
			   const float  (*psi)[4][3][2][4], 
			   const float  (*u)[8][3][3][2][4],
			  float (*out)[4][3][2][4],
			   const int nsites,
			   const bool accumulate[6],
			   const int gauge_prefdist,
			   const int bblock_y_next,
			   const int fblock_y_next,
			   const int bblock_z_next,
			   const int fblock_z_next,
			   const int fblock_t_next,
			   const int spinor_out_prefdist,
			   const bool compress12)
{
  if(compress12) {
#include "sse/dslash_minus_body_sse_float_12"
  }
  else {
#include "sse/dslash_minus_body_sse_float_18" 
  }
}

  
};
