/* 
 * 
 * DSLASH PLUS
 *
 */
#if !defined(SOA) && !defined(COMPRESS12) && !defined(COMPRESS_SUFFIX)

/*
 * DSLASH_PLUS_VEC:
 * Standard template if SOA, COMPRESS12, COMPRESS SUFFIX are all undefined
 *
 */
template<typename FT, int veclen, int soalen, bool compress12> 
inline void
dslash_plus_vec(    
		const typename Geometry<FT,veclen,soalen,compress12>::FourSpinorBlock *xyBase,
		const typename Geometry<FT,veclen,soalen,compress12>::FourSpinorBlock *zbBase,
		const typename Geometry<FT,veclen,soalen,compress12>::FourSpinorBlock *zfBase,
		const typename Geometry<FT,veclen,soalen,compress12>::FourSpinorBlock *tbBase,
		const typename Geometry<FT,veclen,soalen,compress12>::FourSpinorBlock *tfBase,
		typename Geometry<FT,veclen,soalen,compress12>::FourSpinorBlock *oBase,
		const typename Geometry<FT,veclen,soalen,compress12>::SU3MatrixBlock  *gBase,
		const int xbOffs[veclen],
		const int xfOffs[veclen],
		const int ybOffs[veclen],
		const int yfOffs[veclen],
		const int offs[veclen],
		const int gOffs[veclen],
		const int siprefdist1,
		const int siprefdist2,
		const int siprefdist3,
		const int siprefdist4,
		const int gprefdist,
		const int pfyOffs[veclen],
		const typename Geometry<FT,veclen,soalen,compress12>::FourSpinorBlock *pfBase2,
		const typename Geometry<FT,veclen,soalen,compress12>::FourSpinorBlock *pfBase3,
		const typename Geometry<FT,veclen,soalen,compress12>::FourSpinorBlock *pfBase4,
		const unsigned int accumulate[8],
		const FT coeff_s,
		const FT coeff_t_f,
		const FT coeff_t_b
)
{
  // BASE CASE TEMPLATE. Do nothing for now. Define this in dslash_generated_c.h later
  fprintf(stderr, "Generic veclen and soalen not supported yet.\n");
  abort();
}
#else 
/*
 *
 * 
 *  DSLASH_PLUS_VEC Specialized TEMPLATE
 *
 */

#if !defined(FPTYPE)
#error FTYPE not defined
#endif

#if !defined(VEC)
#error VLEN not defined
#endif


template<> 
inline void
dslash_plus_vec<FPTYPE,VEC,SOA,COMPRESS12>(    
		const Geometry<FPTYPE,VEC,SOA,COMPRESS12>::FourSpinorBlock *xyBase,
		const Geometry<FPTYPE,VEC,SOA,COMPRESS12>::FourSpinorBlock *zbBase,
		const Geometry<FPTYPE,VEC,SOA,COMPRESS12>::FourSpinorBlock *zfBase,
		const Geometry<FPTYPE,VEC,SOA,COMPRESS12>::FourSpinorBlock *tbBase,
		const Geometry<FPTYPE,VEC,SOA,COMPRESS12>::FourSpinorBlock *tfBase,
		Geometry<FPTYPE,VEC,SOA,COMPRESS12>::FourSpinorBlock *oBase,
		const Geometry<FPTYPE,VEC,SOA,COMPRESS12>::SU3MatrixBlock  *gBase,
		const int xbOffs[VEC],
		const int xfOffs[VEC],
		const int ybOffs[VEC],
		const int yfOffs[VEC],
		const int offs[VEC],
		const int gOffs[VEC],
		const int siprefdist1,
		const int siprefdist2,
		const int siprefdist3,
		const int siprefdist4,
		const int gprefdist,
		const int pfyOffs[VEC],
		const Geometry<FPTYPE,VEC,SOA,COMPRESS12>::FourSpinorBlock *pfBase2,
		const Geometry<FPTYPE,VEC,SOA,COMPRESS12>::FourSpinorBlock *pfBase3,
		const Geometry<FPTYPE,VEC,SOA,COMPRESS12>::FourSpinorBlock *pfBase4,
		const unsigned int accumulate[8],
		const FPTYPE coeff_s,
		const FPTYPE coeff_t_f,
		const FPTYPE coeff_t_b
)
{
#include INCLUDE_FILE_VAR(qphix/avx512/generated/dslash_plus_body_,FPTYPE,VEC,SOA,COMPRESS_SUFFIX)
}
#endif



/* 
 * 
 * DSLASH MINUS
 *
 */
#if !defined(SOA) && !defined(COMPRESS12) && !defined(COMPRESS_SUFFIX)
/*
 *
 * DSLASH_MINUS_VEC
 * Standard template if SOA, COMPRESS12, COMPRESS SUFFIX are all undefined
 *
 */
template<typename FT, int veclen, int soalen, bool compress12> 
inline void
dslash_minus_vec(    
		const typename Geometry<FT,veclen,soalen,compress12>::FourSpinorBlock *xyBase,
		const typename Geometry<FT,veclen,soalen,compress12>::FourSpinorBlock *zbBase,
		const typename Geometry<FT,veclen,soalen,compress12>::FourSpinorBlock *zfBase,
		const typename Geometry<FT,veclen,soalen,compress12>::FourSpinorBlock *tbBase,
		const typename Geometry<FT,veclen,soalen,compress12>::FourSpinorBlock *tfBase,
		typename Geometry<FT,veclen,soalen,compress12>::FourSpinorBlock *oBase,
		const typename Geometry<FT,veclen,soalen,compress12>::SU3MatrixBlock  *gBase,
		const int xbOffs[veclen],
		const int xfOffs[veclen],
		const int ybOffs[veclen],
		const int yfOffs[veclen],
		const int offs[veclen],
		const int gOffs[veclen],
		const int siprefdist1,
		const int siprefdist2,
		const int siprefdist3,
		const int siprefdist4,
		const int gprefdist,
		const int pfyOffs[veclen],
		const typename Geometry<FT,veclen,soalen,compress12>::FourSpinorBlock *pfBase2,
		const typename Geometry<FT,veclen,soalen,compress12>::FourSpinorBlock *pfBase3,
		const typename Geometry<FT,veclen,soalen,compress12>::FourSpinorBlock *pfBase4,
		const unsigned int accumulate[8],
		const FT coeff_s,
		const FT coeff_t_f,
		const FT coeff_t_b
)
{
  // BASE CASE TEMPLATE. Do nothing for now. Define this in dslash_generated_c.h later
  fprintf(stderr, "Generic veclen and soalen not supported yet.\n");
  abort();
}
#else 
/*
 *
 *
 *  DSLASH_MINUS_VEC Specialized TEMPLATE
 *
 */
template<> 
inline void
dslash_minus_vec<FPTYPE,VEC,SOA,COMPRESS12>(    
		const Geometry<FPTYPE,VEC,SOA,COMPRESS12>::FourSpinorBlock *xyBase,
		const Geometry<FPTYPE,VEC,SOA,COMPRESS12>::FourSpinorBlock *zbBase,
		const Geometry<FPTYPE,VEC,SOA,COMPRESS12>::FourSpinorBlock *zfBase,
		const Geometry<FPTYPE,VEC,SOA,COMPRESS12>::FourSpinorBlock *tbBase,
		const Geometry<FPTYPE,VEC,SOA,COMPRESS12>::FourSpinorBlock *tfBase,
		Geometry<FPTYPE,VEC,SOA,COMPRESS12>::FourSpinorBlock *oBase,
		const Geometry<FPTYPE,VEC,SOA,COMPRESS12>::SU3MatrixBlock  *gBase,
		const int xbOffs[VEC],
		const int xfOffs[VEC],
		const int ybOffs[VEC],
		const int yfOffs[VEC],
		const int offs[VEC],
		const int gOffs[VEC],
		const int siprefdist1,
		const int siprefdist2,
		const int siprefdist3,
		const int siprefdist4,
		const int gprefdist,
		const int pfyOffs[VEC],
		const Geometry<FPTYPE,VEC,SOA,COMPRESS12>::FourSpinorBlock *pfBase2,
		const Geometry<FPTYPE,VEC,SOA,COMPRESS12>::FourSpinorBlock *pfBase3,
		const Geometry<FPTYPE,VEC,SOA,COMPRESS12>::FourSpinorBlock *pfBase4,
		const unsigned int accumulate[8],
		const FPTYPE coeff_s,
		const FPTYPE coeff_t_f,
		const FPTYPE coeff_t_b
)
{
#include INCLUDE_FILE_VAR(qphix/avx512/generated/dslash_minus_body_,FPTYPE,VEC,SOA,COMPRESS_SUFFIX)
}
#endif



 /* 
 * 
 * DSLASH ACHIMBDPSI PLUS
 *
 */
#if !defined(SOA) && !defined(COMPRESS12) && !defined(COMPRESS_SUFFIX)

/*
 * DSLASH ACHIMBDPSI PLUS
 * Standard template if SOA, COMPRESS12, COMPRESS SUFFIX are all undefined
 *
 */
template<typename FT, int veclen, int soalen, bool compress12> 
inline void
dslash_achimbdpsi_plus_vec(    
		 const typename Geometry<FT,veclen,soalen,compress12>::FourSpinorBlock *xyBase,
		 const typename Geometry<FT,veclen,soalen,compress12>::FourSpinorBlock *zbBase,
		 const typename Geometry<FT,veclen,soalen,compress12>::FourSpinorBlock *zfBase,
		 const typename Geometry<FT,veclen,soalen,compress12>::FourSpinorBlock *tbBase,
		 const typename Geometry<FT,veclen,soalen,compress12>::FourSpinorBlock *tfBase,
		 const typename Geometry<FT,veclen,soalen,compress12>::FourSpinorBlock *chiBase,
		 typename Geometry<FT,veclen,soalen,compress12>::FourSpinorBlock *oBase,
		 const typename Geometry<FT,veclen,soalen,compress12>::SU3MatrixBlock  *gBase,
		 const int xbOffs[veclen],
		 const int xfOffs[veclen],
		 const int ybOffs[veclen],
		 const int yfOffs[veclen],
		 const int offs[veclen],
		 const int gOffs[veclen],
		 const int siprefdist1,
		 const int siprefdist2,
		 const int siprefdist3,
		 const int siprefdist4,
		 const int chiprefdist, 
		 const int gprefdist,
		 const int pfyOffs[veclen],
		 const typename Geometry<FT,veclen,soalen,compress12>::FourSpinorBlock *pfBase2,
		 const typename Geometry<FT,veclen,soalen,compress12>::FourSpinorBlock *pfBase3,
		 const typename Geometry<FT,veclen,soalen,compress12>::FourSpinorBlock *pfBase4,
		 const typename Geometry<FT,veclen,soalen,compress12>::FourSpinorBlock *pfBaseChi,
		 const FT alpha, 
		 const FT coeff_s,
		 const FT coeff_t_f,
		 const FT coeff_t_b,
		 const unsigned int accumulate[8])
{
  // BASE CASE TEMPLATE. Do nothing for now. Define this in dslash_generated_c.h later
  fprintf(stderr, "Generic veclen and soalen not supported yet.\n");
  abort();
}
#else 
/*
 *
 *
 *  Dslash AChi Minus BDPsi Plus Body: Specialized TEMPLATE
 *
 */
template<>
inline void
dslash_achimbdpsi_plus_vec<FPTYPE,VEC,SOA,COMPRESS12>(    
		 const Geometry<FPTYPE,VEC,SOA,COMPRESS12>::FourSpinorBlock *xyBase,
		 const Geometry<FPTYPE,VEC,SOA,COMPRESS12>::FourSpinorBlock *zbBase,
		 const Geometry<FPTYPE,VEC,SOA,COMPRESS12>::FourSpinorBlock *zfBase,
		 const Geometry<FPTYPE,VEC,SOA,COMPRESS12>::FourSpinorBlock *tbBase,
		 const Geometry<FPTYPE,VEC,SOA,COMPRESS12>::FourSpinorBlock *tfBase,
		 const Geometry<FPTYPE,VEC,SOA,COMPRESS12>::FourSpinorBlock *chiBase,
		 Geometry<FPTYPE,VEC,SOA,COMPRESS12>::FourSpinorBlock *oBase,
		 const Geometry<FPTYPE,VEC,SOA,COMPRESS12>::SU3MatrixBlock  *gBase,
		 const int xbOffs[VEC],
		 const int xfOffs[VEC],
		 const int ybOffs[VEC],
		 const int yfOffs[VEC],
		 const int offs[VEC],
		 const int gOffs[VEC],
		 const int siprefdist1,
		 const int siprefdist2,
		 const int siprefdist3,
		 const int siprefdist4,
		 const int chiprefdist, 
		 const int gprefdist,
		 const int pfyOffs[VEC],
		 const Geometry<FPTYPE,VEC,SOA,COMPRESS12>::FourSpinorBlock *pfBase2,
		 const Geometry<FPTYPE,VEC,SOA,COMPRESS12>::FourSpinorBlock *pfBase3,
		 const Geometry<FPTYPE,VEC,SOA,COMPRESS12>::FourSpinorBlock *pfBase4,
		 const Geometry<FPTYPE,VEC,SOA,COMPRESS12>::FourSpinorBlock *pfBaseChi,
		 const FPTYPE alpha, 
		 const FPTYPE coeff_s,
		 const FPTYPE coeff_t_f,
		 const FPTYPE coeff_t_b,
		 const unsigned int accumulate[8])
{
#include INCLUDE_FILE_VAR(qphix/avx512/generated/dslash_achimbdpsi_plus_body_,FPTYPE,VEC,SOA,COMPRESS_SUFFIX)
}
#endif

/* 
 * 
 * DSLASH ACHIMBDPSI MINUS
 *
 */
#if !defined(SOA) && !defined(COMPRESS12) && !defined(COMPRESS_SUFFIX)

/*
 * DSLASH ACHIMBDPSI MINUS
 * Standard template if SOA, COMPRESS12, COMPRESS SUFFIX are all undefined
 *
 */
template<typename FT, int veclen, int soalen, bool compress12> 
inline void
dslash_achimbdpsi_minus_vec(    
		 const typename Geometry<FT,veclen,soalen,compress12>::FourSpinorBlock *xyBase,
		 const typename Geometry<FT,veclen,soalen,compress12>::FourSpinorBlock *zbBase,
		 const typename Geometry<FT,veclen,soalen,compress12>::FourSpinorBlock *zfBase,
		 const typename Geometry<FT,veclen,soalen,compress12>::FourSpinorBlock *tbBase,
		 const typename Geometry<FT,veclen,soalen,compress12>::FourSpinorBlock *tfBase,
		 const typename Geometry<FT,veclen,soalen,compress12>::FourSpinorBlock *chiBase,
		 typename Geometry<FT,veclen,soalen,compress12>::FourSpinorBlock *oBase,
		 const typename Geometry<FT,veclen,soalen,compress12>::SU3MatrixBlock  *gBase,
		 const int xbOffs[veclen],
		 const int xfOffs[veclen],
		 const int ybOffs[veclen],
		 const int yfOffs[veclen],
		 const int offs[veclen],
		 const int gOffs[veclen],
		 const int siprefdist1,
		 const int siprefdist2,
		 const int siprefdist3,
		 const int siprefdist4,
		 const int chiprefdist, 
		 const int gprefdist,
		 const int pfyOffs[veclen],
		 const typename Geometry<FT,veclen,soalen,compress12>::FourSpinorBlock *pfBase2,
		 const typename Geometry<FT,veclen,soalen,compress12>::FourSpinorBlock *pfBase3,
		 const typename Geometry<FT,veclen,soalen,compress12>::FourSpinorBlock *pfBase4,
		 const typename Geometry<FT,veclen,soalen,compress12>::FourSpinorBlock *pfBaseChi,
		 const FT alpha, 
		 const FT coeff_s,
		 const FT coeff_t_f,
		 const FT coeff_t_b,
		 const unsigned int accumulate[8])
{
  // BASE CASE TEMPLATE. Do nothing for now. Define this in dslash_generated_c.h later
  fprintf(stderr, "Generic veclen and soalen not supported yet.\n");
  abort();
}
#else 
/*
 *
 *
 *  Dslash AChi Minus BDPsi MINUS Body: Specialized TEMPLATE
 *
 */
template<>
inline void
dslash_achimbdpsi_minus_vec<FPTYPE,VEC,SOA,COMPRESS12>(    
		 const Geometry<FPTYPE,VEC,SOA,COMPRESS12>::FourSpinorBlock *xyBase,
		 const Geometry<FPTYPE,VEC,SOA,COMPRESS12>::FourSpinorBlock *zbBase,
		 const Geometry<FPTYPE,VEC,SOA,COMPRESS12>::FourSpinorBlock *zfBase,
		 const Geometry<FPTYPE,VEC,SOA,COMPRESS12>::FourSpinorBlock *tbBase,
		 const Geometry<FPTYPE,VEC,SOA,COMPRESS12>::FourSpinorBlock *tfBase,
		 const Geometry<FPTYPE,VEC,SOA,COMPRESS12>::FourSpinorBlock *chiBase,
		 Geometry<FPTYPE,VEC,SOA,COMPRESS12>::FourSpinorBlock *oBase,
		 const Geometry<FPTYPE,VEC,SOA,COMPRESS12>::SU3MatrixBlock  *gBase,
		 const int xbOffs[VEC],
		 const int xfOffs[VEC],
		 const int ybOffs[VEC],
		 const int yfOffs[VEC],
		 const int offs[VEC],
		 const int gOffs[VEC],
		 const int siprefdist1,
		 const int siprefdist2,
		 const int siprefdist3,
		 const int siprefdist4,
		 const int chiprefdist, 
		 const int gprefdist,
		 const int pfyOffs[VEC],
		 const Geometry<FPTYPE,VEC,SOA,COMPRESS12>::FourSpinorBlock *pfBase2,
		 const Geometry<FPTYPE,VEC,SOA,COMPRESS12>::FourSpinorBlock *pfBase3,
		 const Geometry<FPTYPE,VEC,SOA,COMPRESS12>::FourSpinorBlock *pfBase4,
		 const Geometry<FPTYPE,VEC,SOA,COMPRESS12>::FourSpinorBlock *pfBaseChi,
		 const FPTYPE alpha, 
		 const FPTYPE coeff_s,
		 const FPTYPE coeff_t_f,
		 const FPTYPE coeff_t_b,
		 const unsigned int accumulate[8])
{
#include INCLUDE_FILE_VAR(qphix/avx512/generated/dslash_achimbdpsi_minus_body_,FPTYPE,VEC,SOA,COMPRESS_SUFFIX)
}
#endif

#ifdef QPHIX_DO_COMMS


#if !defined(SOA) && !defined(COMPRESS12) && !defined(COMPRESS_SUFFIX)
template<typename FT, int veclen, int soalen, bool compress12>
inline void  
face_proj_dir_plus(
		const typename Geometry<FT,veclen,soalen,compress12>::FourSpinorBlock *xyBase,
		const int offs[veclen],
		const int si_prefdist,
		FT *outbuf,
		const int hsprefdist,
		unsigned int mask,
		int dir)
{
  // BASE CASE TEMPLATE. Do nothing for now. Define this in dslash_generated_c.h later
  fprintf(stderr, "Generic veclen and soalen not supported yet.\n");
  abort();
}
#else
template<>
inline void  
face_proj_dir_plus<FPTYPE,VEC,SOA,COMPRESS12>(
		const typename Geometry<FPTYPE,VEC,SOA,COMPRESS12>::FourSpinorBlock *xyBase,
		const int offs[VEC],
		const int si_prefdist,
		FPTYPE *outbuf,
		const int hsprefdist,
		unsigned int mask,
		int dir)
{
  if(dir == 0) {
#include INCLUDE_FILE_VAR(qphix/avx512/generated/dslash_face_pack_to_back_X_plus_,FPTYPE,VEC,SOA,)
  }
  else if(dir == 1) {
#include INCLUDE_FILE_VAR(qphix/avx512/generated/dslash_face_pack_to_forw_X_plus_,FPTYPE,VEC,SOA,)
  }
  else if(dir == 2) {
#include INCLUDE_FILE_VAR(qphix/avx512/generated/dslash_face_pack_to_back_Y_plus_,FPTYPE,VEC,SOA,)
  }
  else if(dir == 3) {
#include INCLUDE_FILE_VAR(qphix/avx512/generated/dslash_face_pack_to_forw_Y_plus_,FPTYPE,VEC,SOA,)
  }
  else if(dir == 4) {
#include INCLUDE_FILE_VAR(qphix/avx512/generated/dslash_face_pack_to_back_Z_plus_,FPTYPE,VEC,SOA,)
  }
  else if(dir == 5) {
#include INCLUDE_FILE_VAR(qphix/avx512/generated/dslash_face_pack_to_forw_Z_plus_,FPTYPE,VEC,SOA,)
  }
  else if(dir == 6) {
#include INCLUDE_FILE_VAR(qphix/avx512/generated/dslash_face_pack_to_back_T_plus_,FPTYPE,VEC,SOA,)
  }
  else if(dir == 7) {
#include INCLUDE_FILE_VAR(qphix/avx512/generated/dslash_face_pack_to_forw_T_plus_,FPTYPE,VEC,SOA,)
  }
  else {
    printf("Invalid dir for pack boundary\n");
    exit(1);
  }
}

#endif

#if !defined(SOA) && !defined(COMPRESS12) && !defined(COMPRESS_SUFFIX)
template<typename FT, int veclen, int soalen, bool compress12>
inline void  
face_proj_dir_minus(
		const typename Geometry<FT,veclen,soalen,compress12>::FourSpinorBlock *xyBase,
		const int offs[veclen],
		const int si_prefdist,
		FT *outbuf,
		const int hsprefdist,
		unsigned int mask,
		int dir)
{
  // BASE CASE TEMPLATE. Do nothing for now. Define this in dslash_generated_c.h later
  fprintf(stderr, "Generic veclen and soalen not supported yet.\n");
  abort();
}
#else 
template<>
inline void  
face_proj_dir_minus<FPTYPE,VEC,SOA,COMPRESS12>(
		const typename Geometry<FPTYPE,VEC,SOA,COMPRESS12>::FourSpinorBlock *xyBase,
		const int offs[VEC],
		const int si_prefdist,
		FPTYPE *outbuf,
		const int hsprefdist,
		unsigned int mask,
		int dir)
{
  if(dir == 0) {
#include INCLUDE_FILE_VAR(qphix/avx512/generated/dslash_face_pack_to_back_X_minus_,FPTYPE,VEC,SOA,)
  }
  else if(dir == 1) {
#include INCLUDE_FILE_VAR(qphix/avx512/generated/dslash_face_pack_to_forw_X_minus_,FPTYPE,VEC,SOA,)
  }
  else if(dir == 2) {
#include INCLUDE_FILE_VAR(qphix/avx512/generated/dslash_face_pack_to_back_Y_minus_,FPTYPE,VEC,SOA,)
  }
  else if(dir == 3) {
#include INCLUDE_FILE_VAR(qphix/avx512/generated/dslash_face_pack_to_forw_Y_minus_,FPTYPE,VEC,SOA,)
  }
  else if(dir == 4) {
#include INCLUDE_FILE_VAR(qphix/avx512/generated/dslash_face_pack_to_back_Z_minus_,FPTYPE,VEC,SOA,)
  }
  else if(dir == 5) {
#include INCLUDE_FILE_VAR(qphix/avx512/generated/dslash_face_pack_to_forw_Z_minus_,FPTYPE,VEC,SOA,)
  }
  else if(dir == 6) {
#include INCLUDE_FILE_VAR(qphix/avx512/generated/dslash_face_pack_to_back_T_minus_,FPTYPE,VEC,SOA,)
  }
  else if(dir == 7) {
#include INCLUDE_FILE_VAR(qphix/avx512/generated/dslash_face_pack_to_forw_T_minus_,FPTYPE,VEC,SOA,)
  }
  else {
    printf("Invalid dir for pack boundary\n");
    exit(1);
  }
}

#endif

#if !defined(SOA) && !defined(COMPRESS12) && !defined(COMPRESS_SUFFIX)
template<typename FT, int veclen, int soalen, bool compress12>
inline void 
  face_finish_dir_plus(
		       const FT *inbuf,
		       const typename Geometry<FT,veclen,soalen,compress12>::SU3MatrixBlock *gBase,
		       typename Geometry<FT,veclen,soalen,compress12>::FourSpinorBlock *oBase,
		       const int gOffs[veclen],
		       const int offs[veclen],
		       const int hsprefdist,
		       const int gprefdist,
		       const int soprefdist,
		       const FT beta,
		       unsigned int mask,
		       int dir) 
{
  // BASE CASE TEMPLATE. Do nothing for now. Define this in dslash_generated_c.h later
  fprintf(stderr, "Generic veclen and soalen not supported yet.\n");
  abort();
}
#else 
template<>
inline void 
face_finish_dir_plus<FPTYPE,VEC,SOA,COMPRESS12>(
		       const FPTYPE *inbuf,
		       const Geometry<FPTYPE,VEC,SOA,COMPRESS12>::SU3MatrixBlock *gBase,
		       Geometry<FPTYPE,VEC,SOA,COMPRESS12>::FourSpinorBlock *oBase,
		       const int gOffs[VEC],
		       const int offs[VEC],
		       const int hsprefdist,
		       const int gprefdist,
		       const int soprefdist,
		       const FPTYPE beta,
		       unsigned int mask,
		       int dir) 
{
  if(dir == 0) {
#include INCLUDE_FILE_VAR(qphix/avx512/generated/dslash_face_unpack_from_back_X_plus_,FPTYPE,VEC,SOA,COMPRESS_SUFFIX)
  }
  else if(dir == 1) {
#include INCLUDE_FILE_VAR(qphix/avx512/generated/dslash_face_unpack_from_forw_X_plus_,FPTYPE,VEC,SOA,COMPRESS_SUFFIX)
  }
  else if(dir == 2) {
#include INCLUDE_FILE_VAR(qphix/avx512/generated/dslash_face_unpack_from_back_Y_plus_,FPTYPE,VEC,SOA,COMPRESS_SUFFIX)
  }
  else if(dir == 3) {
#include INCLUDE_FILE_VAR(qphix/avx512/generated/dslash_face_unpack_from_forw_Y_plus_,FPTYPE,VEC,SOA,COMPRESS_SUFFIX)
  }
  else if(dir == 4) {
#include INCLUDE_FILE_VAR(qphix/avx512/generated/dslash_face_unpack_from_back_Z_plus_,FPTYPE,VEC,SOA,COMPRESS_SUFFIX)
  }
  else if(dir == 5) {
#include INCLUDE_FILE_VAR(qphix/avx512/generated/dslash_face_unpack_from_forw_Z_plus_,FPTYPE,VEC,SOA,COMPRESS_SUFFIX)
  }
  else if(dir == 6) {
#include INCLUDE_FILE_VAR(qphix/avx512/generated/dslash_face_unpack_from_back_T_plus_,FPTYPE,VEC,SOA,COMPRESS_SUFFIX)
  }
  else if(dir == 7) {
#include INCLUDE_FILE_VAR(qphix/avx512/generated/dslash_face_unpack_from_forw_T_plus_,FPTYPE,VEC,SOA,COMPRESS_SUFFIX)
  }
  else {
    printf("Invalid dir for unpack boundary\n");
    exit(1);
  }
}


#endif

#if !defined(SOA) && !defined(COMPRESS12) && !defined(COMPRESS_SUFFIX)
template<typename FT, int veclen, int soalen, bool compress12>
inline void 
  face_finish_dir_minus(
		       const FT *inbuf,
		       const typename Geometry<FT,veclen,soalen,compress12>::SU3MatrixBlock *gBase,
		       typename Geometry<FT,veclen,soalen,compress12>::FourSpinorBlock *oBase,
		       const int gOffs[veclen],
		       const int offs[veclen],
		       const int hsprefdist,
		       const int gprefdist,
		       const int soprefdist,
		       const FT beta,
		       unsigned int mask,
		       int dir) 
{
  // BASE CASE TEMPLATE. Do nothing for now. Define this in dslash_generated_c.h later
  fprintf(stderr, "Generic veclen and soalen not supported yet.\n");
  abort();
}
#else 
template<>
inline void 
face_finish_dir_minus<FPTYPE,VEC,SOA,COMPRESS12>(
						 const FPTYPE *inbuf,
						 const Geometry<FPTYPE,VEC,SOA,COMPRESS12>::SU3MatrixBlock *gBase,
						 Geometry<FPTYPE,VEC,SOA,COMPRESS12>::FourSpinorBlock *oBase,
						 const int gOffs[VEC],
						 const int offs[VEC],
						 const int hsprefdist,
						 const int gprefdist,
						 const int soprefdist,
						 const FPTYPE beta,
						 unsigned int mask,
						 int dir) 
{
  if(dir == 0) {
#include INCLUDE_FILE_VAR(qphix/avx512/generated/dslash_face_unpack_from_back_X_minus_,FPTYPE,VEC,SOA,COMPRESS_SUFFIX)
  }
  else if(dir == 1) {
#include INCLUDE_FILE_VAR(qphix/avx512/generated/dslash_face_unpack_from_forw_X_minus_,FPTYPE,VEC,SOA,COMPRESS_SUFFIX)
  }
  else if(dir == 2) {
#include INCLUDE_FILE_VAR(qphix/avx512/generated/dslash_face_unpack_from_back_Y_minus_,FPTYPE,VEC,SOA,COMPRESS_SUFFIX)
  }
  else if(dir == 3) {
#include INCLUDE_FILE_VAR(qphix/avx512/generated/dslash_face_unpack_from_forw_Y_minus_,FPTYPE,VEC,SOA,COMPRESS_SUFFIX)
  }
  else if(dir == 4) {
#include INCLUDE_FILE_VAR(qphix/avx512/generated/dslash_face_unpack_from_back_Z_minus_,FPTYPE,VEC,SOA,COMPRESS_SUFFIX)
  }
  else if(dir == 5) {
#include INCLUDE_FILE_VAR(qphix/avx512/generated/dslash_face_unpack_from_forw_Z_minus_,FPTYPE,VEC,SOA,COMPRESS_SUFFIX)
  }
  else if(dir == 6) {
#include INCLUDE_FILE_VAR(qphix/avx512/generated/dslash_face_unpack_from_back_T_minus_,FPTYPE,VEC,SOA,COMPRESS_SUFFIX)
  }
  else if(dir == 7) {
#include INCLUDE_FILE_VAR(qphix/avx512/generated/dslash_face_unpack_from_forw_T_minus_,FPTYPE,VEC,SOA,COMPRESS_SUFFIX)
  }
  else {
    printf("Invalid dir for unpack boundary\n");
    exit(1);
  }
}
#endif

#endif // QPHIX_DO_COMMS (outer)



