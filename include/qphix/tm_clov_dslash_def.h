#ifndef QPHIX_TM_CLOVER_DSLASH_DEF_H
#define QPHIX_TM_CLOVER_DSLASH_DEF_H

#include <qphix/dslash_utils.h>
#include <qphix/geometry.h>
#include <qphix/comm.h>

namespace QPhiX
{

template <typename FT,
          int veclen,
          int soalen,
          bool compress12> /* The teplate is the floating point type */
class TMClovDslash
{

 public:
  typedef typename Geometry<FT, veclen, soalen, compress12>::SU3MatrixBlock
      SU3MatrixBlock;
  typedef typename Geometry<FT, veclen, soalen, compress12>::FourSpinorBlock
      FourSpinorBlock;
  typedef typename Geometry<FT, veclen, soalen, compress12>::TwoSpinorBlock
      TwoSpinorBlock;
  typedef typename Geometry<FT, veclen, soalen, compress12>::FullCloverBlock
      FullCloverBlock;

  // Constructor
  TMClovDslash(Geometry<FT, veclen, soalen, compress12> *geom_,
               double t_boundary_,
               double dslash_aniso_s_,
               double dslash_aniso_t_);

  // Destructor
  ~TMClovDslash();

  // Apply the dslash operator,
  // which here will be A^-1 * D * psi,
  // with A^-1 = invclov
  void dslash(FourSpinorBlock *res,
              const FourSpinorBlock *psi,
              const SU3MatrixBlock *u,
              const FullCloverBlock *invclov[2],
              int isign,
              int cb);

  // Apply the A*chi-b*D*psi operator
  // w/ A = alpha 1 +/- i mu gamma_5 + c_SW T = clov,
  // passed (and suitably packed) externally
  void dslashAChiMinusBDPsi(FourSpinorBlock *res,
                            const FourSpinorBlock *psi,
                            const FourSpinorBlock *chi,
                            const SU3MatrixBlock *u,
                            const FullCloverBlock *clov[2],
                            const double beta,
                            int isign,
                            int cb);

  void free(void *p);

  Geometry<FT, veclen, soalen, compress12> &getGeometry(void) { return (*s); }

 private:
  Geometry<FT, veclen, soalen, compress12> *s;
  Comms<FT, veclen, soalen, compress12> *comms;

  int log2veclen;
  int log2soalen;

  const int n_threads_per_core;
  const int By;
  const int Bz;
  const int NCores;
  const int Sy;
  const int Sz;
  const int PadXY;
  const int PadXYZ;
  const int MinCt;

  const double t_boundary;
  const double aniso_coeff_S;
  const double aniso_coeff_T;

  bool amIPtMin;
  bool amIPtMax;

  Barrier *gBar;
  Barrier ***barriers;

  BlockPhase *block_info;
  int *tmpspc_all; // Space for YZ offsets etc

  // Masks for comms in X & Y
  unsigned int xbmask_x0_xodd[2];
  unsigned int xfmask_xn_xodd[2];
  unsigned int ybmask_y0;
  unsigned int yfmask_yn;

  void init();
  TMClovDslash(); // Hide Free Constructor

  void DyzPlus(int tid,
               const FourSpinorBlock *psi,
               FourSpinorBlock *res,
               const SU3MatrixBlock *u,
               const FullCloverBlock *invclov,
               int cb);

  void DyzMinus(int tid,
                const FourSpinorBlock *psi,
                FourSpinorBlock *res,
                const SU3MatrixBlock *u,
                const FullCloverBlock *invclov,
                int cb);

  void DyzPlusAChiMinusBDPsi(int tid,
                             const FourSpinorBlock *psi,
                             const FourSpinorBlock *chi,
                             FourSpinorBlock *res,
                             const SU3MatrixBlock *u,
                             const FullCloverBlock *clov,
                             double beta,
                             int cb);

  void DyzMinusAChiMinusBDPsi(int tid,
                              const FourSpinorBlock *psi,
                              const FourSpinorBlock *chi,
                              FourSpinorBlock *res,
                              const SU3MatrixBlock *u,
                              const FullCloverBlock *clov,
                              double beta,
                              int cb);

  void DPsiPlus(const SU3MatrixBlock *u,
                const FullCloverBlock *invclov,
                const FourSpinorBlock *psi_in,
                FourSpinorBlock *res_out,
                int cb);

  void DPsiMinus(const SU3MatrixBlock *u,
                 const FullCloverBlock *invclov,
                 const FourSpinorBlock *psi_in,
                 FourSpinorBlock *res_out,
                 int cb);

  void DPsiPlusAChiMinusBDPsi(const SU3MatrixBlock *u,
                              const FullCloverBlock *clov,
                              const FourSpinorBlock *psi_in,
                              const FourSpinorBlock *chi,
                              FourSpinorBlock *res_out,
                              double beta,
                              int cb);

  void DPsiMinusAChiMinusBDPsi(const SU3MatrixBlock *u,
                               const FullCloverBlock *clov,
                               const FourSpinorBlock *psi_in,
                               const FourSpinorBlock *chi,
                               FourSpinorBlock *rea_out,
                               double beta,
                               int cb);

// DISABLE COMMS FOR NOW
#ifdef QPHIX_DO_COMMS
  void packFaceDir(int tid,
                   const FourSpinorBlock *psi,
                   FT *res,
                   int cb,
                   int dir,
                   int fb,
                   int isPlus);

  //  RECEIVE AND COMPLETE FACE
  void completeFaceDir(int tid,
                       const FT *psi,
                       FourSpinorBlock *res,
                       const SU3MatrixBlock *u,
                       const FullCloverBlock *invclov,
                       double beta,
                       int cb,
                       int dir,
                       int fb,
                       int isPlus);

  //  RECEIVE AND COMPLETE FACE
  void completeFaceDirAChiMBDPsi(int tid,
                                 const FT *psi,
                                 FourSpinorBlock *res,
                                 const SU3MatrixBlock *u,
                                 const double beta,
                                 int cb,
                                 int dir,
                                 int fb,
                                 int isPlus);
#endif

}; // Class

} // Namespace

#endif // Include guard
