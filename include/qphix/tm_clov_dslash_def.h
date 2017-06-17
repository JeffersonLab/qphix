#pragma once

#include <qphix/dslash_utils.h>
#include <qphix/geometry.h>
#include <qphix/comm.h>

namespace QPhiX
{

template <typename FT, int veclen, int soalen, bool compress12>
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
               double dslash_aniso_t_,
               bool use_tbc_[4] = nullptr,
               double tbc_phases_[4][2] = nullptr);

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

  bool use_tbc[4] = {false, false, false, false};
  FT tbc_phases[4][2] = {{rep<FT, double>(1.0), rep<FT, double>(0.0)},
                         {rep<FT, double>(1.0), rep<FT, double>(0.0)},
                         {rep<FT, double>(1.0), rep<FT, double>(0.0)},
                         {rep<FT, double>(1.0), rep<FT, double>(0.0)}};

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

  TMClovDslash(); // Hide Free Constructor

  void Dyz(int tid,
           const FourSpinorBlock *psi,
           FourSpinorBlock *res,
           const SU3MatrixBlock *u,
           const FullCloverBlock *invclov,
           bool const is_plus,
           int cb);

  void DyzAChiMinusBDPsi(int tid,
                         const FourSpinorBlock *psi,
                         const FourSpinorBlock *chi,
                         FourSpinorBlock *res,
                         const SU3MatrixBlock *u,
                         const FullCloverBlock *clov,
                         double beta,
                         bool const is_plus,
                         int cb);

  void DPsi(const SU3MatrixBlock *u,
            const FullCloverBlock *invclov,
            const FourSpinorBlock *psi_in,
            FourSpinorBlock *res_out,
            bool const is_plus,
            int cb);

  void DPsiAChiMinusBDPsi(const SU3MatrixBlock *u,
                          const FullCloverBlock *clov,
                          const FourSpinorBlock *psi_in,
                          const FourSpinorBlock *chi,
                          FourSpinorBlock *res_out,
                          double beta,
                          bool const is_plus,
                          int cb);

// DISABLE COMMS FOR NOW
#ifdef QPHIX_DO_COMMS
  void packFaceDir(int tid,
                   const FourSpinorBlock *psi,
                   FT *res,
                   int cb,
                   int dir,
                   int fb,
                   bool const is_plus);

  void completeFaceDir(int tid,
                       const FT *psi,
                       FourSpinorBlock *res,
                       const SU3MatrixBlock *u,
                       const FullCloverBlock *invclov,
                       double beta,
                       int cb,
                       int dir,
                       int fb,
                       bool const is_plus);

  void completeFaceDirAChiMBDPsi(int tid,
                                 const FT *psi,
                                 FourSpinorBlock *res,
                                 const SU3MatrixBlock *u,
                                 const double beta,
                                 int cb,
                                 int dir,
                                 int fb,
                                 bool const is_plus);
#endif

}; // Class

} // Namespace

#include "qphix/tm_clov_dslash_body.h"

#ifdef QPHIX_DO_COMMS
#include "qphix/tm_clov_dslash_face.h"
#endif
