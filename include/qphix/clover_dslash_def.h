#pragma once

#include <qphix/dslash_utils.h>
#include <qphix/geometry.h>
#include <qphix/comm.h>

namespace QPhiX
{

/* Specialize - Dslash of float */
/* To note the clover block will already include anisotropy and probably also
 * boundaries */
template <typename FT,
          int veclen,
          int soalen,
          bool compress12> /* The teplate is the floating point type */
class ClovDslash
{
 public:
  // Pack gauges...
  typedef typename Geometry<FT, veclen, soalen, compress12>::SU3MatrixBlock
      SU3MatrixBlock;
  typedef typename Geometry<FT, veclen, soalen, compress12>::FourSpinorBlock
      FourSpinorBlock;
  typedef typename Geometry<FT, veclen, soalen, compress12>::TwoSpinorBlock
      TwoSpinorBlock;
  typedef typename Geometry<FT, veclen, soalen, compress12>::CloverBlock CloverBlock;

  ClovDslash(Geometry<FT, veclen, soalen, compress12> *geom_,
             double t_boundary_,
             double dslash_aniso_s_,
             double dslash_aniso_t_,
             bool use_tbc_[4] = nullptr,
             double tbc_phases_[4][2] = nullptr);

  /* Destructor */
  ~ClovDslash();

  /* Apply the operator */
  void dslash(FourSpinorBlock *res,
              const FourSpinorBlock *psi,
              const SU3MatrixBlock *u,
              const CloverBlock *invclov,
              int isign,
              int cb);

  void dslashAChiMinusBDPsi(FourSpinorBlock *res,
                            const FourSpinorBlock *psi,
                            const FourSpinorBlock *chi,
                            const SU3MatrixBlock *u,
                            const CloverBlock *clov,
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

  // Hide Free Constructor
  ClovDslash();
  void Dyz(int tid,
           const FourSpinorBlock *psi,
           FourSpinorBlock *res,
           const SU3MatrixBlock *u,
           const CloverBlock *invclov,
           bool const is_plus,
           int cb);

  void DyzAChiMinusBDPsi(int tid,
                         const FourSpinorBlock *psi,
                         const FourSpinorBlock *chi,
                         FourSpinorBlock *res,
                         const SU3MatrixBlock *u,
                         const CloverBlock *clov,
                         double beta,
                         bool const is_plus,
                         int cb);

  void DPsi(const SU3MatrixBlock *u,
            const CloverBlock *invclov,
            const FourSpinorBlock *psi_in,
            FourSpinorBlock *res_out,
            bool const is_plus,
            int cb);

  void DPsiAChiMinusBDPsi(const SU3MatrixBlock *u,
                          const CloverBlock *clov,
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

  //  RECEIVE AND COMPLETE FACE
  void completeFaceDir(int tid,
                       const FT *psi,
                       FourSpinorBlock *res,
                       const SU3MatrixBlock *u,
                       const CloverBlock *invclov,
                       double beta,
                       int cb,
                       int dir,
                       int fb,
                       bool const is_plus);

  //  RECEIVE AND COMPLETE FACE
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
#include "qphix/clover_dslash_body.h"
#ifdef QPHIX_DO_COMMS
// Disable comms for nokw
#include "qphix/clov_face.h"
#endif
