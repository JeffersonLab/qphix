#ifndef QPHIX_TM_DSLASH_DEF_H
#define QPHIX_TM_DSLASH_DEF_H

#include "qphix/dslash_utils.h"
#include "qphix/geometry.h"
#include "qphix/comm.h"

namespace QPhiX
{
template <typename FT, int veclen, int soalen, bool compress12>
class TMDslash
{

 public:
  typedef typename Geometry<FT, veclen, soalen, compress12>::SU3MatrixBlock
      SU3MatrixBlock;
  typedef typename Geometry<FT, veclen, soalen, compress12>::FourSpinorBlock
      FourSpinorBlock;
  typedef typename Geometry<FT, veclen, soalen, compress12>::TwoSpinorBlock
      TwoSpinorBlock;

  TMDslash(Geometry<FT, veclen, soalen, compress12> *geom_,
           double t_boundary_,
           double aniso_coeff_S_,
           double aniso_coeff_T_,
           double Mass_,
           double TwistedMass_);

  ~TMDslash();

  /* Apply the operator */
  void dslash(FourSpinorBlock *res,
              const FourSpinorBlock *psi,
              const SU3MatrixBlock *u,
              int isign,
              int cb);

  void dslashAChiMinusBDPsi(FourSpinorBlock *res,
                            const FourSpinorBlock *psi,
                            const FourSpinorBlock *chi,
                            const SU3MatrixBlock *u,
                            double alpha,
                            double beta,
                            int isign,
                            int cb);

  Geometry<FT, veclen, soalen, compress12> &getGeometry(void) { return (*s); }

 private:
  Geometry<FT, veclen, soalen, compress12> *s;
  Comms<FT, veclen, soalen, compress12> *comms;

  int log2veclen;
  int log2soalen;
  void init();

  const int n_threads_per_core;

  const int By;
  const int Bz;
  const int NCores;
  const int Sy;
  const int Sz;
  const int PadXY;
  const int PadXYZ;
  const int MinCt;

  // Antiperiodic Boundary
  const double t_boundary;
  bool amIPtMin;
  bool amIPtMax;

  // Anisotropy
  double aniso_coeff_S;
  double aniso_coeff_T;

  // Wilson & Twisted Mass Parameters
  double Mass;
  double TwistedMass;

  // Derived Twisted Mass Parameters
  double derived_mu; // = \mu / \alpha
  double derived_mu_inv; // = \alpha / ( \mu^2 + \alpha^2 )

  Barrier *gBar;
  Barrier ***barriers;

  BlockPhase *block_info;
  int *tmpspc_all; // Space for YZ offsets etc

  // Masks for comms in X & Y
  unsigned int xbmask_x0_xodd[2];
  unsigned int xfmask_xn_xodd[2];
  unsigned int ybmask_y0;
  unsigned int yfmask_yn;

  TMDslash(); // Hide Free Constructor

  void TMDyzPlus(int tid,
                 const FourSpinorBlock *psi,
                 FourSpinorBlock *res,
                 const SU3MatrixBlock *u,
                 int cb);

  void TMDyzMinus(int tid,
                  const FourSpinorBlock *psi,
                  FourSpinorBlock *res,
                  const SU3MatrixBlock *u,
                  int cb);

  void TMDyzPlusAChiMinusBDPsi(int tid,
                               const FourSpinorBlock *psi,
                               const FourSpinorBlock *chi,
                               FourSpinorBlock *res,
                               const SU3MatrixBlock *u,
                               double alpha,
                               double beta,
                               int cb);

  void TMDyzMinusAChiMinusBDPsi(int tid,
                                const FourSpinorBlock *psi,
                                const FourSpinorBlock *chi,
                                FourSpinorBlock *res,
                                const SU3MatrixBlock *u,
                                double alpha,
                                double beta,
                                int cb);

  void TMDPsiPlus(const SU3MatrixBlock *u,
                  const FourSpinorBlock *psi_in,
                  FourSpinorBlock *res_out,
                  int cb);

  void TMDPsiMinus(const SU3MatrixBlock *u,
                   const FourSpinorBlock *psi_in,
                   FourSpinorBlock *res_out,
                   int cb);

  void TMDPsiPlusAChiMinusBDPsi(const SU3MatrixBlock *u,
                                const FourSpinorBlock *psi_in,
                                const FourSpinorBlock *chi,
                                FourSpinorBlock *res_out,
                                double alpha,
                                double beta,
                                int cb);

  void TMDPsiMinusAChiMinusBDPsi(const SU3MatrixBlock *u,
                                 const FourSpinorBlock *psi_in,
                                 const FourSpinorBlock *chi,
                                 FourSpinorBlock *re_out,
                                 double alpha,
                                 double beta,
                                 int cb);

// packTMFaceDir: same as standard wilson.
#ifdef QPHIX_QMP_COMMS
  // PACK FACE FOR SENDING
  void packTMFaceDir(int tid,
                     const FourSpinorBlock *psi,
                     FT *res,
                     int cb,
                     int dir,
                     int fb,
                     int isPlus);

  //  RECEIVE AND COMPLETE FACE
  void completeTMFaceDir(int tid,
                         const FT *psi,
                         FourSpinorBlock *res,
                         const SU3MatrixBlock *u,
                         const double beta,
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
#endif
