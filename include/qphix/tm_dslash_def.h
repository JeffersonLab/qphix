#pragma once

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
           double TwistedMass_,
           bool use_tbc_[4] = nullptr,
           double tbc_phases_[4][2] = nullptr);

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
                            int cb,
                            const double musign = 1.0);

  void two_flav_AChiMinusBDPsi(FourSpinorBlock * const res[2],
                               const FourSpinorBlock *const psi[2],
                               const FourSpinorBlock *const chi[2],
                               const SU3MatrixBlock *u,
                               double alpha,
                               double beta,
                               double epsilon,
                               int isign,
                               int cb);

#ifdef __INTEL_COMPILER  
  void two_flav_AChiMinusBDPsi(FourSpinorBlock * const res[2],
                               FourSpinorBlock *const psi[2],
                               FourSpinorBlock *const chi[2],
                               const SU3MatrixBlock *u,
                               double alpha,
                               double beta,
                               double epsilon,
                               int isign,
                               int cb)
  {
    this->two_flav_AChiMinusBDPsi(res, 
                                  const_cast<const FourSpinorBlock * const *>(psi),
                                  const_cast<const FourSpinorBlock * const *>(chi),
                                  u, alpha, beta, epsilon, isign, cb); 
  }

  void two_flav_AChiMinusBDPsi(FourSpinorBlock * const res[2],
                               FourSpinorBlock * psi[2],
                               FourSpinorBlock * chi[2],
                               const SU3MatrixBlock *u,
                               double alpha,
                               double beta,
                               double epsilon,
                               int isign,
                               int cb)
  {
    this->two_flav_AChiMinusBDPsi(res, 
                                  const_cast<const FourSpinorBlock * const *>(psi),
                                  const_cast<const FourSpinorBlock * const *>(chi),
                                  u, alpha, beta, epsilon, isign, cb); 
  }
  
  void two_flav_AChiMinusBDPsi(FourSpinorBlock * const res[2],
                               FourSpinorBlock * psi[2],
                               const FourSpinorBlock * const chi[2],
                               const SU3MatrixBlock *u,
                               double alpha,
                               double beta,
                               double epsilon,
                               int isign,
                               int cb)
  {
    this->two_flav_AChiMinusBDPsi(res, 
                                  const_cast<const FourSpinorBlock * const *>(psi),
                                  const_cast<const FourSpinorBlock * const *>(chi),
                                  u, alpha, beta, epsilon, isign, cb); 
  }
#endif

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

  bool use_tbc[4] = {false, false, false, false};
  FT tbc_phases[4][2] = {{rep<FT, double>(1.0), rep<FT, double>(0.0)},
                         {rep<FT, double>(1.0), rep<FT, double>(0.0)},
                         {rep<FT, double>(1.0), rep<FT, double>(0.0)},
                         {rep<FT, double>(1.0), rep<FT, double>(0.0)}};

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

  void TMDyz(int tid,
             const FourSpinorBlock *psi,
             FourSpinorBlock *res,
             const SU3MatrixBlock *u,
             bool const is_plus,
             int cb);

  void TMDyzAChiMinusBDPsi(int tid,
                           const FourSpinorBlock *psi,
                           const FourSpinorBlock *chi,
                           FourSpinorBlock *res,
                           const SU3MatrixBlock *u,
                           double alpha,
                           double beta,
                           bool const is_plus,
                           int cb,
                           const double musign = 1.0);

  void TMDPsi(const SU3MatrixBlock *u,
              const FourSpinorBlock *psi_in,
              FourSpinorBlock *res_out,
              bool const is_plus,
              int cb);

  void TMDPsiAChiMinusBDPsi(const SU3MatrixBlock *u,
                            const FourSpinorBlock *psi_in,
                            const FourSpinorBlock *chi,
                            FourSpinorBlock *res_out,
                            double alpha,
                            double beta,
                            bool const is_plus,
                            int cb,
                            const double musign = 1.0);

// packTMFaceDir: same as standard wilson.
#ifdef QPHIX_QMP_COMMS
  // PACK FACE FOR SENDING
  void packTMFaceDir(int tid,
                     const FourSpinorBlock *psi,
                     FT *res,
                     int cb,
                     int dir,
                     int fb,
                     bool const is_plus);

  //  RECEIVE AND COMPLETE FACE
  void completeTMFaceDir(int tid,
                         const FT *psi,
                         FourSpinorBlock *res,
                         const SU3MatrixBlock *u,
                         const double beta,
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
#include "qphix/tm_dslash_body.h"

#ifdef QPHIX_QMP_COMMS
#include "qphix/tm_dslash_face.h"
#endif
