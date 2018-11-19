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
  typedef typename Geometry<FT, veclen, soalen, compress12>::CloverBlock
      CloverBlock;

  // Constructor
  TMClovDslash(Geometry<FT, veclen, soalen, compress12> *geom_,
               double t_boundary_,
               double dslash_aniso_s_,
               double dslash_aniso_t_,
               bool use_tbc_[4] = nullptr,
               double tbc_phases_[4][2] = nullptr,
               double const prec_mass_rho = 0.0);

  // Destructor
  ~TMClovDslash();

  /**
    Apply the dslash operator, which here will be \f$ A^{-1} D \psi \f$.

    @param[in] invclov \f$ A^{-1} \f$. The two elements of the array shall be
    the original and the hermitian conjugate.
    */
  void dslash(FourSpinorBlock *res,
              const FourSpinorBlock *psi,
              const SU3MatrixBlock *u,
              const FullCloverBlock *const invclov[2],
              int isign,
              int cb,
              int fl = 0);

  // Apply the A*chi-b*D*psi operator
  // w/ A = alpha 1 + c_SW T = clov,
  // passed (and suitably packed) externally
  // The twisted quark mass (plus any preconditioning msss) 
  // is applied via prec_mass_rho
  void dslashAChiMinusBDPsi(FourSpinorBlock *res,
                            const FourSpinorBlock *psi,
                            const FourSpinorBlock *chi,
                            const SU3MatrixBlock *u,
                            const CloverBlock *clov,
                            const double beta,
                            int isign,
                            int cb,
                            int fl = 0);

  /**
    Clover \f$ A \chi - b D \psi \f$ for the non-degenerate twisted mass
    case.

    \param[in] chi Second spinor field \f$ \chi \f$.
    \param[in] clov Clover term not including the twisted mass, that is
                    implemented via the pre-conditioning mass rho
    \param[in] beta The \f$ b \f$ (or \f$ \beta \f$) coefficient.
    \param[in] epsilon Twisted mass splitting \f$ \epsilon \f$ (or \f$
    \mu_\delta \f$).

    Other parameters are exactly like in \ref two_flav_dslash and are not
    explained again here.

    \author Martin Ueding <dev@martin-ueding.de>
    \author Bartosz Kostrzewa <bartosz_kostrzewa@fastmail.com>
    */
  void two_flav_AChiMinusBDPsi(FourSpinorBlock *const res[2],
                               const FourSpinorBlock *const psi[2],
                               const FourSpinorBlock *const chi[2],
                               const SU3MatrixBlock *u,
                               const CloverBlock *clov,
                               const double beta,
                               const double epsilon,
                               const int isign,
                               const int cb);

  /* 
   * Applies the two-flavour inverse clover term from three fields. In the
   * two-flavour EO operator, the hopping matrix is applied to the two
   * flavours individually first, then this function is called.
   * \param[out] res two-flavour output
   * \param[in]  psi two-flavour spinor which should contain D_Wilson \chi
   * \param[in]  fcl inverse clover terms (including twisted mass) on the flavour
   *                 diagonal for the up and down flavours
   * \param[in]  clOffDiag epsilon / ( (alpha+clover)^2 + \mu^2 - \epsilon^2 )
   *                       inverse clover contribution on the flavour off-diagonal                
   */ 
  void two_flav_inverse_clover_term(FourSpinorBlock *const res[2],
                                    const FourSpinorBlock *const psi[2],
                                    const FullCloverBlock *const fcl[2],
                                    const CloverBlock     *const clOffDiag,
                                    int isign);

#ifdef __INTEL_COMPILER  
  void two_flav_AChiMinusBDPsi(FourSpinorBlock *const res[2],
                               FourSpinorBlock *const psi[2],
                               FourSpinorBlock *const chi[2],
                               const SU3MatrixBlock *u,
                               const CloverBlock *clov,
                               double beta,
                               double epsilon,
                               int isign,
                               int cb)
  {
    this->two_flav_AChiMinusBDPsi(res, 
                                  const_cast<const FourSpinorBlock * const *>(psi),
                                  const_cast<const FourSpinorBlock * const *>(chi),
                                  u, clov, beta, epsilon, isign, cb); 
  }

  void two_flav_AChiMinusBDPsi(FourSpinorBlock * const res[2],
                               FourSpinorBlock * psi[2],
                               FourSpinorBlock * chi[2],
                               const SU3MatrixBlock *u,
                               const CloverBlock *clov,
                               double beta,
                               double epsilon,
                               int isign,
                               int cb)
  {
    this->two_flav_AChiMinusBDPsi(res, 
                                  const_cast<const FourSpinorBlock * const *>(psi),
                                  const_cast<const FourSpinorBlock * const *>(chi),
                                  u, clov, beta, epsilon, isign, cb); 
  }
  
  void two_flav_AChiMinusBDPsi(FourSpinorBlock * const res[2],
                               FourSpinorBlock * psi[2],
                               const FourSpinorBlock * const chi[2],
                               const SU3MatrixBlock *u,
                               const CloverBlock *clov,
                               double beta,
                               double epsilon,
                               int isign,
                               int cb)
  {
    this->two_flav_AChiMinusBDPsi(res, 
                                  const_cast<const FourSpinorBlock * const *>(psi),
                                  const_cast<const FourSpinorBlock * const *>(chi),
                                  u, clov, beta, epsilon, isign, cb);
  }
  
  void two_flav_inverse_clover_term(FourSpinorBlock *const res[2],
                                    FourSpinorBlock *const psi[2],
                                    const FullCloverBlock *const fcl[2],
                                    const CloverBlock     *const clOffDiag,
                                    int isign)
  {
    this->two_flav_inverse_clover_term(res,
                                       const_cast<const FourSpinorBlock * const *>(psi),
                                       fcl,
                                       clOffDiag,
                                       isign);
  }
  
  void two_flav_inverse_clover_term(FourSpinorBlock *const res[2],
                                    FourSpinorBlock *const psi[2],
                                    FullCloverBlock *const fcl[2],
                                    const CloverBlock     *const clOffDiag,
                                    int isign)
  {
    this->two_flav_inverse_clover_term(res,
                                       const_cast<const FourSpinorBlock * const *>(psi),
                                       const_cast<const FullCloverBlock * const *>(fcl),
                                       fclDn,
                                       clOffDiag,
                                       isign);
  }
#endif // __INTEL_COMPILER workaround 

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

  double const prec_mass_rho;

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
                         const CloverBlock *clov,
                         double beta,
                         bool const is_plus,
                         int cb,
                         int fl = 0);

  void DPsi(const SU3MatrixBlock *u,
            const FullCloverBlock *invclov,
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
                          int cb,
                          int fl = 0);

  void two_flav_inverse_clover_term_YZ(
      const int tid,
      FourSpinorBlock *const resUp,
      FourSpinorBlock *const resDn,
      const FourSpinorBlock *const psiUp,
      const FourSpinorBlock *const psiDn,
      const FullCloverBlock *const fclUp,
      const FullCloverBlock *const fclDn,
      const CloverBlock *const clOffDiag);

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
