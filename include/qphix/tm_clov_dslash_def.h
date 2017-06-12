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
              const FullCloverBlock *const invclov[2],
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

  /**
    Clover %Dslash for the non-degenerate twisted mass case.

    This implementation of the non-degenerate case uses the degenerate
    twisted clover code. The existing and tested code is called four
    times for each element in the flavor matrix. The flavor structure is
    decomposed into single-flavor %Dslash applications where the existing
    code is used. This way little new code has to be written, the chance
    of introducing errors is reduced.

    The existing code computes \f$ \chi := A^{-1} D \psi \f$, where \f$
    \psi \f$ is a single flavor spinor. \f$ A = 4 + m \pm \mathrm i \mu
    \gamma_5 + c_\mathrm{SW} T_\mathrm{ee} \f$ is the even-even term with
    mass \f$ m \f$, twisted mass \f$ \mu \f$, and clover term \f$ T \f$
    defined on the even lattice sites only. The hopping matrix \f$ D \f$
    is just the Wilson %Dslash operator.

    In the degenerate case, it should actually be \f$ \mathrm i \mu
    \gamma_5 \tau^3 \f$, where this is a two-flavor expression. It is
    diagonal in flavor-space, therefore the determinant of the whole
    expression will eventually be just the product of the determinants of
    the two blocks. Exploding the diagonal structure, one gets away with
    having a single-flavor operator where the sign of the \f$ \mu \f$
    term can be changed by hermitian conjugation.

    The non-degenerate twisted mass case has an additional \f$ - \epsilon
    \tau^1 \f$ in the even-even term \f$ A \f$. Therefore it becomes
    densely populated in flavor-space. The decomposition of the
    determinant is no longer possible, therefore an additional operator
    is needed. Also the inverse of \f$ A \f$ is not the inverse of the
    blocks. The interface to this function accepts \f$ A^{-1} \f$ from
    outside, so it is the caller's responsibility to make the appropriate
    inversion of \f$ A \f$.

    Using the existing one-flavor %Dslash, this function will compute the
    two-flavor %Dslash like this:
    \f[

    \begin{pmatrix}
        \chi_\mathrm u \\ \chi_\mathrm d
    \end{pmatrix}
    :=
    \begin{pmatrix}
        (A^{-1})_\mathrm{uu} &
        (A^{-1})_\mathrm{ud} \\
        (A^{-1})_\mathrm{du} &
        (A^{-1})_\mathrm{dd}
    \end{pmatrix}
    D
    \begin{pmatrix}
        \psi_\mathrm u \\ \psi_\mathrm d
    \end{pmatrix}
    =
    \begin{pmatrix}
        (A^{-1})_\mathrm{uu} D \psi_\mathrm u +
        (A^{-1})_\mathrm{ud} D \psi_\mathrm d \\
        (A^{-1})_\mathrm{du} D \psi_\mathrm u +
        (A^{-1})_\mathrm{dd} D \psi_\mathrm d
    \end{pmatrix}
    \,.
    \f]

    There will be four calls to \ref dslash with different arguments. The
    two resulting up-spinors will be accumulated, same with the two
    down-spinors. The result are two separate spinors.

    \test By setting the flavor-off-diagonal parts of \f$ A^{-1} \f$ to a
    zero matrix, the degenerate case is recovered again. This can then be
    compared to two invocations of \ref dslash. This will make sure that
    the decomposition of the flavor structure is done correctly. Also it
    verifies that the accumulation of the intermediate results within the
    flavor matrix multiplication is done correctly.

    \test Once the real non-degenerate implementation has been finished in
    \ref NDTMClovDslash, the two implementations can be compared to each
    other. The results should be very close numerically, though a
    different order in the summations might lead to small rounding
    deviations.

    \param[out] res Resulting two flavor spinor. The slowest index is the
    flavor.

    \param[in] psi Input spinors. The slowest index is the flavor.

    \param[in] u Gauge field.

    \param[in] invclov The four blocks of \f$ (A^{-1})_{ff'} \f$. The
    first (slowest) index iterates through the flavor index \f$ f \f$, it
    coincides with the resulting spinor flavor index. The second index is
    \f$ f' \f$ which coincides with the input spinor flavor indices. The
    last (fastest) index is the same as for the \ref dslash function: It
    is the hermitian conjugate of the inverse clover block.

    @param[in] isign The value \f$ +1 \f$ gives the operator as-is, the
    value \f$ -1 \f$ applies the hermitian conjugate operator.

    @param[in] cb The checkerbord index of the target spinor, value 0 or
    1.

    \todo Figure out whether the hermitian conjugated FullCloverBlock
    should also conjugate in flavor space.

    \todo We want to have \f$ \gamma_5 D \tau^3 \f$ as the overall
    operator. It does not seem obvious how this translate to this operator
    here, do we want \f$ \gamma_5 A^{-1} D \tau^3 \f$ or rather \f$ A^{-1}
    \gamma_5 D \tau^3 \f$?

    \author Martin Ueding <dev@martin-ueding.de>

    */
  template <typename Spinor1>
  void two_flav_dslash(FourSpinorBlock *res[2],
                       Spinor1 *const psi[2],
                       const SU3MatrixBlock *u,
                       const FullCloverBlock *const invclov[2][2][2],
                       const int isign,
                       const int cb);

  /**
    Clover \f$ A \chi - b D \psi \f$ for the non-degenerate twisted mass
    case.

    Most aspects of this function are similar to \ref two_flav_dslash.
    Please first read the documentation there.

    This function will compute
    \f[
    \begin{pmatrix}
    \Psi_\mathrm u \\ \Psi_\mathrm d
    \end{pmatrix}
    :=
    \begin{pmatrix}
    A_\mathrm{uu} &
    \epsilon \\
    \epsilon &
    A_\mathrm{dd}
    \end{pmatrix}
    \begin{pmatrix}
    \chi_\mathrm u \\ \chi_\mathrm d
    \end{pmatrix}
    - b D
    \begin{pmatrix}
    \psi_\mathrm u \\ \psi_\mathrm d
    \end{pmatrix}
    \,.
    \f]

    In order to leverage the \ref dslashAChiMinusBDPsi implementation,
    this is rewritten as
    \f[
    \begin{pmatrix}
    \Psi_\mathrm u \\ \Psi_\mathrm d
    \end{pmatrix}
    :=
    \underbrace{
    \begin{pmatrix}
    A_\mathrm{uu} \chi_\mathrm u - b D \psi_\mathrm u \\
    A_\mathrm{dd} \chi_\mathrm d - b D \psi_\mathrm d
    \end{pmatrix}
    }_\text{achimdpsi}
    +
    \begin{pmatrix}
    \epsilon \chi_\mathrm d \\
    \epsilon \chi_\mathrm u
    \end{pmatrix}
    \,.
    \f]

    \param[in] chi Second spinor field \f$ \chi \f$.
    \param[in] clov Two flavor parts of the odd-odd term. The first index
    is the flavor index, the second index is for the heritian conjugation,
    just like \p clov of \ref dslashAChiMinusBDPsi or \p invclov of \ref
    two_flav_dslash.
    \param[in] beta The \f$ b \f$ (or \f$ \beta \f$) coefficient.
    \param[in] epsilon Twisted mass splitting \f$ \epsilon \f$ (or \f$
    \mu_\delta \f$).

    Other parameters are exactly like in \ref two_flav_dslash and are not
    explained again here.

    \author Martin Ueding <dev@martin-ueding.de>
    */
  template <typename Spinor1, typename Spinor2>
  void two_flav_achimbdpsi(FourSpinorBlock *res[2],
                           Spinor1 *const chi[2],
                           Spinor2 *const psi[2],
                           const SU3MatrixBlock *u,
                           const FullCloverBlock *clov[2][2],
                           const double beta,
                           const double epsilon,
                           const int isign,
                           const int cb);

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

#include "qphix/tm_clov_dslash_body.h"

#ifdef QPHIX_DO_COMMS
#include "qphix/tm_clov_dslash_face.h"
#endif

#endif // Include guard
