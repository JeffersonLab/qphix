#pragma once

#include "qphix/linearOp.h"
#include "qphix/tm_clov_dslash_def.h"
#include "qphix/print_utils.h"
#include <cstdlib>
#include <memory>

namespace QPhiX
{

/**
  Two flavor twisted mass clover operator, leveraging existing code.

  The \ref TMClovDslash provides all the methods used in \ref
  EvenOddTMCloverOperator. The degenerate case is completely covered. An
  implementation of the two-flavor %Dslash using the one-flavor %Dslash has been
  added to \ref TMClovDslash. This operator now provides a wrapper around those
  functions, _reusing_ most of the existing code.

  \author Martin Ueding <dev@martin-ueding.de>
  */
template <typename FT, int veclen, int soalen, bool compress12>
class EvenOddNDTMCloverReuseOperator
    : public TwoFlavEvenOddLinearOperator<FT, veclen, soalen, compress12>
{
 public:
  typedef typename Geometry<FT, veclen, soalen, compress12>::FullCloverBlock
      FullCloverBlock;
  typedef typename Geometry<FT, veclen, soalen, compress12>::FourSpinorBlock
      FourSpinorBlock;
  typedef typename Geometry<FT, veclen, soalen, compress12>::SU3MatrixBlock
      SU3MatrixBlock;

  /**
    \param[in] u Gauge fields, one for each checkerboard index.

    \param[in] clov Two flavor parts of the odd-odd term. The first index is
    the flavor index, the second index is for the hermitian conjugation, just
    like \p clov of \ref dslashAChiMinusBDPsi or \p invclov of \ref
    two_flav_dslash.

    \param[in] invclov The four blocks of \f$ (A^{-1})_{ff'} \f$. The first
    (slowest) index iterates through the flavor index \f$ f \f$, it coincides
    with the resulting spinor flavor index. The second index is \f$ f' \f$
    which coincides with the input spinor flavor indices. The last (fastest)
    index is the same as for the \ref dslash function: It is the hermitian
    conjugate of the inverse clover block.
    */
  EvenOddNDTMCloverReuseOperator(SU3MatrixBlock *u_[2],
                                 FullCloverBlock *clov_[2][2],
                                 FullCloverBlock *invclov_[2][2][2],
                                 double epsilon,
                                 Geometry<FT, veclen, soalen, compress12> *geom_,
                                 double t_boundary,
                                 double aniso_coeff_s,
                                 double aniso_coeff_t,
                                 bool use_tbc_[4] = nullptr,
                                 double tbc_phases_[4][2] = nullptr)
      : EvenOddNDTMCloverReuseOperator(
            geom_, t_boundary, aniso_coeff_s, aniso_coeff_t, use_tbc_, tbc_phases_)
  {
    this->epsilon = epsilon;
    setFields(u_, clov_, invclov_);
  }

  EvenOddNDTMCloverReuseOperator(Geometry<FT, veclen, soalen, compress12> *geom_,
                                 double t_boundary,
                                 double aniso_coeff_s,
                                 double aniso_coeff_t,
                                 bool use_tbc_[4] = nullptr,
                                 double tbc_phases_[4][2] = nullptr)
      : D(new TMClovDslash<FT, veclen, soalen, compress12>(
            geom_, t_boundary, aniso_coeff_s, aniso_coeff_t, use_tbc_, tbc_phases_)),
        tmp{D->getGeometry().allocCBFourSpinor(),
            D->getGeometry().allocCBFourSpinor()}
  {
  }

  ~EvenOddNDTMCloverReuseOperator()
  {
    Geometry<FT, veclen, soalen, compress12> &geom = D->getGeometry();
    for (int f : {0, 1}) {
      geom.free(tmp[f]);
    }
  }

  void setFields(SU3MatrixBlock *u_[2],
                 FullCloverBlock *clov_[2][2],
                 FullCloverBlock *invclov_[2][2][2])
  {
    for (int pm : {0, 1}) {
      u[pm] = u_[pm];
      for (int f1 : {0, 1}) {
        clov[pm][f1] = clov_[pm][f1];
        for (int f2 : {0, 1}) {
          invclov[pm][f1][f2] = invclov_[pm][f1][f2];
        }
      }
    }
  }

  void operator()(FourSpinorBlock *res[2],
                  const FourSpinorBlock *const in[2],
                  int isign) override
  {
    double beta = 0.25;

    D->two_flav_dslash(tmp, in, u[0], invclov, isign, 0);
    D->two_flav_achimbdpsi(res, tmp, in, u[1], clov, beta, epsilon, isign, 1);
  }

  Geometry<FT, veclen, soalen, compress12> &getGeometry()
  {
    return D->getGeometry();
  }

  void M_offdiag(FourSpinorBlock *res[2],
                 const FourSpinorBlock *const in[2],
                 int isign,
                 int target_cb) const override
  {
    masterPrintf("M_ee_inv not yet implemented for this operator\n");
    std::abort();
  }

  // M_ee_inv is always Hermitian so no need for isign?
  // for wilson it is the identity, for clover it is hermitian, for TWM it is gamma_5?
  void M_ee_inv(FourSpinorBlock *res[2],
                const FourSpinorBlock *const in[2],
                int isign) const override
  {
    masterPrintf("M_ee_inv not yet implemented for this operator\n");
    std::abort();
  }

 private:
  double Mass;
  double epsilon;

  // XXX Use a simple member, then there is no need for a virtual function
  // call within `operator()`.
  std::unique_ptr<TMClovDslash<FT, veclen, soalen, compress12>> D;

  SU3MatrixBlock const *u[2];
  FullCloverBlock const *clov[2][2];
  FullCloverBlock const *invclov[2][2][2];

  FourSpinorBlock *tmp[2];

}; // Class
}; // Namespace
