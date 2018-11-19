#pragma once

#include "qphix/linearOp.h"
#include "qphix/dslash_def.h"
#include "qphix/tm_clov_dslash_def.h"
#include "qphix/print_utils.h"
#include <cstdlib>
#include <memory>

namespace QPhiX
{

/**
  Two flavor twisted mass clover operator, leveraging existing code.

  The \ref TMClovDslash and \ref Dslash provide the funtionality required 
  to implement the two-flavour EO-preconditioned operator with a single
  function of \ref TMClovDslash implementing the two-flavour inverse
  clover term multiplication.

  \author Martin Ueding <dev@martin-ueding.de>
  \author Bartosz Kostrzewa <bartosz_kostrzewa@fastmail.com>
  */
template <typename FT, int veclen, int soalen, bool compress12>
class EvenOddNDTMCloverReuseOperator
    : public TwoFlavEvenOddLinearOperator<FT, veclen, soalen, compress12>
{
 public:
  typedef typename Geometry<FT, veclen, soalen, compress12>::FullCloverBlock
      FullCloverBlock;
  typedef typename Geometry<FT, veclen, soalen, compress12>::CloverBlock
      CloverBlock;
  typedef typename Geometry<FT, veclen, soalen, compress12>::FourSpinorBlock
      FourSpinorBlock;
  typedef typename Geometry<FT, veclen, soalen, compress12>::SU3MatrixBlock
      SU3MatrixBlock;

  /**
    \param[in] u_ Gauge fields, one for each checkerboard index.

    \param[in] clov_ Single Wilson-clover style clover term which will be applied
    in two_flav_AChiMinusBDPsi to all flavour components with the twisted quark
    mass and the mass splitting applied subsequently.

    \param[in] invclov_odiag_ The flavour off-diagonal inverse clover term
    including the mass splitting.

    \param[in] invclov_ FullCloverBlock inverse clover term for the flavour-diagonal
    contribution as applied by the Twisted Clover Dslash.
    */
  EvenOddNDTMCloverReuseOperator(double TwistedMass_,
                                 double Epsilon_,
                                 SU3MatrixBlock *u_[2],
                                 CloverBlock *clov_,
                                 CloverBlock *invclov_odiag_,
                                 FullCloverBlock *invclov_[2],
                                 Geometry<FT, veclen, soalen, compress12> *geom_,
                                 double t_boundary,
                                 double aniso_coeff_s,
                                 double aniso_coeff_t,
                                 bool use_tbc_[4] = nullptr,
                                 double tbc_phases_[4][2] = nullptr)
      : EvenOddNDTMCloverReuseOperator(
            TwistedMass_, Epsilon_, 
            geom_, t_boundary, aniso_coeff_s, aniso_coeff_t, use_tbc_, tbc_phases_)
  {
    setFields(u_, clov_, invclov_odiag_, invclov_);
  }

  EvenOddNDTMCloverReuseOperator(double TwistedMass_,
                                 double Epsilon_,
                                 Geometry<FT, veclen, soalen, compress12> *geom_,
                                 double t_boundary,
                                 double aniso_coeff_s,
                                 double aniso_coeff_t,
                                 bool use_tbc_[4] = nullptr,
                                 double tbc_phases_[4][2] = nullptr)
      : Dtmcl(new TMClovDslash<FT, veclen, soalen, compress12>(
            geom_, t_boundary, aniso_coeff_s, aniso_coeff_t, use_tbc_, tbc_phases_,
            TwistedMass_)),
        Dw(new Dslash<FT, veclen, soalen, compress12>(
             geom_, t_boundary, aniso_coeff_s, aniso_coeff_t, use_tbc_, tbc_phases_)),
        tmp{Dtmcl->getGeometry().allocCBFourSpinor(),
            Dtmcl->getGeometry().allocCBFourSpinor()},
        odiag_tmp{Dtmcl->getGeometry().allocCBFourSpinor(),
                  Dtmcl->getGeometry().allocCBFourSpinor()},
        Epsilon(Epsilon_)
  {
  }

  ~EvenOddNDTMCloverReuseOperator()
  {
    Geometry<FT, veclen, soalen, compress12> &geom = Dtmcl->getGeometry();
    for (int fl : {0, 1}) {
      geom.free(tmp[fl]);
      geom.free(odiag_tmp[fl]);
    }
  }

  void setFields(SU3MatrixBlock *u_[2],
                 CloverBlock *clov_,
                 CloverBlock *invclov_odiag_,
                 FullCloverBlock *invclov_[2])
  {
    for (int cb : {0, 1}) {
      u[cb] = u_[cb];
    }
    clov = clov_;
    invclov_odiag = invclov_odiag_;
    for (int fl : {0, 1}) {
      invclov[fl] = invclov_[fl];
    }
  }

  void operator()(FourSpinorBlock *const res[2],
                  const FourSpinorBlock *const in[2],
                  int isign, 
                  int target_cb = 1) override
  {
    double beta = 0.25;
    const int other_cb = 1-target_cb;
 
    for(int fl : {0, 1}){
      Dw->dslash(odiag_tmp[fl], in[fl], u[other_cb], isign, other_cb);
    }
    Dtmcl->two_flav_inverse_clover_term(tmp, odiag_tmp, invclov, invclov_odiag, isign);
    Dtmcl->two_flav_AChiMinusBDPsi(res, tmp, in, u[target_cb], clov, beta, Epsilon, isign, target_cb);
  }

  Geometry<FT, veclen, soalen, compress12> &getGeometry()
  {
    return Dtmcl->getGeometry();
  }

  void M_offdiag(FourSpinorBlock *res[2],
                 const FourSpinorBlock *const in[2],
                 int isign,
                 int target_cb) const override
  {
    masterPrintf("M_ee_inv not yet implemented for this operator\n");
    std::abort();
  }

  // M_diag_inv is always Hermitian so no need for isign?
  // for wilson it is the identity, for clover it is hermitian, for TWM it is gamma_5?
  void M_diag_inv(FourSpinorBlock *res[2],
                const FourSpinorBlock *const in[2],
                int isign) const override
  {
    masterPrintf("M_ee_inv not yet implemented for this operator\n");
    std::abort();
  }

 private:
  double Epsilon;

  std::unique_ptr<TMClovDslash<FT, veclen, soalen, compress12>> Dtmcl;
  std::unique_ptr<Dslash<FT, veclen, soalen, compress12>> Dw;

  SU3MatrixBlock const *u[2];
  CloverBlock const *clov;
  CloverBlock const *invclov_odiag;
  FullCloverBlock const *invclov[2];

  FourSpinorBlock *tmp[2];
  FourSpinorBlock *odiag_tmp[2];

}; // Class
}; // Namespace
