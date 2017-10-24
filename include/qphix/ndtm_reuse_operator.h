#pragma once

#include "qphix/linearOp.h"
#include "qphix/tm_clov_dslash_def.h"
#include "qphix/print_utils.h"
#include <cstdlib>
#include <memory>

namespace QPhiX
{

template <typename FT, int veclen, int soalen, bool compress12>
class EvenOddNDTMWilsonReuseOperator
    : public TwoFlavEvenOddLinearOperator<FT, veclen, soalen, compress12>
{
 public:
  typedef typename Geometry<FT, veclen, soalen, compress12>::FullCloverBlock
      FullCloverBlock;
  typedef typename Geometry<FT, veclen, soalen, compress12>::FourSpinorBlock
      FourSpinorBlock;
  typedef typename Geometry<FT, veclen, soalen, compress12>::SU3MatrixBlock
      SU3MatrixBlock;

  EvenOddNDTMWilsonReuseOperator(const double mass,
                                 const double mu,
                                 const double epsilon,
                                 SU3MatrixBlock *u[2],
                                 Geometry<FT, veclen, soalen, compress12> *geom,
                                 double t_boundary,
                                 double aniso_coeff_s,
                                 double aniso_coeff_t,
                                 bool mu_plus,
                                 bool use_tbc_[4] = nullptr,
                                 double tbc_phases_[4][2] = nullptr)
      : mass(mass), mass_factor_alpha(4 + mass),
        mass_factor_beta(0.25 / mass_factor_alpha),
        mu(-2.0 * mu * (mu_plus ? 1 : -1) * mass_factor_beta),
        mu_inv(1.0 / (1.0 + Mu * Mu)), epsilon(epsilon),
        D(new TMClovDslash<FT, veclen, soalen, compress12>(
            geom, t_boundary, aniso_coeff_s, aniso_coeff_t, use_tbc_, tbc_phases_)),
        tmp({D->getGeometry().allocCBFourSpinor(),
             D->getGeometry().allocCBFourSpinor()})
  {
    this->u[0] = u[0];
    this->u[1] = u[1];
  }

  ~EvenOddNDTMWilsonReuseOperator()
  {
    Geometry<FT, veclen, soalen, compress12> &geom = D->getGeometry();
    for (int f : {0, 1}) {
      geom.free(tmp[f]);
    }
  }

  void operator()(FourSpinorBlock *res[2],
                  const FourSpinorBlock *const in[2],
                  int isign) override
  {
    D->two_flaw_dslash(tmp, in, u[1], mu, mu_inv, isign, 1);
    D->two_flaw_dslashAChiMinusBDPsi(
        res, tmp, in, u[0], mu, mass_factor_beta, isign, 0);
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

  Geometry<FT, veclen, soalen, compress12> &getGeometry()
  {
    return D->getGeometry();
  }

 private:
  double mass;
  double mass_factor_alpha, mass_factor_beta;
  double mu, mu_inv;
  double epsilon;

  std::unique_ptr<TMClovDslash<FT, veclen, soalen, compress12>> D;

  SU3MatrixBlock *u[2];

  FourSpinorBlock *tmp[2];

}; // Class
}; // Namespace
