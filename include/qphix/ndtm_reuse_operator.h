#pragma once

#include "qphix/linearOp.h"
#include "qphix/tm_dslash_def.h"
#include "qphix/dslash_def.h"
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
  typedef
      typename Geometry<FT, veclen, soalen, compress12>::FourSpinorBlock FourSpinorBlock;
  typedef
      typename Geometry<FT, veclen, soalen, compress12>::SU3MatrixBlock SU3MatrixBlock;

  EvenOddNDTMWilsonReuseOperator(const double mass,
                                 const double mu,
                                 const double epsilon,
                                 SU3MatrixBlock *u[2],
                                 Geometry<FT, veclen, soalen, compress12> *geom,
                                 double t_boundary,
                                 double aniso_coeff_s,
                                 double aniso_coeff_t,
                                 bool use_tbc_[4] = nullptr,
                                 double tbc_phases_[4][2] = nullptr)
      : mass(mass), mass_factor_alpha(4 + mass), mass_factor_beta(0.25), mu(mu),
        epsilon(epsilon),
        mu_inv(1.0 / ((4 + mass) * (4 + mass) + mu * mu - epsilon * epsilon)),
        Dtm(new TMDslash<FT, veclen, soalen, compress12>(geom,
                                                         t_boundary,
                                                         aniso_coeff_s,
                                                         aniso_coeff_t,
                                                         mass,
                                                         mu,
                                                         use_tbc_,
                                                         tbc_phases_)),
        Dw(new Dslash<FT, veclen, soalen, compress12>(
            geom, t_boundary, aniso_coeff_s, aniso_coeff_t, use_tbc_, tbc_phases_)),
        Dw_tmp{Dw->getGeometry().allocCBFourSpinor(),
               Dw->getGeometry().allocCBFourSpinor()},
        mu_p_eps_Dw_tmp{Dtm->getGeometry().allocCBFourSpinor(),
                        Dtm->getGeometry().allocCBFourSpinor()}
  {
    this->u[0] = u[0];
    this->u[1] = u[1];
  }

  ~EvenOddNDTMWilsonReuseOperator()
  {
    Geometry<FT, veclen, soalen, compress12> &geom = Dtm->getGeometry();
    for (int i = 0; i < 2; i++) {
      geom.free(Dw_tmp[i]);
      geom.free(mu_p_eps_Dw_tmp[i]);
    }
  }

  void operator()(FourSpinorBlock *const res[2],
                  FourSpinorBlock const *const in[2],
                  int isign,
                  int target_cb = 1) override
  {
    constexpr int up = 0;
    constexpr int dn = 1;

    const int other_cb = 1 - target_cb;

    /* for the non-degenerate case, we need to apply the
     * hopping matrix ("Wilson dslash" or Dslash::dslash()) rather than TMDslash::dslash()
     * which would also apply the twisted mass term right after */
    Dw->dslash(Dw_tmp[up], in[up], u[other_cb], isign, other_cb);
    Dw->dslash(Dw_tmp[dn], in[dn], u[other_cb], isign, other_cb);

    //// now we apply the inverse twisted mass term and the epsilon flavour cross term
    ////   (alpha - i*mu*gamma_5*tau3 + eps*tau1) / (alpha^2 + mu^2 - eps^2)
    //// hence the change of sign on the twisted mass
    //// ALSO NOTE: the sign of the off-diagonal contribution will be corrected
    //// by AChiMinusBDPsi below
    const double apimu[2] = {mass_factor_alpha * mu_inv, -isign * mu * mu_inv};
    two_flav_twisted_mass(apimu,
                          epsilon * mu_inv,
                          Dw_tmp,
                          mu_p_eps_Dw_tmp,
                          Dw->getGeometry(),
                          Dw->getGeometry().getNSIMT());

    ///* now the two-flavour AChiMinusBDPsi puts it all together
    // *  M_oo^{f=0} = (alpha + i*mu*gamma5 ) \chi^{f=0} - eps*\chi^{f=1} -
    // *               1/(4*(alpha^2 + mu^2 - eps^2)) D_oe *
    // *               ( (alpha - i*mu*gamma5) Deo \chi^{f=0} + eps*Deo \chi^{f=1} )
    // */
    Dtm->two_flav_AChiMinusBDPsi(res,
                                 mu_p_eps_Dw_tmp,
                                 in,
                                 u[target_cb],
                                 mass_factor_alpha,
                                 mass_factor_beta,
                                 epsilon,
                                 isign,
                                 target_cb);
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

  Geometry<FT, veclen, soalen, compress12> &getGeometry() { return Dtm->getGeometry(); }

 private:
  double mass;
  double mass_factor_alpha;
  double mass_factor_beta;
  double epsilon;
  double mu;
  double mu_inv;

  SU3MatrixBlock *u[2];

  // we need the twisted mass Dslash for AChiMinusBDPsi
  std::unique_ptr<TMDslash<FT, veclen, soalen, compress12>> Dtm;
  // and the Wilson Dslash to get the flavour-off-diagonal part right
  std::unique_ptr<Dslash<FT, veclen, soalen, compress12>> Dw;

  FourSpinorBlock *Dw_tmp[2];
  FourSpinorBlock *mu_p_eps_Dw_tmp[2];

}; // Class
}; // Namespace
