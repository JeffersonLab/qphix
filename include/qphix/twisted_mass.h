#pragma once

#include "qphix/linearOp.h"
#include "qphix/tm_dslash_def.h"

#include "qphix/print_utils.h"
#include <cstdlib>
namespace QPhiX
{

template <typename FT, int veclen, int soalen, bool compress12>
class EvenOddTMWilsonOperator
    : public EvenOddLinearOperator<FT, veclen, soalen, compress12>
{

 public:
  typedef typename Geometry<FT, veclen, soalen, compress12>::FourSpinorBlock
      FourSpinorBlock;
  typedef typename Geometry<FT, veclen, soalen, compress12>::SU3MatrixBlock
      SU3MatrixBlock;

  EvenOddTMWilsonOperator(const double Mass_,
                          const double TwistedMass_,
                          SU3MatrixBlock *u_[2],
                          Geometry<FT, veclen, soalen, compress12> *geom_,
                          double t_boundary_,
                          double aniso_fac_s_,
                          double aniso_fac_t_,
                          bool use_tbc_[4] = nullptr,
                          double tbc_phases_[4][2] = nullptr)
      : Mass(Mass_), TwistedMass(TwistedMass_), mass_factor_alpha(4.0 + Mass),
        mass_factor_beta(0.25),
        D(new TMDslash<FT, veclen, soalen, compress12>(geom_,
                                                       t_boundary_,
                                                       aniso_fac_s_,
                                                       aniso_fac_t_,
                                                       Mass,
                                                       TwistedMass,
                                                       use_tbc_,
                                                       tbc_phases_))
  {
    u[0] = u_[0];
    u[1] = u_[1];
    Geometry<FT, veclen, soalen, compress12> &geom = D->getGeometry();
    tmp = (FourSpinorBlock *)geom.allocCBFourSpinor();
  }

  ~EvenOddTMWilsonOperator()
  {
    Geometry<FT, veclen, soalen, compress12> &geom = D->getGeometry();
    geom.free(tmp);
    delete D;
  }

  Geometry<FT, veclen, soalen, compress12> &getGeometry() override
  {
    return D->getGeometry();
  }

  // TODO: Add detailed description of the operator!
  // This is the even-odd preconditioned (odd-odd) fermion matrix
  void operator()(FourSpinorBlock *res, // result spinor field
                  const FourSpinorBlock *in, // input spinor field
                  int isign, // non-conjugate = 1, hermitian conjugate = -1
                  int target_cb = 1) const override
  {
    int source_cb = 1 - target_cb;
    D->dslash(tmp, in, u[source_cb], isign, source_cb);
    D->dslashAChiMinusBDPsi(res,
                            tmp,
                            in,
                            u[target_cb],
                            mass_factor_alpha,
                            mass_factor_beta,
                            isign,
                            target_cb);
  }

  inline void M_offdiag(FourSpinorBlock *res,
                         FourSpinorBlock const *in,
                         int isign,
                         int target_cb) const override {
    masterPrintf("M_offdiag not yet implemented for this operator\n");
    std::abort();
  }

  // M_ee_inv is always Hermitian so no need for isign?
  // for wilson it is the identity, for clover it is hermitian, for TWM it is gamma_5?
  inline void M_diag_inv(FourSpinorBlock *res,
                        FourSpinorBlock const *in,
                        int isign) const override {

    masterPrintf("M_ee_inv not yet implemented for this operator\n");
    std::abort();
  }



 private:
  double Mass;
  double TwistedMass;

  double mass_factor_alpha;
  double mass_factor_beta;

  const SU3MatrixBlock *u[2];
  FourSpinorBlock *tmp;
  TMDslash<FT, veclen, soalen, compress12> *D;

}; // Class
}; // Namespace
