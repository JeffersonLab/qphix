#pragma once

#include "qphix/linearOp.h"
#include "qphix/tm_clov_dslash_def.h"

#include "qphix/print_utils.h"
#include <cstdlib>

namespace QPhiX
{

template <typename FT, int veclen, int soalen, bool compress12>
class EvenOddTMCloverOperator
    : public EvenOddLinearOperator<FT, veclen, soalen, compress12>
{

 public:
  typedef typename Geometry<FT, veclen, soalen, compress12>::FullCloverBlock
      FullCloverBlock;
  typedef typename Geometry<FT, veclen, soalen, compress12>::FourSpinorBlock
      FourSpinorBlock;
  typedef typename Geometry<FT, veclen, soalen, compress12>::SU3MatrixBlock
      SU3MatrixBlock;

  // Constructor
  // No anisotropy, all boundaries periodic for now.
  EvenOddTMCloverOperator(SU3MatrixBlock *u_[2],
                          FullCloverBlock *clov_[2],
                          FullCloverBlock *invclov_[2],
                          Geometry<FT, veclen, soalen, compress12> *geom_,
                          double t_boundary,
                          double aniso_coeff_s,
                          double aniso_coeff_t,
                          bool use_tbc_[4] = nullptr,
                          double tbc_phases_[4][2] = nullptr,
                          double const prec_mass_rho = 0.0)
      : D(new TMClovDslash<FT, veclen, soalen, compress12>(geom_,
                                                           t_boundary,
                                                           aniso_coeff_s,
                                                           aniso_coeff_t,
                                                           use_tbc_,
                                                           tbc_phases_,
                                                           prec_mass_rho))
  {
    Geometry<FT, veclen, soalen, compress12> &geom = D->getGeometry();
    u[0] = u_[0];
    u[1] = u_[1];
    clov[0] = clov_[0];
    clov[1] = clov_[1];
    invclov[0] = invclov_[0];
    invclov[1] = invclov_[1];
    tmp = (FourSpinorBlock *)geom.allocCBFourSpinor();
  }

  EvenOddTMCloverOperator(Geometry<FT, veclen, soalen, compress12> *geom_,
                          double t_boundary,
                          double aniso_coeff_s,
                          double aniso_coeff_t,
                          bool use_tbc_[4] = nullptr,
                          double tbc_phases_[4][2] = nullptr,
                          double const prec_mass_rho = 0.0)
      : D(new TMClovDslash<FT, veclen, soalen, compress12>(geom_,
                                                           t_boundary,
                                                           aniso_coeff_s,
                                                           aniso_coeff_t,
                                                           use_tbc_,
                                                           tbc_phases_,
                                                           prec_mass_rho))
  {
    Geometry<FT, veclen, soalen, compress12> &geom = D->getGeometry();
    tmp = (FourSpinorBlock *)geom.allocCBFourSpinor();
  }

  // Destructor
  ~EvenOddTMCloverOperator()
  {
    Geometry<FT, veclen, soalen, compress12> &geom = D->getGeometry();
    geom.free(tmp);
    delete D;
  }

  void setFields(SU3MatrixBlock *u_[2],
                 FullCloverBlock *clov_[2],
                 FullCloverBlock *invclov_[2])
  {
    u[0] = u_[0];
    u[1] = u_[1];
    clov[0] = clov_[0];
    clov[1] = clov_[1];
    invclov[0] = invclov_[0];
    invclov[1] = invclov_[1];
  }

  void operator()(FourSpinorBlock *res,
                  const FourSpinorBlock *in,
                  int isign,
                  int target_cb = 1) const override
  {
    int source_cb = 1 - target_cb;
    double beta = (double)0.25;

    D->dslash(
        tmp, in, u[source_cb], (const FullCloverBlock **)invclov, isign, source_cb);
    D->dslashAChiMinusBDPsi(res,
                            tmp,
                            in,
                            u[target_cb],
                            (const FullCloverBlock **)clov,
                            beta,
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

  Geometry<FT, veclen, soalen, compress12> &getGeometry()
  {
    return D->getGeometry();
  }

 private:
  double Mass;
  TMClovDslash<FT, veclen, soalen, compress12> *D;
  SU3MatrixBlock *u[2];
  FullCloverBlock *clov[2];
  FullCloverBlock *invclov[2];
  FourSpinorBlock *tmp;

}; // Class
}; // Namespace
