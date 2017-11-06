#pragma once

#include "qphix/linearOp.h"
#include "qphix/print_utils.h"
#include "qphix/blas_new_c.h"
#include "qphix/tsc.h"
#include "qphix/abs_solver.h"

namespace QPhiX
{

/* This is a really silly but useful debug "solver" which checks the hermiticity properties
 * of M and M^dag using the solver interface, which is very useful for calling it from
 * client applications because more likely than not, wrappers exist for calling a solver
 * but not necessarily for simply calling the operator.
 *
 * The idea would be to supply two random fields in 'x' and 'rhs'.
 */ 

// FIXME!!!: passing compress12 as a template to solver breaks 'delegation of
// responsibilities rule.
// Solution may be to take it out of geometry and create a fields class.
// Geometry can then deal
// exclusively with the blocking and stuff...
//
// That will be a second refactor step when everything works
template <typename FT,
          int veclen,
          int soalen,
          bool compress12,
          typename EvenOddLinearOperatorBase =
              EvenOddLinearOperator<FT, veclen, soalen, compress12>>
class InvDummyHermTest : public AbstractSolver<FT,
                                               veclen,
                                               soalen,
                                               compress12,
                                               EvenOddLinearOperatorBase::num_flav>
{
 public:
  typedef typename Geometry<FT, veclen, soalen, compress12>::FourSpinorBlock Spinor;

  static constexpr int num_flav = EvenOddLinearOperatorBase::num_flav;

  InvDummyHermTest(EvenOddLinearOperatorBase &M_, int MaxIters_)
    : M(M_), geom(M_.getGeometry()), MaxIters(MaxIters_)
  {
    // Length of vectors in floats, including PADs
    masterPrintf("Initializing Dummy Hermiticity Test: Nvec=%d Ny=%d Nz=%d Nt=%d\n",
                 geom.nVecs(),
                 geom.Ny(),
                 geom.Nz(),
                 geom.Nt());

    for (uint8_t f = 0; f < num_flav; ++f)
      m_rhs[f] = (Spinor *)geom.allocCBFourSpinor();
    for (uint8_t f = 0; f < num_flav; ++f)
      mdag_m_rhs[f] = (Spinor *)geom.allocCBFourSpinor();
    for (uint8_t f = 0; f < num_flav; ++f)
      m_x[f] = (Spinor *)geom.allocCBFourSpinor();

    copyThreads = geom.getNSIMT();
    aypxThreads = geom.getNSIMT();
    xmyNormThreads = geom.getNSIMT();
    rmammpNorm2rxpapThreads = geom.getNSIMT();
    norm2Threads = geom.getNSIMT();
  }

  ~InvDummyHermTest()
  {
    for (uint8_t f = 0; f < num_flav; ++f) {
      geom.free(m_x[f]);
      geom.free(m_rhs[f]);
      geom.free(mdag_m_rhs[f]);
    }
  }

  // This class overrides the `operator()` from `AbstractSolver`. Due to “name
  // hiding”, the overloads of `operator()` in the base class are no longer
  // visible in this class. Therefore the single-flavor interface is not found
  // when trying to use the solver like it has worked before, namely with an
  // instance of this solver with automatic storage (i.e. no pointers). Here
  // we do want the overload for a single spinor pointer to delegate back to
  // the multi-flavor variant. The overloads need to be included explicitly
  // here. See http://stackoverflow.com/a/42588534/653152 for the full answer.
  using AbstractSolver<FT, veclen, soalen, compress12, num_flav>::operator();

  virtual void operator()(Spinor *const x[num_flav],
                          const Spinor *const rhs[num_flav],
                          const double RsdTarget,
                          int &n_iters,
                          double &rsd_sq_final,
                          unsigned long &site_flops,
                          unsigned long &mv_apps,
                          int isign,
                          bool verboseP,
                          int cb = 1) const override
  {
    if (verboseP) {
      masterPrintf("Entering isign=%d hermiticity test with num_flav=%d and cb=%d\n", isign, num_flav, cb);
    }

    int n_cores;
    int n_simt;
    double cp;

    site_flops = 0;
    mv_apps = 3;
    n_iters = 1;
    n_cores = geom.getNumCores();
    n_simt = geom.getNSIMT();

    double tmp_d;
    M(m_rhs, rhs, isign, cb);
    M(mdag_m_rhs, m_rhs, -isign, cb);
    M(m_x, x, isign, cb);

    double x_inner_mdag_m_rhs[2] = {0.0, 0.0};
    double m_x_inner_m_rhs[2] = {0.0, 0.0};

    innerProduct<FT, veclen, soalen, compress12, num_flav>(x_inner_mdag_m_rhs,
                                                           x,
                                                           mdag_m_rhs,
                                                           geom,
                                                           n_simt);
    
    innerProduct<FT, veclen, soalen, compress12, num_flav>(m_x_inner_m_rhs,
                                                           m_x,
                                                           m_rhs,
                                                           geom,
                                                           n_simt);

    masterPrintf("\n\nHERMITICY TEST\n");
    masterPrintf(" %20s: %lf %lf \n %20s: %lf %lf\n\n\n",
                 "< x | Mdag M y >", x_inner_mdag_m_rhs[0], x_inner_mdag_m_rhs[1],
                 "< M x | M y >", m_x_inner_m_rhs[0], m_x_inner_m_rhs[1]);
  }

  void setCopyThreads(int c) { copyThreads = c; }
  void setAypxThreads(int c) { aypxThreads = c; }
  void setXmyNormThreads(int c) { xmyNormThreads = c; }
  void setRmammpNorm2rxpapThreads(int c) { rmammpNorm2rxpapThreads = c; }
  void setNorm2Threads(int c) { norm2Threads = c; }

  int getCopyThreads(void) { return copyThreads; }
  int getAypxThreads(void) { return aypxThreads; }
  int getXmyNormThreads(void) { return xmyNormThreads; }
  int getRmammpNorm2rxpapThreads(void) { return rmammpNorm2rxpapThreads; }
  int getNorm2Threads(void) { return norm2Threads; }

  Geometry<FT, veclen, soalen, compress12> &getGeometry() { return geom; }

 private:
  EvenOddLinearOperatorBase &M;
  Geometry<FT, veclen, soalen, compress12> &geom;
  int MaxIters;

  Spinor *m_rhs[num_flav];
  Spinor *mdag_m_rhs[num_flav];
  Spinor *m_x[num_flav];

  int copyThreads;
  int aypxThreads;
  int xmyNormThreads;
  int rmammpNorm2rxpapThreads;
  int norm2Threads;

};
};
