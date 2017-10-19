/*
 * invmr.h
 *
 *  Created on: Jul 20, 2017
 *      Author: bjoo
 */

#pragma once
#include "qphix/linearOp.h"
#include "qphix/blas_new_c.h"
#include "qphix/print_utils.h"
#include "qphix/tsc.h"
#include "qphix/abs_solver.h"

namespace QPhiX
{

template <typename FT,
int V,
int S,
bool compress12,
typename EvenOddLinearOperatorBase =
    EvenOddLinearOperator<FT, V, S, compress12>>
    class InvMRBase
    {
    public:
  typedef typename Geometry<FT, V, S, compress12>::FourSpinorBlock Spinor;

  static constexpr int num_flav = EvenOddLinearOperatorBase::num_flav;

  InvMRBase(EvenOddLinearOperatorBase &M_, int MaxIters_, double Omega_)
  : M(M_), geom(M_.getGeometry()), MaxIters(MaxIters_), Omega(Omega_)
  {
    for (uint8_t f = 0; f < num_flav; ++f)
      Mr[f] = geom.allocCBFourSpinor();
    for (uint8_t f = 0; f < num_flav; ++f)
      r[f] = geom.allocCBFourSpinor();

    xmyThreads = geom.getNSIMT();
    norm2Threads = geom.getNSIMT();
    innerProductThreads = geom.getNSIMT();
    rxUpdateThreads = geom.getNSIMT();
  }

  ~InvMRBase()
  {
    for (uint8_t f = 0; f < num_flav; ++f) {
      geom.free(Mr[f]);
      geom.free(r[f]);
    }
  }

  Geometry<FT, V, S, compress12> &getGeometry() { return geom; }
  // This class overrides the `operator()` from `AbstractSolver`. Due to âname
  // hidingâ, the overloads of `operator()` in the base class are no longer
  // visible in this class. Therefore the single-flavor interface is not found
  // when trying to use the solver like it has worked before, namely with an
  // instance of this solver with automatic storage (i.e. no pointers). Here
  // we do want the overload for a single spinor pointer to delegate back to
  // the multi-flavor variant. The overloads need to be included explicitly
  // here. See http://stackoverflow.com/a/42588534/653152 for the full answer.
  void operator()(Spinor *const  x[num_flav],
      const Spinor* const rhs[num_flav],
      const double RsdTarget,
      int &n_iters,
      double &rsd_sq_final,
      unsigned long &site_flops,
      unsigned long &mv_apps,
      int isign,
      bool verbose,
      int cb,
      QPhiX::ResiduumType residType,
      bool terminateOnResiduum) const
  {
    site_flops = 0;
    mv_apps = 0;
    double rsd_sq = 0;
    double norm_chi;
    double cp;
    int k=0;

    if(MaxIters < 0) {
      masterPrintf("ERROR: QPhiX MR invalid value: MaxIters=%d\n", MaxIters);
    }
    if( MaxIters == 0 ) {
      //No work to do
      n_iters=0;
      rsd_sq_final=-1;
      return;
    }

    // Compute r=r0=rhs - A x
    // r <- Mx
    M(r, x, isign, cb);
    mv_apps++;

    // r <- chi - r <- chi - Mx  -- store result in Mr
    bicgstab_xmy<FT, V, S, compress12, num_flav>(rhs, r, geom, xmyThreads);
    site_flops += 24 * num_flav;

    // If we terminate on residuum:
    if( terminateOnResiduum ) {

      // Compute RHS norm
      norm2Spinor<FT, V, S, compress12, num_flav>(norm_chi, rhs, geom, norm2Threads);
      site_flops += 4 * 12 * num_flav;

      // Compute Target RsdSq
      rsd_sq = RsdTarget*RsdTarget;
      if ( residType == RELATIVE ) {
        rsd_sq *= norm_chi;
      }

      
      // Compute Norm of r0
      norm2Spinor<FT, V, S, compress12, num_flav>(cp, r, geom, norm2Threads);
      site_flops += 4 * 12 * num_flav;

      if( verbose ) {
        masterPrintf("QPhiX MR: iter=0: || r ||^2 = %16.8e Target || r ||^2 = %16.8e\n",cp, rsd_sq);
      }

      // Check if converged
      if( cp < rsd_sq ) {
        n_iters = 0;
        rsd_sq_final = cp;

        if( residType == ABSOLUTE) {
          // report absolute residuum
          if(verbose) {
            masterPrintf("QPhiX MR: Final iters=0 || r ||_accum=16.8e || r ||_actual = %16.8e\n",
                sqrt(rsd_sq_final), sqrt(rsd_sq_final));
          }
        }
        else {
          // report relative residuum
          rsd_sq_final /= norm_chi;
          if( verbose ) {
            masterPrintf("QPhiX MR: iters=0 || r ||/|| b ||_accum=16.8e || r ||/|| b ||_actual = %16.8e\n",
                sqrt(rsd_sq_final), sqrt(rsd_sq_final));
          }

        } // resid type

        return;
      } // cp < rsd_Sq


    } // terminateOnResiduum

    bool continueP = true;
    double c[2];
    double a[2];
    double d;

    while( continueP ) {
      ++k;  
      M(Mr,r,isign,cb);
      mv_apps++;

      double c_d_results[3];
      innerProductNorm<FT,V,S,compress12,num_flav>(c_d_results, Mr,r,geom,innerProductThreads);
      site_flops += 12*12*num_flav;

      // c_d_results holds what was { a[0], a[1], d }

      // a[0]/d
      a[0]=c_d_results[0]/c_d_results[2];

      // a[1]/d
      a[1]=c_d_results[1]/c_d_results[2];

      a[0] *= Omega; a[1] *= Omega;


      mr_rxupdate<FT, V, S, compress12, num_flav>(x,r,Mr,a,geom, rxUpdateThreads);
      site_flops += 2*8*12*num_flav;


      if( terminateOnResiduum ) {
        norm2Spinor<FT, V, S, compress12, num_flav>(cp, r, geom, norm2Threads);
        site_flops += 4 * 12 * num_flav;

        if( verbose ) {
          masterPrintf("QPhiX MR: iter=%d: || r ||^2 = %16.8e Target || r ||^2 = %16.8e\n",k, cp, rsd_sq);
        } // verbose
        continueP = ( k < MaxIters ) && ( cp > rsd_sq ) ;
      }
      else {
        if( verbose) {
          masterPrintf("QPhiX MR: iter=%d\n", k);
        }
        continueP =  ( k < MaxIters );
      } // terminate on resiruum
    } // while

    n_iters = k;
    if( terminateOnResiduum ) {

      rsd_sq_final = cp;
      if( residType == ABSOLUTE){
        if (verbose) {
          masterPrintf("QPhiX MR: Final iters=%d || r ||_accum=%16.8e || r ||_actual = %16.8e\n", k,
              sqrt(rsd_sq_final), sqrt(rsd_sq_final));
        }
      }
      else {
        rsd_sq_final /= norm_chi;
        if( verbose ) {
          masterPrintf("QPhiX MR: iters=%d || r ||/|| b ||_accum=%16.8e || r ||/|| b ||_actual = %16.8e\n",k,
              sqrt(rsd_sq_final), sqrt(rsd_sq_final));

        } // ver bose
      } // resid type

      return;
    }
    else {
      rsd_sq_final = -1;
      if(verbose) {
        masterPrintf("QPhiX MR: Smoother finished. %d iterations\n", n_iters);
      }
      return;
    }
  } // function


    private:
  EvenOddLinearOperatorBase &M;
  Geometry<FT, V, S, compress12> &geom;
  int MaxIters;
  double Omega;  
  Spinor* Mr[num_flav];
  Spinor* r[num_flav];
  int norm2Threads;
  int xmyThreads;
  int innerProductThreads;
  int rxUpdateThreads;


    }; // Class

template<typename  FT,
int V,
int S,
bool compress12,
typename EvenOddLinearOperatorBase =
    EvenOddLinearOperator<FT,V,S,compress12>>
    class InvMR : public AbstractSolver<FT,
    V,
    S,
    compress12,
    EvenOddLinearOperatorBase::num_flav>
{
public:
  typedef typename Geometry<FT, V, S, compress12>::FourSpinorBlock Spinor;

  static constexpr int num_flav = EvenOddLinearOperatorBase::num_flav;

  InvMR(EvenOddLinearOperatorBase &M_, int MaxIters_, double Omega)
  : baseSolver(M_, MaxIters_, Omega)
  {}

  Geometry<FT, V, S, compress12> &getGeometry() { return baseSolver.getGeometry(); }
  // This class overrides the `operator()` from `AbstractSolver`. Due to “name
  // hiding”, the overloads of `operator()` in the base class are no longer
  // visible in this class. Therefore the single-flavor interface is not found
  // when trying to use the solver like it has worked before, namely with an
  // instance of this solver with automatic storage (i.e. no pointers). Here
  // we do want the overload for a single spinor pointer to delegate back to
  // the multi-flavor variant. The overloads need to be included explicitly
  // here. See http://stackoverflow.com/a/42588534/653152 for the full answer.
  using AbstractSolver<FT, V, S, compress12, num_flav>::operator();

  void operator()(Spinor *const x[num_flav],
      const Spinor *const rhs[num_flav],
      const double RsdTarget,
      int &n_iters,
      double &rsd_sq_final,
      unsigned long &site_flops,
      unsigned long &mv_apps,
      int isign,
      bool verbose,
      int cb = 1,
      QPhiX::ResiduumType residType=QPhiX::RELATIVE) const override
      {
    baseSolver(x,rhs,RsdTarget,n_iters,rsd_sq_final,site_flops,mv_apps,isign,verbose,cb,residType,true);
      }



private:
  InvMRBase<FT,V,S,compress12,EvenOddLinearOperatorBase> baseSolver;
};

template<typename  FT,
int V,
int S,
bool compress12,
typename EvenOddLinearOperatorBase =
    EvenOddLinearOperator<FT,V,S,compress12>>
    class InvMRSmoother : public AbstractSolver<FT,
    V,
    S,
    compress12,
    EvenOddLinearOperatorBase::num_flav>
{
public:
  typedef typename Geometry<FT, V, S, compress12>::FourSpinorBlock Spinor;

  static constexpr int num_flav = EvenOddLinearOperatorBase::num_flav;

  InvMRSmoother(EvenOddLinearOperatorBase &M_, int MaxIters_, double Omega)
  : baseSolver(M_, MaxIters_, Omega)
  {}

  Geometry<FT, V, S, compress12> &getGeometry() { return baseSolver.getGeometry(); }
  // This class overrides the `operator()` from `AbstractSolver`. Due to “name
      // hiding”, the overloads of `operator()` in the base class are no longer
      // visible in this class. Therefore the single-flavor interface is not found
      // when trying to use the solver like it has worked before, namely with an
  // instance of this solver with automatic storage (i.e. no pointers). Here
  // we do want the overload for a single spinor pointer to delegate back to
  // the multi-flavor variant. The overloads need to be included explicitly
  // here. See http://stackoverflow.com/a/42588534/653152 for the full answer.
  using AbstractSolver<FT, V, S, compress12, num_flav>::operator();

  void operator()(Spinor *const x[num_flav],
      const Spinor *const rhs[num_flav],
      const double RsdTarget,
      int &n_iters,
      double &rsd_sq_final,
      unsigned long &site_flops,
      unsigned long &mv_apps,
      int isign,
      bool verbose,
      int cb = 1,
      QPhiX::ResiduumType residType=QPhiX::RELATIVE) const override
      {
    // Call solver with terminate on REsiduum set to false
    baseSolver(x,rhs,RsdTarget,n_iters,rsd_sq_final,site_flops,mv_apps,isign,verbose,cb,residType,false);
      }



private:
  InvMRBase<FT,V,S,compress12,EvenOddLinearOperatorBase> baseSolver;
};





}; // Namespace
