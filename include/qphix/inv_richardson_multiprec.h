#pragma once

#include "qphix/abs_solver.h"
#include "qphix/linearOp.h"
#include "qphix/blas_new_c.h"
#include "qphix/print_utils.h"

namespace QPhiX
{

template <typename FT,
          int V,
          int S,
          bool Compress,
          typename FTInner,
          int VInner,
          int SInner,
          bool CompressInner,
          bool MdagM = false,
          typename EvenOddLinearOperatorBase =
              EvenOddLinearOperator<FT, V, S, Compress> >
class InvRichardsonMultiPrec
    : public AbstractSolver<FT, V, S, Compress, EvenOddLinearOperatorBase::num_flav>
{
 public:
  typedef typename Geometry<FT, V, S, Compress>::FourSpinorBlock Spinor;
  typedef typename Geometry<FTInner, VInner, SInner, CompressInner>::FourSpinorBlock
      SpinorInner;

  static constexpr int num_flav = EvenOddLinearOperatorBase::num_flav;

  // This class overrides the `operator()` from `AbstractSolver`. Due to “name
  // hiding”, the overloads of `operator()` in the base class are no longer
  // visible in this class. Therefore the single-flavor interface is not found
  // when trying to use the solver like it has worked before, namely with an
  // instance of this solver with automatic storage (i.e. no pointers). Here
  // we do want the overload for a single spinor pointer to delegate back to
  // the multi-flavor variant. The overloads need to be included explicitly
  // here. See http://stackoverflow.com/a/42588534/653152 for the full answer.
  using AbstractSolver<FT, V, S, Compress, num_flav>::operator();

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
    int iter = 0;
    int mv_apps_outer = 0;
    int site_flops_outer = 0;
    int mv_apps_inner_total = 0;
    int site_flops_inner_total = 0;

    double rhs_sq;
    double r_norm_sq;
    norm2Spinor<FT, V, S, Compress, num_flav>(rhs_sq, rhs, geom, norm2Threads);
    site_flops_outer += 24 + 23; // 24 muls, 23 adds

    // This is the target residuum
    double rsd_t = RsdTarget * RsdTarget;
    if( residType == QPhiX::RELATIVE) {
	rsd_t  *= rhs_sq;
        masterPrintf("RICHARDSON: RelativeResiduum requested\n");
    }
    else {
	masterPrintf("RICHARDSON: Absolut Residuum requested\n");
    }

    if( MdagM ){
      m_outer(tmp_MdagM, x, isign, cb);
      m_outer(tmp, tmp_MdagM, -isign, cb);
      mv_apps_outer += 2;
    } else {
      m_outer(tmp, x, isign, cb);
      mv_apps_outer++;
    }

    // r = rhs - M x
    xmyNorm2Spinor<FT, V, S, Compress, num_flav>(
        r, rhs, tmp, r_norm_sq, geom, xmyNormThreads);
    site_flops_outer += 24 + 24 + 23; // (24 -,  24 square, 23 adds)

    masterPrintf(
        "RICHARDSON: Initial || r ||^2 = %16.8e   Target || r ||^2 = %16.8e\n",
        r_norm_sq,
        rsd_t);

    if (r_norm_sq < rsd_t) {

      n_iters = 0;
      rsd_sq_final = r_norm_sq;
      mv_apps = mv_apps_outer + mv_apps_inner_total;
      site_flops = site_flops_outer + site_flops_inner_total * num_flav;
      masterPrintf("RICHARDSON: Solver converged at iter 0\n");
      masterPrintf("RICHARDSON: Inner MV Apps=%lu Outer MV Apps=%lu Inner Site "
                   "Flops=%lu Outer Site Flops=%lu\n",
                   mv_apps_inner_total,
                   mv_apps_outer,
                   site_flops_inner_total,
                   site_flops_outer);

      return;
    }

    for (iter = 1; iter < max_iters; iter++) {

      // Need resid target
      int n_iters_inner;
      double rsd_sq_final_inner;
      unsigned long site_flops_inner;
      unsigned long mv_apps_inner;

      double delta_finish = sqrt(rsd_t / r_norm_sq);
      double delta_use = (delta_finish > delta) ? delta_finish : delta;
      if (verbose) {
        masterPrintf("RICHARDSON: iter=%d Delta = %16.8e  Delta Finish=%16.8e "
                     "Delta_use=%16.8e\n",
                     iter,
                     delta,
                     delta_finish,
                     delta_use);
      }

      // Down Convert r to inner precision. Normalize along the way. Don't count
      // these flops? (Not useful?)
      double r_norm = sqrt(r_norm_sq);
      double scale_factor = ((double)1) / r_norm;
      convert<FTInner, VInner, SInner, CompressInner, FT, V, S, Compress, num_flav>(
          r_inner, scale_factor, r, geom_inner, geom, convertToThreads);

      // Zero out dx_single
      zeroSpinor<FTInner, VInner, SInner, CompressInner, num_flav>(
          dx_inner, geom_inner, zeroThreads);

      solver_inner(dx_inner,
                   r_inner,
                   delta_use,
                   n_iters_inner,
                   rsd_sq_final_inner,
                   site_flops_inner,
                   mv_apps_inner,
                   isign,
                   verbose,
                   cb);

      mv_apps_inner_total += mv_apps_inner;
      site_flops_inner_total += site_flops_inner;

      // Up convert and un-normalize (Again don't count the flops, since its not
      // 'useful'?
      convert<FT, V, S, Compress, FTInner, VInner, SInner, CompressInner, num_flav>(
          delta_x, r_norm, dx_inner, geom, geom_inner, convertFromThreads); 

      if( MdagM ){
        m_outer(tmp_MdagM, delta_x, isign, cb);
        m_outer(tmp, tmp_MdagM, -isign, cb);
        mv_apps_outer += 2;
      }else{
        m_outer(tmp, delta_x, isign, cb);
        mv_apps_outer++;
      }

      richardson_rxupdateNormR<FT, V, S, Compress, num_flav>(
          x, r, delta_x, tmp, r_norm_sq, geom, rxUpdateThreads);
      site_flops_outer += 8 * 12;
      if (verbose)
        masterPrintf("RICHARDSON: Iter %d   r_norm_sq=%e   Target=%e\n",
                     iter,
                     r_norm_sq,
                     rsd_t);

      if (r_norm_sq < rsd_t) {
        // Converged. Compute final residuum.
        if( MdagM ){
          m_outer(tmp_MdagM, x, isign, cb);
          m_outer(tmp, tmp_MdagM, -isign, cb);
          mv_apps_outer += 2;
        } else {
          m_outer(tmp, x, isign, cb);
          mv_apps_outer++;
        }

        xmyNorm2Spinor<FT, V, S, Compress, num_flav>(
            r, rhs, tmp, r_norm_sq, geom, xmyNormThreads);
        site_flops_outer += 24 + 24 + 23; // (24 -,  24 square, 23 adds)

        rsd_sq_final = r_norm_sq / rhs_sq;
        n_iters = iter;
        mv_apps = mv_apps_outer + mv_apps_inner_total;
        site_flops = site_flops_outer + site_flops_inner_total * num_flav;
        masterPrintf("RICHARDSON: Solver converged at iter %d || r || = %16.8e\n",
                     iter,
                     sqrt(rsd_sq_final));
        masterPrintf("RICHARDSON: Inner MV Apps=%lu Outer MV Apps=%lu Inner Site "
                     "Flops=%lu Outer Site Flops=%lu\n",
                     mv_apps_inner_total,
                     mv_apps_outer,
                     site_flops_inner_total,
                     site_flops_outer);
        return;
      }
    }

    // Not Converged. Compute final residuum.
    if( MdagM ){
      m_outer(tmp_MdagM, x, isign, cb);
      m_outer(tmp, tmp_MdagM, -isign, cb);
      mv_apps_outer += 2;
    } else {
      m_outer(tmp, x, isign, cb);
      mv_apps_outer++;
    }

    xmyNorm2Spinor<FT, V, S, Compress, num_flav>(
        r, rhs, tmp, r_norm_sq, geom, xmyNormThreads);
    site_flops_outer += 24 + 24 + 23;

    rsd_sq_final = r_norm_sq / rhs_sq;
    n_iters = iter;

    mv_apps = mv_apps_outer + mv_apps_inner_total;
    site_flops = site_flops_outer + site_flops_inner_total * num_flav;
    masterPrintf("RICHARDSON: Solver NOT converged at iter %d || r || = %16.8e\n",
                 sqrt(rsd_sq_final));
    masterPrintf("RICHARDSON: Inner MV Apps=%lu Outer MV Apps=%lu Inner Site "
                 "Flops=%lu Outer Site Flops=%lu\n",
                 mv_apps_inner_total,
                 mv_apps_outer,
                 site_flops_inner_total,
                 site_flops_outer);

    return;
  }

  InvRichardsonMultiPrec(
      EvenOddLinearOperatorBase &m_outer_,
      AbstractSolver<FTInner, VInner, SInner, CompressInner, EvenOddLinearOperatorBase::num_flav > &solver_inner_,
      const double delta_,
      const int max_iters_)
      : m_outer(m_outer_), solver_inner(solver_inner_), delta(delta_),
        max_iters(max_iters_), geom(m_outer_.getGeometry()),
        geom_inner(solver_inner_.getGeometry())
  {
    for (uint8_t f = 0; f < num_flav; ++f)
      r[f] = (Spinor *)geom.allocCBFourSpinor();
    for (uint8_t f = 0; f < num_flav; ++f)
      tmp[f] = (Spinor *)geom.allocCBFourSpinor();
    for (uint8_t f = 0; f < num_flav; ++f)
      delta_x[f] = (Spinor *)geom.allocCBFourSpinor();
    for (uint8_t f = 0; f < num_flav; ++f)
      r_inner[f] = (SpinorInner *)geom_inner.allocCBFourSpinor();
    for (uint8_t f = 0; f < num_flav; ++f)
      dx_inner[f] = (SpinorInner *)geom_inner.allocCBFourSpinor();
    for (uint8_t f = 0; f < num_flav; ++f){
      if( MdagM ){
        tmp_MdagM[f] = (Spinor *)geom.allocCBFourSpinor();
      } else {
        tmp_MdagM[f] = nullptr;
      }
    }

    // Initial value for norm2threads. Use all threads
    norm2Threads = geom.getNSIMT();
    xmyNormThreads = geom.getNSIMT();
    zeroThreads = geom.getNSIMT();
    convertToThreads = geom.getNSIMT();
    convertFromThreads = geom.getNSIMT();
    rxUpdateThreads = geom.getNSIMT();
  }

  ~InvRichardsonMultiPrec()
  {
    for (uint8_t f = 0; f < num_flav; ++f) {
      geom.free(r[f]);
      geom.free(tmp[f]);
      geom.free(tmp_MdagM[f]);
      geom.free(delta_x[f]);

      geom_inner.free(r_inner[f]);
      geom_inner.free(dx_inner[f]);
    }
  }

  void tune(void) { solver_inner.tune(); }

  Geometry<FT, V, S, Compress> &getGeometry() { return geom; }

 private:
  EvenOddLinearOperatorBase &m_outer;
  AbstractSolver<FTInner, VInner, SInner, CompressInner, num_flav> &solver_inner;
  const double delta;
  const int max_iters;

  // Internal Spinors
  Spinor *r[num_flav];
  Spinor *tmp[num_flav]; // For computing residua  and also Ddelta_x
  Spinor *tmp_MdagM[num_flav]; // when the inner solver is CG, we need to apply MdagM to get the
                              // true residual
  Spinor *delta_x[num_flav];
  SpinorInner *r_inner[num_flav];
  SpinorInner *dx_inner[num_flav];

  int norm2Threads;
  int xmyNormThreads;
  int zeroThreads;
  int convertToThreads;
  int convertFromThreads;
  int rxUpdateThreads;

  // The geometries for allocations/conversions
  Geometry<FT, V, S, Compress> &geom;
  Geometry<FTInner, VInner, SInner, CompressInner> &geom_inner;
};
};
