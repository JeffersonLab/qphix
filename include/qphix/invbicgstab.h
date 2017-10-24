#pragma once

/**
 \TODO user serviceable options -- should be moved to autoconf control
 */
#define QPHIX_VERBOSE_BICGSTAB
#define QPHIX_TIMING_BICGSTAB

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
class InvBiCGStab : public AbstractSolver<FT,
                                          V,
                                          S,
                                          compress12,
                                          EvenOddLinearOperatorBase::num_flav>
{
 public:
  typedef typename Geometry<FT, V, S, compress12>::FourSpinorBlock Spinor;

  static constexpr int num_flav = EvenOddLinearOperatorBase::num_flav;

  InvBiCGStab(EvenOddLinearOperatorBase &M_, int MaxIters_)
      : M(M_), geom(M_.getGeometry()), MaxIters(MaxIters_)
  {
    for (uint8_t f = 0; f < num_flav; ++f)
      r[f] = geom.allocCBFourSpinor();
    for (uint8_t f = 0; f < num_flav; ++f)
      r0[f] = geom.allocCBFourSpinor();
    for (uint8_t f = 0; f < num_flav; ++f)
      p[f] = geom.allocCBFourSpinor();
    for (uint8_t f = 0; f < num_flav; ++f)
      v[f] = geom.allocCBFourSpinor();
    for (uint8_t f = 0; f < num_flav; ++f)
      t[f] = geom.allocCBFourSpinor();

    norm2Threads = geom.getNSIMT();
    xmyThreads = geom.getNSIMT();
    copyThreads = geom.getNSIMT();
    zeroThreads = geom.getNSIMT();
    innerProductThreads = geom.getNSIMT();
    pUpdateThreads = geom.getNSIMT();
    sUpdateThreads = geom.getNSIMT();
    rxUpdateThreads = geom.getNSIMT();
  }

  ~InvBiCGStab()
  {
    for (uint8_t f = 0; f < num_flav; ++f) {
      geom.free(r[f]);
      geom.free(r0[f]);
      geom.free(p[f]);
      geom.free(v[f]);
      geom.free(t[f]);
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
    site_flops = 0;
    mv_apps = 0;
    double rhs_sq = 0;
    double r_norm;
    // Double chi_sq = norm2(chi,s)
    norm2Spinor<FT, V, S, compress12, num_flav>(rhs_sq, rhs, geom, norm2Threads);
    site_flops += 4 * 12 * num_flav;

    double rsd_sq = RsdTarget * RsdTarget;
    if ( residType == QPhiX::RELATIVE) {
      if( verbose ) {
        masterPrintf("BICGSTAB: Relative Residuum requested\n");
      }
      rsd_sq *= rhs_sq;
    }
    else {
      if( verbose ){
          	masterPrintf("BICGSTAB: Absolute Residuum requested\n");
      }
    }

    if (verbose)
    	masterPrintf("BICGSTAB: Target Rsd= %e\n", rsd_sq);


    // Compute r=r0=rhs - A x
    // A(r0,psi,isign)
    M(r0, x, isign, cb);
    mv_apps++;

    // r = chi - r0  -- store result in r0
    bicgstab_xmy<FT, V, S, compress12, num_flav>(rhs, r0, geom, xmyThreads);
    site_flops += 24 * num_flav;

    // Check norm of r0
    norm2Spinor<FT,V,S,compress12, num_flav>(r_norm, r0, geom, norm2Threads);
    if( r_norm < rsd_sq ) {
    	masterPrintf("BICGSTAB converged at iter 0: ||r||=%16.8e target=%16.8e\n",
    			  	  sqrt(r_norm), sqrt(rsd_sq) );
    	rsd_sq_final = r_norm;
    	n_iters = 0;
        return;
    }
    // r = r0
    copySpinor<FT, V, S, compress12, num_flav>(r, r0, geom, copyThreads);

    // Now intitialize p and v to zero.
    zeroSpinor<FT, V, S, compress12, num_flav>(p, geom, zeroThreads);
    zeroSpinor<FT, V, S, compress12, num_flav>(v, geom, zeroThreads);

    // rho_prev = 1
    double rho_prev_c[2] = {(double)1, (double)0};

    // alpha = 1
    double alpha_c[2] = {(double)1, (double)0};

    // omega = 1
    double omega_c[2] = {(double)1, (double)0};

    double rho_c[2] = {(double)0, (double)0};

    bool notConvP = true;
    int k = 1;
    // Now iterate
    for (k = 1; (k <= MaxIters) && notConvP; k++) {

      // rho_{k+1} = < r_0 | r >
      innerProduct<FT, V, S, compress12, num_flav>(
          rho_c, r0, r, geom, innerProductThreads);
      site_flops += 8 * 12 * num_flav;

      if ((rho_c[0] == 0) && (rho_c[1] == 0)) {
        masterPrintf("BICGSTAB: Breakdown in iteration %d. rho = 0\n", k);

        rsd_sq_final = r_norm;
        n_iters = k;
        return;
      }

      double beta_c[2];
      // beta = ( rho_{k+1}/ rho_{k} ) ( alpha / omega )
      double tmp_c[2];
      double tmp2_c[2];
      complex_div(tmp_c, rho_c, rho_prev_c);
      complex_div(tmp2_c, alpha_c, omega_c);
      complex_mul(beta_c, tmp_c, tmp2_c);

      // p = r + beta(p - omega v)
      // y = x + a(y - bz)

      bicgstab_p_update<FT, V, S, compress12, num_flav>(
          r, p, v, beta_c, omega_c, geom, pUpdateThreads);
      site_flops += 16 * 12 * num_flav;

      // v = Ap
      M(v, p, isign, cb);
      mv_apps++;

      innerProduct<FT, V, S, compress12, num_flav>(
          tmp_c, r0, v, geom, innerProductThreads);
      site_flops += 8 * 12 * num_flav;

      if ((tmp_c[0] == 0) && (tmp_c[1] == 0)) {
        masterPrintf("BICGSTAB: Breakdown in iteration %d. <r_0|v> = 0\n", k);
        n_iters = k;
        return;
      }

      // Alpha = rho/<r0,v>
      complex_div(alpha_c, rho_c, tmp_c);

      // Preserve rho as rho_prev
      rho_prev_c[0] = rho_c[0];
      rho_prev_c[1] = rho_c[1];

      // s = r - alpha v
      // I can overlap s with r because I recompute r at the end

      // r = s = r - alpha v:   complex y=ax+y with a=-alpha x=v, y=r
      //	FT alpha_cr[2] = { (FT)(alpha_c[0]), (FT)(alpha_c[1]) };

      bicgstab_s_update<FT, V, S, compress12, num_flav>(
          alpha_c, r, v, geom, sUpdateThreads);
      site_flops += 8 * 12 * num_flav;

      // t = As
      M(t, r, isign, cb);
      mv_apps++;

      double t_norm = 0;
      norm2Spinor<FT, V, S, compress12, num_flav>(t_norm, t, geom, norm2Threads);
      site_flops += 4 * 12 * num_flav;

      if (t_norm == 0) {
        masterPrintf("BICGSTAB: Breakdown in iteration %d. ||t|| = 0\n", k);
        rsd_sq_final = r_norm;
        n_iters = k;
        return;
      }

      innerProduct<FT, V, S, compress12, num_flav>(
          omega_c, t, r, geom, innerProductThreads);
      site_flops += 8 * 12 * num_flav;

      omega_c[0] /= t_norm;
      omega_c[1] /= t_norm;

      // x = omega r + x +  alpha p
      // r = r - omega t
      // r_norm = norm2(r);

      bicgstab_rxupdate<FT, V, S, compress12, num_flav>(
          x, r, t, p, omega_c, alpha_c, r_norm, geom, rxUpdateThreads);
      site_flops += 28 * 12 * num_flav;

      if (verbose)
        masterPrintf(
            "BICGSTAB: iter %d r_norm = %e  target = %e \n", k, r_norm, rsd_sq);
      if (r_norm < rsd_sq) {
        notConvP = false; // Converged
      }

    } // Loop over iterations

    rsd_sq_final = r_norm;
    n_iters = k;

    if (notConvP == true) {
      masterPrintf("Solver did not converge in %d iterations\n", k);
    }
    return;
  }

  Geometry<FT, V, S, compress12> &getGeometry() { return geom; }

 private:
  EvenOddLinearOperatorBase &M;
  Geometry<FT, V, S, compress12> &geom;
  int MaxIters;

  inline void complex_div(double res[2], const double l[2], const double r[2]) const
  {
    double tmp = (double)1 / (r[0] * r[0] + r[1] * r[1]);

    res[0] = (l[0] * r[0] + l[1] * r[1]) * tmp;
    res[1] = (l[1] * r[0] - l[0] * r[1]) * tmp;
  }

  inline void
  complex_mul(double res[2], const double mul1[2], const double mul2[2]) const
  {
    res[0] = mul1[0] * mul2[0] - mul1[1] * mul2[1];
    res[1] = mul1[0] * mul2[1] + mul1[1] * mul2[0];
  }

  Spinor *r[num_flav];
  Spinor *r0[num_flav];
  Spinor *p[num_flav];
  Spinor *v[num_flav];
  Spinor *t[num_flav];

  int norm2Threads;
  int xmyThreads;
  int copyThreads;
  int zeroThreads;
  int innerProductThreads;
  int pUpdateThreads;
  int sUpdateThreads;
  int rxUpdateThreads;

  void tuneNorm2Threads(int iters)
  {
    if (r != 0x0) {
      // Do first with 1 thread
      double rnorm;
      zeroSpinor<FT, V, S, compress12>(r, geom, zeroThreads);
      norm2Threads = 1;
      double start_time = omp_get_wtime();
      for (int i = 0; i < iters; i++) {
        norm2Spinor(rnorm, r, geom, norm2Threads);
      }
      double stop_time = omp_get_wtime();
      double best_time = stop_time - start_time;
      masterPrintf("tuneNorm2Threads: threads = %d, current_time=%g (s)\n",
                   norm2Threads,
                   best_time);
      for (int threads = 2; threads <= geom.getNSIMT(); threads++) {
        start_time = omp_get_wtime();
        for (int i = 0; i < iters; i++) {
          norm2Spinor(rnorm, r, geom, threads);
        }
        stop_time = omp_get_wtime();
        double current_time = stop_time - start_time;

        masterPrintf(
            "tuneNorm2Threads: threads = %d, current_time = %g (s), best=%g (s)\n",
            threads,
            current_time,
            best_time);

        if (current_time < best_time) {
          best_time = current_time;
          norm2Threads = threads;
        }
      }
    }
  }

  void tuneXMYThreads(int iters)
  {
    if (r != 0x0 && v != 0) {
      // Do first with 1 thread

      zeroSpinor<FT, V, S, compress12>(r, geom, zeroThreads);
      zeroSpinor<FT, V, S, compress12>(v, geom, zeroThreads);
      xmyThreads = 1;
      double start_time = omp_get_wtime();
      for (int i = 0; i < iters; i++) {
        bicgstab_xmy<FT, V, S, compress12>(r, v, geom, xmyThreads);
      }
      double stop_time = omp_get_wtime();
      double best_time = stop_time - start_time;
      masterPrintf("tuneXMYThreads: threads = %d, current_time=%g (s)\n",
                   xmyThreads,
                   best_time);
      for (int threads = 2; threads <= geom.getNSIMT(); threads++) {
        start_time = omp_get_wtime();
        for (int i = 0; i < iters; i++) {
          bicgstab_xmy<FT, V, S, compress12>(r, v, geom, threads);
        }
        stop_time = omp_get_wtime();
        double current_time = stop_time - start_time;

        masterPrintf(
            "tuneXMYThreads: threads = %d, current_time = %g (s), best=%g(s)\n",
            threads,
            current_time,
            best_time);
        if (current_time < best_time) {
          best_time = current_time;
          xmyThreads = threads;
        }
      }
    }
  }

  void tuneCopyThreads(int iters)
  {
    if (r != 0x0 && v != 0) {
      // Do first with 1 thread

      zeroSpinor<FT, V, S, compress12>(r, geom, zeroThreads);
      zeroSpinor<FT, V, S, compress12>(v, geom, zeroThreads);
      copyThreads = 1;
      double start_time = omp_get_wtime();
      for (int i = 0; i < iters; i++) {
        copySpinor<FT, V, S, compress12>(r, v, geom, copyThreads);
      }

      double stop_time = omp_get_wtime();
      double best_time = stop_time - start_time;
      masterPrintf("tuneCopyThreads: threads = %d, current_time=%g (s)\n",
                   copyThreads,
                   best_time);

      for (int threads = 2; threads <= geom.getNSIMT(); threads++) {
        start_time = omp_get_wtime();
        for (int i = 0; i < iters; i++) {
          copySpinor<FT, V, S, compress12>(r, v, geom, threads);
        }
        stop_time = omp_get_wtime();
        double current_time = stop_time - start_time;

        masterPrintf(
            "tuneCopyThreads: threads = %d, current_time = %g (s), best=%g(s)\n",
            threads,
            current_time,
            best_time);
        if (current_time < best_time) {
          best_time = current_time;
          copyThreads = threads;
        }
      }
    }
  }

  void tuneZeroThreads(int iters)
  {
    if (r != 0x0) {

      zeroThreads = 1;
      double start_time = omp_get_wtime();
      for (int i = 0; i < iters; i++) {
        zeroSpinor<FT, V, S, compress12>(r, geom, zeroThreads);
      }
      double stop_time = omp_get_wtime();
      double best_time = stop_time - start_time;
      masterPrintf("tuneZeroThreads: threads = %d, current_time=%g (s)\n",
                   zeroThreads,
                   best_time);
      for (int threads = 2; threads <= geom.getNSIMT(); threads++) {
        start_time = omp_get_wtime();
        for (int i = 0; i < iters; i++) {
          zeroSpinor<FT, V, S, compress12>(r, geom, threads);
        }
        stop_time = omp_get_wtime();
        double current_time = stop_time - start_time;

        masterPrintf(
            "tuneZeroThreads: threads = %d, current_time = %g (s), best=%g(s)\n",
            threads,
            current_time,
            best_time);
        if (current_time < best_time) {
          best_time = current_time;
          zeroThreads = threads;
        }
      }
    }
  }

  void tuneInnerProductThreads(int iters)
  {
    if (r != 0x0 && v != 0) {
      // Do first with 1 thread

      zeroSpinor<FT, V, S, compress12>(r, geom, zeroThreads);
      zeroSpinor<FT, V, S, compress12>(v, geom, zeroThreads);
      innerProductThreads = 1;
      double iprod_c[2];
      double start_time = omp_get_wtime();
      for (int i = 0; i < iters; i++) {
        innerProduct<FT, V, S, compress12>(iprod_c, r, v, geom, innerProductThreads);
      }
      double stop_time = omp_get_wtime();
      double best_time = stop_time - start_time;
      masterPrintf("tuneInnerProductThreads: threads = %d, current_time=%g (s)\n",
                   innerProductThreads,
                   best_time);

      for (int threads = 2; threads <= geom.getNSIMT(); threads++) {
        start_time = omp_get_wtime();
        for (int i = 0; i < iters; i++) {
          innerProduct<FT, V, S, compress12>(iprod_c, r, v, geom, threads);
        }
        stop_time = omp_get_wtime();
        double current_time = stop_time - start_time;

        masterPrintf("tuneInnerProductThreads: threads = %d, current_time = %g (s), "
                     "best=%g(s)\n",
                     threads,
                     current_time,
                     best_time);
        if (current_time < best_time) {
          best_time = current_time;
          innerProductThreads = threads;
        }
      }
    }
  }

  void tunePUpdateThreads(int iters)
  {
    if (r != 0x0 && v != 0 && p != 0) {
      zeroSpinor<FT, V, S, compress12>(r, geom, zeroThreads);
      zeroSpinor<FT, V, S, compress12>(p, geom, zeroThreads);
      zeroSpinor<FT, V, S, compress12>(v, geom, zeroThreads);
      pUpdateThreads = 1;
      double beta_cr[2] = {(double)1.0, (double)0.5};
      double omega_cr[2] = {(double)2.0, (double)-0.5};
      double start_time = omp_get_wtime();
      for (int i = 0; i < iters; i++) {
        bicgstab_p_update<FT, V, S, compress12>(
            r, p, v, beta_cr, omega_cr, geom, pUpdateThreads);
      }
      double stop_time = omp_get_wtime();
      double best_time = stop_time - start_time;

      masterPrintf("tunePUpdateThreads: threads = %d, current_time=%g (s)\n",
                   pUpdateThreads,
                   best_time);

      for (int threads = 2; threads <= geom.getNSIMT(); threads++) {
        start_time = omp_get_wtime();
        for (int i = 0; i < iters; i++) {
          bicgstab_p_update<FT, V, S, compress12>(
              r, p, v, beta_cr, omega_cr, geom, threads);
        }
        stop_time = omp_get_wtime();
        double current_time = stop_time - start_time;

        masterPrintf(
            "tunePUpdateThreads: threads = %d, current_time = %g (s), best=%g(s)\n",
            threads,
            current_time,
            best_time);
        if (current_time < best_time) {
          best_time = current_time;
          pUpdateThreads = threads;
        }
      }
    }
  }

  void tuneSUpdateThreads(int iters)
  {
    if (r != 0x0 && v != 0) {
      zeroSpinor<FT, V, S, compress12>(r, geom, zeroThreads);
      zeroSpinor<FT, V, S, compress12>(v, geom, zeroThreads);
      sUpdateThreads = 1;
      double alpha_cr[2] = {(double)1.0, (double)0.5};
      double start_time = omp_get_wtime();
      for (int i = 0; i < iters; i++) {
        bicgstab_s_update<FT, V, S, compress12>(
            alpha_cr, r, v, geom, sUpdateThreads);
      }
      double stop_time = omp_get_wtime();
      double best_time = stop_time - start_time;

      masterPrintf("tuneSUpdateThreads: threads = %d, current_time=%g (s)\n",
                   sUpdateThreads,
                   best_time);
      for (int threads = 2; threads <= geom.getNSIMT(); threads++) {
        start_time = omp_get_wtime();
        for (int i = 0; i < iters; i++) {
          bicgstab_s_update<FT, V, S, compress12>(alpha_cr, r, v, geom, threads);
        }
        stop_time = omp_get_wtime();
        double current_time = stop_time - start_time;

        masterPrintf(
            "tuneSUpdateThreads: threads = %d, current_time = %g (s), best=%g(s)\n",
            threads,
            current_time,
            best_time);
        if (current_time < best_time) {
          best_time = current_time;
          sUpdateThreads = threads;
        }
      }
    }
  }

  void tuneRXUpdateThreads(int iters)
  {
    if (r != 0x0 && v != 0 && p != 0 && t != 0) {
      zeroSpinor<FT, V, S, compress12>(r, geom, zeroThreads);
      zeroSpinor<FT, V, S, compress12>(v, geom, zeroThreads);
      zeroSpinor<FT, V, S, compress12>(t, geom, zeroThreads);
      zeroSpinor<FT, V, S, compress12>(p, geom, zeroThreads);
      double r_norm = 0;
      rxUpdateThreads = 1;
      double alpha_cr[2] = {(double)1.0, (double)0.5};
      double omega_cr[2] = {(double)-1.0, (double)-0.25};
      double start_time = omp_get_wtime();
      for (int i = 0; i < iters; i++) {
        bicgstab_rxupdate<FT, V, S, compress12>(
            r, v, t, p, omega_cr, alpha_cr, r_norm, geom, rxUpdateThreads);
      }
      double stop_time = omp_get_wtime();
      double best_time = stop_time - start_time;

      masterPrintf("tuneRXUpdateThreads: threads = %d, current_time=%g (s)\n",
                   rxUpdateThreads,
                   best_time);
      for (int threads = 2; threads <= geom.getNSIMT(); threads++) {
        start_time = omp_get_wtime();
        for (int i = 0; i < iters; i++) {
          bicgstab_rxupdate<FT, V, S, compress12>(
              r, v, t, p, omega_cr, alpha_cr, r_norm, geom, threads);
        }
        stop_time = omp_get_wtime();
        double current_time = stop_time - start_time;

        masterPrintf(
            "tuneRXUpdateThreads: threads = %d, current_time = %g (s), best=%g(s)\n",
            threads,
            current_time,
            best_time);
        if (current_time < best_time) {
          best_time = current_time;
          rxUpdateThreads = threads;
        }
      }
    }
  }
};
};
