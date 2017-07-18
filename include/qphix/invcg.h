#pragma once

/**
  \TODO user serviceable options -- should be moved to autoconf control
  */
#define QPHIX_TIMING_CG
// #define CGDEBUG

#include "qphix/linearOp.h"
#include "qphix/print_utils.h"
#include "qphix/blas_new_c.h"
#include "qphix/tsc.h"
#include "qphix/abs_solver.h"

namespace QPhiX
{

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
class InvCG : public AbstractSolver<FT,
                                    veclen,
                                    soalen,
                                    compress12,
                                    EvenOddLinearOperatorBase::num_flav>
{
 public:
  typedef typename Geometry<FT, veclen, soalen, compress12>::FourSpinorBlock Spinor;

  static constexpr int num_flav = EvenOddLinearOperatorBase::num_flav;

  InvCG(EvenOddLinearOperatorBase &M_, int MaxIters_)
      : M(M_), geom(M_.getGeometry()), MaxIters(MaxIters_)
  {
    // Length of vectors in floats, including PADs
    masterPrintf("Initializing CG Solver: Nvec=%d Ny=%d Nz=%d Nt=%d\n",
                 geom.nVecs(),
                 geom.Ny(),
                 geom.Nz(),
                 geom.Nt());

    for (uint8_t f = 0; f < num_flav; ++f)
      mp[f] = (Spinor *)geom.allocCBFourSpinor();
    for (uint8_t f = 0; f < num_flav; ++f)
      mmp[f] = (Spinor *)geom.allocCBFourSpinor();
    for (uint8_t f = 0; f < num_flav; ++f)
      p[f] = (Spinor *)geom.allocCBFourSpinor();
    for (uint8_t f = 0; f < num_flav; ++f)
      r[f] = (Spinor *)geom.allocCBFourSpinor();

    copyThreads = geom.getNSIMT();
    aypxThreads = geom.getNSIMT();
    xmyNormThreads = geom.getNSIMT();
    rmammpNorm2rxpapThreads = geom.getNSIMT();
    norm2Threads = geom.getNSIMT();
  }

  ~InvCG()
  {
    for (uint8_t f = 0; f < num_flav; ++f) {
      geom.free(mp[f]);
      geom.free(mmp[f]);
      geom.free(p[f]);
      geom.free(r[f]);
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

  void operator()(Spinor *const x[num_flav],
                          const Spinor *const rhs[num_flav],
                          const double RsdTarget,
                          int &n_iters,
                          double &rsd_sq_final,
                          unsigned long &site_flops,
                          unsigned long &mv_apps,
                          int isign,
                          bool verboseP,
                          int cb = 1,
			  QPhiX::ResiduumType residType=QPhiX::RELATIVE) const override
  {
    if (verboseP) {
      masterPrintf("Entering the CG inverter with num_flav=%d and cb=%d\n", num_flav, cb);
    }

#ifdef QPHIX_TIMING_CG
    TSC_tick cg_start;
    TSC_tick cg_end;
    TSC_tick cg_time = 0;

    // LINOP_TIMINGS
    TSC_tick lstart;
    TSC_tick lend;
    TSC_tick ltime = 0;

    // BLAS TIMINGS: norm2Spinor
    TSC_tick norm2_start;
    TSC_tick norm2_end;
    TSC_tick norm2_time = 0;

    TSC_tick blas_time = 0;
    // BLAS TIMINGS: xmyNorm
    TSC_tick xmyNorm_start;
    TSC_tick xmyNorm_end;
    TSC_tick xmyNorm_time = 0;

    // BLAS_TIMINGS: aypx
    TSC_tick aypx_start;
    TSC_tick aypx_end;
    TSC_tick aypx_time = 0;

    // BLAS_TIMINGS: copy
    TSC_tick copy_start;
    TSC_tick copy_end;
    TSC_tick copy_time = 0;

    // BLAS_TIMINGS: rMammpNorm2Rpmax
    TSC_tick rmammp_start;
    TSC_tick rmammp_end;
    TSC_tick rmammp_time = 0;

    TSC_tick missing_time;
#endif

    int n_cores;
    int n_simt;
    double cp;

#ifdef TIMING_CG
    CLOCK_NOW(cg_start);
#endif

    site_flops = 0;
    mv_apps = 0;
    n_cores = geom.getNumCores();
    n_simt = geom.getNSIMT();

#ifdef TIMING_CG
    CLOCK_NOW(norm2_start);
#endif
    double chi_sq;

    norm2Spinor<FT, veclen, soalen, compress12, num_flav>(
        chi_sq, rhs, geom, norm2Threads);
    masterPrintf("Chi_sq=%g RsdTarget=%g\n", chi_sq, RsdTarget);

#ifdef CGDEBUG
    norm2Spinor<FT, veclen, soalen, compress12, num_flav>(
        chi_sq, x, geom, norm2Threads);
    masterPrintf("||x||=%lf\n", chi_sq);

    norm2Spinor<FT, veclen, soalen, compress12, num_flav>(
        chi_sq, rhs, geom, norm2Threads);
    masterPrintf("Chi_sq=%g RsdTarget=%g\n", chi_sq, RsdTarget);
#endif

#ifdef TIMING_CG
    CLOCK_NOW(norm2_end);
    norm2_time += (norm2_end - norm2_start);
#endif

    site_flops += 4 * 12 * num_flav;
    double rsd_sq = (RsdTarget * RsdTarget);
    if( residType == QPhiX::RELATIVE )
    	rsd_sq *= chi_sq;

    double tmp_d;
// M^\dagger M psi
#ifdef TIMING_CG
    CLOCK_NOW(lstart);
#endif
    M(mp, x, isign, cb);

#ifdef CGDEBUG
    norm2Spinor<FT, veclen, soalen, compress12, num_flav>(
        tmp_d, mp, geom, norm2Threads);
    masterPrintf("M p = %lf \n", tmp_d);
#endif

    M(mmp, mp, -isign, cb);

#ifdef CGDEBUG
    norm2Spinor<FT, veclen, soalen, compress12, num_flav>(
        tmp_d, mmp, geom, norm2Threads);
    masterPrintf("MM p = %lf \n", tmp_d);
#endif

#ifdef TIMING_CG
    CLOCK_NOW(lend);
    ltime += (lend - lstart);
#endif
    mv_apps += 2;

// r = chi_internal - mmp,   cp=norm2(r)
#ifdef TIMING_CG
    CLOCK_NOW(xmyNorm_start);
#endif
    xmyNorm2Spinor<FT, veclen, soalen, compress12, num_flav>(
        r, rhs, mmp, cp, geom, xmyNormThreads);

#ifdef TIMING_CG
    CLOCK_NOW(xmyNorm_end);
    xmyNorm_time += (xmyNorm_end - xmyNorm_start);
#endif

    site_flops += 6 * 12 * num_flav;

    if (verboseP)
      masterPrintf("CG: r0 = %g   target = %g \n", cp, rsd_sq);

    if (cp <= rsd_sq) {
      n_iters = 0;
      rsd_sq_final = cp / chi_sq;

#ifdef TIMING_CG
      CLOCK_NOW(cg_end);
      cg_time = cg_end - cg_start;

      masterPrintf("CG Total time: %ld ticks", cg_time);
      masterPrintf("Linop Total time: %ld ticks %u applications, %g percent of CG\n",
                   ltime,
                   mv_apps,
                   100.0 * ((double)ltime) / ((double)cg_time));
      blas_time = norm2_time + xmyNorm_time + aypx_time + copy_time + rmammp_time;
      masterPrintf("BLAS Total time: %ld ticks, %g percent of CG\n",
                   blas_time,
                   100.0 * ((double)blas_time) / ((double)cg_time));
      masterPrintf("\t norm2_time: %ld ticks,   %g percent of CG\n",
                   norm2_time,
                   100.0 * ((double)norm2_time) / ((double)cg_time));
      masterPrintf("\t xmyNorm_time: %ld ticks,   %g percent of CG\n",
                   xmyNorm_time,
                   100.0 * ((double)xmyNorm_time) / ((double)cg_time));
      masterPrintf("\t aypx_time: %ld ticks,   %g percent of CG\n",
                   aypx_time,
                   100.0 * ((double)aypx_time) / ((double)cg_time));
      masterPrintf("\t rmammp_norm2r_xpap_time: %ld ticks,   %g percent of CG\n",
                   rmammp_time,
                   100.0 * ((double)rmammp_time) / ((double)cg_time));
      masterPrintf("\t copy_time: %ld ticks,   %g percent of CG\n",
                   copy_time,
                   100.0 * ((double)copy_time) / ((double)cg_time));
      missing_time = cg_time - ltime - blas_time;
      masterPrintf("Missing time = %ld ticks = %g percent of CG\n",
                   missing_time,
                   100.0 * (double)(missing_time) / ((double)cg_time));
#endif
      return;
    }

#ifdef TIMING_CG
    CLOCK_NOW(copy_start);
#endif

    copySpinor<FT, veclen, soalen, compress12, num_flav>(p, r, geom, copyThreads);

#ifdef TIMING_CG
    CLOCK_NOW(copy_end);
    copy_time += (copy_end - copy_start);
#endif

    double a, b, c, d;
    for (int k = 0; k <= MaxIters; ++k) {
      c = cp;

#ifdef TIMING_CG
      CLOCK_NOW(lstart);
#endif
      M(mp, p, isign, cb);
      M(mmp, mp, -isign, cb);

#ifdef TIMING_CG
      CLOCK_NOW(lend);
      ltime += (lend - lstart);
#endif
      mv_apps += 2;

#ifdef TIMING_CG
      CLOCK_NOW(norm2_start);
#endif

      norm2Spinor<FT, veclen, soalen, compress12, num_flav>(
          d, mp, geom, norm2Threads);

#ifdef TIMING_CG
      CLOCK_NOW(norm2_end);
      norm2_time += (norm2_end - norm2_start);
#endif

      site_flops += 4 * 12 * num_flav;

      a = c / d;

#ifdef CGDEBUG
      masterPrintf("CG:   iter %d c=%lf d=%lf a=%lf \n", k, c, d, a);
#endif

#ifdef TIMING_CG
      CLOCK_NOW(rmammp_start);
#endif

      rmammpNorm2rxpap<FT, veclen, soalen, compress12, num_flav>(
          r, a, mmp, cp, x, p, geom, rmammpNorm2rxpapThreads);

#ifdef CGDEBUG
      norm2Spinor<FT, veclen, soalen, compress12, num_flav>(
          tmp_d, r, geom, norm2Threads);
      masterPrintf("CG:   iter %d: r2 = %lf \n", k, tmp_d);

      norm2Spinor<FT, veclen, soalen, compress12, num_flav>(
          tmp_d, x, geom, norm2Threads);
      masterPrintf("CG:   iter %d: x2 = %lf \n", k, tmp_d);
#endif

#ifdef TIMING_CG
      CLOCK_NOW(rmammp_end);
      rmammp_time += (rmammp_end - rmammp_start);
#endif
      site_flops += 12 * 12 * num_flav;

      if (verboseP)
        masterPrintf("CG: iter %d:  r2 = %g   target = %g \n", k, cp, rsd_sq);

      if (cp < rsd_sq) {
        n_iters = k;

#ifdef TIMING_CG
        CLOCK_NOW(lstart);
        M(mp, x, isign, cb);
        M(mmp, mp, -isign, cb);
        CLOCK_NOW(lend);
        ltime += (lend - lstart);
        mv_apps += 2;

        CLOCK_NOW(xmyNorm_start);
        xmyNorm2Spinor<FT, veclen, soalen, compress12, num_flav>(
            r, rhs, mmp, cp, geom, xmyNormThreads);
        CLOCK_NOW(xmyNorm_end);
        xmyNorm_time += (xmyNorm_end - xmyNorm_start);

        site_flops += 6 * 12 * num_flav;
        rsd_sq_final = cp;
        CLOCK_NOW(cg_end);
        cg_time = cg_end - cg_start;

        masterPrintf("Iters=%d Final Rsd = %e Final RsdSq=%e Final ResidRelSq=%e "
                     "Final ResidRel=%e\n",
                     n_iters,
                     sqrt(rsd_sq_final),
                     rsd_sq_final,
                     rsd_sq_final / chi_sq,
                     sqrt(rsd_sq_final / chi_sq));
        masterPrintf("CG Total time: %ld ticks\n", cg_time);
        masterPrintf(
            "Linop Total time: %ld ticks %u applications, %g percent of CG\n",
            ltime,
            mv_apps,
            100.0 * ((double)ltime) / ((double)cg_time));
        blas_time = norm2_time + xmyNorm_time + aypx_time + copy_time + rmammp_time;
        masterPrintf("BLAS Total time: %ld ticks, %g percent of CG\n",
                     blas_time,
                     100.0 * ((double)blas_time) / ((double)cg_time));
        masterPrintf("\t norm2_time: %ld ticks,   %g percent of CG\n",
                     norm2_time,
                     100.0 * ((double)norm2_time) / ((double)cg_time));
        masterPrintf("\t xmyNorm_time: %ld ticks,   %g percent of CG\n",
                     xmyNorm_time,
                     100.0 * ((double)xmyNorm_time) / ((double)cg_time));
        masterPrintf("\t aypx_time: %ld ticks,   %g percent of CG\n",
                     aypx_time,
                     100.0 * ((double)aypx_time) / ((double)cg_time));
        masterPrintf("\t rmammp_norm2r_xpap_time: %ld ticks,   %g percent of CG\n",
                     rmammp_time,
                     100.0 * ((double)rmammp_time) / ((double)cg_time));
        masterPrintf("\t copy_time: %ld ticks,   %g percent of CG\n",
                     copy_time,
                     100.0 * ((double)copy_time) / ((double)cg_time));
        missing_time = cg_time - ltime - blas_time;
        masterPrintf("Missing time = %ld ticks = %g percent of CG \n",
                     missing_time,
                     100.0 * (double)(missing_time) / ((double)cg_time));
#endif

        return;
      }

      b = cp / c;

#ifdef TIMING_CG
      CLOCK_NOW(aypx_start);
#endif
      aypx<FT, veclen, soalen, compress12, num_flav>(b, r, p, geom, aypxThreads);

#ifdef TIMING_CG
      CLOCK_NOW(aypx_end);
      aypx_time += (aypx_end - aypx_start);
#endif
      site_flops += 4 * 12 * num_flav;
    }

    n_iters = MaxIters;
    rsd_sq_final = cp;

#ifdef TIMING_CG
    CLOCK_NOW(cg_end);
    cg_time = cg_end - cg_start;

    masterPrintf("CG Total time: %ld ticks\n", cg_time);
    masterPrintf("Linop Total time: %ld ticks %u applications, %g percent of CG\n",
                 ltime,
                 mv_apps,
                 100.0 * ((double)ltime) / ((double)cg_time));
    blas_time = norm2_time + xmyNorm_time + aypx_time + copy_time + rmammp_time;
    masterPrintf("BLAS Total time: %ld ticks, %g percent of CG\n",
                 blas_time,
                 100.0 * ((double)blas_time) / ((double)cg_time));
    masterPrintf("\t norm2_time: %ld ticks,   %g percent of CG\n",
                 norm2_time,
                 100.0 * ((double)norm2_time) / ((double)cg_time));
    masterPrintf("\t xmyNorm_time: %ld ticks,   %g percent of CG\n",
                 xmyNorm_time,
                 100.0 * ((double)xmyNorm_time) / ((double)cg_time));
    masterPrintf("\t aypx_time: %ld ticks,   %g percent of CG\n",
                 aypx_time,
                 100.0 * ((double)aypx_time) / ((double)cg_time));
    masterPrintf("\t rmammp_norm2r_xpap_time: %ld ticks,   %g percent of CG\n",
                 rmammp_time,
                 100.0 * ((double)rmammp_time) / ((double)cg_time));
    masterPrintf("\t copy_time: %ld ticks,   %g percent of CG\n",
                 copy_time,
                 100.0 * ((double)copy_time) / ((double)cg_time));
    missing_time = cg_time - ltime - blas_time;
    masterPrintf("Missing time = %ld ticks = %g percent of CG\n",
                 missing_time,
                 100.0 * (double)(missing_time) / ((double)cg_time));
#endif
    return;
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

  Spinor *mp[num_flav];
  Spinor *mmp[num_flav];
  Spinor *p[num_flav];
  Spinor *r[num_flav];

  int copyThreads;
  int aypxThreads;
  int xmyNormThreads;
  int rmammpNorm2rxpapThreads;
  int norm2Threads;

  void tuneCopyThreads(int iters)
  {
    if (r != 0x0 && p != 0) {
      // Do first with 1 thread

      zeroSpinor<FT, veclen, soalen, compress12>(r, geom, 1);
      zeroSpinor<FT, veclen, soalen, compress12>(p, geom, 1);
      copyThreads = 1;
      double start_time = omp_get_wtime();
      for (int i = 0; i < iters; i++) {
        copySpinor<FT, veclen, soalen, compress12>(r, p, geom, copyThreads);
      }

      double stop_time = omp_get_wtime();
      double best_time = stop_time - start_time;
      masterPrintf("tuneCopyThreads: threads = %d, current_time=%g (s)\n",
                   copyThreads,
                   best_time);

      for (int threads = 2; threads <= geom.getNSIMT(); threads++) {
        start_time = omp_get_wtime();
        for (int i = 0; i < iters; i++) {
          copySpinor<FT, veclen, soalen, compress12>(r, p, geom, threads);
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

  void tuneAypxThreads(int iters)
  {
    if (r != 0x0 && p != 0) {
      // Do first with 1 thread

      zeroSpinor<FT, veclen, soalen, compress12>(r, geom, 1);
      zeroSpinor<FT, veclen, soalen, compress12>(p, geom, 1);
      aypxThreads = 1;
      double br = (double)0.5;

      double start_time = omp_get_wtime();
      for (int i = 0; i < iters; i++) {
        aypx<FT, veclen, soalen, compress12>(br, r, p, geom, aypxThreads);
      }

      double stop_time = omp_get_wtime();
      double best_time = stop_time - start_time;
      masterPrintf("tuneAypxThreads: threads = %d, current_time=%g (s)\n",
                   aypxThreads,
                   best_time);

      for (int threads = 2; threads <= geom.getNSIMT(); threads++) {
        start_time = omp_get_wtime();
        for (int i = 0; i < iters; i++) {
          aypx<FT, veclen, soalen, compress12>(br, r, p, geom, threads);
        }
        stop_time = omp_get_wtime();
        double current_time = stop_time - start_time;

        masterPrintf(
            "tuneAypxThreads: threads = %d, current_time = %g (s), best=%g(s)\n",
            threads,
            current_time,
            best_time);
        if (current_time < best_time) {
          best_time = current_time;
          aypxThreads = threads;
        }
      }
    }
  }

  void tuneNorm2Threads(int iters)
  {
    if (r != 0x0) {
      // Do first with 1 thread
      double rnorm;
      zeroSpinor<FT, veclen, soalen, compress12>(r, geom, 1);
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

  void tuneXMYNorm2Threads(int iters)
  {
    if (r != 0x0 && p != 0 && mmp != 0) {
      // Do first with 1 thread

      zeroSpinor<FT, veclen, soalen, compress12>(r, geom, 1);
      zeroSpinor<FT, veclen, soalen, compress12>(p, geom, 1);
      zeroSpinor<FT, veclen, soalen, compress12>(mmp, geom, 1);
      xmyNormThreads = 1;
      double reduction = 0;

      double start_time = omp_get_wtime();
      for (int i = 0; i < iters; i++) {
        xmyNorm2Spinor<FT, veclen, soalen, compress12>(
            r, p, mmp, reduction, geom, xmyNormThreads);
      }
      double stop_time = omp_get_wtime();
      double best_time = stop_time - start_time;
      masterPrintf("tuneXmyNorm2Threads: threads = %d, current_time=%g (s)\n",
                   xmyNormThreads,
                   best_time);
      for (int threads = 2; threads <= geom.getNSIMT(); threads++) {
        start_time = omp_get_wtime();
        for (int i = 0; i < iters; i++) {
          xmyNorm2Spinor<FT, veclen, soalen, compress12>(
              r, p, mmp, reduction, geom, threads);
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
          xmyNormThreads = threads;
        }
      }
    }
  }

  void tuneRXUpdateThreads(int iters)
  {
    if (r != 0x0 && p != 0 && mp != 0 && mmp != 0) {
      zeroSpinor<FT, veclen, soalen, compress12>(r, geom, 1);
      zeroSpinor<FT, veclen, soalen, compress12>(p, geom, 1);
      zeroSpinor<FT, veclen, soalen, compress12>(mp, geom, 1);
      zeroSpinor<FT, veclen, soalen, compress12>(mmp, geom, 1);
      double r_norm = 0;
      double a = (double)0.5;
      int rmammpNorm2rxpapThreads = 1;
      double start_time = omp_get_wtime();
      for (int i = 0; i < iters; i++) {
        rmammpNorm2rxpap<FT, veclen, soalen, compress12>(
            r, a, mmp, r_norm, mp, p, geom, rmammpNorm2rxpapThreads);
      }
      double stop_time = omp_get_wtime();
      double best_time = stop_time - start_time;

      masterPrintf("tuneRXUpdateThreads: threads = %d, current_time=%g (s)\n",
                   rmammpNorm2rxpapThreads,
                   best_time);
      for (int threads = 2; threads <= geom.getNSIMT(); threads++) {
        start_time = omp_get_wtime();
        for (int i = 0; i < iters; i++) {
          rmammpNorm2rxpap<FT, veclen, soalen, compress12>(
              r, a, mmp, r_norm, mp, p, geom, rmammpNorm2rxpapThreads);
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
          rmammpNorm2rxpapThreads = threads;
        }
      }
    }
  }
};
};
