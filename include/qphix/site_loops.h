#pragma once

#include <qphix/geometry.h>
#include <qphix/comm.h>

#include <omp.h>
#include "qphix/print_utils.h"
#include "qphix/thread_limits.h"
#include <array>
#include <numeric>

namespace QPhiX
{

/*#if defined (__GNUG__) && !defined (__INTEL_COMPILER)
static double new_norm_array[MAX_THREADS][MAXV] __attribute__
((aligned(QPHIX_LLC_CACHE_ALIGN)));
static double new_iprod_array[MAX_THREADS][2][MAXV]  __attribute__
((aligned(QPHIX_LLC_CACHE_ALIGN)));
#else
__declspec(align(QPHIX_LLC_CACHE_ALIGN)) static double
new_norm_array[MAX_THREADS][MAXV];
__declspec(align(QPHIX_LLC_CACHE_ALIGN)) static double
new_iprod_array[MAX_THREADS][2][MAXV];
#endif*/

template <typename FT, int V, int S, bool compress, typename SpinorFunctor>
void siteLoopNoReduction(const SpinorFunctor theFunctor,
                         const Geometry<FT, V, S, compress> &geom,
                         int n_blas_simt)
{

  const int n_simt = geom.getNSIMT();
  const int n_cores = geom.getNumCores();
  const int Nxh = geom.Nxh();
  const int Ny = geom.Ny();
  const int Nz = geom.Nz();
  const int Nt = geom.Nt();
  const int Pxy = geom.getPxy();
  const int Pxyz = geom.getPxyz();
  const int Nvecs = geom.nVecs();

  // This is the total number of spinor vectors
  const int n_soavec = (Nxh * Ny * Nz * Nt) / S;

#pragma omp parallel
  {
    // Self ID
    int tid = omp_get_thread_num();
    int cid = tid / n_simt;
    int smtid = tid - n_simt * cid;
    int num_bthreads = n_cores * n_blas_simt;

    if (smtid < n_blas_simt) {
      int btid = smtid + n_blas_simt * cid;
      int n_soavec_per_core = n_soavec / n_cores;
      if (n_soavec % n_cores != 0)
        n_soavec_per_core++;

      int low = cid * n_soavec_per_core;
      int next = (cid + 1) * n_soavec_per_core;
      int hi = n_soavec < next ? n_soavec : next;

      // Loop over spinors. Local
      for (int spinor_i = low + smtid; spinor_i < hi; spinor_i += n_blas_simt) {

        // Decompose spinor_i into vec, ybase, z and t coords
        int t1 = spinor_i / Nvecs;
        int vec = spinor_i - Nvecs * t1;

        int t2 = t1 / Ny;
        int y = t1 - Ny * t2;

        int t = t2 / Nz;
        int z = t2 - Nz * t;

        // ztbase for padding
        int ztbase = Pxy * z + Pxyz * t;
        int block = vec + Nvecs * y + ztbase;
        theFunctor.func(block);
      }
    }
  }
}

template <typename FT, int V, int S, bool compress, typename Reduce1Functor>
void siteLoop1Reduction(const Reduce1Functor theFunctor,
                        double &reduction,
                        const Geometry<FT, V, S, compress> &geom,
                        int n_blas_simt, bool globalSum=true)
{

  const int n_simt = geom.getNSIMT();
  const int n_cores = geom.getNumCores();
  const int Nxh = geom.Nxh();
  const int Ny = geom.Ny();
  const int Nz = geom.Nz();
  const int Nt = geom.Nt();
  const int Pxy = geom.getPxy();
  const int Pxyz = geom.getPxyz();
  const int Nvecs = geom.nVecs();

  // This is the total number of spinor vectors
  const int n_soavec = (Nxh * Ny * Nz * Nt) / S;

  double reduction0 = 0.0;
#pragma omp parallel reduction(+: reduction0)
  {
    // Self ID
    int tid = omp_get_thread_num();
    int cid = tid / n_simt;
    int smtid = tid - n_simt * cid;
    int num_bthreads = n_cores * n_blas_simt;

    if (smtid < n_blas_simt) {
      int btid = smtid + n_blas_simt * cid;

      int n_soavec_per_core = n_soavec / n_cores;
      if (n_soavec % n_cores != 0)
        n_soavec_per_core++;

      int low = cid * n_soavec_per_core;
      int next = (cid + 1) * n_soavec_per_core;
      int hi = n_soavec < next ? n_soavec : next;

      // Loop over spinors. Local
      std::array<double, S> r0; r0.fill(0.0);
      for (int spinor_i = low + smtid; spinor_i < hi; spinor_i += n_blas_simt) {

        // Decompose spinor_i into vec, ybase, z and t coords
        int t1 = spinor_i / Nvecs;
        int vec = spinor_i - Nvecs * t1;

        int t2 = t1 / Ny;
        int y = t1 - Ny * t2;

        int t = t2 / Nz;
        int z = t2 - Nz * t;

        // ztbase for padding
        int ztbase = Pxy * z + Pxyz * t;
        int block = vec + Nvecs * y + ztbase;
        theFunctor.func(block, r0);
  }
      reduction0 = std::accumulate(r0.begin(), r0.end(), 0.0);
    }
  }

  // If globalSum enabled by default is false each node will have its local reduction
  if( globalSum ) {
    CommsUtils::sumDouble(&reduction0);
  }

  reduction = reduction0;
}

template <typename FT, int V, int S, bool compress, typename Reduce2Functor>
void siteLoop2Reductions(const Reduce2Functor theFunctor,
                         double reduction[2],
                         const Geometry<FT, V, S, compress> &geom,
                         int n_blas_simt,
                         bool globalSum=true)
{

  const int n_simt = geom.getNSIMT();
  const int n_cores = geom.getNumCores();
  const int Nxh = geom.Nxh();
  const int Ny = geom.Ny();
  const int Nz = geom.Nz();
  const int Nt = geom.Nt();
  const int Pxy = geom.getPxy();
  const int Pxyz = geom.getPxyz();
  const int Nvecs = geom.nVecs();

  // This is the total number of spinor vectors
  const int n_soavec = (Nxh * Ny * Nz * Nt) / S;

  double reduction0 = 0.0, reduction1 = 0.0;
#pragma omp parallel reduction(+: reduction0, reduction1)
  {
    // Self ID
    int tid = omp_get_thread_num();
    int cid = tid / n_simt;
    int smtid = tid - n_simt * cid;
    int num_bthreads = n_cores * n_blas_simt;

    if (smtid < n_blas_simt) {
      int btid = smtid + n_blas_simt * cid;

      int n_soavec_per_core = n_soavec / n_cores;
      if (n_soavec % n_cores != 0)
        n_soavec_per_core++;

      int low = cid * n_soavec_per_core;
      int next = (cid + 1) * n_soavec_per_core;
      int hi = n_soavec < next ? n_soavec : next;

      // Loop over spinors. Local
      std::array<double, S> r0; r0.fill(0.0);
      std::array<double, S> r1; r1.fill(0.0);
      for (int spinor_i = low + smtid; spinor_i < hi; spinor_i += n_blas_simt) {

        // Decompose spinor_i into vec, ybase, z and t coords
        int t1 = spinor_i / Nvecs;
        int vec = spinor_i - Nvecs * t1;

        int t2 = t1 / Ny;
        int y = t1 - Ny * t2;

        int t = t2 / Nz;
        int z = t2 - Nz * t;

        // ztbase for padding
        int ztbase = Pxy * z + Pxyz * t;
        int block = vec + Nvecs * y + ztbase;
        theFunctor.func(block, r0, r1);
      }
      reduction0 = std::accumulate(r0.begin(), r0.end(), 0.0);
      reduction1 = std::accumulate(r1.begin(), r1.end(), 0.0);
    }
  }

  reduction[0] = reduction0;
  reduction[1] = reduction1;

  // DO A GLOBAL SUM HERE
  if( globalSum ) {
    CommsUtils::sumDoubleArray(reduction, 2);
  }
}

template <typename FT, int V, int S, bool compress, typename Reduce3Functor>
void siteLoop3Reductions(const Reduce3Functor theFunctor,
                         double reduction[3],
                         const Geometry<FT, V, S, compress> &geom,
                         int n_blas_simt, bool globalSum=true)
{

  const int n_simt = geom.getNSIMT();
  const int n_cores = geom.getNumCores();
  const int Nxh = geom.Nxh();
  const int Ny = geom.Ny();
  const int Nz = geom.Nz();
  const int Nt = geom.Nt();
  const int Pxy = geom.getPxy();
  const int Pxyz = geom.getPxyz();
  const int Nvecs = geom.nVecs();

  // This is the total number of spinor vectors
  const int n_soavec = (Nxh * Ny * Nz * Nt) / S;

  double reduction0 = 0.0, reduction1 = 0.0, reduction2 = 0.0;
#pragma omp parallel reduction(+: reduction0, reduction1, reduction2)
  {
    // Self ID
    int tid = omp_get_thread_num();
    int cid = tid / n_simt;
    int smtid = tid - n_simt * cid;
    int num_bthreads = n_cores * n_blas_simt;

    if (smtid < n_blas_simt) {
      int btid = smtid + n_blas_simt * cid;

      int n_soavec_per_core = n_soavec / n_cores;
      if (n_soavec % n_cores != 0)
        n_soavec_per_core++;

      int low = cid * n_soavec_per_core;
      int next = (cid + 1) * n_soavec_per_core;
      int hi = n_soavec < next ? n_soavec : next;

      // Loop over spinors. Local
      std::array<double, S> r0; r0.fill(0.0);
      std::array<double, S> r1; r1.fill(0.0);
      std::array<double, S> r2; r2.fill(0.0);
      for (int spinor_i = low + smtid; spinor_i < hi; spinor_i += n_blas_simt) {

        // Decompose spinor_i into vec, ybase, z and t coords
        int t1 = spinor_i / Nvecs;
        int vec = spinor_i - Nvecs * t1;

        int t2 = t1 / Ny;
        int y = t1 - Ny * t2;

        int t = t2 / Nz;
        int z = t2 - Nz * t;

        // ztbase for padding
        int ztbase = Pxy * z + Pxyz * t;
        int block = vec + Nvecs * y + ztbase;
        theFunctor.func(block, r0, r1, r2);
      }
      reduction0 = std::accumulate(r0.begin(), r0.end(), 0.0);
      reduction1 = std::accumulate(r1.begin(), r1.end(), 0.0);
      reduction2 = std::accumulate(r2.begin(), r2.end(), 0.0);
    }
  }

  // Horizontally sum the btid=0 result.
  reduction[0] = reduction0;
  reduction[1] = reduction1;
  reduction[2] = reduction2;

  // DO A GLOBAL SUM HERE
  if( globalSum )  {
    CommsUtils::sumDoubleArray(reduction, 3);
  }
}

}; // Namespace
