#pragma once

#include <iostream>
#include "qphix/qphix_config.h"
#include "qphix/dslash_utils.h"
#include "qphix/kernel_selector.h"
#include "qphix/print_utils.h"
#include <immintrin.h>
#include <omp.h>

#include "qphix_codegen/dslash_generated.h"
#include "qphix_codegen/tm_dslash_generated.h"

namespace QPhiX
{
template <typename FT, int veclen, int soalen, bool compress12>
TMDslash<FT, veclen, soalen, compress12>::TMDslash(
    Geometry<FT, veclen, soalen, compress12> *geom_,
    double t_boundary_,
    double aniso_coeff_S_,
    double aniso_coeff_T_,
    double Mass_,
    double TwistedMass_,
    bool use_tbc_[4],
    double tbc_phases_[4][2])
    : s(geom_), comms(new Comms<FT, veclen, soalen, compress12>(geom_)),
      By(geom_->getBy()), Bz(geom_->getBz()), NCores(geom_->getNumCores()),
      Sy(geom_->getSy()), Sz(geom_->getSz()), PadXY(geom_->getPadXY()),
      PadXYZ(geom_->getPadXYZ()), MinCt(geom_->getMinCt()),
      n_threads_per_core(geom_->getSy() * geom_->getSz()), t_boundary(t_boundary_),
      aniso_coeff_S(aniso_coeff_S_), aniso_coeff_T(aniso_coeff_T_), Mass(Mass_),
      TwistedMass(TwistedMass_), derived_mu(TwistedMass / (4.0 + Mass)),
      derived_mu_inv((4.0 + Mass) /
                     (TwistedMass * TwistedMass + (4.0 + Mass) * (4.0 + Mass)))
{
  if (use_tbc_ != nullptr && tbc_phases_ != nullptr) {
    for (int dim = 0; dim < 4; ++dim) {
      use_tbc[dim] = use_tbc_[dim];

      for (int dir : {0, 1}) {
        tbc_phases[dim][dir] = rep<FT, double>(tbc_phases_[dim][dir]);
      }
    }
  }

  // OK we need to set up log of veclen
  log2veclen = 0;
  int veclen_bits = veclen;
  while (veclen_bits > 1) {
    log2veclen++;
    veclen_bits = (veclen_bits >> 1);
  }

  log2soalen = 0;
  int soalen_bits = soalen;
  while (soalen_bits > 1) {
    log2soalen++;
    soalen_bits = (soalen_bits >> 1);
  }

  int Nxh = s->Nxh();
  int Ny = s->Ny();
  int Nz = s->Nz();
  int Nt = s->Nt();

  // Get info about our location in the time dimension
  amIPtMin = comms->amIPtMin();
  amIPtMax = comms->amIPtMax();

  // We must have Nxh be divisible by soalen
  if (Nxh % soalen != 0) {
    printf("X length after checkerboarding (%d) must be divisible by soalen (%d)\n",
           Nxh,
           soalen);
    std::abort();
  }

  // We must have Ny be divisible by nGY (ratio of VECLEN to SOALEN)
  int ngy = s->nGY();
  if (Ny % ngy != 0) {
    printf(
        "Y length (%d) must be divisible by ratio of VECLEN/SOALEN=%d\n", Ny, ngy);
    std::abort();
  }

  if (Sy > By / ngy) {
    printf("Warning Sy > By/nyg. Some threads may be idle\n");
  }

  // sanity: we ought to have at least 'expected_no of threads available'

  int expected_threads = NCores * n_threads_per_core;
  if (expected_threads != omp_get_max_threads()) {
    std::cout << "Expected (Cores per Socket x Threads per Core)="
              << expected_threads << " but found " << omp_get_max_threads() << "..."
              << std::endl;
    std::cout
        << "Check your OMP_NUM_THREADS or QMT_NUM_THREADS env variable, or the "
           "CORES_PER_SOCKET and THREADS_PER_CORE env variables"
        << std::endl;
    std::abort();
  }

  int nvecs = s->nVecs();
  int num_cores = s->getNumCores();
  int num_phases = s->getNumPhases();

  // Allocate barriers
  gBar = new Barrier(num_cores * Sy * Sz, n_threads_per_core);
  barriers = new Barrier **[num_phases];
  for (int ph = 0; ph < num_phases; ph++) {
    const CorePhase &phase = s->getCorePhase(ph);
    barriers[ph] = new Barrier *[phase.Ct];
    for (int i = 0; i < phase.Ct; i++) {
      barriers[ph][i] = new Barrier(phase.Cyz * Sy * Sz, n_threads_per_core);
    }
  }

// Set up barriers
#pragma omp parallel
  {
    int tid = omp_get_thread_num();
    gBar->init(tid);
    int cid = tid / n_threads_per_core;
    for (int ph = 0; ph < num_phases; ph++) {
      const CorePhase &phase = s->getCorePhase(ph);
      int nActiveCores = phase.Cyz * phase.Ct;
      if (cid >= nActiveCores)
        continue;
      int cid_t = cid / phase.Cyz;
      int ngroup = phase.Cyz * Sy * Sz;
      int group_tid = tid % ngroup;
      barriers[ph][cid_t]->init(group_tid);
    }
  }

  // Set up block info array. These are thread local
  // Indices run as (phases fastest and threads slowest)
  int num_blockinfo = num_phases * s->getNumCores() * n_threads_per_core;
  block_info = (BlockPhase *)ALIGNED_MALLOC(num_blockinfo * sizeof(BlockPhase), 64);
  if (block_info == 0x0) {
    fprintf(stderr, "Could not allocate Block Info array\n");
    std::abort();
  }

  // Set up blockinfo
  masterPrintf("Setting Up Blockinfo array: num_phases=%d\n", num_phases);
#pragma omp parallel shared(num_phases)
  {
    int tid = omp_get_thread_num();
    int cid = tid / n_threads_per_core;
    int smtid = tid - n_threads_per_core * cid;
    int ly = Ny / By;

    for (int ph = 0; ph < num_phases; ph++) {
      const CorePhase &phase = s->getCorePhase(ph);
      BlockPhase &binfo = block_info[num_phases * tid + ph];

      int nActiveCores = phase.Cyz * phase.Ct;
      if (cid > nActiveCores)
        continue;
      binfo.cid_t = cid / phase.Cyz;
      binfo.cid_yz = cid - binfo.cid_t * phase.Cyz;
      int syz = phase.startBlock + binfo.cid_yz;
      binfo.bz = syz / ly;
      binfo.by = syz - binfo.bz * ly;
      binfo.bt = (Nt * binfo.cid_t) / phase.Ct;
      binfo.nt = (Nt * (binfo.cid_t + 1)) / phase.Ct - binfo.bt;
      binfo.by *= By;
      binfo.bz *= Bz;
      int ngroup = phase.Cyz * Sy * Sz;
      binfo.group_tid = tid % ngroup;
    }
  } // OMP parallel
  masterPrintf("Phase info set up\n");

  masterPrintf("Precomputing offsets\n");
  // Alloc tmpspc. It is thread local so we need one for
  // every thread
  size_t tmpspc_size = num_cores * n_threads_per_core * veclen * 16 * sizeof(int);
  tmpspc_all = (int *)ALIGNED_MALLOC(tmpspc_size, 64);
  if (tmpspc_all == 0x0) {
    printf("Failed to allocate xy offset tmpspc\n");
    std::abort();
  }

#pragma omp parallel
  {
    int tid = omp_get_thread_num();
    int *tmpspc = &(tmpspc_all[veclen * 16 * tid]);

    int *offs, *xbOffs, *xfOffs, *ybOffs, *yfOffs, *gOffs, *pfyOffs;
    int *xbOffs_xodd[2], *xbOffs_x0_xodd[2];
    int *xfOffs_xodd[2], *xfOffs_xn_xodd[2];
    int *ybOffs_yn0, *ybOffs_y0, *yfOffs_ynn, *yfOffs_yn;

    // Why is this needed. tmpspc should already be aligned.
    int *atmp = (int *)((((unsigned long long)tmpspc) + 0x3F) & ~0x3F);
    offs = &atmp[0];
    xbOffs_xodd[0] = &atmp[veclen * 1];
    xbOffs_xodd[1] = &atmp[veclen * 2];
    xbOffs_x0_xodd[0] = &atmp[veclen * 3];
    xbOffs_x0_xodd[1] = &atmp[veclen * 4];
    xfOffs_xodd[0] = &atmp[veclen * 5];
    xfOffs_xodd[1] = &atmp[veclen * 6];
    xfOffs_xn_xodd[0] = &atmp[veclen * 7];
    xfOffs_xn_xodd[1] = &atmp[veclen * 8];
    ybOffs_yn0 = &atmp[veclen * 9];
    ybOffs_y0 = &atmp[veclen * 10];
    yfOffs_ynn = &atmp[veclen * 11];
    yfOffs_yn = &atmp[veclen * 12];
    gOffs = &atmp[veclen * 13];
    pfyOffs = &atmp[veclen * 14];

    int nvec = s->nVecs();
    int nyg = s->nGY();
    const int gauge_line_in_floats =
        sizeof(SU3MatrixBlock) / sizeof(FT); // One gauge scanline, in floats
    const int spinor_line_in_floats =
        sizeof(FourSpinorBlock) / sizeof(FT); //  One spinor scanline, in floats

    if (tid == 0) {
      // Initialize masks
      xbmask_x0_xodd[0] = -1;
      xbmask_x0_xodd[1] = -1;
      xfmask_xn_xodd[0] = -1;
      xfmask_xn_xodd[1] = -1;
      ybmask_y0 = -1;
      yfmask_yn = -1;
    }

    for (int y = 0; y < nyg; y++) {
      // Various indexing things
      int ind = y * soalen;
      int X = nvecs * y * spinor_line_in_floats;
      int y1 = y & 1;
      int y2 = 1 - y1;
      for (int x = 0; x < soalen; x++) {
        xbOffs_x0_xodd[y1][ind] = X + x - 1;
        xbOffs_xodd[y1][ind] = X + x - 1;
        if (x == 0) {
          if (comms->localX()) {
            xbOffs_x0_xodd[y1][ind] -=
                (spinor_line_in_floats - soalen - nvecs * spinor_line_in_floats);
          } else {
            xbOffs_x0_xodd[y1][ind] +=
                soalen; // This lane is disabled, just set it within same cache line
            if (tid == 0)
              xbmask_x0_xodd[y1] &= ~(1 << ind); // reset a bit in the mask
          }
          xbOffs_xodd[y1][ind] -= (spinor_line_in_floats - soalen);
        }
        xfOffs_xodd[y1][ind] = X + x;
        xfOffs_xn_xodd[y1][ind] = X + x;

        xbOffs_x0_xodd[y2][ind] = X + x;
        xbOffs_xodd[y2][ind] = X + x;
        xfOffs_xodd[y2][ind] = X + x + 1;
        xfOffs_xn_xodd[y2][ind] = X + x + 1;
        if (x == soalen - 1) {
          xfOffs_xodd[y2][ind] += (spinor_line_in_floats - soalen);
          if (comms->localX()) {
            xfOffs_xn_xodd[y2][ind] +=
                (spinor_line_in_floats - soalen - nvecs * spinor_line_in_floats);
          } else {
            xfOffs_xn_xodd[y2][ind] -=
                soalen; // This lane is disabled, just set it within same cache line
            if (tid == 0)
              xfmask_xn_xodd[y2] &= ~(1 << ind); // reset the ind bit in the mask
          }
        }

        ybOffs_y0[ind] = X - nvecs * spinor_line_in_floats +
                         x; // previous y-neighbor site offsets
        if (y == 0) {
          if (comms->localY()) {
            ybOffs_y0[ind] += Ny * nvecs * spinor_line_in_floats;
          } else {
            ybOffs_y0[ind] =
                X + x; // This lane is disabled, just set it within same cache line
            if (tid == 0)
              ybmask_y0 &= ~(1 << ind); // reset the ind bit in the mask
          }
        }
        ybOffs_yn0[ind] = X - nvecs * spinor_line_in_floats +
                          x; // previous y-neighbor site offsets
        yfOffs_yn[ind] =
            X + nvecs * spinor_line_in_floats + x; // next y-neighbor site offsets
        if (y == nyg - 1) {
          if (comms->localY()) {
            yfOffs_yn[ind] -= Ny * nvecs * spinor_line_in_floats;
          } else {
            yfOffs_yn[ind] =
                X + x; // This lane is disabled, just set it within same cache line
            if (tid == 0)
              yfmask_yn &= ~(1 << ind); // reset the ind bit in the mask
          }
        }
        yfOffs_ynn[ind] =
            X + nvecs * spinor_line_in_floats + x; // next y-neighbor site offsets
        offs[ind] = X + x; // site offsets for z & t neighbors

        gOffs[ind] = ind; // this not used really

        ind++;
      }
    }
  } // OMP parallel
  // Info
  if (compress12) {
    masterPrintf("WILL Use 12 compression\n");
  }
}

template <typename FT, int veclen, int soalen, bool compress12>
TMDslash<FT, veclen, soalen, compress12>::~TMDslash()
{
  masterPrintf("Freeing BlockInfo\n");
  ALIGNED_FREE(block_info);

  masterPrintf("Freeing tmpspc_all\n");
  ALIGNED_FREE(tmpspc_all);

  masterPrintf("Freeing gbar\n");
  delete gBar;

  masterPrintf("Deleting Barriers\n");
  for (int ph = 0; ph < s->getNumPhases(); ph++) {
    const CorePhase &phase = s->getCorePhase(ph);
    for (int i = 0; i < phase.Ct; i++) {
      delete barriers[ph][i];
    }
    delete[] barriers[ph];
  }
  delete[] barriers;

  masterPrintf("Deleting Comms\n");
  delete comms;

  masterPrintf("All Destructed\n");
}

template <typename FT, int veclen, int soalen, bool compress12>
void TMDslash<FT, veclen, soalen, compress12>::dslash(FourSpinorBlock *res,
                                                      const FourSpinorBlock *psi,
                                                      const SU3MatrixBlock *u,
                                                      int isign,
                                                      int cb)
{
  TMDPsi(u, psi, res, isign == 1, cb);
}

template <typename FT, int veclen, int soalen, bool compress12>
void TMDslash<FT, veclen, soalen, compress12>::dslashAChiMinusBDPsi(
    FourSpinorBlock *res,
    const FourSpinorBlock *psi,
    const FourSpinorBlock *chi,
    const SU3MatrixBlock *u,
    double alpha,
    double beta,
    int isign,
    int cb,
    const double musign)
{
  TMDPsiAChiMinusBDPsi(u, psi, chi, res, alpha, beta, isign == 1, cb, musign);
}

template <typename FT, int veclen, int soalen, bool compress12>
void TMDslash<FT, veclen, soalen, compress12>::TMDyz(int tid,
                                                     const FourSpinorBlock *psi,
                                                     FourSpinorBlock *res,
                                                     const SU3MatrixBlock *u,
                                                     bool const is_plus,
                                                     int cb)
{
  auto kernel =
      (is_plus
           ? QPHIX_KERNEL_SELECT(tm_dslash_plus_vec, FT, veclen, soalen, compress12)
           : QPHIX_KERNEL_SELECT(
                 tm_dslash_minus_vec, FT, veclen, soalen, compress12));

  const int Nxh = s->Nxh();
  const int Nx = s->Nx();
  const int Ny = s->Ny();
  const int Nz = s->Nz();
  const int Nt = s->Nt();
  const int By = s->getBy();
  const int Bz = s->getBz();
  const int Sy = s->getSy();
  const int Sz = s->getSz();
  const int ngy = s->nGY();
  const int Pxy = s->getPxy();
  const int Pxyz = s->getPxyz();

  // Get Core ID and SIMT ID
  int cid = tid / n_threads_per_core;
  int smtid = tid - n_threads_per_core * cid;

  // Compute smt ID Y and Z indices
  int smtid_z = smtid / Sy;
  int smtid_y = smtid - Sy * smtid_z;

  // unsigned int accumulate[8] = { -1, -1, -1, -1, -1, -1, -1, -1 };
  unsigned int accumulate[8] = {~0U, ~0U, ~0U, ~0U, ~0U, ~0U, ~0U, ~0U};

  int nvecs = s->nVecs();

  const int gauge_line_in_floats =
      sizeof(SU3MatrixBlock) / sizeof(FT); // One gauge soavector
  const int spinor_line_in_floats =
      sizeof(FourSpinorBlock) / sizeof(FT); //  One spinor soavecto

  // Indexing constants
  const int V1 = 2 * nvecs; // No of vectors in x (without checkerboarding)
  const int NyV1 = Ny * V1;
  const int NzNyV1 = Nz * Ny * V1;

  const int Nxm1 = 2 * Nxh - 1;
  const int Nym1 = Ny - 1;
  const int Nzm1 = Nz - 1;
  const int Ntm1 = Nt - 1;

  const int NyV1mV1 = V1 * (Ny - 1);
  const int NzNyV1mNyV1 = V1 * Ny * (Nz - 1);
  const int NtNzNyV1mNzNyV1 = V1 * Nz * Ny * (Nt - 1);

  const int nyg = s->nGY();
  // Get the number of checkerboarded sites and various indexing constants
  int gprefdist = 0;
  int soprefdist = 0;

#if defined(__GNUG__) && !defined(__INTEL_COMPILER)
  int *tmpspc __attribute__((aligned(64))) = &(tmpspc_all[veclen * 16 * tid]);
#else
  __declspec(align(64)) int *tmpspc = &(tmpspc_all[veclen * 16 * tid]);
#endif

  int *offs, *xbOffs, *xfOffs, *ybOffs, *yfOffs, *gOffs, *pfyOffs;
  int *xbOffs_xodd[2], *xbOffs_x0_xodd[2];
  int *xfOffs_xodd[2], *xfOffs_xn_xodd[2];
  int *ybOffs_yn0, *ybOffs_y0, *yfOffs_ynn, *yfOffs_yn;
  int *atmp = (int *)((((unsigned long long)tmpspc) + 0x3F) & ~0x3F);
  offs = &atmp[0];
  xbOffs_xodd[0] = &atmp[veclen * 1];
  xbOffs_xodd[1] = &atmp[veclen * 2];
  xbOffs_x0_xodd[0] = &atmp[veclen * 3];
  xbOffs_x0_xodd[1] = &atmp[veclen * 4];
  xfOffs_xodd[0] = &atmp[veclen * 5];
  xfOffs_xodd[1] = &atmp[veclen * 6];
  xfOffs_xn_xodd[0] = &atmp[veclen * 7];
  xfOffs_xn_xodd[1] = &atmp[veclen * 8];
  ybOffs_yn0 = &atmp[veclen * 9];
  ybOffs_y0 = &atmp[veclen * 10];
  yfOffs_ynn = &atmp[veclen * 11];
  yfOffs_yn = &atmp[veclen * 12];
  gOffs = &atmp[veclen * 13];
  pfyOffs = &atmp[veclen * 14];

  int num_phases = s->getNumPhases();

  for (int ph = 0; ph < num_phases; ph++) {
    const CorePhase &phase = s->getCorePhase(ph);
    const BlockPhase &binfo = block_info[tid * num_phases + ph];

    int nActiveCores = phase.Cyz * phase.Ct;
    if (cid >= nActiveCores)
      continue;

    int ph_next = ph;
    int Nct = binfo.nt;

    // Loop over timeslices
    for (int ct = 0; ct < Nct; ct++) {
      int t = ct + binfo.bt;
      int ct_next = ct;

      double forw_t_coeff = aniso_coeff_T;
      double back_t_coeff = aniso_coeff_T;

      accumulate[6] = -1;
      accumulate[7] = -1;

      // If we are on timeslice 0, we should set back t_coeff and accumulate
      if (t == 0) {
        if (!comms->localT()) {
          accumulate[6] = 0;
        } else {
          if (amIPtMin) {
            // This can be in half precision now

            back_t_coeff *= t_boundary;
          }
        }
      }

      // If we are on timeslice 1, we should set forw t_coeff and accumulate flags
      if (t == Nt - 1) {
        if (!comms->localT()) {
          accumulate[7] = 0;
        } else {
          if (amIPtMax) {
            forw_t_coeff *= t_boundary;
          }
        }
      }

      // Loop over z. Start at smtid_z and work up to Ncz
      // (Ncz truncated for the last block so should be OK)
      for (int cz = smtid_z; cz < Bz; cz += Sz) {

        int z = cz + binfo.bz; // Add on origin of block
        int cz_next = cz;
        if (!comms->localZ()) {
          if (z == 0) {
            accumulate[4] = 0;
          } else {
            accumulate[4] = -1;
          }

          if (z == Nz - 1) {
            accumulate[5] = 0;
          } else {
            accumulate[5] = -1;
          }
        }

        const FourSpinorBlock *xyBase =
            &psi[t * Pxyz + z * Pxy]; // base address for x & y neighbours
        const FourSpinorBlock *zbBase =
            &psi[t * Pxyz] +
            (z == 0 ? (Nz - 1) * Pxy
                    : (z - 1) * Pxy); // base address for prev z neighbour
        const FourSpinorBlock *zfBase =
            &psi[t * Pxyz] +
            (z == Nz - 1 ? 0 : (z + 1) * Pxy); // base address for next z neighbour
        const FourSpinorBlock *tbBase =
            &psi[z * Pxy] +
            (t == 0 ? (Nt - 1) * Pxyz
                    : (t - 1) * Pxyz); // base address for prev t neighbour
        const FourSpinorBlock *tfBase =
            &psi[z * Pxy] +
            (t == Nt - 1 ? 0 : (t + 1) * Pxyz); // base address for next t neighbour
        FourSpinorBlock *oBase = &res[t * Pxyz + z * Pxy];

        // Loop over y. Start at smtid_y and work up to Ncy
        // (Ncy truncated for the last block so should be OK)
        for (int cy = nyg * smtid_y; cy < By; cy += nyg * Sy) {
          int yi = cy + binfo.by;
          int cy_next = cy;
          const int xodd = (yi + z + t + cb) & 1;

          // cx loops over the soalen partial vectors
          for (int cx = 0; cx < nvecs; cx++) {
            int ind = 0;
            int cx_next = cx + 1;

            if (cx_next == nvecs) {
              cx_next = 0;
              cy_next += nyg * Sy;
              if (cy_next >= By) {
                cy_next = nyg * smtid_y;
                cz_next += Sz;
                if (cz_next >= Bz) {
                  cz_next = smtid_z;
                  ct_next++;
                  if (ct_next == Nct) {
                    ct_next = 0;
                    ph_next++;
                    if (ph_next == num_phases) {
                      ph_next = 0;
                    }
                  }
                }
              }
            }

            const BlockPhase &binfo_next = block_info[tid * num_phases + ph_next];
            int yi_next = cy_next + binfo_next.by;
            int z_next = cz_next + binfo_next.bz;
            int t_next = ct_next + binfo_next.bt;

            int off_next = (t_next - t) * Pxyz + (z_next - z) * Pxy +
                           (yi_next - yi) * nvecs + (cx_next - cx);
            int si_off_next = off_next * spinor_line_in_floats;

            const SU3MatrixBlock *gBase =
                &u[(t * Pxyz + z * Pxy + yi * nvecs) / nyg + cx];
            int g_off_next = (((t_next - t) * Pxyz + (z_next - z) * Pxy +
                               (yi_next - yi) * nvecs) /
                                  nyg +
                              (cx_next - cx)) *
                             gauge_line_in_floats;

            int X = nvecs * yi + cx;
            xbOffs = (cx == 0 ? xbOffs_x0_xodd[xodd] : xbOffs_xodd[xodd]);

#if 1
            accumulate[0] = (cx == 0 ? xbmask_x0_xodd[xodd] : -1);
#endif

            xfOffs = (cx == nvecs - 1 ? xfOffs_xn_xodd[xodd] : xfOffs_xodd[xodd]);

#if 1
            accumulate[1] = (cx == nvecs - 1 ? xfmask_xn_xodd[xodd] : -1);
#endif

            ybOffs = (yi == 0 ? ybOffs_y0 : ybOffs_yn0);
#if 1
            accumulate[2] = (yi == 0 ? ybmask_y0 : -1);
#endif
            yfOffs = (yi == Ny - nyg ? yfOffs_yn : yfOffs_ynn);
#if 1
            accumulate[3] = (yi == Ny - nyg ? yfmask_yn : -1);
#endif

// using Cilk array notations:
#ifdef QPHIX_USE_CEAN
            pfyOffs [0:(veclen + 1) / 2] = ybOffs [0:(veclen + 1) / 2];
            pfyOffs [(veclen + 1) / 2:veclen / 2] =
            yfOffs [(veclen + 1) / 2:veclen / 2];
#else
            for (int it = 0; it < (veclen + 1) / 2; it++) {
              pfyOffs[it] = ybOffs[it];
              pfyOffs[it + (veclen + 1) / 2] = yfOffs[it + (veclen + 1) / 2];
            }
#endif

            if (soalen == veclen) {
              if (!comms->localY()) {
                accumulate[2] = (yi == 0 ? 0 : -1);
                accumulate[3] = (yi == Ny - 1 ? 0 : -1);
              }
            }

            FT aniso_coeff_s_T = rep<FT, double>(aniso_coeff_S);
            FT forw_t_coeff_T = rep<FT, double>(forw_t_coeff);
            FT back_t_coeff_T = rep<FT, double>(back_t_coeff);

            kernel(xyBase + X,
                   zbBase + X,
                   zfBase + X,
                   tbBase + X,
                   tfBase + X,
                   oBase + X,
                   gBase,
                   xbOffs,
                   xfOffs,
                   ybOffs,
                   yfOffs,
                   offs,
                   gOffs,
                   si_off_next,
                   si_off_next,
                   si_off_next,
                   si_off_next,
                   g_off_next,
                   pfyOffs,
                   zbBase + X,
                   zfBase + X,
                   tfBase + X,
                   accumulate,
                   aniso_coeff_s_T,
                   forw_t_coeff_T,
                   back_t_coeff_T,
                   tbc_phases,
                   derived_mu,
                   derived_mu_inv);
          }
        } // End for over scanlines y
      } // End for over scalines z

      if (ct % BARRIER_TSLICES == 0)
        barriers[ph][binfo.cid_t]->wait(binfo.group_tid);

    } // end for over t
  } // phases
}

template <typename FT, int veclen, int soalen, bool compress12>
void TMDslash<FT, veclen, soalen, compress12>::TMDyzAChiMinusBDPsi(
    int tid,
    const FourSpinorBlock *psi,
    const FourSpinorBlock *chi,
    FourSpinorBlock *res,
    const SU3MatrixBlock *u,
    double alpha,
    double beta,
    bool const is_plus,
    int cb,
    const double musign)
{
  auto kernel =
      (is_plus
           ? QPHIX_KERNEL_SELECT(
                 tm_dslash_achimbdpsi_plus_vec, FT, veclen, soalen, compress12)
           : QPHIX_KERNEL_SELECT(
                 tm_dslash_achimbdpsi_minus_vec, FT, veclen, soalen, compress12));

  const int Nxh = s->Nxh();
  const int Nx = s->Nx();
  const int Ny = s->Ny();
  const int Nz = s->Nz();
  const int Nt = s->Nt();
  const int By = s->getBy();
  const int Bz = s->getBz();
  const int Sy = s->getSy();
  const int Sz = s->getSz();
  const int ngy = s->nGY();
  const int Pxy = s->getPxy();
  const int Pxyz = s->getPxyz();

  double beta_s = beta * aniso_coeff_S;
  double beta_t = beta * aniso_coeff_T;

  // Get Core ID and SIMT ID
  int cid = tid / n_threads_per_core;
  int smtid = tid - n_threads_per_core * cid;

  // Compute smt ID Y and Z indices
  int smtid_z = smtid / Sy;
  int smtid_y = smtid - Sy * smtid_z;

  // unsigned int accumulate[8] = { -1, -1, -1, -1, -1, -1, -1, -1 };
  unsigned int accumulate[8] = {~0U, ~0U, ~0U, ~0U, ~0U, ~0U, ~0U, ~0U};
  int nvecs = s->nVecs();

  const int gauge_line_in_floats =
      sizeof(SU3MatrixBlock) / sizeof(FT); // One gauge soavector
  const int spinor_line_in_floats =
      sizeof(FourSpinorBlock) / sizeof(FT); //  One spinor soavecto

  // Indexing constants
  const int V1 = 2 * nvecs; // No of vectors in x (without checkerboarding)
  const int NyV1 = Ny * V1;
  const int NzNyV1 = Nz * Ny * V1;

  const int Nxm1 = 2 * Nxh - 1;
  const int Nym1 = Ny - 1;
  const int Nzm1 = Nz - 1;
  const int Ntm1 = Nt - 1;

  const int NyV1mV1 = V1 * (Ny - 1);
  const int NzNyV1mNyV1 = V1 * Ny * (Nz - 1);
  const int NtNzNyV1mNzNyV1 = V1 * Nz * Ny * (Nt - 1);

  const int nyg = s->nGY();
  // Get the number of checkerboarded sites and various indexing constants
  int gprefdist = 0;
  int soprefdist = 0;

#if defined(__GNUG__) && !defined(__INTEL_COMPILER)
  int *tmpspc __attribute__((aligned(64))) = &(tmpspc_all[veclen * 16 * tid]);
#else
  __declspec(align(64)) int *tmpspc = &(tmpspc_all[veclen * 16 * tid]);
#endif

  int *offs, *xbOffs, *xfOffs, *ybOffs, *yfOffs, *gOffs, *pfyOffs;
  int *xbOffs_xodd[2], *xbOffs_x0_xodd[2];
  int *xfOffs_xodd[2], *xfOffs_xn_xodd[2];
  int *ybOffs_yn0, *ybOffs_y0, *yfOffs_ynn, *yfOffs_yn;
  int *atmp = (int *)((((unsigned long long)tmpspc) + 0x3F) & ~0x3F);
  offs = &atmp[0];
  xbOffs_xodd[0] = &atmp[veclen * 1];
  xbOffs_xodd[1] = &atmp[veclen * 2];
  xbOffs_x0_xodd[0] = &atmp[veclen * 3];
  xbOffs_x0_xodd[1] = &atmp[veclen * 4];
  xfOffs_xodd[0] = &atmp[veclen * 5];
  xfOffs_xodd[1] = &atmp[veclen * 6];
  xfOffs_xn_xodd[0] = &atmp[veclen * 7];
  xfOffs_xn_xodd[1] = &atmp[veclen * 8];
  ybOffs_yn0 = &atmp[veclen * 9];
  ybOffs_y0 = &atmp[veclen * 10];
  yfOffs_ynn = &atmp[veclen * 11];
  yfOffs_yn = &atmp[veclen * 12];
  gOffs = &atmp[veclen * 13];
  pfyOffs = &atmp[veclen * 14];

  int num_phases = s->getNumPhases();

  for (int ph = 0; ph < num_phases; ph++) {
    const CorePhase &phase = s->getCorePhase(ph);
    const BlockPhase &binfo = block_info[tid * num_phases + ph];

    int nActiveCores = phase.Cyz * phase.Ct;
    if (cid >= nActiveCores)
      continue;

    int ph_next = ph;
    int Nct = binfo.nt;

    // Loop over timeslices
    for (int ct = 0; ct < Nct; ct++) {
      int t = ct + binfo.bt;
      double beta_t_f = beta_t;
      double beta_t_b = beta_t;
      accumulate[6] = -1;
      accumulate[7] = -1;

      if (t == 0) {
        if (!comms->localT()) {
          accumulate[6] = 0;
        } else {
          if (amIPtMin) {
            beta_t_b *= t_boundary;
          }
        }
      }

      if (t == Nt - 1) {
        if (!comms->localT()) {
          accumulate[7] = 0;
        } else {
          if (amIPtMax) {
            beta_t_f *= t_boundary;
          }
        }
      }

      int ct_next = ct;

      // Loop over z. Start at smtid_z and work up to Ncz
      // (Ncz truncated for the last block so should be OK)
      for (int cz = smtid_z; cz < Bz; cz += Sz) {

        int z = cz + binfo.bz; // Add on origin of block
        int cz_next = cz;
        if (!comms->localZ()) {
          if (z == 0) {
            accumulate[4] = 0;
          } else {
            accumulate[4] = -1;
          }

          if (z == Nz - 1) {
            accumulate[5] = 0;
          } else {
            accumulate[5] = -1;
          }
        }

        const FourSpinorBlock *xyBase =
            &psi[t * Pxyz + z * Pxy]; // base address for x & y neighbours
        const FourSpinorBlock *zbBase =
            &psi[t * Pxyz] +
            (z == 0 ? (Nz - 1) * Pxy
                    : (z - 1) * Pxy); // base address for prev z neighbour
        const FourSpinorBlock *zfBase =
            &psi[t * Pxyz] +
            (z == Nz - 1 ? 0 : (z + 1) * Pxy); // base address for next z neighbour
        const FourSpinorBlock *tbBase =
            &psi[z * Pxy] +
            (t == 0 ? (Nt - 1) * Pxyz
                    : (t - 1) * Pxyz); // base address for prev t neighbour
        const FourSpinorBlock *tfBase =
            &psi[z * Pxy] +
            (t == Nt - 1 ? 0 : (t + 1) * Pxyz); // base address for next t neighbour
        const FourSpinorBlock *chiBase = &chi[t * Pxyz + z * Pxy];
        FourSpinorBlock *oBase = &res[t * Pxyz + z * Pxy];

        // Loop over y. Start at smtid_y and work up to Ncy
        // (Ncy truncated for the last block so should be OK)
        for (int cy = nyg * smtid_y; cy < By; cy += nyg * Sy) {
          int yi = cy + binfo.by;
          int cy_next = cy;
          const int xodd = (yi + z + t + cb) & 1;

          // cx loops over the soalen partial vectors
          for (int cx = 0; cx < nvecs; cx++) {
            int ind = 0;
            int cx_next = cx + 1;

            if (cx_next == nvecs) {
              cx_next = 0;
              cy_next += nyg * Sy;
              if (cy_next >= By) {
                cy_next = nyg * smtid_y;
                cz_next += Sz;
                if (cz_next >= Bz) {
                  cz_next = smtid_z;
                  ct_next++;
                  if (ct_next == Nct) {
                    ct_next = 0;
                    ph_next++;
                    if (ph_next == num_phases) {
                      ph_next = 0;
                    }
                  }
                }
              }
            }

            const BlockPhase &binfo_next = block_info[tid * num_phases + ph_next];
            int yi_next = cy_next + binfo_next.by;
            int z_next = cz_next + binfo_next.bz;
            int t_next = ct_next + binfo_next.bt;

            int off_next = (t_next - t) * Pxyz + (z_next - z) * Pxy +
                           (yi_next - yi) * nvecs + (cx_next - cx);
            int si_off_next = off_next * spinor_line_in_floats;

            const SU3MatrixBlock *gBase =
                &u[(t * Pxyz + z * Pxy + yi * nvecs) / nyg + cx];
            int g_off_next = (((t_next - t) * Pxyz + (z_next - z) * Pxy +
                               (yi_next - yi) * nvecs) /
                                  nyg +
                              (cx_next - cx)) *
                             gauge_line_in_floats;

            int X = nvecs * yi + cx;
            xbOffs = (cx == 0 ? xbOffs_x0_xodd[xodd] : xbOffs_xodd[xodd]);
#if 1
            accumulate[0] = (cx == 0 ? xbmask_x0_xodd[xodd] : -1);
#endif
            xfOffs = (cx == nvecs - 1 ? xfOffs_xn_xodd[xodd] : xfOffs_xodd[xodd]);

#if 1
            accumulate[1] = (cx == nvecs - 1 ? xfmask_xn_xodd[xodd] : -1);
#endif

            ybOffs = (yi == 0 ? ybOffs_y0 : ybOffs_yn0);

#if 1
            accumulate[2] = (yi == 0 ? ybmask_y0 : -1);

#endif
            yfOffs = (yi == Ny - nyg ? yfOffs_yn : yfOffs_ynn);

#if 1
            accumulate[3] = (yi == Ny - nyg ? yfmask_yn : -1);
#endif

#ifdef QPHIX_USE_CEAN
            pfyOffs [0:(veclen + 1) / 2] = ybOffs [0:(veclen + 1) / 2];
            pfyOffs [(veclen + 1) / 2:veclen / 2] =
            yfOffs [(veclen + 1) / 2:veclen / 2];

#else
            for (int it = 0; it < (veclen + 1) / 2; it++) {
              pfyOffs[it] = ybOffs[it];
              pfyOffs[it + (veclen + 1) / 2] = yfOffs[it + (veclen + 1) / 2];
            }
#endif

            if (soalen == veclen) {
              if (!comms->localY()) {
                accumulate[2] = (yi == 0 ? 0 : -1);
                accumulate[3] = (yi == Ny - 1 ? 0 : -1);
              }
            }

            FT alpha_T = rep<FT, double>(alpha);
            FT beta_s_T = rep<FT, double>(beta_s);
            FT beta_t_f_T = rep<FT, double>(beta_t_f);
            FT beta_t_b_T = rep<FT, double>(beta_t_b);

            kernel(xyBase + X,
                   zbBase + X,
                   zfBase + X,
                   tbBase + X,
                   tfBase + X,
                   chiBase + X,
                   oBase + X,
                   gBase,
                   xbOffs,
                   xfOffs,
                   ybOffs,
                   yfOffs,
                   offs,
                   gOffs,
                   si_off_next,
                   si_off_next,
                   si_off_next,
                   si_off_next,
                   si_off_next,
                   g_off_next,
                   pfyOffs,
                   zbBase + X,
                   zfBase + X,
                   tfBase + X,
                   chiBase + X,
                   alpha_T,
                   beta_s_T,
                   beta_t_f_T,
                   beta_t_b_T,
                   tbc_phases,
                   musign*derived_mu,
                   accumulate);
          }
        } // End for over scanlines y
      } // End for over scalines z

      if (ct % BARRIER_TSLICES == 0)
        barriers[ph][binfo.cid_t]->wait(binfo.group_tid);

    } // end for over t
  } // phases
}

template <typename FT, int veclen, int soalen, bool compress12>
void TMDslash<FT, veclen, soalen, compress12>::TMDPsi(const SU3MatrixBlock *u,
                                                      const FourSpinorBlock *psi_in,
                                                      FourSpinorBlock *res_out,
                                                      bool const is_plus,
                                                      int cb)
{
  double beta_s = -aniso_coeff_S;
  double beta_t_f = -aniso_coeff_T;
  double beta_t_b = -aniso_coeff_T;

  // Antiperiodic BCs on back links
  if (amIPtMin) {
    beta_t_b *= t_boundary;
  }

  // Antiperiodic BCs on forw links
  if (amIPtMax) {
    beta_t_f *= t_boundary;
  }

  TSC_tick t_start, t_end;
#ifdef QPHIX_DO_COMMS
  // Pre-initiate all receives

  for (int d = 3; d >= 0; d--) {
    if (!comms->localDir(d)) {
      comms->startRecvFromDir(2 * d + 0);
      comms->startRecvFromDir(2 * d + 1);

#pragma omp parallel
      {
        int tid = omp_get_thread_num();

        packTMFaceDir(tid, psi_in, comms->sendToDir[2 * d + 1], cb, d, 1, is_plus);
        packTMFaceDir(tid, psi_in, comms->sendToDir[2 * d + 0], cb, d, 0, is_plus);
      }
      comms->startSendDir(2 * d + 1);
      comms->startSendDir(2 * d + 0);
    }
  }
#endif // QPHIX_DO_COMMS

// DO BODY DON"T ACCUMULATE BOUNDARY
#pragma omp parallel
  {
    int tid = omp_get_thread_num();
    TMDyz(tid, psi_in, res_out, u, is_plus, cb);
  }

#ifdef QPHIX_DO_COMMS
  for (int d = 3; d >= 0; d--) {
    if (!comms->localDir(d)) {
      comms->finishSendDir(2 * d + 1);
      comms->finishRecvFromDir(2 * d + 0);
      comms->finishSendDir(2 * d + 0);
      comms->finishRecvFromDir(2 * d + 1);

#pragma omp parallel
      {
        int tid = omp_get_thread_num();
        completeTMFaceDir(tid,
                          comms->recvFromDir[2 * d + 0],
                          res_out,
                          u,
                          (d == 3 ? beta_t_b : beta_s),
                          cb,
                          d,
                          0,
                          is_plus);
        completeTMFaceDir(tid,
                          comms->recvFromDir[2 * d + 1],
                          res_out,
                          u,
                          (d == 3 ? beta_t_f : beta_s),
                          cb,
                          d,
                          1,
                          is_plus);
      }
    }
  }
#endif // QPHIX_DO_COMMS
}

template <typename FT, int veclen, int soalen, bool compress12>
void TMDslash<FT, veclen, soalen, compress12>::TMDPsiAChiMinusBDPsi(
    const SU3MatrixBlock *u,
    const FourSpinorBlock *psi_in,
    const FourSpinorBlock *chi_in,
    FourSpinorBlock *res_out,
    double alpha,
    double beta,
    bool const is_plus,
    int cb,
    const double musign)
{
  double beta_s = beta * aniso_coeff_S;
  double beta_t_f = beta * aniso_coeff_T;
  double beta_t_b = beta * aniso_coeff_T;

  // Antiperiodic BCs on back links
  if (amIPtMin) {
    beta_t_b *= t_boundary;
  }

  // Antiperiodic BCs on forw links
  if (amIPtMax) {
    beta_t_f *= t_boundary;
  }

#ifdef QPHIX_DO_COMMS
  // Pre-initiate all receives

  for (int d = 3; d >= 0; d--) {
    if (!comms->localDir(d)) {
      comms->startRecvFromDir(2 * d + 0);
      comms->startRecvFromDir(2 * d + 1);

#pragma omp parallel
      {
        int tid = omp_get_thread_num();

        packTMFaceDir(tid, psi_in, comms->sendToDir[2 * d + 1], cb, d, 1, is_plus);
        packTMFaceDir(tid, psi_in, comms->sendToDir[2 * d + 0], cb, d, 0, is_plus);
      }
      comms->startSendDir(2 * d + 1);
      comms->startSendDir(2 * d + 0);
    }
  }
#endif // QPHIX_DO_COMMS

#pragma omp parallel
  {
    int tid = omp_get_thread_num();
    TMDyzAChiMinusBDPsi(tid, psi_in, chi_in, res_out, u, alpha, beta, is_plus, cb, musign);
  }

#ifdef QPHIX_DO_COMMS
  for (int d = 3; d >= 0; d--) {
    if (!comms->localDir(d)) {
      comms->finishSendDir(2 * d + 1);
      comms->finishRecvFromDir(2 * d + 0);
      comms->finishSendDir(2 * d + 0);
      comms->finishRecvFromDir(2 * d + 1);

#pragma omp parallel
      {
        int tid = omp_get_thread_num();

        completeFaceDirAChiMBDPsi(tid,
                                  comms->recvFromDir[2 * d + 0],
                                  res_out,
                                  u,
                                  (d == 3 ? beta_t_b : beta_s),
                                  cb,
                                  d,
                                  0,
                                  is_plus);
        completeFaceDirAChiMBDPsi(tid,
                                  comms->recvFromDir[2 * d + 1],
                                  res_out,
                                  u,
                                  (d == 3 ? beta_t_f : beta_s),
                                  cb,
                                  d,
                                  1,
                                  is_plus);
      }
    } // end if
  } // end for
#endif // QPHIX_DO_COMMS
}

template <typename FT, int veclen, int soalen, bool compress12>
void TMDslash<FT, veclen, soalen, compress12>::two_flav_AChiMinusBDPsi(
    FourSpinorBlock * const res[2],
    const FourSpinorBlock *const psi[2],
    const FourSpinorBlock *const chi[2],
    const SU3MatrixBlock *u,
    double alpha,
    double beta,
    double epsilon,
    int isign,
    int cb)
{
  const int n_blas_simt = s->getNSIMT();

  // Iterate over the two result flavors …
  for (int f : {0, 1}) {
    // The following does:
    //
    //   \alpha*{(1 + i(\mu/\alpha)\tau3\gamma5)} \chi - \beta*Doe \psi
    // = {\alpha + i\mu\tau3\gamma5} \chi - \beta*Doe\psi
    //
    // where \psi at this stage is expected to be
    //
    //        (\alpha - i\mu\tau3\gamma5 + \eps\tau1) Deo \chi
    // \psi = ------------------------------------------------
    //                ( \alpha^2 + \mu^2 - \eps^2 )
    //
    // such that all signs come out correct and \beta should be 1/4
    // we use "musign" to implement tau3
    double musign = 1.0-2.0*f;
    dslashAChiMinusBDPsi(res[f], psi[f], chi[f], u, alpha, beta, isign, cb, musign);

    // now 
    //   \res += -\eps\tau1 \chi
    // (note the minus sign) to get the full \tilde{M}_{cb,cb}
    axpy(-epsilon, chi[1 - f], res[f], *s, n_blas_simt);
  }
}

} // Namespace
