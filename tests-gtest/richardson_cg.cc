/**
  \file Tests for the Richardson multiprecision solver using CG as inner and
  outer solver.
  */

#include <qphix/blas_new_c.h>
#include <qdp.h>
#include <qphix/qphix_cli_args.h>
#include <qphix/clover.h>
#include <qphix/invcg.h>
#include <qphix/inv_richardson_multiprec.h>

#include "../tests/compare_qdp_spinors_gtest.h"
#include "../tests/RandomGauge.h"
#include "../tests/veclen.h"

extern int nrow_in[];
extern QPhiX::QPhiXCLIArgs some_user_args;

using namespace QPhiX;

TEST(RichardsonCG, TMClover)
{
  bool constexpr compress12 = true;
  int const cb = 1;
  int const other_cb = 1 - cb;


  Geometry<double, VECLEN_DP, VECLEN_DP, compress12> geom_outer(
      Layout::subgridLattSize().slice(),
      some_user_args.getBy(),
      some_user_args.getBz(),
      some_user_args.getNCores(),
      some_user_args.getSy(),
      some_user_args.getSz(),
      some_user_args.getPxy(),
      some_user_args.getPxyz(),
      some_user_args.getMinCt());

  Geometry<float, VECLEN_SP, VECLEN_SP / 2, compress12> geom_inner(
      Layout::subgridLattSize().slice(),
      some_user_args.getBy(),
      some_user_args.getBz(),
      some_user_args.getNCores(),
      some_user_args.getSy(),
      some_user_args.getSz(),
      some_user_args.getPxy(),
      some_user_args.getPxyz(),
      some_user_args.getMinCt());

  RandomGauge<double,
              VECLEN_DP,
              VECLEN_DP,
              compress12,
              QDP::LatticeColorMatrixD,
              QDP::LatticeDiracFermionD>
      gauge_outer(geom_outer, 1.0);
  RandomGauge<float,
              VECLEN_SP,
              VECLEN_SP / 2,
              compress12,
              QDP::LatticeColorMatrixF,
              QDP::LatticeDiracFermionF>
      gauge_inner(geom_inner, gauge_outer);

  EvenOddCloverOperator<double, VECLEN_DP, VECLEN_DP, compress12> linop_outer(
      gauge_outer.u_packed,
      gauge_outer.clov_packed[cb],
      gauge_outer.invclov_packed[other_cb],
      &geom_outer,
      gauge_outer.t_boundary,
      gauge_outer.aniso_fac_s,
      gauge_outer.aniso_fac_t);
  EvenOddCloverOperator<float, VECLEN_SP, VECLEN_SP / 2, compress12> linop_inner(
      gauge_inner.u_packed,
      gauge_inner.clov_packed[cb],
      gauge_inner.invclov_packed[other_cb],
      &geom_inner,
      gauge_inner.t_boundary,
      gauge_inner.aniso_fac_s,
      gauge_inner.aniso_fac_t);

  int const max_iters = 10000;
  // The value is copied from https://github.com/JeffersonLab/qphix/issues/98
  double const delta = 1e-2;

  InvCG<float, VECLEN_SP, VECLEN_SP / 2, compress12> inv_cg(linop_inner, max_iters);
  InvRichardsonMultiPrec<double,
                         VECLEN_DP,
                         VECLEN_DP,
                         compress12,
                         float,
                         VECLEN_SP,
                         VECLEN_SP / 2,
                         compress12,
                         true>
      inv_richardson(linop_outer, inv_cg, delta, max_iters);

  HybridSpinor<double, VECLEN_DP, VECLEN_DP, compress12, QDP::LatticeDiracFermionD>
      hs_source(geom_outer), hs_qphix1(geom_outer), hs_qphix2(geom_outer);

  // Create a source vector that we are going to invert. You may switch to a
  // point source here easily.
  bool const point_source = false;
  if (point_source) {
    hs_source.zero();
    hs_source.qdp()
        .elem(0) // Lattice
        .elem(0) // Spin
        .elem(0) // Color
        .real() = 1.0;
    hs_source.qdp()
        .elem(0) // Lattice
        .elem(0) // Spin
        .elem(0) // Color
        .imag() = 0.0;
  } else {
    gaussian(hs_source.qdp());
  }
  hs_source.pack();

  double const rsd_target = 1e-7;
  int n_iters;
  double rsd_sq_final;
  unsigned long site_flops;
  unsigned long mv_apps;
  int const isign = 1;
  bool const verbose = true;

  // Invert once.
  inv_richardson(hs_qphix1[cb],
                 hs_source[cb],
                 rsd_target,
                 n_iters,
                 rsd_sq_final,
                 site_flops,
                 mv_apps,
                 isign,
                 verbose,
                 cb);

  // Apply the operator on the solution, should give the original spinor.
  linop_outer(hs_qphix2[cb], hs_qphix1[cb], isign, cb);

  hs_qphix2.unpack();
  expect_near(hs_source, hs_qphix2, rsd_target, geom_outer, cb);
}
