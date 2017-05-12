#pragma once

#include <qphix/geometry.h>

#include <qdp.h>

namespace QPhiX
{
template <typename FT, int veclen, int soalen, bool compress12>
class RandomGauge
{
 public:
  typedef typename Geometry<FT, veclen, soalen, compress12>::FourSpinorBlock Spinor;
  typedef typename Geometry<FT, veclen, soalen, compress12>::SU3MatrixBlock Gauge;
  typedef typename Geometry<FT, veclen, soalen, compress12>::CloverBlock Clover;

  //static Real const gauge_random_factor(0.08);
  static double const gauge_random_factor(0.08);

  static double const xi_0_f = 0.3;
  static double const nu_f = 1.4;
  //static QDP::Real const clover_mass(0.1);
  //static QDP::Real const clover_coeff_R(1.2);
  //static QDP::Real const clover_coeff_T(0.9);
  static double const clover_mass(0.1);
  static double const clover_coeff_R(1.2);
  static double const clover_coeff_T(0.9);
  static double const t_boundary = 1.0;

  RandomGauge(Geometry<FT, veclen, soalen, compress12> &geom);

 private:
  void init_random_gauge();
  void init_clover();

  multi1d<U> u;
  multi1d<U> u_aniso;

  GaugeHandle<FT, veclen, soalen, compress12> gauge_even, gauge_odd;
  Gauge *u_packed[2];

  CloverHandle<FT, veclen, soalen, compress12> A_even, A_odd, A_inv_even, A_inv_odd;
  Clover *invclov_packed[2];
  Clover *clov_packed[2];

  AnisoParam_t aniso;
  CloverFermActParams clparam;
  double aniso_fac_s;
  double aniso_fac_t;
};

template <typename FT, int veclen, int soalen, bool compress12>
RandomGauge<FT, veclen, soalen, compress12>::RandomGauge(
    Geometry<FT, veclen, soalen, compress12> &geom)
    : u(4), u_aniso(4), gauge_even(geom), gauge_odd(geom), A_even(geom), A_odd(geom),
      A_inv_even(geom), A_inv_odd(geom),
      aniso_fac_s(static_cast<double>(nu_f) / xi_0_f), aniso_fac_t(1.0)
{
  clov_packed[0] = A_even.get();
  clov_packed[1] = A_odd.get();
  invclov_packed[0] = A_inv_even.get();
  invclov_packed[1] = A_inv_odd.get();
  u_packed[0] = packed_gauge_cb0;
  u_packed[1] = packed_gauge_cb1;

  init_random_gauge();
  init_clover();
}

template <typename FT, int veclen, int soalen, bool compress12>
RandomGauge<FT, veclen, soalen, compress12>::init_random_gauge() {
  aniso.anisoP = true;
  aniso.xi_0 = xi_0_f;
  aniso.nu = nu_f;
  aniso.t_dir = 3;

  U g;
  U uf;
  for (int mu = 0; mu < 4; mu++) {
    uf = 1; // Unit gauge
    gaussian(g);
    u[mu] = uf + gauge_random_factor * g;
    reunit(u[mu]);
  }
  qdp_pack_gauge<>(u, packed_gauge_cb0, packed_gauge_cb1, geom);

  for (int mu = 0; mu < Nd; mu++) {
    Real factor;
    if (mu == 3) {
      factor = Real(aniso_fac_t);
    } else {
      factor = Real(aniso_fac_s);
    }
    u_aniso[mu] = factor * u[mu];
  }
}

template <typename FT, int veclen, int soalen, bool compress12>
RandomGauge<FT, veclen, soalen, compress12>::init_clover() {
  clparam.anisoParam = aniso;
  clparam.Mass = clover_mass;
  clparam.clovCoeffR = clover_coeff_R;
  clparam.clovCoeffT = clover_coeff_T;

  CloverTermT<Phi, U> clov_qdp;
  clov_qdp.create(u, clparam);

  CloverTermT<Phi, U> invclov_qdp(clov_qdp);

  // Test spinors.
  QDP::LatticeDiracFermionD qdp_spinor_1;
  QDP::LatticeDiracFermionD qdp_spinor_2;
  QDP::LatticeDiracFermionD qdp_spinor_3;
  gaussian(qdp_spinor_1);

  // Apply the clover term and its copy (not inverted yet!) to the spinor as a sanity
  // check.
  for (int cb = 0; cb < 2; ++cb) {
    clov_qdp.apply(qdp_spinor_2, qdp_spinor_1, 1, cb);
    invclov_qdp.apply(qdp_spinor_3, qdp_spinor_1, 1, cb);
    expect_near(qdp_spinor_1, qdp_spinor_3, tolerance<FT>::small, geom, cb);
  }

  // Invert the clover term.
  for (int cb = 0; cb < 2; cb++) {
    invclov_qdp.choles(cb);
  }

  // Apply the clover term and inverse, result should be the original spinor.
  for (int cb = 0; cb < 2; ++cb) {
    clov_qdp.apply(qdp_spinor_2, qdp_spinor_1, 1, cb);
    invclov_qdp.apply(qdp_spinor_3, qdp_spinor_2, 1, cb);
    expect_near(qdp_spinor_1, qdp_spinor_3, tolerance<FT>::small, geom, cb);
  }

  for (int cb = 0; cb < 2; cb++) {
    qdp_pack_clover<>(clov_qdp, clov_packed[cb], geom, cb);
    qdp_pack_clover<>(invclov_qdp, invclov_packed[cb], geom, cb);
  }
}
}
