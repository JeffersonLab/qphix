#pragma once

#include "clover_term.h"
#include "reunit.h"
#include "tolerance.h"

#include <qphix/geometry.h>
#include <qphix/qdp_packer.h>

#include <qdp.h>

namespace QPhiX
{
template <typename FT,
          int veclen,
          int soalen,
          bool compress12,
          typename QdpGauge_ = QDP::LatticeColorMatrixD,
          typename QdpSpinor_ = QDP::LatticeDiracFermionD>
class RandomGauge
{
 public:
  typedef typename Geometry<FT, veclen, soalen, compress12>::FourSpinorBlock Spinor;
  typedef typename Geometry<FT, veclen, soalen, compress12>::SU3MatrixBlock Gauge;
  typedef typename Geometry<FT, veclen, soalen, compress12>::CloverBlock Clover;

  typedef QdpGauge_ QdpGauge;
  typedef QdpSpinor_ QdpSpinor;

  RandomGauge(Geometry<FT, veclen, soalen, compress12> &geom,
              double const t_boundary = 1.0,
              double const gauge_random_factor = 0.08);

  double const gauge_random_factor;
#if 0
  double const xi_0_f = 0.3;
  double const nu_f = 1.4;
#else
  double const xi_0_f = 1.0;
  double const nu_f = 1.0;
#endif
  double const clover_mass = 0.1;
  double const clover_coeff_R = 1.2;
  double const clover_coeff_T = 0.9;
  double const t_boundary;

  double const aniso_fac_s;
  double const aniso_fac_t;

  Gauge *u_packed[2];
  Clover *invclov_packed[2];
  Clover *clov_packed[2];

  ::QDP::multi1d<QdpGauge> u_aniso;

  ::QDP::multi1d<QdpGauge> const &get_u_aniso() const { return u_aniso; }

  CloverTermT<QdpSpinor, QdpGauge> clov_qdp, invclov_qdp;

 private:
  Geometry<FT, veclen, soalen, compress12> &geom;

  void init_random_gauge();
  void init_clover();

  ::QDP::multi1d<QdpGauge> u;

  GaugeHandle<FT, veclen, soalen, compress12> gauge_even, gauge_odd;
  CloverHandle<FT, veclen, soalen, compress12> A_even, A_odd, A_inv_even, A_inv_odd;

  AnisoParam_t aniso;
  CloverFermActParams clparam;
};

template <typename FT,
          int veclen,
          int soalen,
          bool compress12,
          typename QdpGauge,
          typename QdpSpinor>
RandomGauge<FT, veclen, soalen, compress12, QdpGauge, QdpSpinor>::RandomGauge(
    Geometry<FT, veclen, soalen, compress12> &geom,
    double const t_boundary,
    double const gauge_random_factor)
    : gauge_random_factor(gauge_random_factor), geom(geom), u(4), u_aniso(4),
      gauge_even(geom), gauge_odd(geom), A_even(geom), A_odd(geom), A_inv_even(geom),
      A_inv_odd(geom), aniso_fac_s(static_cast<double>(nu_f) / xi_0_f),
      aniso_fac_t(1.0), t_boundary(t_boundary)
{
  clov_packed[0] = A_even.get();
  clov_packed[1] = A_odd.get();
  invclov_packed[0] = A_inv_even.get();
  invclov_packed[1] = A_inv_odd.get();
  u_packed[0] = gauge_even.get();
  u_packed[1] = gauge_odd.get();

  init_random_gauge();
  init_clover();
}

template <typename FT,
          int veclen,
          int soalen,
          bool compress12,
          typename QdpGauge,
          typename QdpSpinor>
void RandomGauge<FT, veclen, soalen, compress12, QdpGauge, QdpSpinor>::
    init_random_gauge()
{
  aniso.anisoP = true;
  aniso.xi_0 = xi_0_f;
  aniso.nu = nu_f;
  aniso.t_dir = 3;

  QdpGauge g;
  QdpGauge uf;
  for (int mu = 0; mu < 4; mu++) {
    uf = 1; // Unit gauge
    gaussian(g);
    if (std::fabs(gauge_random_factor) > std::numeric_limits<double>::epsilon()) {
      u[mu] = uf + gauge_random_factor * g;
      reunit(u[mu]);
    }
  }
  qdp_pack_gauge<>(u, gauge_even.get(), gauge_odd.get(), geom);

  for (int mu = 0; mu < Nd; mu++) {
    Real factor;
    if (mu == 3) {
      factor = Real(aniso_fac_t);
    } else {
      factor = Real(aniso_fac_s);
    }
    u_aniso[mu] = factor * u[mu];
  }

  int const mu_t = QDP::Nd - 1;
  u_aniso[mu_t] *= QDP::where(QDP::Layout::latticeCoordinate(mu_t) ==
                                  (QDP::Layout::lattSize()[mu_t] - 1),
                              QDP::Real(t_boundary),
                              QDP::Real(1));
}

template <typename FT,
          int veclen,
          int soalen,
          bool compress12,
          typename QdpGauge,
          typename QdpSpinor>
void RandomGauge<FT, veclen, soalen, compress12, QdpGauge, QdpSpinor>::init_clover()
{
  clparam.anisoParam = aniso;
  clparam.Mass = clover_mass;
  clparam.clovCoeffR = clover_coeff_R;
  clparam.clovCoeffT = clover_coeff_T;

  clov_qdp.create(u, clparam);

  invclov_qdp = clov_qdp;


  // Invert the clover term.
  for (int cb = 0; cb < 2; cb++) {
    invclov_qdp.choles(cb);
  }

  for (int cb = 0; cb < 2; cb++) {
    qdp_pack_clover<>(clov_qdp, clov_packed[cb], geom, cb);
    qdp_pack_clover<>(invclov_qdp, invclov_packed[cb], geom, cb);
  }
}
}
