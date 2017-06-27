/**
  \file Tests for the \f$ \alpha + \mathrm i \mu \gamma_5$ multiplication BLAS
  routine.
  */

#include <gtest/gtest.h>
#include <qdp.h>
#include <qphix/blas_new_c.h>

#include "../tests/compare_qdp_spinors.h"

using namespace QPhiX;

/**
  Verifies that multiplying with the inverse returns the original spinor.

  We want to have
  \f[
    \psi = (\alpha + \mathrm i \mu \gamma_5)^{-1}
    (\alpha + \mathrm i \mu \gamma_5) \psi \,.
  \f]
  The inverse is given by \f$
    (\alpha + \mathrm i \mu \gamma_5)^{-1}
    = \frac{\alpha}{\alpha^2 + \mu^2}
    + \frac{- \mathrm i \mu \gamma_5}{\alpha^2 + \mu^2} \,.
  \f$
  */
TEST(apimugamma5, inverse) {
    int latt_size[] = {16, 16, 16, 16};
    typedef double FT;
    int constexpr veclen = 8;
    int constexpr soalen = 4;
    bool constexpr compress12 = true;

    Geometry<FT, veclen, soalen, compress12> geom(
        latt_size, 4, 4, 1, 1, 1, 1, 0, 1, true);

    HybridSpinor<FT, veclen, soalen, compress12, QDP::LatticeFermionD> source(geom);
    HybridSpinor<FT, veclen, soalen, compress12, QDP::LatticeFermionD> step1(geom);
    HybridSpinor<FT, veclen, soalen, compress12, QDP::LatticeFermionD> step2(geom);

    QDP::gaussian(source.qdp());
    source.pack();

    double const alpha = 3.4;
    double const mu = 1.2;
    double const denominator = alpha * alpha + mu * mu;
    double const apimu[] = {alpha, mu};
    double const apimu_inv[] = {alpha / denominator, -mu / denominator};

    twisted_mass(apimu, source[0], step1[0], geom, 1);
    twisted_mass(apimu_inv, step1[0], step2[0], geom, 1);

    step2.unpack();
    expect_near(source, step2, 1e-10, geom, 0, "α+iμγ_5 normal and inverse");
}
