/**
  \file Tests for the \f$ \alpha + \mathrm i \mu \gamma_5$ multiplication BLAS
  routine.
  */

#include <gtest/gtest.h>
#include <qdp.h>
#include <qphix/blas_new_c.h>

#include "../tests/compare_qdp_spinors.h"

using namespace QPhiX;

class BlasTest : public ::testing::Test
{
 public:
  typedef double FT;

  BlasTest()
      : latt_size_{16, 16, 16, 16}, geom_(latt_size_, 4, 4, 1, 1, 1, 1, 0, 1, true),
        source_(geom_), step1_(geom_), step2_(geom_)
  {
    QDP::gaussian(source_.qdp());
    source_.pack();
  }

  static int constexpr veclen = 8;
  static int constexpr soalen = 4;
  static bool constexpr compress12 = true;

  int latt_size_[4] = {16, 16, 16, 16};

  Geometry<FT, veclen, soalen, compress12> geom_;
  HybridSpinor<FT, veclen, soalen, compress12, QDP::LatticeFermionD> source_;
  HybridSpinor<FT, veclen, soalen, compress12, QDP::LatticeFermionD> step1_;
  HybridSpinor<FT, veclen, soalen, compress12, QDP::LatticeFermionD> step2_;
};

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
TEST_F(BlasTest, apimugamma5)
{
  double const alpha = 3.4;
  double const mu = 1.2;
  double const denominator = alpha * alpha + mu * mu;
  double const apimu[] = {alpha, mu};
  double const apimu_inv[] = {alpha / denominator, -mu / denominator};

  twisted_mass(apimu, source_[0], step1_[0], geom_, 1);
  twisted_mass(apimu_inv, step1_[0], step2_[0], geom_, 1);

  step2_.unpack();
  expect_near(source_, step2_, 1e-10, geom_, 0, "α+iμγ_5 normal and inverse");
}

TEST(blas, cm)
{
}
