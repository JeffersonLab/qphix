/**
  \file Tests for the \f$ \alpha + \mathrm i \mu \gamma_5$ multiplication BLAS
  routine.
  */


#include <gtest/gtest.h>
#include <qdp.h>
#include <qphix/blas_new_c.h>

#include "../tests/compare_qdp_spinors_gtest.h"

#include <cmath>

using namespace QPhiX;

class BlasTest : public ::testing::Test
{
 public:
  typedef double FT;

  BlasTest()
      : latt_size_{16, 16, 16, 16}, geom_(latt_size_, 4, 4, 1, 1, 1, 1, 0, 1, false),
        source_(geom_), source_dn_(geom_), step1_(geom_), step1_dn_(geom_),
        step2_(geom_), step2_dn_(geom_), two_flav_source_{source_[0], source_dn_[0]},
        two_flav_step1_{step1_[0], step1_dn_[0]},
        two_flav_step2_{step2_[0], step2_dn_[0]}
  {
    source_.qdp()
        .elem(0) // Lattice
        .elem(0) // Spin
        .elem(0) // Color
        .real() = 1.0;

    //QDP::gaussian(source_.qdp());
    //QDP::gaussian(source_dn_.qdp());

    source_.pack();
    source_dn_.pack();
  }

  static int constexpr veclen = 8;
  static int constexpr soalen = 4;
  static bool constexpr compress12 = true;

  int latt_size_[4] = {16, 16, 16, 16};

  Geometry<FT, veclen, soalen, compress12> geom_;
  HybridSpinor<FT, veclen, soalen, compress12, QDP::LatticeFermionD> source_;
  HybridSpinor<FT, veclen, soalen, compress12, QDP::LatticeFermionD> source_dn_;
  HybridSpinor<FT, veclen, soalen, compress12, QDP::LatticeFermionD> step1_;
  HybridSpinor<FT, veclen, soalen, compress12, QDP::LatticeFermionD> step1_dn_;
  HybridSpinor<FT, veclen, soalen, compress12, QDP::LatticeFermionD> step2_;
  HybridSpinor<FT, veclen, soalen, compress12, QDP::LatticeFermionD> step2_dn_;

  typename Geometry<FT, veclen, soalen, compress12>::FourSpinorBlock *two_flav_source_[2];
  typename Geometry<FT, veclen, soalen, compress12>::FourSpinorBlock *two_flav_step1_[2];
  typename Geometry<FT, veclen, soalen, compress12>::FourSpinorBlock *two_flav_step2_[2];
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
TEST_F(BlasTest, apimugamma5Inverse)
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

TEST_F(BlasTest, twoFlavTwistedMassSignChange) {
    double const alpha = 2.3;
    double const mu = 4.2;
    double const epsilon = 0.37;

    double const apimu[2] = {alpha, mu};
    double const apimu_minus[2] = {-alpha, -mu};

    two_flav_twisted_mass(
        apimu, epsilon, two_flav_source_, two_flav_step1_, geom_, 1);
    two_flav_twisted_mass(
        apimu_minus, -epsilon, two_flav_source_, two_flav_step2_, geom_, 1);

    axpy<FT, veclen, soalen, compress12, 2>(
        1.0, two_flav_step1_, two_flav_step2_, geom_, 1);

    step2_.unpack();
    step2_dn_.unpack();
    expect_near(step1_, step2_, 1e-10, geom_, 0, "α+iμγ_5 plus and minus (up)");
    expect_near(step2_dn_, step2_dn_, 1e-10, geom_, 0, "α+iμγ_5 plus and minus (down)");
}

TEST_F(BlasTest, twoFlavTwistedMassInverse) {
    double const alpha = 3.0;
    double const mu = 0.0;
    double const epsilon = 5.0;
    double const denominator = alpha * alpha + mu * mu - epsilon * epsilon;

    ASSERT_GT(std::fabs(denominator), 0.0);

    double const apimu[2] = {alpha, mu};
    double const apimu_inv[2] = {alpha / denominator, -mu / denominator};
    double const epsilon_inv = -epsilon / denominator;

    two_flav_twisted_mass(
        apimu, epsilon, two_flav_source_, two_flav_step1_, geom_, 1);
    two_flav_twisted_mass(
        apimu_inv, epsilon_inv, two_flav_step1_, two_flav_step2_, geom_, 1);

    step1_.unpack();
    step1_dn_.unpack();
    step2_.unpack();
    step2_dn_.unpack();

    QDPIO::cout << "Source:\n  " << source_.qdp().elem(0).elem(0).elem(0).real()
                << " + i " << source_.qdp().elem(0).elem(0).elem(0).imag() << "\n  "
                << source_dn_.qdp().elem(0).elem(0).elem(0).real() << "  "
                << source_dn_.qdp().elem(0).elem(0).elem(0).imag() << "\n"
                << "Step 1:\n  " << step1_.qdp().elem(0).elem(0).elem(0).real() << " + i "
                << step1_.qdp().elem(0).elem(0).elem(0).imag() << "\n  "
                << step1_dn_.qdp().elem(0).elem(0).elem(0).real() << " + i "
                << step1_dn_.qdp().elem(0).elem(0).elem(0).imag() << "\n"
                << "Step 2:\n  " << step2_.qdp().elem(0).elem(0).elem(0).real() << " + i "
                << step2_.qdp().elem(0).elem(0).elem(0).imag() << "\n  "
                << step2_dn_.qdp().elem(0).elem(0).elem(0).real() << " + i "
                << step2_dn_.qdp().elem(0).elem(0).elem(0).imag() << "\n";

    expect_near(source_, step2_, 1e-10, geom_, 0, "α+iμγ_5 normal and inverse (up)");
    expect_near(
        source_dn_, step2_dn_, 1e-10, geom_, 0, "α+iμγ_5 normal and inverse (down)");
}
