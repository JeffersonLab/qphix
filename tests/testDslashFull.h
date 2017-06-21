#pragma once

#include "unittest.h"

#include "prec.h"

class TestDslash : public TestFixture
{
 public:
  TestDslash(CliArgs &args) : args_(args) {}

  void run() override;

  template <typename FT,
            int veclen,
            int soalen,
            bool compress12,
            typename QdpGauge,
            typename QdpSpinor>
  void operator()();

 private:
  CliArgs &args_;
  Seed rng_seed;

  template <typename T, int V, int S, bool compress, typename U, typename Phi>
  void testDslash(const multi1d<U> &u, int t_bc);

  template <typename T, int V, int S, bool compress, typename U, typename Phi>
  void testDslashAChiMBDPsi(const multi1d<U> &u, int t_bc);

  template <typename T, int V, int S, bool compress, typename U, typename Phi>
  void testM(const multi1d<U> &u, int t_bc);

  template <typename T, int V, int S, bool compress, typename U, typename Phi>
  void testCG(const multi1d<U> &u, int t_bc);

  template <typename T, int V, int S, bool compress, typename U, typename Phi>
  void testBiCGStab(const multi1d<U> &u, int t_bc);

  template <typename T1,
            int VEC1,
            int SOA1,
            bool compress,
            typename T2,
            int VEC2,
            int SOA2,
            typename U,
            typename Phi>
  void testRichardson(const multi1d<U> &u, int t_bc);
};
