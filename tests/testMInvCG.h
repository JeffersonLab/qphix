#pragma once

#include "unittest.h"
#include "cli_args.h"

class TestMultishift : public TestFixture
{
 public:
  TestMultishift(CliArgs &args) : args_(args) {}

  void run(void);

  template <typename FT,
            int veclen,
            int soalen,
            bool compress12,
            typename QdpGauge,
            typename QdpSpinor>
  void operator()();

 private:
  template <typename T, int V, int S, bool compress, typename U, typename Phi>
  void testMInvCG(int t_bc);

  CliArgs &args_;
  Seed rng_seed;
};
