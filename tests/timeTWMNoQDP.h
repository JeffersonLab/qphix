#pragma once

#include "cli_args.h"

class TimeTMDslash
{
 public:
  TimeTMDslash(CliArgs &args) : args_(args) {}

  void run();

  template <typename FT,
            int veclen,
            int soalen,
            bool compress12,
            typename QdpGauge,
            typename QdpSpinor>
  void operator()();

 private:
  template <typename FT, int V, int S, bool compress12>
  void runTest();

  CliArgs &args_;
};
