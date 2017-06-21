#pragma once

#include "cli_args.h"
#include "unittest.h"

class TestClovFile : public TestFixture
{
 public:
  TestClovFile(CliArgs &args) : args_(args) {}

  void run(void);

 private:
  template <typename FT,
            int V,
            int S,
            bool compress,
            typename FT2,
            int V2,
            int SOA2,
            typename U,
            typename Phi>
  void runTest(double mass, double clov_coeff, const std::string &filename);

  CliArgs &args_;
};
