#pragma once

#include "cli_args.h"
#include "prec.h"
#include "tparam_selector.h"
#include "unittest.h"

#include <qphix/qphix_cli_args.h>

class TestClover : public TestFixture
{
 public:
  TestClover(CliArgs &args) : args_(args) {}

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
};

