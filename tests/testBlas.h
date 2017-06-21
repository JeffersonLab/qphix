#pragma once

#include "cli_args.h"

class TestBlas
{
 public:
  TestBlas(CliArgs &args) : args_(args) {}

  void run();

 private:
  CliArgs &args_;
};
