#pragma once

#include "qphix/qphix_cli_args.h"
#include "prec.h"

#include <string>

struct CliArgs {
  int nrow_in[4];
  int iters = 1;

  bool compress12;
  int qmp_geometry[4];

  Prec prec;
  int soalen;
  bool thread_bind;

  bool do_dslash;
  bool do_m;
  bool do_cg;
  bool do_bicgstab;

  QPhiX::QPhiXCLIArgs some_user_args;

  // Arguments contained in `QPhiXCLIArgs`.
  int By;
  int Bz;
  int PadXY;
  int PadXYZ;
  int NCores;
  int Sy;
  int Sz;
  int MinCt;
};

void printHelp(bool const is_timing);
void printArgHelp(bool const is_timing);
void printTimingCaseHelp();

CliArgs processArgs(int &argc, char **&argv, bool const is_timing = false);
