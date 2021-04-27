#pragma once

#include "qphix/qphix_cli_args.h"
#include "prec.h"

#include <string>

struct CliArgs {
	CliArgs() {
		// Defaults
		for(int mu=0; mu < 4; ++mu) {
			nrow_in[mu]=8;
			qmp_geometry[mu]=1;
		}

		iters=50;
		compress12=false;
		prec=FLOAT_PREC;
		soalen=4; // Safe
		thread_bind=false;
		do_dslash=false;
		do_m = false;
		do_cg = false;
		do_bicgstab = false;
		// Leave some uninitialized

	}

  int nrow_in[4];
  int iters = 1;

  bool compress12 = false;
  int qmp_geometry[4];

  Prec prec;
  int soalen;
  bool thread_bind = false;

  bool do_dslash = false;
  bool do_m = false;
  bool do_cg = false;
  bool do_bicgstab = false;

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
