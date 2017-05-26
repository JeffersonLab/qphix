#pragma once

#include "qphix/qphix_cli_args.h"
#include "prec.h"

#include <string>

extern int nrow_in[4];
extern int iters;

extern bool compress12;
extern int qmp_geometry[4];

extern Prec prec_user;
extern bool thread_bind;

extern bool do_dslash;
extern bool do_m;
extern bool do_cg;
extern bool do_bicgstab;

extern QPhiX::QPhiXCLIArgs some_user_args;

// Arguments contained in `QPhiXCLIArgs`.
extern int By_user;
extern int Bz_user;
extern int PadXY_user;
extern int PadXYZ_user;
extern int NCores_user;
extern int Sy_user;
extern int Sz_user;
extern int MinCt_user;

void printHelp();

void processArgs(int &argc, char **&argv);
