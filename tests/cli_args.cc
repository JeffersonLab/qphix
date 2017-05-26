#include "cli_args.h"

#include <qphix/print_utils.h>

using QPhiX::masterPrintf;

int nrow_in[4] = {4, 4, 4, 4};
int iters = 10000;

bool compress12 = false;
int qmp_geometry[4] = {1, 1, 1, 1};

Prec prec_user = FLOAT_PREC;
bool thread_bind = false;

bool do_dslash = false;
bool do_m = false;
bool do_cg = false;
bool do_bicgstab = false;

QPhiX::QPhiXCLIArgs some_user_args;

// Arguments contained in `QPhiXCLIArgs`.
int By_user = -1;
int Bz_user = -1;
int PadXY_user = 1;
int PadXYZ_user = 0;
int NCores_user = -1;
int Sy_user = -1;
int Sz_user = -1;
int MinCt_user = 1;

void printArgHelp() {
  masterPrintf("Lattice size:\n"
               "  -x Lx              lattice size in X [default: %i]\n"
               "  -y Ly              lattice size in Y [default: %i]\n"
               "  -z Lz              lattice size in Z [default: %i]\n"
               "  -t Lt              lattice size in T [default: %i]\n"
               "\n",
               nrow_in[0],
               nrow_in[1],
               nrow_in[2],
               nrow_in[3]);
               "\n",
  masterPrintf("MPI options (optional but recommended):\n"
               "  -geom Px Py Pz Pt  4D grid of MPI processes\n"
               "\n");
  masterPrintf("Timing parameters:\n"
               "  -i iters           number of iterations [default: %i]\n"
               "  -prec f|h|d        precision (float, half, double) [default: %s]\n"
               "  -compress12        enable gauge compression [default: %s]\n"
               "\n",
               iters,
               prec_user == FLOAT_PREC
                   ? "float"
                   : (prec_user == DOUBLE_PREC ? "double" : "half"),
               compress12 ? "given" : "not given");
  masterPrintf("Timing case selection (select at least one):\n"
               "  -dslash            Dslash\n"
               "  -mmat              Even-odd linear operator\n"
               "  -cg                Conjugate gradient\n"
               "  -bicgstab          BiCGStab\n"
               "\n");
}

void printHelp()
{
  masterPrintf("This program measures the performance of the fundamental operations "
               "in QPhiX.\n"
               "\n"
               "The following options are available. If no default is listed, that "
               "option is a mandatory one.\n"
               "\n");

  some_user_args.printArgHelp();
  printArgHelp();

  exit(1);
}

void processArgs(int &argc, char **&argv)
{
  if (argc == 1) {
    printHelp();
  }

  masterPrintf("argc = %i\n", argc);
  for (int i = 0; i < argc; ++i) {
      masterPrintf("argv[%02d] = %s\n", i, argv[i]);
  }

  // Let the class parse the arguments first. It will consume whatever it can parse.
  some_user_args.init(argc, argv);

  By_user = some_user_args.getBy();
  Bz_user = some_user_args.getBz();
  PadXY_user = some_user_args.getPxy();
  PadXYZ_user = some_user_args.getPxyz();
  NCores_user = some_user_args.getNCores();
  Sy_user = some_user_args.getSy();
  Sz_user = some_user_args.getSz();
  MinCt_user = some_user_args.getMinCt();

  masterPrintf("argc = %i\n", argc);
  for (int i = 0; i < argc; ++i) {
    masterPrintf("argv[%02d] = %s\n", i, argv[i]);
  }

  for (int i = 1; i < argc;) {
    std::string const arg(argv[i]);

    if (arg == "-x") {
      nrow_in[0] = atoi(argv[i + 1]);
      i += 2;
    } else if (arg == "-y") {
      nrow_in[1] = atoi(argv[i + 1]);
      i += 2;
    } else if (arg == "-z") {
      nrow_in[2] = atoi(argv[i + 1]);
      i += 2;
    } else if (arg == "-t") {
      nrow_in[3] = atoi(argv[i + 1]);
      i += 2;
    } else if (arg == "-i") {
      iters = atoi(argv[i + 1]);
      i += 2;
    } else if (arg == "-compress12") {
      compress12 = true;
      i++;
    } else if (arg == "-geom") {
      qmp_geometry[0] = atoi(argv[i + 1]);
      qmp_geometry[1] = atoi(argv[i + 2]);
      qmp_geometry[2] = atoi(argv[i + 3]);
      qmp_geometry[3] = atoi(argv[i + 4]);
      i += 4;

    } else if (arg == "-dslash") {
      do_dslash = true;
      i++;

    } else if (arg == "-mmat") {
      do_m = true;
      i++;
    } else if (arg == "-cg") {
      do_cg = true;
      i++;

    } else if (arg == "-bicgstab") {
      do_bicgstab = true;
      i++;
    }

    else if (arg == "-prec") {
      string user_arg(argv[i + 1]);
      if (user_arg.compare("f") == 0) {
        prec_user = FLOAT_PREC;
      }
      if (user_arg.compare("h") == 0) {
        prec_user = HALF_PREC;
      }

      if (user_arg.compare("d") == 0) {
        prec_user = DOUBLE_PREC;
      }
      i += 2;
    } else if (arg == "-bind") {
      thread_bind = true;
      i++;
    } else {
      QPhiX::masterPrintf("ERROR: The option “%s” could not be parsed. Please have "
                          "a look at the supported options below.\n\n",
                          arg.c_str());
      printArgHelp();
      exit(1);
    }
  }


  if (!do_dslash && !do_m && !do_cg && !do_bicgstab) {
    QPhiX::masterPrintf("No timing case has been selected. Please select at least "
                        "one of them such that this program can do some actual "
                        "work.\n\n");
    printArgHelp();
    exit(1);
  }
}
