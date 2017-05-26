#include "cli_args.h"

#include <qphix/print_utils.h>

int nrow_in[4] = {4, 4, 4, 4};
int iters = 10000;
int By_user = -1;
int Bz_user = -1;
int NCores_user = -1;
int Sy_user = -1;
int Sz_user = -1;

// Hardwire these for now.
int PadXY_user = 1;
int PadXYZ_user = 0;
int MinCt_user = 1;

bool compress12 = false;
int qmp_geometry[4] = {1, 1, 1, 1};

Prec prec_user = FLOAT_PREC;
bool thread_bind = false;

bool do_dslash = false;
bool do_m = false;
bool do_cg = false;
bool do_bicgstab = false;

void printHelp()
{
  using QPhiX::masterPrintf;

  masterPrintf("This program measures the performance of the fundamental operations "
               "in QPhiX.\n"
               "\n"
               "The following options are available. If no default is listed, that "
               "option is a mandatory one.\n"
               "\n");
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
  masterPrintf("Internal data layout:\n"
               "  -by BY             block size in Y [defaults to Ly]\n"
               "  -bz BZ             block size in Z [defaults to Lz]\n"
               "  -pxy Pxy           extra pad in the XY plane [default: %i]\n"
               "  -pxyz Pxyz         extra pad in the XYZ plane [default: %i]\n"
               "  -minct MinCt       MinCt (default %i)\n"
               "\n",
               PadXY_user,
               PadXYZ_user,
               MinCt_user);
  masterPrintf("Threading options (all required):\n"
               "  -c NCores          number of cores\n"
               "  -sy Sy             number of SMT threads in Y\n"
               "  -sz Sz             number of SMT threads in Z\n"
               "\n");
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
               "");

  exit(1);
}

void processArgs(int argc, char *argv[])
{
  if (argc == 1) {
    printHelp();
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
    } else if (arg == "-by") {
      By_user = atoi(argv[i + 1]);
      i += 2;
    } else if (arg == "-bz") {
      Bz_user = atoi(argv[i + 1]);
      i += 2;
    } else if (arg == "-c") {
      NCores_user = atoi(argv[i + 1]);
      i += 2;
    } else if (arg == "-sy") {
      Sy_user = atoi(argv[i + 1]);
      i += 2;
    } else if (arg == "-sz") {
      Sz_user = atoi(argv[i + 1]);
      i += 2;
    } else if (arg == "-pxy") {
      PadXY_user = atoi(argv[i + 1]);
      i += 2;
    } else if (arg == "-pxyz") {
      PadXYZ_user = atoi(argv[i + 1]);
      i += 2;
    } else if (arg == "-minct") {
      MinCt_user = atoi(argv[i + 1]);
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
                          "a look at the supported options below.\n\n----\n\n",
                          arg.c_str());
      printHelp();
      exit(1);
    }
  }

  if (NCores_user < 0) {
    printHelp();
  }
  if (Sy_user < 0) {
    printHelp();
  }
  if (Sz_user < 0) {
    printHelp();
  }

  // Ct does not have to divide t, we can pick that up.
  if (By_user < 0) {
    By_user = nrow_in[1];
  }
  if (Bz_user < 0) {
    Bz_user = nrow_in[2];
  }

  if (!do_dslash && !do_m && !do_cg && !do_bicgstab) {
    QPhiX::masterPrintf("No timing case has been selected. Please select at least "
                        "one of them such that this program can do some actual "
                        "work.\n\n");
    printHelp();
  }
}
