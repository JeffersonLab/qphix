#include "cli_args.h"

#include "veclen.h"

#include <qdp.h>
#include <qphix/print_utils.h>

using QPhiX::masterPrintf;

void printArgHelp(CliArgs &args, bool const is_timing)
{
  masterPrintf("Lattice size:\n"
               "  -x Lx              lattice size in X [default: %i]\n"
               "  -y Ly              lattice size in Y [default: %i]\n"
               "  -z Lz              lattice size in Z [default: %i]\n"
               "  -t Lt              lattice size in T [default: %i]\n"
               "\n",
               args.nrow_in[0],
               args.nrow_in[1],
               args.nrow_in[2],
               args.nrow_in[3]);
  masterPrintf("MPI options (optional but recommended):\n"
               "  -geom Px Py Pz Pt  4D grid of MPI processes\n"
               "\n");
  masterPrintf("Timing parameters:\n"
               "  -i iters           number of iterations [default: %i]\n"
               "  -prec f|h|d        precision (float, half, double) [default: %s]\n"
               "  -soalen soalen     structure of array length [default: %d]\n"
               "  -compress12        enable gauge compression [default: %s]\n"
               "\n",
               args.iters,
               args.prec == FLOAT_PREC
                   ? "float"
                   : (args.prec == DOUBLE_PREC ? "double" : "half"),
               args.soalen,
               args.compress12 ? "given" : "not given");

  if (is_timing) {
    printTimingCaseHelp();
  }
}

void printTimingCaseHelp()
{
  masterPrintf("Timing case selection (select at least one):\n"
               "  -dslash            Dslash\n"
               "  -mmat              Even-odd linear operator\n"
               "  -cg                Conjugate gradient\n"
               "  -bicgstab          BiCGStab\n"
               "\n");
}

void printHelp(CliArgs &args, bool const is_timing)
{
  if (is_timing) {
    masterPrintf("This program measures the performance of the fundamental "
                 "operations in QPhiX.\n\n");
  } else {
    masterPrintf("This program verifies the correctness of the fundamental "
                 "operations in QPhiX by comparing them to QDP++.\n\n");
  }
  masterPrintf("The following options are available. If no default is listed, that "
               "option is a mandatory one.\n"
               "\n");

  args.some_user_args.printArgHelp();
  printArgHelp(args, is_timing);
}

CliArgs processArgs(int &argc, char **&argv, bool const is_timing)
{
  CliArgs args;

  if (argc == 1) {
    printHelp(args, is_timing);
    exit(1);
  }

  /*
  masterPrintf("argc = %i\n", argc);
  for (int i = 0; i < argc; ++i) {
      masterPrintf("argv[% 2d] = %s\n", i, argv[i]);
  }
  */


  // Let the class parse the arguments first. It will consume whatever it can parse.
  args.some_user_args.init(argc, argv);

  args.By = args.some_user_args.getBy();
  args.Bz = args.some_user_args.getBz();
  args.PadXY = args.some_user_args.getPxy();
  args.PadXYZ = args.some_user_args.getPxyz();
  args.NCores = args.some_user_args.getNCores();
  args.Sy = args.some_user_args.getSy();
  args.Sz = args.some_user_args.getSz();
  args.MinCt = args.some_user_args.getMinCt();

  /*
  masterPrintf("argc = %i\n", argc);
  for (int i = 0; i < argc; ++i) {
    masterPrintf("argv[% 2d] = %s\n", i, argv[i]);
  }
  */

  for (int i = 1; i < argc;) {
    std::string const arg(argv[i]);

    if (arg == "-x") {
      args.nrow_in[0] = atoi(argv[i + 1]);
      i += 2;
    } else if (arg == "-y") {
      args.nrow_in[1] = atoi(argv[i + 1]);
      i += 2;
    } else if (arg == "-z") {
      args.nrow_in[2] = atoi(argv[i + 1]);
      i += 2;
    } else if (arg == "-t") {
      args.nrow_in[3] = atoi(argv[i + 1]);
      i += 2;
    } else if (arg == "-i") {
      args.iters = atoi(argv[i + 1]);
      i += 2;
    } else if (arg == "-compress12") {
      args.compress12 = true;
      i++;
    } else if (arg == "-geom") {
      args.qmp_geometry[0] = atoi(argv[i + 1]);
      args.qmp_geometry[1] = atoi(argv[i + 2]);
      args.qmp_geometry[2] = atoi(argv[i + 3]);
      args.qmp_geometry[3] = atoi(argv[i + 4]);
      i += 5;

    } else if (is_timing && arg == "-dslash") {
      args.do_dslash = true;
      i++;

    } else if (is_timing && arg == "-mmat") {
      args.do_m = true;
      i++;
    } else if (is_timing && arg == "-cg") {
      args.do_cg = true;
      i++;

    } else if (is_timing && arg == "-bicgstab") {
      args.do_bicgstab = true;
      i++;
    }

    else if (arg == "-prec") {
      std::string user_arg(argv[i + 1]);
      if (user_arg.compare("f") == 0) {
        args.prec = FLOAT_PREC;
      }
      if (user_arg.compare("h") == 0) {
        args.prec = HALF_PREC;
      }
      if (user_arg.compare("d") == 0) {
        args.prec = DOUBLE_PREC;
      }
      i += 2;
    } else if (arg == "-soalen") {
      args.soalen = atoi(argv[i + 1]);
      i += 2;
    } else if (arg == "-bind") {
      args.thread_bind = true;
      i++;
    } else {
#if 0
      QPhiX::masterPrintf("ERROR: The option “%s” could not be parsed. Please have "
                          "a look at the supported options below.\n\n",
                          arg.c_str());
      args.some_user_args.printArgHelp();
      printArgHelp(args, is_timing);
      exit(1);
#endif
      /* Ignore unrecognized arguments -- they may be needed by other components */
      i++;
    }
  }

  if (is_timing &&
      (!args.do_dslash && !args.do_m && !args.do_cg && !args.do_bicgstab)) {
    QPhiX::masterPrintf("No timing case has been selected. Please select at least "
                        "one of them such that this program can do some actual "
                        "work.\n\n");
    printTimingCaseHelp();
    exit(1);
  }

  return args;
}
