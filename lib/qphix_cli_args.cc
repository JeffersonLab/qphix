/*
 * qphix_cli_args.h
 *
 *  Created on: Aug 1, 2016
 *      Author: bjoo
 */

#include "qphix/qphix_cli_args.h"
#include "qphix/print_utils.h"

#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

namespace QPhiX
{

void QPhiXCLIArgs::init(int &argc, char **&argv)
{
  bool FoundBy = false;
  bool FoundBz = false;
  bool FoundPxy = false;
  bool FoundPxyz = false;
  bool FoundNCores = false;
  bool FoundSy = false;
  bool FoundSz = false;
  bool FoundMinCt = false;

  std::vector<char *> new_argv;

  for (int i = 0; i < argc;) {
    std::string const arg(argv[i]);

    if (arg == "-by") {
      By = atoi(argv[i + 1]);
      i += 2;
      FoundBy = true;
    } else if (arg == "-bz") {
      Bz = atoi(argv[i + 1]);
      i += 2;
      FoundBz = true;
    } else if (arg == "-c") {
      NCores = atoi(argv[i + 1]);
      i += 2;
      FoundNCores = true;
    } else if (arg == "-sy") {
      Sy = atoi(argv[i + 1]);
      i += 2;
      FoundSy = true;
    } else if (arg == "-sz") {
      Sz = atoi(argv[i + 1]);
      i += 2;
      FoundSz = true;
    } else if (arg == "-pxy") {
      Pxy = atoi(argv[i + 1]);
      i += 2;
      FoundPxy = true;
    } else if (arg == "-pxyz") {
      Pxyz = atoi(argv[i + 1]);
      i += 2;
      FoundPxyz = true;
    } else if (arg == "-minct") {
      MinCt = atoi(argv[i + 1]);
      i += 2;
      FoundMinCt = true;
    } else {
      new_argv.push_back(argv[i]);
      i++;
    }

  } // While loop.

  if (!FoundPxy) {
    Pxy = 0;
    FoundPxy = true;
  }
  if (!FoundPxyz) {
    Pxyz = 0;
    FoundPxyz = true;
  }
  if (!FoundMinCt) {
    MinCt = 1;
    FoundMinCt = true;
  }

  bool FoundAll = FoundBy && FoundBz && FoundPxy && FoundPxyz && FoundNCores &&
                  FoundSy && FoundSz & FoundMinCt;

  if (!FoundAll) {
    printHelp();
    std::abort();
  }
  initedP = true;

#if 0
  // Copy the remaining arguments back to the argument list. This function will have
  // consumed the arguments that it has understood.
  for (int i = 0; i < argc; ++i) {
    argv[i] = nullptr;
  }
  argc = new_argv.size();
  for (int i = 0; i < new_argv.size(); ++i) {
    argv[i] = new_argv[i];
  }
#endif
}

void QPhiXCLIArgs::printHelp() const
{
  masterPrintf("There has been an issue with parsing the QPhiX options. The "
               "following is a list of all QPhiX options. All of them are "
               "required!\n"
               "\n");
  printArgHelp();
}

void QPhiXCLIArgs::printArgHelp() const
{
  masterPrintf("Internal data layout:\n"
               "  -by BY             block size in Y \n"
               "  -bz BZ             block size in Z \n"
               "  -pxy Pxy           extra pad in the XY plane \n"
               "  -pxyz Pxyz         extra pad in the XYZ plane \n"
               "  -minct MinCt       MinCt\n"
               "\n"
               "Threading options:\n"
               "  -c NCores          number of cores\n"
               "  -sy Sy             number of SMT threads in Y\n"
               "  -sz Sz             number of SMT threads in Z\n"
               "\n");
}
}
