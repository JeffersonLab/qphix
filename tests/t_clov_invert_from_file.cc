#include "qdp.h"
#include "unittest.h"
#include "qphix/print_utils.h"

#include "testClovInvertFromFile.h"
#include <omp.h>
#include <iostream>
#include <cstdio>

#include "cli_args.h"

using namespace QDP;
using namespace std;

int main(int argc, char **argv)
{
  // Initialize UnitTest jig
  processArgs(argc, argv);

  omp_set_num_threads(NCores_user * Sy_user * Sz_user);
  TestRunner tests(&argc, &argv, nrow_in);
  const multi1d<int> &localLattSize = Layout::subgridLattSize();

  if (By_user < 0) {
    By_user = localLattSize[1];
  }
  if (Bz_user < 0) {
    Bz_user = localLattSize[2];
  }

  tests.addTest(new testClovInvertFromFile(By_user,
                                           Bz_user,
                                           NCores_user,
                                           Sy_user,
                                           Sz_user,
                                           PadXY_user,
                                           PadXYZ_user,
                                           MinCt_user,
                                           compress12,
                                           prec_user),
                "testClovInvertFromFile");

  // Run all tests
  tests.run();

  // Testjig is destroyed
  tests.summary();
}
