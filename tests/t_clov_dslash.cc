#include "qdp.h"
#include "unittest.h"

#include "testClovDslashFull.h"
#include <iostream>
#include <cstdio>
#include <omp.h>

#include "cli_args.h"

using namespace QDP;
using namespace std;

int main(int argc, char **argv)
{
  // Initialize UnitTest jig
  processArgs(argc, argv);

  omp_set_num_threads(some_user_args.getNCores() * some_user_args.getSy() *
                      some_user_args.getSz());

  TestRunner tests(&argc, &argv, nrow_in);

  const multi1d<int> &localLattSize = Layout::subgridLattSize();

  tests.addTest(new testClovDslashFull(some_user_args, compress12, prec_user),
                "testClovDslashFull");

  // Run all tests
  tests.run();

  // Testjig is destroyed
  tests.summary();
}
