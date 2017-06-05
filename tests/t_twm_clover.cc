#include "qdp.h"
#include "unittest.h"

#include "testTWMCloverFull.h"
#include "qphix/qphix_cli_args.h"
#include <iostream>
#include <cstdio>
#include <omp.h>

#include "cli_args.h"

using namespace QDP;
using namespace std;

int main(int argc, char **argv)
{
  TestRunner tests(&argc, &argv, nrow_in);

  tests.addTest(new testTWMCloverFull(tests.args()), "testTWMCloverFull");

  tests.run();
  tests.summary();
}
