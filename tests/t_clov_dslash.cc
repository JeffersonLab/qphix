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
  TestRunner tests(&argc, &argv, nrow_in);

  tests.addTest(new testClovDslashFull(tests.args()), "testClovDslashFull");

  tests.run();
  tests.summary();
}
