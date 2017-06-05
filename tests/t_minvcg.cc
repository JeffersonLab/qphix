#include "qdp.h"
#include "unittest.h"

#include "testMInvCG.h"

#include <omp.h>
#include <iostream>
#include <cstdio>

#include "cli_args.h"

using namespace QDP;
using namespace std;

int main(int argc, char **argv)
{
  TestRunner tests(&argc, &argv, nrow_in);

  for (int i = 0; i < iters; i++) {
    tests.addTest(new MInvCGTester(tests.args()), "MinvCGTester\n");
  }

  tests.run();
  tests.summary();
}
