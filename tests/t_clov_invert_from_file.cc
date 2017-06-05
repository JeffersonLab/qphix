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
  TestRunner tests(&argc, &argv, nrow_in);

  tests.addTest(new testClovInvertFromFile(tests.args()),
                "testClovInvertFromFile");

  tests.run();
  tests.summary();
}
