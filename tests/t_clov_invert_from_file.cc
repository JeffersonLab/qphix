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

  tests.run();
  tests.summary();
}
