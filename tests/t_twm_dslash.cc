// $Id: t_dslash.cc,v 1.1 2008-08-26 13:24:29 bjoo Exp $

#include "qdp.h"
#include "unittest.h"

#include "testTWMDslashFull.h"

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
    tests.addTest(new testTWMDslashFull(By_user,
                                        Bz_user,
                                        NCores_user,
                                        Sy_user,
                                        Sz_user,
                                        PadXY_user,
                                        PadXYZ_user,
                                        MinCt_user,
                                        compress12,
                                        prec_user,
                                        g_soalen),
                  "testDslashFull");
  }

  tests.run();
  tests.summary();
}
