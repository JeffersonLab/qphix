#include "qdp.h"
#include "unittest.h"

#include "testDslashFull.h"

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

  // QDP++ gets set up here
  TestRunner tests(&argc, &argv, nrow_in);

  const multi1d<int> &localLattSize = Layout::subgridLattSize();

  // If user doesn't specify block size, use local volume dimensions.
  if (By_user < 0) {
    By_user = localLattSize[1];
  }
  if (Bz_user < 0) {
    Bz_user = localLattSize[2];
  }

  for (int i = 0; i < iters; i++) {
#if defined(QPHIX_SCALAR_SOURCE)
    tests.addTest(new testDslashFull(By_user,
                                     Bz_user,
                                     NCores_user,
                                     Sy_user,
                                     Sz_user,
                                     PadXY_user,
                                     PadXYZ_user,
                                     MinCt_user,
                                     compress12,
                                     prec_user,
                                     1),
                  "testDslashFull_S1\n");

#elif defined(QPHIX_QPX_SOURCE)
    tests.addTest(new testDslashFull(By_user,
                                     Bz_user,
                                     NCores_user,
                                     Sy_user,
                                     Sz_user,
                                     PadXY_user,
                                     PadXYZ_user,
                                     MinCt_user,
                                     compress12,
                                     prec_user,
                                     4),
                  "testDslashFull_S4\n");
#else
    tests.addTest(new testDslashFull(By_user,
                                     Bz_user,
                                     NCores_user,
                                     Sy_user,
                                     Sz_user,
                                     PadXY_user,
                                     PadXYZ_user,
                                     MinCt_user,
                                     compress12,
                                     prec_user,
                                     4),
                  "testDslashFull_S4\n");

    tests.addTest(new testDslashFull(By_user,
                                     Bz_user,
                                     NCores_user,
                                     Sy_user,
                                     Sz_user,
                                     PadXY_user,
                                     PadXYZ_user,
                                     MinCt_user,
                                     compress12,
                                     prec_user,
                                     8),
                  "testDslashFull_S8\n");
    tests.addTest(new testDslashFull(By_user,
                                     Bz_user,
                                     NCores_user,
                                     Sy_user,
                                     Sz_user,
                                     PadXY_user,
                                     PadXYZ_user,
                                     MinCt_user,
                                     compress12,
                                     prec_user,
                                     16),
                  "testDslashFull_S16\n");

#endif
  }

  tests.run();

  // Testjig is destroyed
  tests.summary();
}
