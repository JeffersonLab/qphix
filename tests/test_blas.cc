#include <iostream>
#include <cstdio>

using namespace std;
#include "qphix/qphix_config.h"
#include "qphix/print_utils.h"
#include "qphix/threadbind.h"

#ifdef QPHIX_QMP_COMMS
#include "qmp.h"
#endif

#include <omp.h>
#include "testBlas.h"

#include "cli_args.h"

int qmp_geom[4] = {1, 1, 1, 1};

using namespace QPhiX;

int main(int argc, char **argv)
{

  // Initialize UnitTest jig
  processArgs(argc, argv);

  for (int i = 0; i < 4; ++i) {
    qmp_geom[i] = nrow_in[i];
  }

  // Max no of threads = NCores * Sy * Sz
  omp_set_num_threads(NCores_user * Sy_user * Sz_user);

#ifdef QPHIX_QMP_COMMS

  // Initialize QMP
  if (QMP_is_initialized() == QMP_FALSE) {
    QMP_thread_level_t prv;
    if (QMP_init_msg_passing(&argc, &argv, QMP_THREAD_SINGLE, &prv) != QMP_SUCCESS) {
      QMP_error("Failed to initialize QMP\n");
      abort();
    }
  }

  // Declare the logical topology
  if (QMP_declare_logical_topology(qmp_geom, 4) != QMP_SUCCESS) {
    QMP_error("Failed to declare QMP Logical Topology\n");
    abort();
  }

#endif

  masterPrintf("Declared QMP Topology: %d %d %d %d\n",
               qmp_geom[0],
               qmp_geom[1],
               qmp_geom[2],
               qmp_geom[3]);
  masterPrintf("Launching TestCase\n");

  // Launch the test case.
  testBlas test(NCores_user, Sy_user, Sz_user, PadXY_user, PadXYZ_user, iters);

  test.run(nrow_in, qmp_geom);

  masterPrintf("Test Case Done\n");

#ifdef QPHIX_QMP_COMMS
  QMP_finalize_msg_passing();
#endif
}
