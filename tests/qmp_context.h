#pragma once

#include "cli_args.h"

#include <qphix/memory_usage.h>
#include <qphix/print_utils.h>
#include <qphix/threadbind.h>

#ifdef QPHIX_QMP_COMMS
#include <qmp.h>
#endif

#include <omp.h>

class QmpContext
{
 public:
  QmpContext(int &argc, char **&argv)
  {
#ifdef QPHIX_QMP_COMMS
    QMP_thread_level_t prv;
    if (QMP_init_msg_passing(&argc, &argv, QMP_THREAD_SINGLE, &prv) != QMP_SUCCESS) {
      QMP_error("Failed to initialize QMP\n");
      abort();
    }
    if (QMP_is_primary_node()) {
      printf("QMP IS INITIALIZED\n");
    }
#endif

    args_ = processArgs(argc, argv, true);
    omp_set_num_threads( (args_.NCores + args_.NCommCores) * args_.Sy * args_.Sz);

#ifdef QPHIX_QMP_COMMS
    if (QMP_declare_logical_topology(args_.qmp_geometry, 4) != QMP_SUCCESS) {
      QMP_error("Failed to declare QMP Logical Topology\n");
      abort();
    }

    if (QMP_is_primary_node()) {
      QPhiX::masterPrintf("Declared QMP Topology: %d %d %d %d\n",
                          args_.qmp_geometry[0],
                          args_.qmp_geometry[1],
                          args_.qmp_geometry[2],
                          args_.qmp_geometry[3]);
    }
#endif

#ifdef QPHIX_QPX_SOURCE
    if (thread_bind) {
      QPhiX::setThreadAffinity(args_.NCores, args_.Sy * args_.Sz);
    }

    QPhiX::reportAffinity();
#endif
  }

  ~QmpContext()
  {
    QPhiX::masterPrintf("Maximum resident memory size: %g MiB\n",
                        QPhiX::get_max_resident_mib());

#ifdef QPHIX_QMP_COMMS
    QMP_finalize_msg_passing();
#endif
  }

  CliArgs &args() { return args_; };

 private:
  CliArgs args_;
};
