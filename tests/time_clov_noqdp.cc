// $Id: t_dslash.cc,v 1.1 2008-08-26 13:24:29 bjoo Exp $

#include <iostream>
#include <cstdio>

using namespace std;
#include "qphix/qphix_config.h"
#include "qphix/print_utils.h"
#include "qphix/threadbind.h"
#include "qphix/memory_usage.h"

#ifdef QPHIX_QMP_COMMS
#include "qmp.h"
#endif

#include <omp.h>
#include "timeClovNoQDP.h"

#include "cli_args.h"


#if 0
/******************************************************************************/
// Bind the openmp threads to hardware threads
//#define _GNU_SOURCE
#include <sys/types.h>
#include <unistd.h>
#include <sys/syscall.h>
#include <sched.h>

void setThreadAffinity(int nCores, int threadsPerCore)
{
#pragma omp parallel 
    {

        // Get the OpenMP thread number
        int tid = omp_get_thread_num();

        // Split into core and SIMT ID (assuming SIMT runs fastest) 
        int core  = tid/threadsPerCore;
        int smtid = tid - core*threadsPerCore;

        // Convert to hardware processor ID, basically using same scheme as 
        // Table 3-2 of the IBM redbook: 
        // http://www.redbooks.ibm.com/redbooks/pdfs/sg247948.pdf            

        // NB: 4 is hardwired here for BG/Q.  (NB: This would let us run with 
        // 'gaps' i.e. even 2 threads per core but still get the binding right.)       

        int hw_procid = smtid + 4*core;   
	
	cpu_set_t set;

        CPU_ZERO(&set);
        CPU_SET(hw_procid, &set);

        pid_t pid = (pid_t) syscall(SYS_gettid);
        // Bind the OMP threads to hardware threads
        if((sched_setaffinity(pid, sizeof(set), &set)) == -1)
	    std::cerr << "WARN: Cound not do sched_setaffinity\n" << std::endl;
  }
}

#include <spi/include/kernel/location.h>

void reportAffinity()
{

  uint32_t cids[64], htids[64];

#pragma omp parallel
  {
    htids[omp_get_thread_num()] = Kernel_ProcessorThreadID();
    cids[omp_get_thread_num()] = Kernel_ProcessorCoreID();
  }
  
  QPhiX::masterPrintf("ThreadBindings\n");
  for (int i = 0; i < omp_get_max_threads(); ++i)
    QPhiX::masterPrintf("OMP thread %d: core = %d, hardware thread = %d\n",
		 i, cids[i], htids[i]);
}

/******************************************************************************/
#endif

int main(int argc, char **argv)
{
  // Initialize UnitTest jig
  processArgs(argc, argv);
  omp_set_num_threads(NCores_user * Sy_user * Sz_user);

#ifdef QPHIX_QMP_COMMS
  // Initialize QMP
  QMP_thread_level_t prv;
  if (QMP_init_msg_passing(&argc, &argv, QMP_THREAD_SINGLE, &prv) != QMP_SUCCESS) {
    QMP_error("Failed to initialize QMP\n");
    abort();
  }
  if (QMP_is_primary_node()) {
    printf("QMP IS INITIALIZED\n");
  }

  // Declare the logical topology
  if (QMP_declare_logical_topology(qmp_geometry, 4) != QMP_SUCCESS) {
    QMP_error("Failed to declare QMP Logical Topology\n");
    abort();
  }

  if (QMP_is_primary_node()) {
    printf("Declared QMP Topology: %d %d %d %d\n",
           qmp_geometry[0],
           qmp_geometry[1],
           qmp_geometry[2],
           qmp_geometry[3]);
  }
#endif

#ifdef QPHIX_QPX_SOURCE
  if (thread_bind) {
    QPhiX::setThreadAffinity(NCores_user, Sy_user * Sz_user);
  }

  QPhiX::reportAffinity();
#endif

  QPhiX::masterPrintf("Launching TestCase\n");

  // Launch the test case.
  timeClovNoQDP test(By_user,
                     Bz_user,
                     NCores_user,
                     Sy_user,
                     Sz_user,
                     PadXY_user,
                     PadXYZ_user,
                     MinCt_user,
                     iters,
                     compress12,
                     prec_user,
                     do_dslash,
                     do_m,
                     do_cg,
                     do_bicgstab);

  test.run(nrow_in, qmp_geometry);

  QPhiX::masterPrintf("Maximum resident memory size: %g MiB\n",
                      QPhiX::get_max_resident_mib());

#ifdef QPHIX_QMP_COMMS
  QMP_finalize_msg_passing();
#endif
}
