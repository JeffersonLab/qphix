// $Id: t_dslash.cc,v 1.1 2008-08-26 13:24:29 bjoo Exp $

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
#include "timeDslashNoQDP.h"

int nrow_in[4]={4,4,4,4};
int iters=10000;
int By_user = -1;
int Bz_user = -1;
int NCores_user = -1;
int Sy_user = -1;
int Sz_user = -1;

// Hardwire these for now.
int PadXY_user = 1;
int PadXYZ_user = 0;
int MinCt_user = 1;

bool compress12=false;
int qmp_geometry[4]={1,1,1,1};
Prec prec_user = FLOAT_PREC;
bool thread_bind = false;

bool do_dslash = false;
bool do_m = false;
bool do_cg = false;
bool do_bicgstab = false;

void printHelp() 
{ 
       cout << "t_dslash -x Lx -y Ly -z Lz -t Lt -i iters -by BY -bz BZ -c NCores -sy SY -sz SZ  -pxy Pxy -pxyz Pxyz -minct MinCt -compress12 -geom Px Py Pz Pt -prec Prec -bind" << endl;
       cout << "   Lx is the lattice size in X" << endl;
       cout << "   Ly is the lattice size in Y" << endl;
       cout << "   Lz is the lattice size in Z" << endl;
       cout << "   iters is the number of iterations " << endl;
       cout << "   BY is the block size in Y " << endl;
       cout << "   BZ is the block size in Z " << endl;
       cout << "   NCores is the no of cores" << endl;
       cout << "   Sy is the no of simt threads in y" << endl;
       cout << "   Sz is the no of simt threads in z" << endl;
       cout << "   Pxy is the extra pad in the XY plane" << endl;
       cout << "   Pxyz is the extra pad in the XYZ plane" << endl;
       cout << "   MinCt is the MinCt parameter in the blocking scheme" << endl;
       cout << "   Px Py Pz Pt define a 4D grid of MPI tasks" << endl;
       cout << "   Prec for precision " << endl;	
	cout << "  -dslash Do dslash timing " << endl;
	cout << "  -mmat Do M testing " << endl;
	cout << "  -cg Do CG Timing " << endl;
	cout <<"   -bicgstab Do BiCGStab timing " << endl;
}

void processArgs(int argc, char *argv[]) 
{
  int i=1;
  if (argc == 1) {
    printHelp();
  }
  
  while( i < argc)  {
    if( string(argv[i]).compare("-x") == 0 ) {
      nrow_in[0]=atoi(argv[i+1]);
      i+=2;
    }
    else if ( string(argv[i]).compare("-y") == 0 ) {
      nrow_in[1]=atoi(argv[i+1]);
      i+=2;
    }
    else if ( string(argv[i]).compare("-z") == 0) {
      nrow_in[2]=atoi(argv[i+1]);
      i+=2;
    }
    else if ( string(argv[i]).compare("-t") == 0) {
      nrow_in[3]=atoi(argv[i+1]);
      i+=2;
    }
    else if ( string(argv[i]).compare("-i") == 0) {
      iters=atoi(argv[i+1]);
      i+=2;
    }
    else if ( string(argv[i]).compare("-by") == 0 ) {
      By_user=atoi(argv[i+1]);
      i+=2;
    }
    else if (string(argv[i]).compare("-bz") == 0 ) {
      Bz_user=atoi(argv[i+1]);
      i+=2;
    }
    else if ( string(argv[i]).compare("-c") == 0 ) {
      NCores_user=atoi(argv[i+1]);
      i+=2;
    }
    else if ( string(argv[i]).compare("-sy") == 0 ) {
      Sy_user=atoi(argv[i+1]);
      i+=2;
    }
    else if (string(argv[i]).compare("-sz") == 0 ) {
      Sz_user=atoi(argv[i+1]);
      i+=2;
    }
    else if ( string(argv[i]).compare("-pxy") == 0 ) {
      PadXY_user=atoi(argv[i+1]);
      i+=2;
    }
    else if (string(argv[i]).compare("-pxyz") == 0 ) {
      PadXYZ_user=atoi(argv[i+1]);
      i+=2;
    }
    else if (string(argv[i]).compare("-minct") == 0 ) {
      MinCt_user=atoi(argv[i+1]);
      i+=2;
    }
    else if (string(argv[i]).compare("-compress12") == 0 ) {
      compress12 =true;
      i++;
    }
    else if (string(argv[i]).compare("-prec") == 0 ) { 
	string user_arg(argv[i+1]);
	if( user_arg.compare("f") == 0 ) { 
	  prec_user = FLOAT_PREC;
	}

	if( user_arg.compare("h") == 0 ) { 
	  prec_user = HALF_PREC;
	}

	if( user_arg.compare("d") == 0 ) { 
	  prec_user = DOUBLE_PREC;
	}
	i+=2 ;
    }
    else if (string(argv[i]).compare("-geom") == 0 ) {
      qmp_geometry[0] = atoi(argv[i+1]);
      qmp_geometry[1] = atoi(argv[i+2]);
      qmp_geometry[2] = atoi(argv[i+3]);
      qmp_geometry[3] = atoi(argv[i+4]);
      i+=4;
      
    }
    else if (string(argv[i]).compare("-dslash") == 0 ) {
	do_dslash = true;
 	i++;

    } 
   else if (string(argv[i]).compare("-mmat") == 0 ) {
        do_m = true;
	i++;
    } 
   else if (string(argv[i]).compare("-cg") == 0 ) {
        do_cg = true;
	i++;

    } 
  else if (string(argv[i]).compare("-bicgstab") == 0 ) {
        do_bicgstab = true;
	i++;
    } 
    else {
      i++;
    }
    
  }

  if( NCores_user < 0 ) { printHelp(); exit(1); }
  if( Sy_user < 0 ) { printHelp(); exit(1); }
  if( Sz_user < 0 ) { printHelp(); exit(1); }
  
  
  // Ct does not have to divide t, we can pick that up.
  if( By_user < 0 ) { By_user = nrow_in[1]; }
  if( Bz_user < 0 ) { Bz_user = nrow_in[2]; }
  

}


int main(int argc, char **argv)
{
  // Initialize UnitTest jig
  processArgs(argc,argv);
  omp_set_num_threads(NCores_user*Sy_user*Sz_user);
  
#ifdef QPHIX_QMP_COMMS
  // Initialize QMP
    QMP_thread_level_t prv;
    if( QMP_init_msg_passing(&argc, &argv, QMP_THREAD_SINGLE, &prv) != QMP_SUCCESS ) { 
      QMP_error("Failed to initialize QMP\n");
      abort();
      
    }
  if ( QMP_is_primary_node() ) { 
	printf("QMP IS INITIALIZED\n");
  }
  
  // Declare the logical topology
  if ( QMP_declare_logical_topology(qmp_geometry, 4)!= QMP_SUCCESS ) { 
    QMP_error("Failed to declare QMP Logical Topology\n");
    abort();
  }
 #endif

  QPhiX::masterPrintf("Declared QMP Topology: %d %d %d %d\n", 
		qmp_geometry[0], qmp_geometry[1], qmp_geometry[2], qmp_geometry[3]);


#ifdef QPHIX_QPX_SOURCE
  if( thread_bind ) { 
    QPhiX::setThreadAffinity(NCores_user, Sy_user*Sz_user);
  }
  QPhiX::reportAffinity();
#endif

  QPhiX::masterPrintf("Launching TestCase\n");

  // Launch the test case. 
  timeDslashNoQDP test(By_user, Bz_user, NCores_user, Sy_user, Sz_user, PadXY_user, PadXYZ_user, MinCt_user,  iters, compress12, prec_user,do_dslash, do_m, do_cg, do_bicgstab);
  
  test.run(nrow_in, qmp_geometry);
#ifdef QPHIX_QMP_COMMS
  QMP_finalize_msg_passing();
#endif

}
