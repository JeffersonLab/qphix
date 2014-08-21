// $Id: t_dslash.cc,v 1.1 2008-08-26 13:24:29 bjoo Exp $

#include <iostream>
#include <cstdio>

using namespace std;
#include <dslash_config.h>
#ifdef CPP_DSLASH_PARSCALAR
#include "qmp.h"
#endif

#include <omp.h>
#include "testDslashNoQDP.h"

int nrow_in[4]={4,4,4,4};
int iters=10000;
int By_user = -1;
int Bz_user = -1;
int Cy_user = -1;
int Cz_user = -1;
int Ct_user = -1;
int N_simt_user=-1;
bool compress12=false;
int qmp_geom[4]={-1,-1,-1,-1};

void printHelp() 
{ 
       cout << "t_dslash -x Lx -y Ly -z Lz -t Lt -i iters -by BY -bz BZ -cy CY -cz CZ -ct CT -nsimt N -compress12 -geom Px Py Pz Pt" << endl;
       cout << "   Lx is the lattice size in X" << endl;
       cout << "   Ly is the lattice size in Y" << endl;
       cout << "   Lz is the lattice size in Z" << endl;
       cout << "   iters is the number of iterations " << endl;
       cout << "   BY is the block size in Y " << endl;
       cout << "   BZ is the block size in Z " << endl;
       cout << "   Cy is the no of cores in Y" << endl;
       cout << "   Cz is the no of cores in Z" << endl;
       cout << "   Ct is the no of cores in T" << endl;
       cout << "   N is the number of SIMT threads per core " << endl;
       cout << "   Px Py Pz Pt define a 4D grid of MPI tasks" << endl;
}

void processArgs(int argc, char *argv[]) {
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
        else if ( string(argv[i]).compare("-cy") == 0 ) {
	  Cy_user=atoi(argv[i+1]);
            i+=2;
        }
        else if (string(argv[i]).compare("-cz") == 0 ) {
          Cz_user=atoi(argv[i+1]);
            i+=2;
        }
	else if (string(argv[i]).compare("-ct") == 0 ) {
          Ct_user=atoi(argv[i+1]);
	  i+=2;
        }
        else if (string(argv[i]).compare("-nsimt") == 0 ) {
            N_simt_user =atoi(argv[i+1]);
            i+=2;
        }
	else if (string(argv[i]).compare("-compress12") == 0 ) {
	  compress12 =true;
	  i++;
        }
	else if (string(argv[i]).compare("-geom") == 0 ) {
	  qmp_geom[0] = atoi(argv[i+1]);
	  qmp_geom[1] = atoi(argv[i+2]);
	  qmp_geom[2] = atoi(argv[i+3]);
	  qmp_geom[3] = atoi(argv[i+4]);
	  i+=4;

        }
        else {
           i++;
        }

    }

    if( Cy_user < 0 ) { printHelp(); exit(1); }
    if( Cz_user < 0 ) { printHelp(); exit(1); }
    if( Ct_user < 0 ) { 
      Ct_user = 1;   // Default (only YZ)
    }

    if( nrow_in[1] % Cy_user != 0  ) { 
      printf("Error. Cy_user=%d, does not divide Y_dimension = %d\n",
	     Cy_user, nrow_in[1]);
    }

    if( nrow_in[2] % Cz_user != 0  ) { 
      printf("Error. Cz_user=%d, does not divide Z_dimension = %d\n",
	     Cz_user, nrow_in[2]);
    }
    
    // Ct does not have to divide t, we can pick that up.
    if( By_user < 0 ) { By_user = nrow_in[1]/Cy_user; }
    if( Bz_user < 0 ) { Bz_user = nrow_in[2]/Cz_user; }
    if( N_simt_user < 0 ) { printHelp() ; exit(1); }


}

int main(int argc, char **argv)
{
  // Initialize UnitTest jig
  processArgs(argc,argv);
  omp_set_num_threads(Cy_user*Cz_user*Ct_user*N_simt_user);

  // Initialize QMP
  if( QMP_is_initialized() == QMP_FALSE ) { 
    QMP_thread_level_t prv;
    if( QMP_init_msg_passing(&argc, &argv, QMP_THREAD_SINGLE, &prv) != QMP_SUCCESS ) { 
      QMP_error("Failed to initialize QMP\n");
      abort();
  
    }
  }
  
  // Declare the logical topology
  if ( QMP_declare_logical_topology(qmp_geom, 4)!= QMP_SUCCESS ) { 
    QMP_error("Failed to declare QMP Logical Topology\n");
    abort();
  }
 
  if ( QMP_is_primary_node() ) { 
    printf("Declared QMP Topology: %d %d %d %d\n", 
	   qmp_geom[0], qmp_geom[1], qmp_geom[2], qmp_geom[3]);
  }

  if (QMP_is_primary_node()) {
    printf("Launching TestCase\n");
  }

  // Launch the test case. 
  testDslashNoQDP test(By_user, Bz_user, Cy_user, Cz_user, Ct_user, 
		       N_simt_user, compress12);

  test.run(nrow_in, qmp_geom);

  if(QMP_master_io_node()) { 
    printf("Test Case Done\n");
  }

  QMP_finalize_msg_passing();

}

