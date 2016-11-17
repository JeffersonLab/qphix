#include "qdp.h"
#include "unittest.h"

#include "testTWMCloverFull.h"
#include "qphix/qphix_cli_args.h"
#include <iostream>
#include <cstdio>
#include <omp.h>

using namespace QDP;
using namespace std;

int nrow_in[4]={4,4,4,4};
int iters=10000;
QPhiXCLIArgs some_user_args;
bool compress12=false;
Prec prec_user = FLOAT_PREC;

void printHelp()
{
  cout << "t_twm_clover Test Specific Args: -x Lx -y Ly -z Lz -t Lt -i iters -prec Prec" << endl;
  cout << "   Lx is the lattice size in X" << endl;
  cout << "   Ly is the lattice size in Y" << endl;
  cout << "   Lz is the lattice size in Z" << endl;
  cout << "   iters is the number of iterations " << endl;
  cout << "   Prec for precision " << endl;
  some_user_args.printHelp();
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
    else if (string(argv[i]).compare("-compress12") == 0 ) {
      compress12 =true;
      i++;
    }
    else if (string(argv[i]).compare("-prec") == 0 ) {
      string user_arg(argv[i+1]);
      if( user_arg.compare("h") == 0 ) {
        prec_user = HALF_PREC;
      }
      if( user_arg.compare("f") == 0 ) {
        prec_user = FLOAT_PREC;
      }
      if( user_arg.compare("d") == 0 ) {
        prec_user = DOUBLE_PREC;
      }
      i+=2 ;
    }
    else {
      i++;
    }

  }

  some_user_args.init(argc,argv);
}

int main(int argc, char **argv)
{
  // Initialize UnitTest jig
  processArgs(argc,argv);

  omp_set_num_threads(
      some_user_args.getNCores() *
      some_user_args.getSy() *
      some_user_args.getSz()
      );

  TestRunner tests(&argc, &argv, nrow_in);

  const multi1d<int>& localLattSize = Layout::subgridLattSize();

  tests.addTest(new testTWMCloverFull(some_user_args,compress12,prec_user), "testTWMCloverFull" );

  // Run all tests
  tests.run();

  // Testjig is destroyed
  tests.summary();
}

