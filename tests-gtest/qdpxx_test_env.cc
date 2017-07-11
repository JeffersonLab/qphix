/*
 * qdpxx_test_env.cpp
 *
 *  Created on: Jul 10, 2017
 *      Author: traveluser
 */


#include "qdpxx_test_env.h"
#include "process_args.h"
#include "omp.h"

#ifndef QDP_INCLUDE
#include "qdp.h"
#endif

using namespace QDP;

namespace QPhiXTesting {

	/** The Constructor to set up a test environment.
	 *   Its job is essentially to set up QMP
	 */
	QDPXXTestEnv::QDPXXTestEnv(int  *argc, char ***argv)
	{
		_qphix_cli_args.init(*argc, *argv);
		// Set Up the Lattice from _nrow_in
		 omp_set_num_threads(_qphix_cli_args.getNCores() * _qphix_cli_args.getSy() *
		                      _qphix_cli_args.getSz());

		QDP_initialize(argc,argv);
		_nrow_in[0]=8; _nrow_in[1]=4; _nrow_in[2]=4; _nrow_in[3]=4;

		processArgs(*argc, *argv, _nrow_in);



		  multi1d<int> nrow(Nd);
		  nrow = static_cast<const int *>(_nrow_in);
		  Layout::setLattSize(nrow);
		  Layout::create();

		  const multi1d<int> &localLattSize = Layout::subgridLattSize();
		  QDPIO::cout << "Local Lattice Size = [ ";
		  for (int i = 0; i < 4; ++i) {
		    QDPIO::cout << localLattSize[i] << " ";
		  }
		  QDPIO::cout << "]\n";
	}

	QDPXXTestEnv::~QDPXXTestEnv()
	{
		// tear down QDPXX here.
		QDP_finalize();
	}

	static QDPXXTestEnv* _the_env = nullptr;
	const QDPXXTestEnv& getQDPXXTestEnv()
	{
		if( _the_env == nullptr ) {
			std::cerr << "QDPXXTestEnv Accessed while uninitialized" << std::endl;
			std::abort();
		}
		return (*_the_env);
	}

	/* This is a convenience routine to setup the test environment for GTest and its layered test environments */
	int QDPXXTestMain(int *argc, char **argv)
	{
		  ::testing::InitGoogleTest(argc, argv);
		  _the_env = new  QPhiXTesting::QDPXXTestEnv(argc,&argv);
		  ::testing::AddGlobalTestEnvironment(_the_env);
		  auto test_result = RUN_ALL_TESTS();
		  return test_result;
	}

} // Namespace
