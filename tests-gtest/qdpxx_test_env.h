/*
 * qdpxx_test_env.h
 *
 *  Created on: Jul 10, 2017
 *      Author: traveluser
 */

#ifndef TESTS_GTEST_QDPXX_TEST_ENV_H_
#define TESTS_GTEST_QDPXX_TEST_ENV_H_

#include "qphix/qphix_cli_args.h"
#include "gtest/gtest.h"

/** A Namespace for testing utilities */
namespace QPhiXTesting {

/** A Test Environment to set up QMP */
class QDPXXTestEnv : public ::testing::Environment {
public:
	QDPXXTestEnv(int *argc, char ***argv);
	~QDPXXTestEnv();

	const int* getNrowIn() const {
		return _nrow_in;
	}

	const QPhiX::QPhiXCLIArgs& getCLIArgs() const {
		return _qphix_cli_args;
	}
private:
	int _nrow_in[4];
	QPhiX::QPhiXCLIArgs _qphix_cli_args;
};

const QDPXXTestEnv& getQDPXXTestEnv();

int QDPXXTestMain(int *argc, char **argv);


} // Namespace QPhiXTesting



#endif /* TESTS_GTEST_QDPXX_TEST_ENV_H_ */
