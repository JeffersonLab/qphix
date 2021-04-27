#include "qphix/qphix_cli_args.h"

#include <gtest/gtest.h>
#include <omp.h>
#include <qdp.h>
#include <string>
#include "process_args.h"

QPhiX::QPhiXCLIArgs some_user_args;
int nrow_in[4] = {4, 4, 4, 4};



int main(int argc, char **argv)
{
  ::testing::InitGoogleTest(&argc, argv);
  QDP::QDP_initialize(&argc,&argv);
  QPhiXTesting::processArgs(argc, argv, nrow_in);
  some_user_args.init(argc, argv);

  omp_set_num_threads(some_user_args.getNCores() * some_user_args.getSy() *
                      some_user_args.getSz());

  QDP::multi1d<int> nrow(QDP::Nd);
  nrow = static_cast<const int *>(nrow_in);
  QDP::Layout::setLattSize(nrow);
  QDP::Layout::create();

  const QDP::multi1d<int> &localLattSize = QDP::Layout::subgridLattSize();
  for (int i = 0; i < 4; ++i) {
    std::cout << localLattSize[i] << std::endl;
  }

  auto result =  RUN_ALL_TESTS();

  QDP::QDP_finalize();

  return result;
}
