#include "qphix/qphix_cli_args.h"

#include <gtest/gtest.h>
#include <omp.h>
#include <qdp.h>

#include <string>

QPhiX::QPhiXCLIArgs some_user_args;
int nrow_in[4] = {4, 4, 4, 4};

void processArgs(int argc, char *argv[])
{

    int i = 1;
    while (i < argc) {
        if (std::string(argv[i]).compare("-x") == 0) {
            nrow_in[0] = atoi(argv[i + 1]);
            i += 2;
        } else if (std::string(argv[i]).compare("-y") == 0) {
            nrow_in[1] = atoi(argv[i + 1]);
            i += 2;
        } else if (std::string(argv[i]).compare("-z") == 0) {
            nrow_in[2] = atoi(argv[i + 1]);
            i += 2;
        } else if (std::string(argv[i]).compare("-t") == 0) {
            nrow_in[3] = atoi(argv[i + 1]);
            i += 2;
        }
#if 0 // Not implemented yet.
        else if (std::string(argv[i]).compare("-i") == 0) {
            iters = atoi(argv[i + 1]);
            i += 2;
        }
        else if (std::string(argv[i]).compare("-compress12") == 0) {
            compress12 = true;
            i++;
        } else if (std::string(argv[i]).compare("-prec") == 0) {
            std::string user_arg(argv[i + 1]);
            if (user_arg.compare("h") == 0) {
                prec_user = HALF_PREC;
            }
            if (user_arg.compare("f") == 0) {
                prec_user = FLOAT_PREC;
            }
            if (user_arg.compare("d") == 0) {
                prec_user = DOUBLE_PREC;
            }
            i += 2;
        }
#endif
        else {
            i++;
        }
    }
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);

    processArgs(argc, argv);
    some_user_args.init(argc, argv);

    omp_set_num_threads(some_user_args.getNCores() * some_user_args.getSy() *
                        some_user_args.getSz());

    QDP::QDP_initialize(&argc, &argv);
    QDP::multi1d<int> nrow(Nd);
    nrow = static_cast<const int *>(nrow_in);
    QDP::Layout::setLattSize(nrow);
    QDP::Layout::create();

    const multi1d<int>& localLattSize = Layout::subgridLattSize();
    for (int i = 0; i < 4; ++i) {
        std::cout << localLattSize[i] << std::endl;
    }

    return RUN_ALL_TESTS();

    QDP::QDP_finalize();
}
