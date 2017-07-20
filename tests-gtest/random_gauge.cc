#include <gtest/gtest.h>

#include "../tests/RandomGauge.h"
#include "compare_spinors.h"
#include <qphix/qphix_cli_args.h>

using namespace QPhiX;

extern int nrow_in[];
extern QPhiXCLIArgs some_user_args;

TEST(RandomGauge, clover) {
  Geometry<double, 1, 1, false> geom(Layout::subgridLattSize().slice(),
                                     some_user_args.getBy(),
                                     some_user_args.getBz(),
                                     some_user_args.getNCores(),
                                     some_user_args.getSy(),
                                     some_user_args.getSz(),
                                     some_user_args.getPxy(),
                                     some_user_args.getPxyz(),
                                     some_user_args.getMinCt());

  int t_bc = -1;
  RandomGauge<double, 1, 1, false> gauge(geom, t_bc);

  typedef typename RandomGauge<double, 1, 1, false>::QdpSpinor QdpSpinor;

  // Test spinors.
  QdpSpinor qdp_spinor_1;
  QdpSpinor qdp_spinor_2;
  QdpSpinor qdp_spinor_3;
  gaussian(qdp_spinor_1);

  // Apply the clover term and inverse, result should be the original spinor.
  for (int cb = 0; cb < 2; ++cb) {
    gauge.clov_qdp.apply(qdp_spinor_2, qdp_spinor_1, 1, cb);
    gauge.invclov_qdp.apply(qdp_spinor_3, qdp_spinor_2, 1, cb);
    expect_near(qdp_spinor_1, qdp_spinor_3, 1e-6, geom, cb);
  }
}
