#ifndef TEST_CLOVDSLASH_FULL
#define TEST_CLOVDSLASH_FULL

#include "unittest.h"
#include "prec.h"
#include "tparam_selector.h"
#include "cli_args.h"

#include <qphix/qphix_cli_args.h>
using namespace QPhiX;

class testClovDslashFull : public TestFixture
{
 public:
  testClovDslashFull(const QPhiXCLIArgs &GeomArgs_, bool c12, Prec precision_)
      : By(GeomArgs_.getBy()), Bz(GeomArgs_.getBz()), NCores(GeomArgs_.getNCores()),
        Sy(GeomArgs_.getSy()), Sz(GeomArgs_.getSz()), PadXY(GeomArgs_.getPxy()),
        PadXYZ(GeomArgs_.getPxyz()), MinCt(GeomArgs_.getMinCt()), N_simt(Sy * Sz),
        compress12(c12), precision(precision_), GeomArgs(GeomArgs_)
  {
  }

  void run() override;

  template <typename FT,
            int veclen,
            int soalen,
            bool compress12,
            typename QdpGauge,
            typename QdpSpinor>
  void operator()();

 private:

  const int By;
  const int Bz;
  const int NCores;
  const int Sy;
  const int Sz;
  const int PadXY;
  const int PadXYZ;
  const int MinCt;
  const int N_simt;
  const bool compress12;
  const Prec precision;
  QPhiXCLIArgs GeomArgs;
};

#endif
