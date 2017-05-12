#ifndef TEST_TWM_CLOVER_FULL
#define TEST_TWM_CLOVER_FULL

#ifndef UNITTEST_H
#include "unittest.h"
#endif

enum Prec { FLOAT_PREC = 0, HALF_PREC, DOUBLE_PREC };

#include <qphix/qphix_cli_args.h>
using namespace QPhiX;

class testTWMCloverFull : public TestFixture
{

 public:
  testTWMCloverFull(const QPhiXCLIArgs &GeomArgs_, bool c12, Prec precision_)
      : By(GeomArgs_.getBy()), Bz(GeomArgs_.getBz()), NCores(GeomArgs_.getNCores()),
        Sy(GeomArgs_.getSy()), Sz(GeomArgs_.getSz()), PadXY(GeomArgs_.getPxy()),
        PadXYZ(GeomArgs_.getPxyz()), MinCt(GeomArgs_.getMinCt()), N_simt(Sy * Sz),
        compress12(c12), precision(precision_), GeomArgs(GeomArgs_)
  {
  }

  // Toplevel test function wrapper
  void run(void);

 private:
  // Templated test function
  template <typename FT, int V, int S, bool compress, typename U, typename Phi>
  void runTest(void);

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
