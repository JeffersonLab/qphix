#pragma once

#include "unittest.h"

enum Prec { FLOAT_PREC = 0, HALF_PREC, DOUBLE_PREC };

#include <qphix/qphix_cli_args.h>
using namespace QPhiX;

/**
  Compare Martin's NDTM clover implementation with Peter's TM clover
  implementation.

  The non-degenerate case for zero mass splitting should give comparable
  results to the degenerate case. One has to be a bit careful because NDTM is a
  two-flavor operator whereas TM is a one-flavor operator.
  */
class testNDTMvsTMCloverFull : public TestFixture
{
public:
    testNDTMvsTMCloverFull(const QPhiXCLIArgs &GeomArgs_, bool c12, Prec precision_)
        : By(GeomArgs_.getBy()), Bz(GeomArgs_.getBz()),
          NCores(GeomArgs_.getNCores()), Sy(GeomArgs_.getSy()),
          Sz(GeomArgs_.getSz()), PadXY(GeomArgs_.getPxy()),
          PadXYZ(GeomArgs_.getPxyz()), MinCt(GeomArgs_.getMinCt()),
          N_simt(Sy * Sz), compress12(c12), precision(precision_),
          GeomArgs(GeomArgs_)
    {
    }

    /// Toplevel test function wrapper.
    void run(void);

private:
    /// Templated test function.
    template <typename FT,
              int V,
              int S,
              bool compress,
              typename U,
              typename Phi>
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
