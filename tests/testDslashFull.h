#ifndef TEST_DSLASH_FULL
#define TEST_DSLASH_FULL

#include "unittest.h"

#include "prec.h"

class testDslashFull : public TestFixture
{
 public:
  testDslashFull(int By_,
                 int Bz_,
                 int NCores_,
                 int Sy_,
                 int Sz_,
                 int PadXY_,
                 int PadXYZ_,
                 int MinCt_,
                 bool c12,
                 Prec precision_,
                 const int soalen_)
      : By(By_), Bz(Bz_), NCores(NCores_), Sy(Sy_), Sz(Sz_), PadXY(PadXY_),
        PadXYZ(PadXYZ_), MinCt(MinCt_), N_simt(Sy_ * Sz_), compress12(c12),
        precision(precision_), soalen(soalen_)
  {
  }
  void run(void);

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
  const int soalen;
  Seed rng_seed;

  template <typename T, int V, int S, bool compress, typename U, typename Phi>
  void testDslash(const multi1d<U> &u, int t_bc);

  template <typename T, int V, int S, bool compress, typename U, typename Phi>
  void testDslashAChiMBDPsi(const multi1d<U> &u, int t_bc);

  template <typename T, int V, int S, bool compress, typename U, typename Phi>
  void testM(const multi1d<U> &u, int t_bc);

  template <typename T, int V, int S, bool compress, typename U, typename Phi>
  void testCG(const multi1d<U> &u, int t_bc);

  template <typename T, int V, int S, bool compress, typename U, typename Phi>
  void testBiCGStab(const multi1d<U> &u, int t_bc);

  template <typename T1,
            int VEC1,
            int SOA1,
            bool compress,
            typename T2,
            int VEC2,
            int SOA2,
            typename U,
            typename Phi>
  void testRichardson(const multi1d<U> &u, int t_bc);
};

#endif
