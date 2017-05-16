#ifndef TEST_TWM_DSLASH_FULL
#define TEST_TWM_DSLASH_FULL

#ifndef UNITTEST_H
#include "unittest.h"
#endif

enum Prec { FLOAT_PREC = 0, HALF_PREC, DOUBLE_PREC };

class testTWMDslashFull : public TestFixture
{
 public:
  testTWMDslashFull(int By_,
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

  template <typename T, int V, int S, typename U, typename Phi>
  void testTWMDslashWrapper()
  {
    for (int t_bc = +1; t_bc >= -1; t_bc -= 2) {
      if (compress12) {
        testTWMDslash<T, V, S, true, U, Phi>(t_bc);
      } else {
        testTWMDslash<T, V, S, false, U, Phi>(t_bc);
      }
    }
  }

  template <typename T, int V, int S, typename U, typename Phi>
  void testTWMDslashAChiMBDPsiWrapper()
  {
    for (int t_bc = 1; t_bc >= -1; t_bc -= 2) {
      if (compress12) {
        testTWMDslashAChiMBDPsi<T, V, S, true, U, Phi>(t_bc);
      } else {
        testTWMDslashAChiMBDPsi<T, V, S, false, U, Phi>(t_bc);
      }
    }
  }

#if 1
  template <typename T, int V, int S, typename U, typename Phi>
  void testTWMMWrapper()
  {
    for (int t_bc = -1; t_bc <= +1; t_bc += 2) {
      if (compress12) {
        testTWMM<T, V, S, true, U, Phi>(t_bc);
      } else {
        testTWMM<T, V, S, false, U, Phi>(t_bc);
      }
    }
  }
#endif

  template <typename T, int V, int S, typename U, typename Phi>
  void testTWMCGWrapper()
  {
    for (int t_bc = -1; t_bc <= +1; t_bc += 2) {
      if (compress12) {
        testTWMCG<T, V, S, true, U, Phi>(t_bc);
      } else {
        testTWMCG<T, V, S, false, U, Phi>(t_bc);
      }
    }
  }

#if 1
  template <typename T, int V, int S, typename U, typename Phi>
  void testTWMBiCGStabWrapper(const U &u)
  {
    // for(int t_bc=-1; t_bc <= +1; t_bc+=2) {
    int t_bc = -1;
    if (compress12) {
      testTWMBiCGStab<T, V, S, true, U, Phi>(u, t_bc);
    } else {
      testTWMBiCGStab<T, V, S, false, U, Phi>(u, t_bc);
    }
    //}
  }
#endif

  template <typename T1,
            int VEC1,
            int SOA1,
            typename T2,
            int VEC2,
            int SOA2,
            typename U,
            typename Phi>
  void testTWMRichardsonWrapper(const U &u)
  {
    for (int t_bc = -1; t_bc <= +1; t_bc += 2) {
      if (compress12) {
        testTWMRichardson<T1, VEC1, SOA1, true, T2, VEC2, SOA2, U, Phi>(u, t_bc);
      } else {
        testTWMRichardson<T1, VEC1, SOA1, false, T2, VEC2, SOA2, U, Phi>(u, t_bc);
      }
    }
  }

  template <typename T, int V, int S, bool compress, typename U, typename Phi>
  void testTWMDslash(int t_bc);

  template <typename T, int V, int S, bool compress, typename U, typename Phi>
  void testTWMDslashAChiMBDPsi(int t_bc);

  template <typename T, int V, int S, bool compress, typename U, typename Phi>
  void testTWMM(int t_bc);

  template <typename T, int V, int S, bool compress, typename U, typename Phi>
  void testTWMCG(int t_bc);

  template <typename T, int V, int S, bool compress, typename U, typename Phi>
  void testTWMBiCGStab(const U &u, int t_bc);

  template <typename T1,
            int VEC1,
            int SOA1,
            bool compress,
            typename T2,
            int VEC2,
            int SOA2,
            typename U,
            typename Phi>
  void testTWMRichardson(const U &u, int t_bc);

  // FIXME There is no need for those functions to be members. It is
  // convenient, but passing the Geometry object into them makes much more
  // sense.

  template <typename QDPSpinor>
  void applyTwist(QDPSpinor &psi, double Mu, double Alpha, int isign, int target_cb);

  template <typename QDPSpinor>
  void
  applyInvTwist(QDPSpinor &psi, double Mu, double MuInv, int isign, int target_cb);

  template <typename QdpGauge, typename QdpSpinor>
  void qdp_apply_operator(QdpSpinor &out,
                          QdpSpinor const &chi,
                          QdpSpinor const &psi,
                          QDP::multi1d<QdpGauge> const &u_aniso,
                          double const Mu,
                          double const MuInv,
                          double const alpha,
                          double const beta,
                          int const isign,
                          int const target_cb);

  template <typename QdpGauge, typename QdpSpinor>
  void qdp_achimbdpsi(QdpSpinor &out,
                      QdpSpinor const &chi,
                      QdpSpinor const &psi,
                      QDP::multi1d<QdpGauge> const &u_aniso,
                      double const Mu,
                      double const MuInv,
                      double const alpha,
                      double const beta,
                      int const isign,
                      int const target_cb);
};

#endif
