#ifndef TEST_MINVCG
#define TEST_MINVCG


#ifndef UNITTEST_H
#include "unittest.h"
#endif

enum Prec { FLOAT_PREC=0, HALF_PREC, DOUBLE_PREC };

class MInvCGTester : public TestFixture { 
public: 
 MInvCGTester(int By_, int Bz_, int NCores_, int Sy_, int Sz_, int PadXY_, int PadXYZ_, int MinCt_, bool c12, Prec precision_, const int soalen_) : By(By_), Bz(Bz_), NCores(NCores_), Sy(Sy_), Sz(Sz_), PadXY(PadXY_), PadXYZ(PadXYZ_), MinCt(MinCt_), N_simt(Sy_*Sz_), compress12(c12), precision(precision_), soalen(soalen_) {}
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

  template<typename T, int V, int S, typename U, typename Phi>
  void testMInvCGWrapper(const U& u) 
  {
    for(int t_bc=-1; t_bc <= +1; t_bc+=2) {
      if( compress12  ) { 
	testMInvCG<T,V,S,true,U,Phi>(u,t_bc);
      }
      else {
	testMInvCG<T,V,S,false,U,Phi>(u,t_bc);
      }
    }
  }


  template<typename T, int V, int S, bool compress, typename U, typename Phi>
    void testMInvCG(const U& u, int t_bc);
};

#endif
