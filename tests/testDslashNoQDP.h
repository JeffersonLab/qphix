#ifndef TEST_DSLASH_NOQDP
#define TEST_DSLASH_NOQDP


class testDslashNoQDP { 
public: 
 testDslashNoQDP(int By_, int Bz_, int Cy_, int Cz_, int Ct_, int N_simt_, bool c12) : By(By_), Bz(Bz_), Cy(Cy_), Cz(Cz_), Ct(Ct_), N_simt(N_simt_), compress12(c12) {}
  void run(const int lattSize[], const int qmp_geom[]); 
 private:
  const int By;
  const int Bz;
  const int Cy;
  const int Cz;
  const int Ct;
  const int N_simt;
  const bool compress12;
};

#endif
