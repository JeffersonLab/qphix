#ifndef TIME_CLOV_NOQDP
#define TIME_CLOV_NOQDP

enum Prec { FLOAT_PREC=0, HALF_PREC, DOUBLE_PREC };

class timeClovNoQDP { 
public: 
 timeClovNoQDP(int By_, int Bz_, int NCores_, int Sy_, int Sz_, int PadXY_, int PadXYZ_, int MinCt_, int niters_, bool c12, Prec precision_, bool do_dslash_, bool do_m_, bool do_cg_, bool do_bicgstab_) : By(By_), Bz(Bz_), NCores(NCores_), Sy(Sy_), Sz(Sz_), PadXY(PadXY_), PadXYZ(PadXYZ_), MinCt(MinCt_), N_simt(Sy_*Sz_),  compress12(c12), iters(niters_), precision(precision_) , do_dslash(do_dslash_), do_m(do_m_), do_cg(do_cg_), do_bicgstab(do_bicgstab_) {}

  void run(const int lattSize[], const int qmp_geom[]); 
 private:

  template<typename FT, int V, int S, bool compress12>
  void runTest(const int lattSize[], const int qmp_geom[]); 

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
  const int iters;
  const Prec precision;
  const bool do_dslash;
  const bool do_m;
  const bool do_cg;
  const bool do_bicgstab;
};

#endif
