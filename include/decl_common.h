#include "immintrin.h"

#ifndef CODEGEN_GEOMETRY
#define CODEGEN_GEOMETRY

template <typename T, int V, int S, bool compressP>
struct CodegenGeometry {
  typedef T FourSpinorBlock[3][4][2][S];
  typedef T TwoSpinorBlock[3][2][2][V];
  typedef T SU3MatrixBlock[8][(compressP ? 2 : 3)][3][2][V];

  struct CloverBlock {
    T diag1[6][V];
    T off_diag1[15][2][V];
    T diag2[6][V];
    T off_diag2[15][2][V];
  };

  struct FullCloverBlock {
    T block1[6][6][2][V];
    T block2[6][6][2][V];
  };
};

#endif
