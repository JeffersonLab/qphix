#pragma once

#include "immintrin.h"

namespace QPhiX
{
typedef unsigned short half;

template <typename T, int V, int S, bool compressP>
struct Types {
  typedef T FourSpinorBlock[3][4][2][S];
  typedef T TwoSpinorBlock[3][2][2][V];
  typedef T SU3MatrixBlock[8][(compressP ? 2 : 3)][3][2][V];

  struct CloverBlock {
    /**
      Diagonal part of upper spin-block, real.

      The indices are:

      - Spin (slow) and spin (fast).
      - SIMD length.
      */
    T diag1[6][V];

    /**
      Off-diagonal part of upper spin-block, complex.
      */
    T off_diag1[15][2][V];
    T diag2[6][V]; // Real Diagonal part of block 2
    T off_diag2[15][2][V]; // Complex, off diagonal part of block 2
  };

  /**
    Clover structure suitable for twisted mass.

    It has the following indices:

    - The fields `block1` and `block2` refer to the half-spinors, the blocks
    are in spin-space. Each of those, then has the following indices.

    - Combined spin and color index, length 6. The color part is faster, so let
    the combined index be `sc`. Then the (two-spinor) spin index is computed
    with `s2 = sc / 3` and the color index with `c = sc % 3`. A four-spinor
    spin index would be `s4 = 2 * block + s2`, where `block` is either 0 or 1,
    depending on whether `block1` or `block2` has been chosen.

    - Same as the previous index, this is a block-diagonal matrix in spin-color
    space and therefore has a second index of that kind.

    - Real and imaginary part, length 2.

    - SIMD vector, iterates through the \$f x \f$ coordinate, length
    `V`.
    */
  struct FullCloverBlock {
    /// Full complex, non-hermitian clover block 1
    T block1[6][6][2][V];
    /// Full complex, non-hermitian clover block 2
    T block2[6][6][2][V];
  };
};
}
