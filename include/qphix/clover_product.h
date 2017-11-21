#pragma once

#include <qphix/blas_new_c.h>
#include <qphix/geometry.h>

#include <cassert>
#include <initializer_list>
#include <iomanip>
#include <iostream>
#include <vector>

namespace
{
size_t constexpr re = 0;
size_t constexpr im = 1;
int const n_blas_simt = 1;

// The even checkerboard is given by ( (x + y + z + t ) & 1 == 0 ) -> cb0 is even
int constexpr cb_even = 0;
int constexpr cb_odd = 1;
}

/**
  Complex multiplication accumulate.

  Computes \f$ (r + \mathrm i i) += (a + \mathrm i b) * (c + \mathrm i d) \f$.
  */
template <typename FT>
void cplx_mul_acc(
    FT &r_out, FT &i_out, FT const &a, FT const &b, FT const &c, FT const &d)
{
  r_out += a * c - b * d;
  i_out += a * d + b * c;
}

/**
  Wrapper for the clover multiplication function.

  The `struct` is needed in order to allow for partial template specialization in the
  `Clover`
  parameter.

  \tparam Clover Type of clover block to use, must be a type from Geometry such that
  there exists a
  specialization for it.
  */
template <typename FT, int veclen, int soalen, bool compress12, typename Clover>
struct InnerCloverProduct {
  /**
  Multiplies the clover term for a single lattice size to a spinor.

  This function is intended to be used in a loop over all lattice sites. It is
  expected from the
  caller to have figured out all the correct indices. There are template
  specializations for the two
  different types of clover term that are used in QPhiX.

  \param[out] out Output spinor block. It is assumed to be zeroed properly, the
  function will just
  accumulate values into that output variable. Use \ref QPhiX::zeroSpinor for that.
  \param[in] in Input spinor block.
  \param[in] clover Single clover block that contains the lattice site of the spinor.
  \param[in] xi SIMD index for the arrays with length `soalen`, as in the spinors.
  \param[in] veclen_idx SIMD index for the arrays with length `veclen`, as in the
  clover term.
  */
  static void multiply(
      typename ::QPhiX::Geometry<FT, veclen, soalen, compress12>::FourSpinorBlock &out,
      typename ::QPhiX::Geometry<FT, veclen, soalen, compress12>::FourSpinorBlock const
          &in,
      Clover const &clover,
      int const xi,
      int const veclen_idx);
};

template <typename FT, int veclen, int soalen, bool compress12>
struct InnerCloverProduct<
    FT,
    veclen,
    soalen,
    compress12,
    typename ::QPhiX::Geometry<FT, veclen, soalen, compress12>::CloverBlock> {
  static void multiply(
      typename ::QPhiX::Geometry<FT, veclen, soalen, compress12>::FourSpinorBlock
          &spinor_out,
      typename ::QPhiX::Geometry<FT, veclen, soalen, compress12>::FourSpinorBlock const
          &spinor_in,
      typename ::QPhiX::Geometry<FT, veclen, soalen, compress12>::CloverBlock const
          &clov_block,
      int const xi,
      int const veclen_idx)
  {
    // The clover term is block-diagonal in spin. Therefore we need
    // to iterate over the two blocks of spin.
    for (auto s_block : {0, 1}) {
      // Extract the diagonal and triangular parts.
      auto const &diag_in = s_block == 0 ? clov_block.diag1 : clov_block.diag2;
      auto const &off_diag_in =
          s_block == 0 ? clov_block.off_diag1 : clov_block.off_diag2;
      // Input two-spinor component.
      for (auto two_s_in : {0, 1}) {
        // Reconstruct four spinor index.
        auto const four_s_in = 2 * s_block + two_s_in;
        // Output two-spinor component.
        for (auto two_s_out : {0, 1}) {
          // Reconstruct four spinor index.
          auto const four_s_out = 2 * s_block + two_s_out;
          // Input color.
          for (auto c_in : {0, 1, 2}) {
            // Spin-color index (0, ..., 5).
            auto const sc_in = 3 * two_s_in + c_in;
            // Output color.
            for (auto c_out : {0, 1, 2}) {
              // Spin-color index (0, ..., 5).
              auto const sc_out = 3 * two_s_out + c_out;

              // See `qphix-codegen` file `dslash_common.cc`
              // function
              // `clover_term` for the index manipulations done
              // here.

              // Using separate loops over the actual indices is
              // probably
              // faster than the branching in the innermost loop.

              if (sc_out == sc_in) {
                cplx_mul_acc(spinor_out[c_out][four_s_out][re][xi],
                             spinor_out[c_out][four_s_out][im][xi],
                             diag_in[sc_in][veclen_idx],
                             QPhiX::rep<FT, double>(0.0),
                             spinor_in[c_in][four_s_in][re][xi],
                             spinor_in[c_in][four_s_in][im][xi]);
              } else if (sc_out < sc_in) {
                auto const idx15 = sc_in * (sc_in - 1) / 2 + sc_out;
                cplx_mul_acc(
                    spinor_out[c_out][four_s_out][re][xi],
                    spinor_out[c_out][four_s_out][im][xi],
                    off_diag_in[idx15][re][veclen_idx],
                    // aww hell, maybe one should just add negation to QPhiX::half ?
                    QPhiX::rep<FT, double>(
                        -QPhiX::rep<double, FT>(off_diag_in[idx15][im][veclen_idx])),
                    spinor_in[c_in][four_s_in][re][xi],
                    spinor_in[c_in][four_s_in][im][xi]);
              } else {
                auto const idx15 = sc_out * (sc_out - 1) / 2 + sc_in;
                cplx_mul_acc(spinor_out[c_out][four_s_out][re][xi],
                             spinor_out[c_out][four_s_out][im][xi],
                             off_diag_in[idx15][re][veclen_idx],
                             off_diag_in[idx15][im][veclen_idx],
                             spinor_in[c_in][four_s_in][re][xi],
                             spinor_in[c_in][four_s_in][im][xi]);
              }
            }
          }
        }
      }
    }
  }
};

template <typename FT, int veclen, int soalen, bool compress12>
struct InnerCloverProduct<
    FT,
    veclen,
    soalen,
    compress12,
    typename ::QPhiX::Geometry<FT, veclen, soalen, compress12>::FullCloverBlock> {
  static void multiply(
      typename ::QPhiX::Geometry<FT, veclen, soalen, compress12>::FourSpinorBlock
          &spinor_out,
      typename ::QPhiX::Geometry<FT, veclen, soalen, compress12>::FourSpinorBlock const
          &spinor_in,
      typename ::QPhiX::Geometry<FT, veclen, soalen, compress12>::FullCloverBlock const
          &clov_block,
      int const xi,
      int const veclen_idx)
  {
    // The clover term is block-diagonal in spin. Therefore we need
    // to iterate over the two blocks of spin.
    for (auto s_block : {0, 1}) {
      // handy reference to half-spinor block
      auto const &block_in = s_block == 0 ? clov_block.block1 : clov_block.block2;
      // Input two-spinor component.
      for (auto two_s_in : {0, 1}) {
        // Reconstruct four spinor index.
        auto const four_s_in = 2 * s_block + two_s_in;
        // Output two-spinor component.
        for (auto two_s_out : {0, 1}) {
          // Reconstruct four spinor index.
          auto const four_s_out = 2 * s_block + two_s_out;
          // Input color.
          for (auto c_in : {0, 1, 2}) {
            // Spin-color index (0, ..., 5).
            auto const sc_in = 3 * two_s_in + c_in;
            // Output color.
            for (auto c_out : {0, 1, 2}) {
              // Spin-color index (0, ..., 5).
              auto const sc_out = 3 * two_s_out + c_out;

              cplx_mul_acc(spinor_out[c_out][four_s_out][re][xi],
                           spinor_out[c_out][four_s_out][im][xi],
                           block_in[sc_out][sc_in][re][veclen_idx],
                           block_in[sc_out][sc_in][im][veclen_idx],
                           spinor_in[c_in][four_s_in][re][xi],
                           spinor_in[c_in][four_s_in][im][xi]);
            }
          }
        }
      }
    }
  }
};

/**
  Multiplies a checkerboarded QPhiX Clover term with a checkerboarded QPhiX spinor.

  Padding is taken care of. A test case for (a copy of) this function exists in
  QPhiX.

  If the preprocessor macro `PRINT_MAPPING` is defined, it will print out the mapping
  of `(x, y, z,
  t)` coordinates to block indices. Also it will check that each block is accessed
  the proper number
  of times, that is `soalen` for spinors and `veclen` for clover blocks.

  \param[out] out Output spinor
  \param[in] in Input spinor
  \param[in] clover Clover block
  \param[in] geom Geometry object holding the dimension of clover and spinor
  */
template <typename FT, int veclen, int soalen, bool compress12, typename Clover>
void clover_product(
    typename ::QPhiX::Geometry<FT, veclen, soalen, compress12>::FourSpinorBlock
        *const out,
    typename ::QPhiX::Geometry<FT, veclen, soalen, compress12>::FourSpinorBlock const
        *const in,
    Clover *clover,
    ::QPhiX::Geometry<FT, veclen, soalen, compress12> &geom)
{
  ::QPhiX::zeroSpinor<FT, veclen, soalen, compress12>(out, geom, n_blas_simt);

#ifdef PRINT_MAPPING
  std::vector<int> spin_touches(geom.getPxyz() * geom.Nt(), 0);
  std::vector<int> clover_touches(geom.getPxyz() * geom.Nt() * soalen / veclen, 0);

  std::cout << std::setw(3) << "x" << std::setw(3) << "y" << std::setw(3) << "z"
            << std::setw(3) << "t"
            << ":" << std::setw(5) << "spin" << std::setw(5) << "clov"
            << "\n";
#endif

  // Iterate through all the block.
  for (int t = 0; t < geom.Nt(); ++t) {
    for (int z = 0; z < geom.Nz(); ++z) {
      for (int y = 0; y < geom.Ny(); ++y) {
        for (int x = 0; x < geom.Nxh(); ++x) {
          // First element in the current XY plane at desired Z and T.
          auto const xyBase = t * geom.getPxyz() + z * geom.getPxy();
          // Index of the SoA along the X direction.
          auto const xb = x / soalen;
          // Index within the SoA.
          auto const xi = x % soalen;
          // Global spin block index.
          auto const spin_block_idx = xb + geom.Nxh() / soalen * y + xyBase;
          // Global clover/gauge block index.
          auto const clov_block_idx =
              xb + (y / geom.nGY()) * geom.Nxh() / soalen + xyBase / geom.nGY();
          // Index of the SoA structure within the current tile.
          // auto const tile = (geom.Nxh() / soalen * y + xyBase) % geom.nGY();
          auto const tile = y % geom.nGY();
          // Vector index for clover/gauge. The SoA index only runs to
          // `soalen`, this index needs to run to `veclen`, that is across the
          // various SoA within the tile.
          auto const veclen_idx = soalen * tile + xi;

#ifdef PRINT_MAPPING
          ++spin_touches[spin_block_idx];
          ++clover_touches[clov_block_idx];

          std::cout << std::setw(3) << x << std::setw(3) << y << std::setw(3) << z
                    << std::setw(3) << t << ":" << std::setw(5) << spin_block_idx
                    << std::setw(5) << clov_block_idx << "\n";
#endif

          assert(xi + xb * soalen == x);

          // References to the objects at desired block.
          auto const &clov_block = clover[clov_block_idx];
          auto const &spinor_in = in[spin_block_idx];
          auto &spinor_out = out[spin_block_idx];

          InnerCloverProduct<FT, veclen, soalen, compress12, Clover>::multiply(
              spinor_out, spinor_in, clov_block, xi, veclen_idx);
        }
      }
    }
  }

#ifdef PRINT_MAPPING
  std::cout << std::flush;

  // Make sure that each block got touched the correct number of times.
  for (int i = 0; i != spin_touches.size(); ++i) {
    if (spin_touches[i] != soalen) {
      std::cout << "Spin missmatch: Block " << std::setw(4) << i << " accessed "
                << std::setw(4) << spin_touches[i] << " times instead of " << soalen
                << "\n";
    }
  }

  for (int i = 0; i != clover_touches.size(); ++i) {
    if (clover_touches[i] != veclen) {
      std::cout << "Clover missmatch: Block " << std::setw(4) << i << " accessed "
                << std::setw(4) << clover_touches[i] << " times instead of " << veclen
                << "\n";
    }
  }

  std::cout << std::flush;
#endif
}
