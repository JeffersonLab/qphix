#pragma once

// FIXME This header file is not stand-alone.

#include <cassert>

/**
  A wrapper for a spinor and its checkerboard packings.

  In the test code, there are a lot of conversions between packed and unpacked
  spinors. These three variables have different names in the bulk test code.
  Then there is a fourth variable (an array) which holds both checkerboard parts
  such that it can be indexed with the checkerboard index. Also replicated is
  allocation and deallocation code, as well as random initialization and an
  initial packing. This class bundles these as one handy class.

  The three member fields are deliberately public. Consistency between the
  packed and unpacked fields is _not_ an invariant of this class. The user has
  to call \ref pack and \ref unpack when something has changed.
  */
template <typename FT, int V, int S, bool compress, typename Phi>
class PackedSpinor
{
 public:
  typedef typename Geometry<FT, V, S, compress>::FourSpinorBlock Spinor;

  PackedSpinor(Geometry<FT, V, S, compress> &geom);
  ~PackedSpinor();

  /// Updates the members `even` and `odd`.
  void pack() { qdp_pack_spinor<>(all, even, odd, geom); }

  /// Updates the member `all`.
  void unpack() { qdp_unpack_spinor<>(even, odd, all, geom); }

  /**
    Returns the checkerboarded part.

    \param[in] cb Checkerboarding index, must be 0 or 1.
    */
  Spinor *operator[](const int cb)
  {
    assert(cb == 0 || cb == 1);
    return cb == 0 ? even : odd;
  }

  /// Full field.
  Phi all;

  /// Fields on one checkerboard.
  Spinor *even, *odd;

 private:
  Geometry<FT, V, S, compress> &geom;
};

template <typename FT, int V, int S, bool compress, typename Phi>
PackedSpinor<FT, V, S, compress, Phi>::PackedSpinor(
    Geometry<FT, V, S, compress> &geom)
    : even(geom.allocCBFourSpinor()), odd(geom.allocCBFourSpinor()), geom(geom)
{
  gaussian(all);
  pack();
}

template <typename FT, int V, int S, bool compress, typename Phi>
PackedSpinor<FT, V, S, compress, Phi>::~PackedSpinor()
{
  geom.free(even);
  geom.free(odd);
}
