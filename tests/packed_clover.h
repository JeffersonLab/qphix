#pragma once

// FIXME This header file is not stand-alone.

#include <cassert>

/**
  The member fields are deliberately public. Consistency between the
  packed and unpacked fields is _not_ an invariant of this class. The user has
  to call \ref pack and \ref unpack when something has changed.
  */
template <typename FT, int V, int S, bool compress>
struct PackedClover {
    typedef typename Geometry<FT, V, S, compress>::SU3MatrixBlock Gauge;
    typedef typename Geometry<FT, V, S, compress>::FullCloverBlock FullClover;

    PackedClover(Geometry<FT, V, S, compress> &geom);
    ~PackedClover();

    /**
      Initializes the clover term from a random gauge field.

      This function is a bit tricky to call. One needs to supply an additional
      `template` like this:

      ```{.cpp}
      MyClover zero_clover(geom);
      zero_clover.template init<U, Phi>(clparam);
      ```

      \tparam U Gauge field type
      \tparam Phi Spinor field type
      */
    template <class U, class Phi>
    void init(const CloverFermActParams &clparam);

    /**
      Returns the checkerboarded part.

      \param[in] cb Checkerboarding index, must be 0 or 1.
      */
    FullClover *operator[](const int cb)
    {
        assert(cb == 0 || cb == 1);
        return cb == 0 ? even : odd;
    }

    FullClover *clov_packed[2][2];
    FullClover *invclov_packed[2][2];

    Geometry<FT, V, S, compress> &geom;
};

template <typename FT, int V, int S, bool compress>
PackedClover<FT, V, S, compress>::PackedClover(
    Geometry<FT, V, S, compress> &geom):geom(geom)
{
    for (int cb : {0, 1}) {
        for (int pm : {0, 1}) {
            clov_packed[cb][pm] = geom.allocCBFullClov();
            invclov_packed[cb][pm] = geom.allocCBFullClov();
        }
    }
}

template <typename FT, int V, int S, bool compress>
PackedClover<FT, V, S, compress>::~PackedClover()
{
    for (int cb : {0, 1}) {
        for (int pm : {0, 1}) {
            geom.free(clov_packed[cb][pm]);
            geom.free(invclov_packed[cb][pm]);
        }
    }
}

template <typename FT, int V, int S, bool compress>
template <class U, class Phi>
void PackedClover<FT, V, S, compress>::init(const CloverFermActParams &clparam)
{
    multi1d<U> u(4);
    U g;
    U uf;
    for (int mu = 0; mu < 4; mu++) {
        uf = 1; // Unit gauge
        Real factor = Real(0.08);
        gaussian(g);
        u[mu] = uf + factor * g;
        reunit(u[mu]);
    }

    Gauge *u_packed[2];
    u_packed[0] = geom.allocCBGauge();
    u_packed[1] = geom.allocCBGauge();

    qdp_pack_gauge<>(u, u_packed[0], u_packed[1], geom);

    CloverTermT<Phi, U> clov_qdp;
    clov_qdp.create(u, clparam);

    CloverTermT<Phi, U> invclov_qdp(clov_qdp);
    for (int cb : {0, 1}) {
        invclov_qdp.choles(cb);
    }

    for (int cb : {0, 1}) {
        geom.free(u_packed[cb]);
    }

    for (int cb : {0, 1}) {
        for (int pm : {0, 1}) {
            qdp_pack_full_clover<>(clov_qdp, clov_packed[cb][pm], geom, cb);
            qdp_pack_full_clover<>(
                invclov_qdp, invclov_packed[cb][pm], geom, cb);
        }
    }
}

template <typename FT, int V, int S, bool compress>
void zero_packed_clover(PackedClover<FT, V, S, compress> &packed_clover)
{
    typedef typename Geometry<FT, V, S, compress>::FullCloverBlock FullClover;

    const auto &geom = packed_clover.geom;

    // XXX This expression has been taken out of the `Geometry` class. It does
    // not make sense to replicate it. In `Geometry`, it is private, for some
    // reason. Either make that public and const there, or write a getter
    // function.
    const size_t full_clover_bytes =
        ((geom.getPadXYZ() * geom.Nt() * S) / V) * sizeof(FullClover);

    for (int cb : {0, 1}) {
        for (int pm : {0, 1}) {
            for (auto field : {packed_clover.clov_packed[cb][pm],
                               packed_clover.invclov_packed[cb][pm]}) {
                memset(field, 0, full_clover_bytes);
            }
        }
    }
}
