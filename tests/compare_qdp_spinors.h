#pragma once

#if !defined(COMPARE_QDP_SPINORS_GTEST) && !defined(COMPARE_QDP_SPINORS_CUSTOM)
#error "This header must not be included directly, use the versions with _gtest or _custom in their filename."
#endif

#include <qphix/geometry.h>
#include <qphix/print_utils.h>
#include <qphix/qdp_packer.h>

#include <qdp.h>

#include <iomanip>

template <typename FT,
          int veclen,
          int soalen,
          bool compress12,
          typename QdpSpinor = QDP::LatticeDiracFermionD>
class HybridSpinor
{
 public:
  typedef typename QPhiX::Geometry<FT, veclen, soalen, compress12>::FourSpinorBlock
      Spinor;

  HybridSpinor(QPhiX::Geometry<FT, veclen, soalen, compress12> &geom)
      : geom_(geom), even_(geom), odd_(geom)
  {
  }

  Spinor *operator[](int const cb) { return cb == 0 ? even_.get() : odd_.get(); }

  void pack() { QPhiX::qdp_pack_spinor<>(qdp_, even_.get(), odd_.get(), geom_); }
  void unpack() { QPhiX::qdp_unpack_spinor<>(even_.get(), odd_.get(), qdp_, geom_); }
  void zero()
  {
    qdp_ = QDP::zero;
    pack();
  }

  Spinor *even() { return even_.get(); }
  Spinor *odd() { return odd_.get(); }
  QdpSpinor &qdp() { return qdp_; }

 private:
  QPhiX::Geometry<FT, veclen, soalen, compress12> &geom_;

  QPhiX::FourSpinorHandle<FT, veclen, soalen, compress12> even_, odd_;
  QdpSpinor qdp_;
};

template <typename FT,
          int veclen,
          int soalen,
          bool compress12,
          typename QdpSpinor = QDP::LatticeDiracFermionD>
void expect_near(QdpSpinor const &spinor_a,
                 QdpSpinor const &spinor_b,
                 double const abs_err,
                 QPhiX::Geometry<FT, veclen, soalen, compress12> const &geom,
                 int const target_cb,
                 char const *const message = nullptr)
{
  QdpSpinor const diff = spinor_b - spinor_a;

  QDP::Double const diff_norm =
      sqrt(QDP::norm2(diff, QDP::rb[target_cb])) /
      (QDP::Real(4 * 3 * 2 * QDP::Layout::vol()) / QDP::Real(2));

  if (message != nullptr) {
    QDPIO::cout << "Spinor comparison: " << message << ": ";
  }
  QDPIO::cout << "diff/volume = " << diff_norm << ", limit = " << abs_err << std::endl;

  if (QDP::toBool(diff_norm < abs_err)) {
    return;
  }

  QDPIO::cout << "A = control, B = candidate" << std::endl;

  uint64_t printed_out = 0;

  for (int t = 0; t < geom.Nt(); t++) {
    for (int z = 0; z < geom.Nz(); z++) {
      for (int y = 0; y < geom.Ny(); y++) {
        for (int x = 0; x < geom.Nxh(); x++) {

          // These are unpadded QDP++ indices...
          int const ind = x + geom.Nxh() * (y + geom.Ny() * (z + geom.Nz() * t));
          for (int s = 0; s < QDP::Ns; s++) {
            for (int c = 0; c < QDP::Nc; c++) {
              auto const &a =
                  spinor_a.elem(QDP::rb[target_cb].start() + ind).elem(s).elem(c);
              auto const &b =
                  spinor_b.elem(QDP::rb[target_cb].start() + ind).elem(s).elem(c);
              double const diff_real = a.real() - b.real();
              double const diff_imag = a.imag() - b.imag();

              if (std::fabs(diff_real) > abs_err || std::fabs(diff_imag) > abs_err) {
                QPhiX::masterPrintf("(xyzt)=(%2d,%2d,%2d,%2d) site=%5d s=%d c=%d "
                                    "A=(% 12.5e,% 12.5e) B=(% 12.5e,% 12.5e) "
                                    "A-B=(% 12.5e,% 12.5e)\n",
                                    x,
                                    y,
                                    z,
                                    t,
                                    ind,
                                    s,
                                    c,
                                    a.real(),
                                    a.imag(),
                                    b.real(),
                                    b.imag(),
                                    diff_real,
                                    diff_imag);

                ++printed_out;

                if (printed_out > 20) {
                  QPhiX::masterPrintf("More elements are not printed in order to "
                                      "make the output readable.\n");
#ifdef FAIL
                  FAIL();
#else
                  assertion(false);
#endif
                  break;
                }
              }
            }
          }
        } // x
      } // y
    } // z
  } // t

#ifdef FAIL
                  FAIL();
#else
                  assertion(false);
#endif
}
