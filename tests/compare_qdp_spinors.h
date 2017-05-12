#pragma once

#include <qphix/geometry.h>
#include <qphix/print_utils.h>

#include <qdp.h>

#include <iomanip>

template <typename FT, int veclen, int soalen, bool compress12, typename QdpSpinor>
void expect_near(QdpSpinor &spinor_a,
                 QdpSpinor &spinor_b,
                 double const abs_err,
                 QPhiX::Geometry<FT, veclen, soalen, compress12> &geom,
                 int const target_cb)
{
  QdpSpinor diff = spinor_b - spinor_a;

  QDP::Double diff_norm = sqrt(QDP::norm2(diff, QDP::rb[target_cb])) /
                          (QDP::Real(4 * 3 * 2 * QDP::Layout::vol()) / QDP::Real(2));

  if (QDP::toBool(diff_norm < abs_err)) {
    return;
  }

  QDPIO::cout << "Spinors are not near! diff/volume = " << diff_norm << std::endl;

  for (int t = 0; t < geom.Nt(); t++) {
    for (int z = 0; z < geom.Nz(); z++) {
      for (int y = 0; y < geom.Ny(); y++) {
        for (int x = 0; x < geom.Nxh(); x++) {

          // These are unpadded QDP++ indices...
          int ind = x + geom.Nxh() * (y + geom.Ny() * (z + geom.Nz() * t));
          for (int s = 0; s < QDP::Ns; s++) {
            for (int c = 0; c < QDP::Nc; c++) {
              auto &a =
                  spinor_a.elem(QDP::rb[target_cb].start() + ind).elem(s).elem(c);
              auto &b =
                  spinor_b.elem(QDP::rb[target_cb].start() + ind).elem(s).elem(c);
              double const diff_real = a.real() - b.real();
              double const diff_imag = a.imag() - b.imag();

              if (std::fabs(diff_real) > abs_err || std::fabs(diff_imag) > abs_err) {
                masterPrintf("(xyzt)=(%2d,%2d,%2d,%2d) site=%5d s=%d c=%d "
                             "A=(% 14.7e,% 14.7e) B=(% 14.7e,% 14.7e) "
                             "A-B=(% 14.7e,% 14.7e)\n",
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
              }
            }
          }
        } // x
      } // y
    } // z
  } // t

  assertion(false);
}
