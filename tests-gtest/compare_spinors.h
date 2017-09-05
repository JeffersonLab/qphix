#pragma once

#include <qphix/geometry.h>
#include <qphix/print_utils.h>

#include <gtest/gtest.h>
#include <qdp.h>

template <typename FT, int veclen, int soalen, bool compress12>
void expect_near(QDP::LatticeDiracFermionD &spinor_a,
                 QDP::LatticeDiracFermionD &spinor_b,
                 double const abs_err,
                 QPhiX::Geometry<FT, veclen, soalen, compress12> &geom,
                 int const target_cb)
{
  typedef QDP::LatticeDiracFermionD Phi;

  Phi diff = spinor_b - spinor_a;

  QDP::Double diff_norm = sqrt(QDP::norm2(diff, QDP::rb[target_cb])) /
                          (QDP::Real(4 * 3 * 2 * QDP::Layout::vol()) / QDP::Real(2));

  if (QDP::toBool(diff_norm < abs_err)) {
    return;
  }

  int printed_out = 0;
  bool failed = false;

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

              // EXPECT_NEAR(a.real(), b.real(), abs_err)
              //     << "(x,y,z,t)=(" << x << "," << y << "," << z << "," << t
              //     << ") site=" << ind << " spin=" << s << " color=" << c;
              // EXPECT_NEAR(a.imag(), b.imag(), abs_err)
              //     << "(x,y,z,t)=(" << x << "," << y << "," << z << "," << t
              //     << ") site=" << ind << " spin=" << s << " color=" << c;

              // if (std::fabs(a.real() - b.real()) > abs_err ||
              //     std::fabs(a.imag() - b.imag()) > abs_err) {
              //   ++printed_out;
              // }

              // if (printed_out > 10) {
              //   masterPrintf("More elements are not printed in order to make the "
              //                "output readable.\n");
              //   FAIL();
              // }

              if (std::fabs(diff_real) > abs_err || std::fabs(diff_imag) > abs_err) {
                failed = true;

                QPhiX::masterPrintf("(xyzt)=(%2d,%2d,%2d,%2d) site=%5d s=%d c=%d "
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

                ++printed_out;

                if (printed_out > 50) {
                  QPhiX::masterPrintf(
                      "More elements are not printed in order to make the "
                      "output readable.\n");
                  FAIL();
                  return;
                }
              }
            }
          }
        } // x
      } // y
    } // z
  } // t

  if (failed) {
      FAIL();
  }
}
