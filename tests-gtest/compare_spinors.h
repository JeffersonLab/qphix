#pragma once
#include <qphix/geometry.h>
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

              EXPECT_NEAR(a.real(), b.real(), abs_err)
                  << "(x,y,z,t)=(" << x << "," << y << "," << z << "," << t
                  << ") site=" << ind << " spin=" << s << " color=" << c;
              EXPECT_NEAR(a.imag(), b.imag(), abs_err)
                  << "(x,y,z,t)=(" << x << "," << y << "," << z << "," << t
                  << ") site=" << ind << " spin=" << s << " color=" << c;
            }
          }
        } // x
      } // y
    } // z
  } // t
}
