#pragma once

#include <qphix/geometry.h>

#include <qdp.h>

#include <iomanip>

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
              double const diff_real = a.real() - b.real();
              double const diff_imag = a.imag() - b.imag();

              if (std::fabs(diff_real) > abs_err || std::fabs(diff_imag) > abs_err) {
                QDPIO::cout << "(x,y,z,t)=(" << x << "," << y << "," << z << "," << t
                            //<< ") site=" << std::setw(5) << ind << " spin=" << s
                            //<< " color=" << c << "a=(" << std::scientific
                            //<< std::setw(15) << std::showpos << a.real() << ","
                            //<< std::scientific << std::setw(15) << std::showpos
                            //<< a.imag() << ") b=(" << std::scientific
                            //<< std::setw(15) << std::showpos << b.real() << ","
                            //<< std::scientific << std::setw(15) << std::showpos
                            //<< b.imag() << ") diff=(" << std::scientific
                            //<< std::setw(15) << std::showpos << diff_real << ","
                            //<< std::scientific << std::setw(15) << std::showpos
                            //<< diff_imag << ")" << std::endl;
                            << ") site=" << ind << " spin=" << s << " color=" << c
                            << "a=(" << a.real() << "," << a.imag() << ") b=("
                            << b.real() << "," << b.imag() << ") diff=(" << diff_real
                            << "," << diff_imag << ")" << std::endl;
              }
            }
          }
        } // x
      } // y
    } // z
  } // t

  assertion(false);
}
