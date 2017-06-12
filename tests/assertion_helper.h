#pragma once

// FIXME This header file is not stand-alone.

/**
  \file Helpers for various checks.
  */

/**
  \todo Update the function signature such that `target_cb` and `isign` are
  removed. They have nothing to do with this function, really.
  */
template <typename FT, int V, int soalen, bool compress, typename Phi>
void print_spinor_differences(const Phi &chi1,
                              const Phi &chi2,
                              const Geometry<FT, V, soalen, compress> &geometry,
                              const int target_cb,
                              const int isign,
                              const uint_fast64_t max_errors = 0)
{
    Phi diff = chi1 - chi2;

    Double diff_norm = sqrt(norm2(diff, rb[target_cb])) /
                       (Real(4 * 3 * 2 * Layout::vol()) / Real(2));
    QDPIO::cout << "\t cb = " << target_cb << "  isign = " << isign
                << "  diff_norm = " << diff_norm << endl;

    const auto Nt = geometry.Nt();
    const auto Nx = geometry.Nx();
    const auto Ny = geometry.Ny();
    const auto Nz = geometry.Nz();
    const auto Pxy = geometry.getPadXY();
    const auto Pxyz = geometry.getPadXYZ();
    const auto nvecs = geometry.nVecs();
    const auto nyg = geometry.nGY();

    const auto Nxh = Nx / 2;

    if (toBool(diff_norm < tolerance<FT>::small)) {
        assertion(toBool(diff_norm < tolerance<FT>::small));
        return;
    }

    uint_fast64_t errors = 0;

    for (int t = 0; t < Nt; t++) {
        for (int z = 0; z < Nz; z++) {
            for (int y = 0; y < Ny; y++) {
                for (int x = 0; x < Nxh; x++) {

                    // These are unpadded QDP++ indices...
                    int ind = x + Nxh * (y + Ny * (z + Nz * t));
                    for (int s = 0; s < Ns; s++) {
                        for (int c = 0; c < Nc; c++) {
                            auto d = diff.elem(rb[target_cb].start() + ind)
                                         .elem(s)
                                         .elem(c);
                            auto c1 = chi1.elem(rb[target_cb].start() + ind)
                                          .elem(s)
                                          .elem(c);
                            auto c2 = chi2.elem(rb[target_cb].start() + ind)
                                          .elem(s)
                                          .elem(c);

                            REAL dr = d.real();
                            REAL di = d.imag();
                            if (toBool(fabs(dr) > tolerance<FT>::small) ||
                                toBool(fabs(di) > tolerance<FT>::small)) {
                                QDPIO::cout << "(x,y,z,t)=(" << x << "," << y
                                            << "," << z << "," << t
                                            << ") site=" << ind << " spin=" << s
                                            << " color=" << c << " Diff=" << d
                                            << " χ₁=" << c1 << " χ₂=" << c2
                                            << "\n";
                                ++errors;

                                if (max_errors > 0 && errors >= max_errors) {
                                    QDPIO::cout << "Too many errors emitted, "
                                                   "aborting.\n";
                                    assert(false);
                                    return;
                                }
                            }
                        }
                    }
                } // x
            } // y
        } // z
    } // t
}
