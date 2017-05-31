#pragma once

/**
  Selects the correct kernel based on the selected twisted boundary conditions.

  The parameter names are so short such that this macro can be formatted by
  clang-format in a reasonable and readable way.

  \param k Name of the kernel function
  \param F Floating point type (usually `FT`)
  \param v Vector length (usually `veclen`)
  \param s SoA length (usually `soalen`)
  \param c Gauge compression (usually `compress12`)
  */
#define QPHIX_KERNEL_SELECT(k, F, v, s, c)                                          \
  (use_tbc[0]                                                                       \
       ? (use_tbc[1]                                                                \
              ? (use_tbc[2]                                                         \
                     ? (use_tbc[3] ? k<F, v, s, c, true, true, true, true>          \
                                   : k<F, v, s, c, true, true, true, false>)        \
                     : (use_tbc[3] ? k<F, v, s, c, true, true, false, true>         \
                                   : k<F, v, s, c, true, true, false, false>))      \
              : (use_tbc[2]                                                         \
                     ? (use_tbc[3] ? k<F, v, s, c, true, false, true, true>         \
                                   : k<F, v, s, c, true, false, true, false>)       \
                     : (use_tbc[3] ? k<F, v, s, c, true, false, false, true>        \
                                   : k<F, v, s, c, true, false, false, false>)))    \
       : (use_tbc[1]                                                                \
              ? (use_tbc[2]                                                         \
                     ? (use_tbc[3] ? k<F, v, s, c, false, true, true, true>         \
                                   : k<F, v, s, c, false, true, true, false>)       \
                     : (use_tbc[3] ? k<F, v, s, c, false, true, false, true>        \
                                   : k<F, v, s, c, false, true, false, false>))     \
              : (use_tbc[2]                                                         \
                     ? (use_tbc[3] ? k<F, v, s, c, false, false, true, true>        \
                                   : k<F, v, s, c, false, false, true, false>)      \
                     : (use_tbc[3] ? k<F, v, s, c, false, false, false, true>       \
                                   : k<F, v, s, c, false, false, false, false>))))
