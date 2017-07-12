#pragma once

#include "qphix/geometry.h"

namespace QPhiX
{

/**
  Base class for solvers.

  This type of (single-shift) solvers solves the equation \f$ M x = b \f$.
  */
template <typename FT, int V, int S, bool compress12, int num_flav = 1>
class AbstractSolver
{
 public:
  virtual ~AbstractSolver(){};

  typedef typename Geometry<FT, V, S, compress12>::FourSpinorBlock Spinor;

  /**
    Solves \f$ M x = b \f$ for \f$ x \f$.

    This is the original interface for the solver. Use this for working with a
    single flavor. It will delegate to the case of arbitrary number of
    flavors.

    \param[in,out] x Solution vector \f$ x \f$. The value it has when entering
    this function is used as an initial guess.
    \param[in] rhs Right hand side of the equation, \f$ b \f$.
    \param[in] RsdTarget Target residual, not squared.
    \param[out] niters Number of iterations used for the solve.
    \param[out] rsd_sq_final Squared final normalized residual, \f$ \| r \|^2
    / \| b \|^2 \f$.
    \param[out] site_flops Number of floating point operations per lattice
    site.
    \param[out] mv_apps Number of applications of \f$ M \f$ onto \f$ x \f$.
    \param[in] isign `+1` for the \f$ M \f$, `-1` for \f$ M^\dagger \f$.
    \param[in] verboseP Print additional information to `stdout` during solve.
    */
  virtual void operator()(Spinor *x,
                          const Spinor *rhs,
                          const double RsdTarget,
                          int &niters,
                          double &rsd_sq_final,
                          unsigned long &site_flops,
                          unsigned long &mv_apps,
                          int isign,
                          bool verboseP,
                          int target_cb = 1) const
  {
    Spinor *x_array[1] = {x};
    const Spinor *rhs_array[1] = {rhs};
    (*this)(x_array,
            rhs_array,
            RsdTarget,
            niters,
            rsd_sq_final,
            site_flops,
            mv_apps,
            isign,
            verboseP,
            target_cb);
  }

  /**
    Solves \f$ M x = b \f$ for \f$ x \f$ for arbitrary many flavors.

    \param[in,out] x Solution vector. The explicit array index is the flavor
    index. It is assumed that the linear operator will take this array of
    pointers as its argument.
    \param[in] rhs Right hand side, flavor index is explicit.

    The other parameters are exactly like the other `operator()`.
    */
  virtual void operator()(Spinor *const x[num_flav],
                          const Spinor *const rhs[num_flav],
                          const double RsdTarget,
                          int &niters,
                          double &rsd_sq_final,
                          unsigned long &site_flops,
                          unsigned long &mv_apps,
                          int isign,
                          bool verboseP,
                          int target_cb = 1) const = 0;

#ifdef __INTEL_COMPILER
  /**
    Workaround for the Intel C++ 17 compiler.

    \see The article \ref intel-cpp-compiler-workaround contains the
    motivation for this extra overload.
    */
  virtual void operator()(Spinor *const x[num_flav],
                          Spinor *const rhs[num_flav],
                          const double RsdTarget,
                          int &niters,
                          double &rsd_sq_final,
                          unsigned long &site_flops,
                          unsigned long &mv_apps,
                          int isign,
                          bool verboseP,
                          int target_cb = 1) const
  {
    (*this)(x,
            const_cast<Spinor const *const *>(rhs),
            RsdTarget,
            niters,
            rsd_sq_final,
            site_flops,
            mv_apps,
            isign,
            verboseP,
            target_cb);
  }
#endif

  virtual Geometry<FT, V, S, compress12> &getGeometry() = 0;
};

template <typename FT, int V, int S, bool compress12, int num_flav = 1>
class AbstractMultiSolver
{
 public:
  virtual ~AbstractMultiSolver(){};

  typedef typename Geometry<FT, V, S, compress12>::FourSpinorBlock Spinor;

  /**
    Solves \f$ (M + \mu_i \mathbb 1) x_i = b \f$ for \f$ x_i \f$.

    \param[out] x Solutions for the shifts. First index must be \p n_shifts
    long.
    \param[in] n_shift Number of shifts to compute.
    \param[in] shifts Array of length \p n_shifts with the shifted masses.

    For the remaining parameters, see \ref AbstractSolver::operator()().
    */
  virtual void operator()(Spinor **x,
                          const Spinor *rhs,
                          const int n_shift,
                          const double *shifts,
                          const double *RsdTarget,
                          int &niters,
                          double *rsd_sq_final,
                          unsigned long &site_flops,
                          unsigned long &mv_apps,
                          int isign,
                          bool verboseP,
                          int cb = 1) const
  {
    masterPrintf("Given is the following:\n");
    for (int s = 0; s < n_shift; ++s) {
      masterPrintf("  %p: x[%d] = %p\n", &x[s], s, x[s]);
    }

    // Create the a data structure for the `Spinor *` and populate it.
    std::vector<std::array<Spinor *, 1>> x_vec(n_shift);
    for (int s = 0; s < n_shift; ++s) {
      x_vec[s][0] = x[s];
    }

    std::vector<Spinor **> x_vec2(n_shift);
    for (int s = 0; s < n_shift; ++s) {
        x_vec2[s] = x_vec[s].data();
    }
    Spinor ***x_array = x_vec2.data();

    masterPrintf("Result is this::\n");
    masterPrintf("  %p: x_array = %p\n", &x_array, x_array);
    for (int s = 0; s < n_shift; ++s) {
      masterPrintf("  %p: x_array[%d] = %p\n", &x_array[0], s, x_array[s]);
      masterPrintf("  %p: x_array[%d][0] = %p\n", &x_array[s][0], s, x_array[s][0]);
    }

    const Spinor *rhs_array[1] = {rhs};

    (*this)(x_array,
            rhs_array,
            n_shift,
            shifts,
            RsdTarget,
            niters,
            rsd_sq_final,
            site_flops,
            mv_apps,
            isign,
            verboseP,
            cb);
  }

  /**
    Solves \f$ M x = b \f$ for \f$ x \f$ for arbitrary many flavors.

    \param[in,out] x Solution vector. The explicit array index is the flavor
    index. It is assumed that the linear operator will take this array of
    pointers as its argument.
    \param[in] rhs Right hand side, flavor index is explicit.

    The other parameters are exactly like the other `operator()`.
    */
  virtual void operator()(Spinor **const *x,
                          const Spinor *const *rhs,
                          const int n_shift,
                          const double *shifts,
                          const double *RsdTarget,
                          int &niters,
                          double *rsd_sq_final,
                          unsigned long &site_flops,
                          unsigned long &mv_apps,
                          int isign,
                          bool verboseP,
                          int cb = 1) const = 0;

#ifdef __INTEL_COMPILER
  /**
    Workaround for the Intel C++ 17 compiler.

    \see The article \ref intel-cpp-compiler-workaround contains the
    motivation for this extra overload.
    */
  virtual void operator()(Spinor **const *x,
                          Spinor *const *rhs,
                          const int n_shift,
                          const double *shifts,
                          const double *RsdTarget,
                          int &niters,
                          double *rsd_sq_final,
                          unsigned long &site_flops,
                          unsigned long &mv_apps,
                          int isign,
                          bool verboseP,
                          int cb = 1) const
  {
    (*this)(x,
            const_cast<Spinor const *const *>(rhs),
            n_shift,
            shifts,
            RsdTarget,
            niters,
            rsd_sq_final,
            site_flops,
            mv_apps,
            isign,
            verboseP,
            cb);
  }
#endif

  virtual Geometry<FT, V, S, compress12> &getGeometry() = 0;
};
}
