#pragma once

#include "qphix/geometry.h"
#include "qphix/full_spinor.h"

namespace QPhiX
{
	enum ResiduumType { ABSOLUTE, RELATIVE, INVALID };

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
                          int target_cb = 1,
						  QPhiX::ResiduumType residType=QPhiX::RELATIVE) const
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
            target_cb,
			residType);
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
                          int target_cb = 1,
						  QPhiX::ResiduumType residType=QPhiX::RELATIVE) const = 0;

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
                          int target_cb = 1,
						  QPhiX::ResiduumType residType=QPhiX::RELATIVE) const
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
            target_cb,
			residType);
  }
#endif

  virtual Geometry<FT, V, S, compress12> &getGeometry() = 0;
};

template <typename FT, int V, int S, bool compress12>
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
                          int cb = 1) const = 0;

  virtual Geometry<FT, V, S, compress12> &getGeometry() = 0;
};

// Not sure yet how to make this for multiple num_flav
template<typename FT, int V, int S, bool compress12, int num_flav=1>
class AbstractUnprecSolver {
public:
    using Spinor =  FullSpinor<FT,V,S,compress12>;
    using CBSpinor = typename Spinor::CBSpinor;

    virtual
    AbstractSolver<FT,V,S,compress12>& getEOSolver() const = 0;

    virtual
    void SourcePrepare(CBSpinor *const out[num_flav],
        CBSpinor const *const in_cb[num_flav],
        CBSpinor const *const in_other_cb[num_flav],
        int isign,
        int solve_cb) const = 0;

    virtual
    void SolutionReconstruct(CBSpinor *const out_other_cb[num_flav],
          CBSpinor const *const in_other_cb[num_flav],
          CBSpinor const *const out_cb[num_flav],
          int isign,
          int solve_cb) const = 0;


    // Virtual destructor
    virtual ~AbstractUnprecSolver() {}

    virtual void operator()(Spinor* x,
                             const Spinor* rhs,
                             const double RsdTarget,
                             int &niters,
                             double &rsd_sq_final,
                             unsigned long &site_flops,
                             unsigned long &mv_apps,
                             int isign,
                             bool verboseP,
                             int solve_cb = 1,
                             QPhiX::ResiduumType residType=QPhiX::RELATIVE) const {
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
          solve_cb,
          residType);

    }

    virtual void operator()(Spinor *const x[num_flav],
                            const Spinor* const rhs[num_flav],
                            const double RsdTarget,
                            int &niters,
                            double &rsd_sq_final,
                            unsigned long &site_flops,
                            unsigned long &mv_apps,
                            int isign,
                            bool verboseP,
                            int solve_cb = 1,
                            QPhiX::ResiduumType residType=QPhiX::RELATIVE) const {


      AbstractSolver<FT,V,S,compress12>& theEOSolver = getEOSolver();
      Geometry<FT,V,S,compress12>& geom = theEOSolver.getGeometry();

      CBSpinor* in_cb[num_flav];
      CBSpinor* in_other_cb[num_flav];

      CBSpinor* out_cb[num_flav];
      CBSpinor* out_other_cb[num_flav];

      // Th cb_solutiosn will be the actual solutions.
      for(int soln=0; soln < num_flav; ++soln) {
        in_cb[soln] = rhs[soln]->getCBData(solve_cb);
        in_other_cb[soln] = rhs[soln]->getCBData(1-solve_cb);

        // The actual solutions (on the target CB)
        out_cb[soln] = x[soln]->getCBData(solve_cb);
        out_other_cb[soln]= x[soln]->getCBData(1-solve_cb);
      }

      // Reuse space for out_other_cb to store the prepared sources.
      //
      SourcePrepare(out_other_cb, in_cb, in_other_cb, isign, solve_cb);

      // CB Solve
      theEOSolver(out_cb,
                  out_other_cb,
                  RsdTarget,
                  niters,
                  rsd_sq_final,
                  site_flops,
                  mv_apps,
                  isign,
                  verboseP,
                  solve_cb,
                  residType);

      // Recompute out_other_cb -- this writes directly into the x-s.
      SolutionReconstruct(out_other_cb,in_other_cb,out_cb,isign, solve_cb);

    }

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
                          int solve_cb = 1,
              QPhiX::ResiduumType residType=QPhiX::RELATIVE) const
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
            solve_cb,
            residType);
  }

  virtual
    void SourcePrepare(CBSpinor *const out[num_flav],
        CBSpinor  *const in_cb[num_flav],
        CBSpinor  *const in_other_cb[num_flav],
        int isign,
        int solve_cb) const
  {
    this->SourcePrepare(out,
                        const_cast<CBSpinor const *const *>(in_cb),
                        const_cast<CBSpinor const *const *>(in_other_cb),
                        isign,
                        solve_cb);
  }

   virtual
   void SolutionReconstruct(CBSpinor *const out_other_cb[num_flav],
         CBSpinor *const in_other_cb[num_flav],
         CBSpinor *const out_cb[num_flav],
         int isign,
         int solve_cb) const
   {
     this->SolutionReconstruct(out_other_cb,
           const_cast<CBSpinor const *const *>(in_other_cb),
           const_cast<CBSpinor const *const *>(out_cb),
           isign,
           solve_cb);
   }
#endif
};
}
