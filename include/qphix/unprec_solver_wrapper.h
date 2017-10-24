/*
 * unprec_solver_wrapper.h
 *
 *  Created on: Oct 10, 2017
 *      Author: bjoo
 */

#ifndef INCLUDE_QPHIX_UNPREC_SOLVER_WRAPPER_H_
#define I
#include "qphix/abs_solver.h"
#include "qphix/linearOp.h"
#include "qphix/blas.h"
#include "qphix/full_spinor.h"
#include "qphix/print_utils.h"

namespace QPhiX {
template<typename FT, int V, int S, bool compress12, typename LinearOpType, int num_flav=1>
class UnprecSolverWrapper : public AbstractUnprecSolver<FT,V,S,compress12,num_flav> {
public:

  using Spinor = FullSpinor<FT,V,S,compress12>;
  using CBSpinor = typename Spinor::CBSpinor;

  UnprecSolverWrapper(AbstractSolver<FT,V,S,compress12>& theEOSolver_,
      LinearOpType& M_) : theEOSolver(theEOSolver_), M(M_), geom(M_.getGeometry()) {

    Geometry<FT,V,S,compress12>& geom = M.getGeometry();
    for(int flav=0; flav < num_flav; ++flav) {
      tmp_space1[flav] = geom.allocCBFourSpinor();
      tmp_space2[flav] = geom.allocCBFourSpinor();
    }

  }

 ~UnprecSolverWrapper() {
    for(int flav=0; flav < num_flav; ++flav) {
      geom.free(tmp_space1[flav]);
      geom.free(tmp_space2[flav]);
      tmp_space1[flav]=nullptr;
      tmp_space2[flav]=nullptr;
    }
 }

 void SourcePrepare(CBSpinor *const out[num_flav],
     CBSpinor const *const in_cb[num_flav],
     CBSpinor const *const in_other_cb[num_flav],
     int isign,
     int solve_cb) const override {

     for(int flav=0; flav < num_flav;++flav) {
       copySpinor(out[flav],in_cb[flav],geom, geom.getNSIMT());
     }

     M.M_diag_inv(tmp_space1,in_other_cb,isign);
     M.M_offdiag(tmp_space2,tmp_space1,isign,solve_cb);

     // out = out - tmp_space2 =>  y = ax + y with  y=out, x = tmp_space2 and a = -1
     for(int flav=0; flav < num_flav; ++flav) {
       axpy((double)(-1),tmp_space2[flav],out[flav],geom, geom.getNSIMT());
     }
  }

 void SolutionReconstruct(CBSpinor *const out_other_cb[num_flav],
     CBSpinor const *const in_other_cb[num_flav],
     CBSpinor const *const out_cb[num_flav],
     int isign,
     int solve_cb) const override {

     // tmp_space2 = in_other_cb
     for(int flav=0; flav < num_flav; ++flav) {
       copySpinor(tmp_space2[flav],in_other_cb[flav],geom, geom.getNSIMT());
     }

     // tmp_space_1 = M_(1-cb,cb) out_cb
     M.M_offdiag(tmp_space1, out_cb, isign, 1-solve_cb);

     // tmp_space2 =  -tmp_space1 + tmp_space2 =  in_other_cb - M_(1-cb,cb) out_cb

     for(int flav=0; flav < num_flav; ++flav) {
       axpy((double)(-1),tmp_space1[flav],tmp_space2[flav],geom, geom.getNSIMT());
     }

     // out_other_cb = A^{-1} ( in_other_cb - M_(1-cb,cb) out_cb )
     M.M_diag_inv(out_other_cb,tmp_space2,isign);
 }

 AbstractSolver<FT,V,S,compress12>& getEOSolver() const
 {
   return theEOSolver;
 }

private:
 AbstractSolver<FT,V,S,compress12>& theEOSolver;
 // Maybe I can do this nicer
 LinearOpType& M;
 Geometry<FT,V,S,compress12>& geom;
 CBSpinor* tmp_space1[num_flav];
 CBSpinor* tmp_space2[num_flav];
};

}



#endif /* INCLUDE_QPHIX_UNPREC_SOLVER_WRAPPER_H_ */
