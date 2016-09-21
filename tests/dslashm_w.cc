// $Id: dslashm_w.cc,v 1.16 2007-02-24 01:00:29 bjoo Exp $
/*! \file
 *  \brief Wilson-Dirac operator
 */

#include "dslashm_w.h"

using namespace QDP;
namespace QPhiX
{
//! General Wilson-Dirac dslash
/*!
 * DSLASH
 *
 * This routine is specific to Wilson fermions!
 *
 * Description:
 *
 * This routine applies the operator D' to Psi, putting the result in Chi.
 *
 *	       Nd-1
 *	       ---
 *	       \
 *   chi(x)  :=  >  U  (x) (1 - isign gamma  ) psi(x+mu)
 *	       /    mu			  mu
 *	       ---
 *	       mu=0
 *
 *	             Nd-1
 *	             ---
 *	             \    +
 *                +    >  U  (x-mu) (1 + isign gamma  ) psi(x-mu)
 *	             /    mu			   mu
 *	             ---
 *	             mu=0
 *
 * Arguments:
 *
 *  \param chi	      Pseudofermion field				(Write)
 *  \param u	      Gauge field					(Read)
 *  \param psi	      Pseudofermion field				(Read)
 *  \param isign      D'^dag or D'  ( +1 | -1 ) respectively		(Read)
 *  \param cb	      Checkerboard of output vector			(Read) 
 */

void dslash(LatticeFermionF& chi,
	    const multi1d<LatticeColorMatrixF>& u,
	    const LatticeFermionF& psi,
	    int isign, int cb)
{
  /*     F 
   *   a2  (x)  :=  U  (x) (1 - isign gamma  ) psi(x)
   *     mu          mu                    mu
   */
  /*     B           +
   *   a2  (x)  :=  U  (x-mu) (1 + isign gamma  ) psi(x-mu)
   *     mu          mu                       mu
   */
  // Recontruct the bottom two spinor components from the top two
  /*                        F           B
   *   chi(x) :=  sum_mu  a2  (x)  +  a2  (x)
   *                        mu          mu
   */
#if 1
  // NOTE: this is unrolled for 2 dimensions tests or some preproc hooks needed
  // for other Nd
  if (isign > 0)
  {
    chi[rb[cb]] = spinReconstructDir0Minus(u[0] * shift(spinProjectDir0Minus(psi), FORWARD, 0))
                + spinReconstructDir0Plus(shift(adj(u[0]) * spinProjectDir0Plus(psi), BACKWARD, 0))
#if QDP_ND >= 2
                + spinReconstructDir1Minus(u[1] * shift(spinProjectDir1Minus(psi), FORWARD, 1))
                + spinReconstructDir1Plus(shift(adj(u[1]) * spinProjectDir1Plus(psi), BACKWARD, 1))
#endif
#if QDP_ND >= 3
                + spinReconstructDir2Minus(u[2] * shift(spinProjectDir2Minus(psi), FORWARD, 2))
                + spinReconstructDir2Plus(shift(adj(u[2]) * spinProjectDir2Plus(psi), BACKWARD, 2))
#endif
#if QDP_ND >= 4
                + spinReconstructDir3Minus(u[3] * shift(spinProjectDir3Minus(psi), FORWARD, 3))
                + spinReconstructDir3Plus(shift(adj(u[3]) * spinProjectDir3Plus(psi), BACKWARD, 3))
#endif
#if QDP_ND >= 5
#error "Unsupported number of dimensions"
#endif
    ;
  }
  else
  {
    chi[rb[cb]] = spinReconstructDir0Plus(u[0] * shift(spinProjectDir0Plus(psi), FORWARD, 0))
                + spinReconstructDir0Minus(shift(adj(u[0]) * spinProjectDir0Minus(psi), BACKWARD, 0))
#if QDP_ND >= 2
                + spinReconstructDir1Plus(u[1] * shift(spinProjectDir1Plus(psi), FORWARD, 1))
                + spinReconstructDir1Minus(shift(adj(u[1]) * spinProjectDir1Minus(psi), BACKWARD, 1))
#endif
#if QDP_ND >= 3
                + spinReconstructDir2Plus(u[2] * shift(spinProjectDir2Plus(psi), FORWARD, 2))
                + spinReconstructDir2Minus(shift(adj(u[2]) * spinProjectDir2Minus(psi), BACKWARD, 2))
#endif
#if QDP_ND >= 4
                + spinReconstructDir3Plus(u[3] * shift(spinProjectDir3Plus(psi), FORWARD, 3))
                + spinReconstructDir3Minus(shift(adj(u[3]) * spinProjectDir3Minus(psi), BACKWARD, 3))
#endif
#if QDP_ND >= 5
#error "Unsupported number of dimensions"
#endif
    ;
  }
#else

  // NOTE: the loop is not unrolled - it should be all in a single line for
  // optimal performance
  chi[rb[cb]] = zero;

  // NOTE: temporarily has conversion call of LatticeHalfFermion - will be removed
  for(int mu = 0; mu < Nd; ++mu)
  {
    chi[rb[cb]] += spinReconstruct(LatticeHalfFermionF(u[mu] * shift(spinProject(psi,mu,-isign), FORWARD, mu)),mu,-isign)
      + spinReconstruct(LatticeHalfFermionF(shift(adj(u[mu]) * spinProject(psi,mu,+isign), BACKWARD, mu)),mu,+isign);
  }
#endif

}


void dslash(LatticeFermionD& chi,
	    const multi1d<LatticeColorMatrixD>& u,
	    const LatticeFermionD& psi,
	    int isign, int cb)
{
  /*     F
   *   a2  (x)  :=  U  (x) (1 - isign gamma  ) psi(x)
   *     mu          mu                    mu
   */
  /*     B           +
   *   a2  (x)  :=  U  (x-mu) (1 + isign gamma  ) psi(x-mu)
   *     mu          mu                       mu
   */
  // Recontruct the bottom two spinor components from the top two
  /*                        F           B
   *   chi(x) :=  sum_mu  a2  (x)  +  a2  (x)
   *                        mu          mu
   */
#if 1
  // NOTE: this is unrolled for 2 dimensions tests or some preproc hooks needed
  // for other Nd
  if (isign > 0)
  {
    chi[rb[cb]] = spinReconstructDir0Minus(u[0] * shift(spinProjectDir0Minus(psi), FORWARD, 0))
                + spinReconstructDir0Plus(shift(adj(u[0]) * spinProjectDir0Plus(psi), BACKWARD, 0))
#if QDP_ND >= 2
                + spinReconstructDir1Minus(u[1] * shift(spinProjectDir1Minus(psi), FORWARD, 1))
                + spinReconstructDir1Plus(shift(adj(u[1]) * spinProjectDir1Plus(psi), BACKWARD, 1))
#endif
#if QDP_ND >= 3
                + spinReconstructDir2Minus(u[2] * shift(spinProjectDir2Minus(psi), FORWARD, 2))
                + spinReconstructDir2Plus(shift(adj(u[2]) * spinProjectDir2Plus(psi), BACKWARD, 2))
#endif
#if QDP_ND >= 4
                + spinReconstructDir3Minus(u[3] * shift(spinProjectDir3Minus(psi), FORWARD, 3))
                + spinReconstructDir3Plus(shift(adj(u[3]) * spinProjectDir3Plus(psi), BACKWARD, 3))
#endif
#if QDP_ND >= 5
#error "Unsupported number of dimensions"
#endif
    ;
  }
  else
  {
    chi[rb[cb]] = spinReconstructDir0Plus(u[0] * shift(spinProjectDir0Plus(psi), FORWARD, 0))
                + spinReconstructDir0Minus(shift(adj(u[0]) * spinProjectDir0Minus(psi), BACKWARD, 0))
#if QDP_ND >= 2
                + spinReconstructDir1Plus(u[1] * shift(spinProjectDir1Plus(psi), FORWARD, 1))
                + spinReconstructDir1Minus(shift(adj(u[1]) * spinProjectDir1Minus(psi), BACKWARD, 1))
#endif
#if QDP_ND >= 3
                + spinReconstructDir2Plus(u[2] * shift(spinProjectDir2Plus(psi), FORWARD, 2))
                + spinReconstructDir2Minus(shift(adj(u[2]) * spinProjectDir2Minus(psi), BACKWARD, 2))
#endif
#if QDP_ND >= 4
                + spinReconstructDir3Plus(u[3] * shift(spinProjectDir3Plus(psi), FORWARD, 3))
                + spinReconstructDir3Minus(shift(adj(u[3]) * spinProjectDir3Minus(psi), BACKWARD, 3))
#endif
#if QDP_ND >= 5
#error "Unsupported number of dimensions"
#endif
    ;
  }
#else

  // NOTE: the loop is not unrolled - it should be all in a single line for
  // optimal performance
  chi[rb[cb]] = zero;

  // NOTE: temporarily has conversion call of LatticeHalfFermion - will be removed
  for(int mu = 0; mu < Nd; ++mu)
  {
    chi[rb[cb]] += spinReconstruct(LatticeHalfFermionD(u[mu] * shift(spinProject(psi,mu,-isign), FORWARD, mu)),mu,-isign)
      + spinReconstruct(LatticeHalfFermionD(shift(adj(u[mu]) * spinProject(psi,mu,+isign), BACKWARD, mu)),mu,+isign);
  }
#endif

}

}
