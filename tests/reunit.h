#ifndef REUNIT_H
#define REUNIT_H

#ifndef QDP_INCLUDE
#include "qdp.h"
#endif

using namespace QDP;

// Reunitarize a Lattice Color Matrix
void reunit(LatticeColorMatrixF &a);
void reunit(LatticeColorMatrixD &a);

#endif
