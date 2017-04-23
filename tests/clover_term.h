#ifndef QPHIX_TEST_CLOVER_TERM_H
#define QPHIX_TEST_CLOVER_TERM_H

#include "qphix/qphix_config.h"

#ifdef QPHIX_BUILD_CLOVER
#ifdef QPHIX_BUILD_QDPJIT
#warning using LLVM Clover Term

#include "./clover_term_llvm_w.h"
namespace QPhiX
{
template<typename Phi, typename U>
using CloverTermT = LLVMCloverTermT<Phi,U>;
};
#else
// Using regular cloverTerm
#warning using QDP Clover Term

#include "./clover_term_qdp_w.h"

namespace QPhiX 
{
template<typename Phi,typename U>
using CloverTermT = QDPCloverTermT<Phi,U>;
};

#endif // if else defined QDP-JIT
#endif // if else defined BUILD-CLOVER

#endif
