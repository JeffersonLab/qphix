#ifndef QPHIX_TEST_CLOVER_TERM_H
#define QPHIX_TEST_CLOVER_TERM_H

#include <type_traits>

#include "qphix/diagnostics.h"

#include "qphix/qphix_config.h"

#if defined(QPHIX_BUILD_CLOVER) || defined(QPHIX_BUILD_TWISTED_MASS_WITH_CLOVER)
#ifdef QPHIX_BUILD_QDPJIT
QPHIX_MESSAGE("using LLVM Clover Term")

#include "./clover_term_llvm_w.h"
namespace QPhiX
{
template <typename Phi, typename U>
using CloverTermT = LLVMCloverTermT<Phi, U>;
};
#else
// Using regular cloverTerm
QPHIX_MESSAGE("using QDP Clover Term")

#include "./clover_term_qdp_w.h"

namespace QPhiX
{
template <typename Phi, typename U>
using CloverTermT =
    typename std::enable_if<std::is_same<U, QDP::LatticeColorMatrixF>::value ||
                                std::is_same<U, QDP::LatticeColorMatrixD>::value,
                            QDPCloverTermT<Phi, U>>::type;
};

#endif // if else defined QDP-JIT
#endif // if else defined WILSON/TWISTED-MASS CLOVER

#endif
