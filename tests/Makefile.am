INCFLAGS=-I@top_srcdir@/include -I@top_builddir@/include

AM_CXXFLAGS=$(INCFLAGS) @QPHIX_AUX_CFLAGS@
AM_LDFLAGS=-L@top_builddir@/lib @LDFLAGS@ @QPHIX_AUX_LDFLAGS@
LDADD = -lqphix_solver @LIBS@ @QPHIX_AUX_LIBS@

check_PROGRAMS = time_dslash_noqdp
if QPHIX_BUILD_WITH_QDP
check_PROGRAMS += t_dslash t_minvcg
check_PROGRAMS += test_blas
endif

if QPHIX_BUILD_CLOVER
check_PROGRAMS += time_clov_noqdp
if QPHIX_BUILD_WITH_QDP
check_PROGRAMS +=  t_clov_dslash t_clov_invert_from_file
endif
endif

if QPHIX_BUILD_TWISTED_MASS
if QPHIX_BUILD_WITH_QDP
check_PROGRAMS += t_twm_dslash
endif
check_PROGRAMS += time_twm_noqdp
endif

if QPHIX_BUILD_TWISTED_MASS_WITH_CLOVER
if QPHIX_BUILD_WITH_QDP
check_PROGRAMS += t_twm_clover
endif
check_PROGRAMS += time_twm_clover
endif

if QPHIX_BUILD_ND_TWISTED_MASS
if QPHIX_BUILD_WITH_QDP
check_PROGRAMS += t_ndtm_dslash
endif
check_PROGRAMS += time_ndtm_noqdp
endif

if QPHIX_BUILD_ND_TWISTED_MASS_WITH_CLOVER
if QPHIX_BUILD_WITH_QDP
check_PROGRAMS += t_ndtm_clover t_ndtm_vs_tm_clover t_ndtm_invert_from_file
endif
check_PROGRAMS += time_ndtm_clover
endif

bin_PROGRAMS = $(check_PROGRAMS)

FRAMEWORK_SRCS =  unittest.h \
	dslashm_w.h \
	dslashm_w.cc \
	reunit.h \
	invbicgstab_test.h \
	reunit.cc

###############################################################################
#                                Plain Wilson                                 #
###############################################################################

t_dslash_SOURCES = $(FRAMEWORK_SRCS) \
	testDslashFull.h \
	testDslashFull.cc \
	t_dslash.cc



time_dslash_noqdp_SOURCES= \
	time_dslash_noqdp.cc \
	timeDslashNoQDP.h \
	timeDslashNoQDP.cc

t_minvcg_SOURCES = $(FRAMEWORK_SRCS) \
	testMInvCG.h \
	testMInvCG.cc \
	t_minvcg.cc

t_dslash_DEPENDENCIES = build_lib
time_dslash_noqdp_DEPENDENCIES=build_lib
t_minvcg_DEPENDENCIES=build_lib

if QPHIX_BUILD_CLOVER

if QPHIX_BUILD_QDPJIT_CLOVER
CLOVER_SRCS=clover_term_llvm_w.h
endif

if QPHIX_BUILD_QDP_CLOVER
CLOVER_SRCS=clover_term_qdp_w.h
endif

t_clov_dslash_SOURCES = $(FRAMEWORK_SRCS) \
	testClovDslashFull.h \
	testClovDslashFull.cc \
	t_clov_dslash.cc \
	mesfield.h \
	mesfield.cc \
	$(CLOVER_SRCS) \
	clover_fermact_params_w.h

t_clov_invert_from_file_SOURCES = $(FRAMEWORK_SRCS) \
	testClovInvertFromFile.h \
	testClovInvertFromFile.cc \
	t_clov_invert_from_file.cc \
	mesfield.h \
	mesfield.cc \
	$(CLOVER_SRCS) \
	clover_fermact_params_w.h

time_clov_noqdp_SOURCES= \
	time_clov_noqdp.cc \
	timeClovNoQDP.h \
	timeClovNoQDP.cc

t_clov_dslash_DEPENDENCIES = build_lib
t_clov_invert_from_file_DEPENDENCIES = build_lib
time_clov_noqdp_DEPENDENCIES = build_lib

endif

###############################################################################
#                                Twisted Mass                                 #
###############################################################################

if QPHIX_BUILD_TWISTED_MASS
t_twm_dslash_SOURCES = $(FRAMEWORK_SRCS) \
	testTWMDslashFull.h \
	testTWMDslashFull.cc \
	t_twm_dslash.cc

time_twm_noqdp_SOURCES=$(FRAMEWORK_SRCS) \
	time_twm_noqdp.cc \
	timeTWMNoQDP.h \
	timeTWMNoQDP.cc

t_twm_dslash_DEPENDENCIES = build_lib
time_twm_noqdp_DEPENDENCIES = build_lib
endif


if QPHIX_BUILD_TWISTED_MASS_WITH_CLOVER

if QPHIX_BUILD_QDPJIT_CLOVER
TM_CLOVER_SRCS = clover_term_llvm_w.h
endif

if QPHIX_BUILD_QDP_CLOVER
TM_CLOVER_SRCS = clover_term_qdp_w.h
endif

t_twm_clover_SOURCES = $(FRAMEWORK_SRCS) \
	testTWMCloverFull.h \
	testTWMCloverFull.cc \
	t_twm_clover.cc \
	mesfield.h \
	mesfield.cc \
	$(TM_CLOVER_SRCS) \
	clover_fermact_params_w.h

time_twm_clover_SOURCES=$(FRAMEWORK_SRCS) \
	time_twm_clover.cc \
	timeTWMClover.h \
	timeTWMClover.cc

t_twm_clover_DEPENDENCIES = build_lib
time_twm_clover_DEPENDENCIES = build_lib
endif

###############################################################################
#                         Non-Degenerate Twisted Mass                         #
###############################################################################

if QPHIX_BUILD_ND_TWISTED_MASS
t_ndtm_dslash_SOURCES = $(FRAMEWORK_SRCS) \
	testNDTMDslashFull.h \
	testNDTMDslashFull.cc \
	t_ndtm_dslash.cc

time_ndtm_noqdp_SOURCES=$(FRAMEWORK_SRCS) \
	time_ndtm_noqdp.cc \
	timeNDTMNoQDP.h \
	timeNDTMNoQDP.cc

t_ndtm_dslash_DEPENDENCIES = build_lib
time_ndtm_noqdp_DEPENDENCIES = build_lib
endif


if QPHIX_BUILD_ND_TWISTED_MASS_WITH_CLOVER

if QPHIX_BUILD_QDPJIT_CLOVER
TM_CLOVER_SRCS = clover_term_llvm_w.h
endif

if QPHIX_BUILD_QDP_CLOVER
TM_CLOVER_SRCS = clover_term_qdp_w.h
endif

t_ndtm_clover_SOURCES = $(FRAMEWORK_SRCS) \
	testNDTMCloverFull.h \
	testNDTMCloverFull.cc \
	t_ndtm_clover.cc \
	mesfield.h \
	mesfield.cc \
	$(TM_CLOVER_SRCS) \
	clover_fermact_params_w.h

t_ndtm_vs_tm_clover_SOURCES = $(FRAMEWORK_SRCS) \
	testNDTMvsTMCloverFull.h \
	testNDTMvsTMCloverFull.cc \
	t_ndtm_vs_tm_clover.cc \
	mesfield.h \
	mesfield.cc \
	$(TM_CLOVER_SRCS) \
	clover_fermact_params_w.h \
	packed_spinor.h

t_ndtm_invert_from_file_SOURCES = $(FRAMEWORK_SRCS) \
	testNDTMInvertFromFile.h \
	testNDTMInvertFromFile.cc \
	t_ndtm_invert_from_file.cc \
	mesfield.h \
	mesfield.cc \
	$(TM_CLOVER_SRCS) \
	clover_fermact_params_w.h \
	packed_spinor.h

time_ndtm_clover_SOURCES=$(FRAMEWORK_SRCS) \
	time_ndtm_clover.cc \
	timeNDTMClover.h \
	timeNDTMClover.cc

t_ndtm_clover_DEPENDENCIES = build_lib
t_ndtm_vs_tm_clover_DEPENDENCIES = build_lib
t_ndtm_invert_from_file_DEPENDENCIES = build_lib
time_ndtm_clover_DEPENDENCIES = build_lib
endif

###############################################################################

test_blas_SOURCES=$(FRAMEWORK_SRCS) \
	test_blas.cc \
	testBlas.h \
	testBlas.cc

test_blas_DEPENDENCIES= build_lib

DEPENDENCIES = build_lib
${bin_PROGRAMS}: ${DEPENDENCIES}
${check_PROGRAMS}: ${DEPENDENCIES}
${EXTRA_PROGRAMS}: ${DEPENDENCIES}

.PHONY: build_lib
	cd @top_builddir@/lib; $(MAKE)
