TOPDIR=@top_srcdir@
BUILDDIR=@top_builddir@

AM_CXXFLAGS= -I$(TOPDIR)/include -I$(BUILDDIR)/include @QPHIX_AUX_CFLAGS@

lib_LIBRARIES = libqphix_solver.a

libqphix_solver_a_SOURCES = \
	print_utils.cc \
	qphix_cli_args.cc

if QPHIX_BUILD_QPX
libqphix_solver_a_SOURCES += bgq_threadbind.cc
else
libqphix_solver_a_SOURCES += generic_threadbind.cc
endif


###############################################################################
#                      Explicit Template Instantiations                       #
###############################################################################

libqphix_solver_a_SOURCES += Dslash_instantiation.cpp

if QPHIX_BUILD_CLOVER
libqphix_solver_a_SOURCES += ClovDslash_instantiation.cpp
endif

if QPHIX_BUILD_TWISTED_MASS
libqphix_solver_a_SOURCES += TMDslash_instantiation.cpp
endif

if QPHIX_BUILD_TWISTED_MASS_WITH_CLOVER
libqphix_solver_a_SOURCES += TMClovDslash_instantiation.cpp
endif
