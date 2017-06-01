TOPDIR=@top_srcdir@
BUILDDIR=@top_builddir@

AM_CXXFLAGS= -I$(TOPDIR)/include -I$(BUILDDIR)/include @QPHIX_AUX_CFLAGS@ -I@QPHIX_CODEGEN_HOME@/include

lib_LIBRARIES = libqphix_solver.a

libqphix_solver_a_SOURCES = \
	print_utils.cc \
	memory_usage.cc \
	qphix_cli_args.cc

if QPHIX_BUILD_QPX
libqphix_solver_a_SOURCES += bgq_threadbind.cc
else
libqphix_solver_a_SOURCES += generic_threadbind.cc
endif
