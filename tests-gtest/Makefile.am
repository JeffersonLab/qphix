SUBDIRS = gtest

INCFLAGS=-I@top_srcdir@/include -I@top_builddir@/include -I@top_srcdir@/tests-gtest/gtest/include
AM_CXXFLAGS=$(INCFLAGS) @QPHIX_AUX_CFLAGS@
AM_LDFLAGS=-L@top_builddir@/lib @LDFLAGS@ @QPHIX_AUX_LDFLAGS@ -L@top_builddir@/tests-gtest/gtest/lib

check_PROGRAMS = tests-g

if QPHIX_BUILD_TESTS
bin_PROGRAMS = tests-g
endif

tests_g_SOURCES = main.cc \
		  clover_product.cc \
		  ../tests/reunit.cc \
		  ../tests/mesfield.cc
tests_g_CXXFLAGS = $(AM_CXXFLAGS)

tests_g_LDADD = -lqphix_solver -lqphix_codegen @LIBS@ @QPHIX_AUX_LIBS@ -lgtest
