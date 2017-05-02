SUBDIRS = gtest

INCFLAGS=-I@top_srcdir@/include -I@top_builddir@/include
AM_CXXFLAGS=$(INCFLAGS) @QPHIX_AUX_CFLAGS@
AM_LDFLAGS=-L@top_builddir@/lib @LDFLAGS@ @QPHIX_AUX_LDFLAGS@ -L@top_builddir@/tests-gtest/gtest/lib

check_PROGRAMS = tests-g
bin_PROGRAMS = tests-g

tests_g_SOURCES = main.cc \
		  clover_product.cc \
		  ../tests/reunit.cc \
		  ../tests/mesfield.cc
tests_g_CXXFLAGS = -fsanitize=address -g $(AM_CXXFLAGS)

tests_g_LDADD = -lgtest -lqphix_solver @LIBS@ @QPHIX_AUX_LIBS@
