ACLOCAL_AMFLAGS = -I m4

# gtest source files that we don't compile directly.  They are
# #included by gtest-all.cc.
GTEST_SRC = \
  src/gtest-death-test.cc \
  src/gtest-filepath.cc \
  src/gtest-internal-inl.h \
  src/gtest-port.cc \
  src/gtest-printers.cc \
  src/gtest-test-part.cc \
  src/gtest-typed-test.cc \
  src/gtest.cc

# Distribute and install M4 macro
#m4datadir = $(datadir)/aclocal
#m4data_DATA = m4/gtest.m4
#EXTRA_DIST += $(m4data_DATA)

# We define the global AM_CPPFLAGS as everything we compile includes from these
# directories.
AM_CPPFLAGS = -I$(srcdir) -I$(srcdir)/include

# Build rules for libraries.
lib_LIBRARIES = lib/libgtest.a lib/libgtest_main.a

lib_libgtest_a_SOURCES = src/gtest-all.cc

lib_libgtest_main_a_SOURCES = src/gtest_main.cc
lib_libgtest_main_a_LIBADD = lib/libgtest.a
