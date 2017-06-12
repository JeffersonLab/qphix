CMAKE build system
==================

This is an introduction to the CMake way of building QPhiX. There are essentially two ways:

 i) Use an installed qphix_codegen
 ii) Recursively build your own qphix_codegen

The following flags are needed

 -    -Disa=<ARCH>             Select processor arch: arch = avx, avx2, avx512, sse, scalar (Autoconf: --enable-proc=)
 
 -    -Dclover=TRUE/FALSE            Build Clover Components
 -    -Dtwisted_mass=TRUE/FALSE	     Build Twisted Mass Components
 -    -Dtm_clover=TRUE/FALSE         Build Clover with Twisted Mass
 -    -DQPHIX_CODEGEN=<path>         Path to installed QPhiX-Codegen (if left out, package will try to build one)


 -    -Dhost_cxx=<compiler>          Host Compiler for building QPhiX-Codegen generators
 -    -Dhost_cxxflags=" flags "      Flags to host compiler 
 -    -Drecursive-jN=<number>        Number to pass to make -j in sub-builds (NB this is needed on Linux, but seemingly not on Mac)

 -    -DQDPXX_DIR=path               Path to QDPXX if building with QDP++
 -    -Dqdpjit=TRUE/FALSE            Signal that QDPXX is QDP-JIT or not (pick up appropriate import functions)
 -    -Dqdpalloc=TRUE/FALSE          Signal that QPhiX should use QDP++'s allocator

 -    -Dmm_malloc=TRUE/FALSE         If not using QDP++ allocator, use mm_malloc (if TRUE) or posix_memalign (if FALSE)
 -    -Dtesting=TRUE/FALSE           build tests (timer programs always built regardsless)
 -    -Dparallel_arch=arch           the parallel arch for which we are building. Either scalar or parscalar 
                                     if building with QDP++, inherited from QDP++
                                     if parscalar and not building with QDP++, one must specify QMP_DIR for a QMP location
				     
 -    -Dfake_comms=TRUE/FALSE        if in scalar parallel arch, specify whether we should fake comms
 -    -Dextra_messages=TRUE/FALSE    should we enable #warning style messages
 -    -Dcean=TRUE/FALSE		     enable Cilk++ Extended Array Notation

THe best approach is to build out of source.

Please note:

a)  recusrive building of the code-generator requires python3 be available on the 
command line, and the Jinja2 Python Module be installed for generating header file prototypes.

b)  C++-11 is required for both the codegen and qphix parts.

c)  host_cxx and host_cxx refers to the machine which is carrying out the build. This may be different from the target architecture.
    Use CXX, CXXFLAGS and CMAKE_EXE_LINKER_FLAGS for adding options to the target architecture. For example, curreently when TBB is used
    QDP++ does not signal this in its qdpxx-config.sh file. So one has to add this on manually. 

d)  g++ and clang do not respect the restrict keyword and need an explicit -Drestrict=__restrict__ cause potentially

e) mic (KNC) target has been obsoleted, as has qpx (or at least it is untested)

An example configuration command to build the code-generator as a sub-module is below

CXX=mpiicpc \
CXXFLAGS="-xCORE-AVX2 -std=c++11 -O3 -fopenmp -g"  \
cmake -Disa=avx2  \
      -Dqdpalloc=TRUE -DQDPXX_DIR=/home/bjoo/package-5-17-17/avx/install/qdp++-double \
      -Dtwisted_mass=TRUE \
      -Dtm_clover=TRUE \
      -Dclover=TRUE \
      -Dhost_cxx=g++ \
      -Dhost_cxxflags="-std=c++11 -O3" \
      -Drecursive_jN=24 \
      -Dtesting=TRUE  \
      -DCMAKE_EXE_LINKER_FLAGS="-L${TBBLIBDIR} -ltbbmalloc -ltbb"   \
      -DCMAKE_INSTALL_PREFIX=${HOME}/install/qphix/avx2 ../qphi
    
make 
make install
 
