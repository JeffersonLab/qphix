# QPhiX Dslash and Solver Library: version 1.0.0

Build status:

- `devel` ![](https://api.travis-ci.org/JeffersonLab/qphix.svg?branch=devel)
- `master` ![](https://api.travis-ci.org/JeffersonLab/qphix.svg?branch=master)

## Licensing Copying and Distribution:

Please see the `LICENSE` file for the License. Jefferson Science Associates
Copyright notice and Licens and Intel Corporation `Copyright` notices are also
reproduced in the file `COPYING`.

## Disclaimers

This is research software and while every effort has been made to ensure it
performs correctly it is not certified or warranted to be free of errors, bugs.
Nor is it certified or warranted that it will not crash, hang or in other way
malfunction. We likewise cannot and will not guarantee any specific performance
on any specific system and use of this software is entirely at your own risk,
as is disclaimed in the `COPYRIGHT` and `LICENSE` files

## Getting Building and Installing this library

The library is a C++ library and is mostly a collection of header files. This package also
contains some test programs.

| Program | Description |
| --- | --- |
| `./t_dslash` | Test Wilson Dslash operator (against QDP++) |
| `./t_clov_dslash` | Test the Clover Dslash operators (against implementation outlined from Chroma) |
| `./time_dslash_noqdp` | Time the Dslash operator, Wilson operator, and solvers,  without linkage to QDP++ |
| `./time_clov_noqdp` | Time the "Clover Dslash", Clover operator, and solver without linkage to QDP++ |

The installation and building of test programs is done using the GNU automake
system. Configuration is performed with `autoconf`. The `Makefile.in`-s are not
checked into the repository and need to be regenerated.

### Getting the library

The library can be downloaded from the [GitHub repository](https://github.com/JeffersonLab/qphix.git)

### Dependencies

- The [QDP++ library](http://usqcd-software.github.io/qdpxx/) -- to build the test programs
- Gnu Autoconf and Automake -- to regenerate the build system
- `make` -- to build and install
- Intel C++ Compiler (Code tested with icpc (ICC) 14.0.1 20131008)
- MPI for multi-node builds
- Intel MPI and CML Proxy for multi Xeon Phi running. Please contact Intel Parallel Computing Labs for CML Proxy source code.

NB: The code was tested with ICC 14.0.1 from Intel Composer XE 2013 SP1.1.106. 
Some funnies were observed with the BLAS routines using ICC fom Composer XE 2015

Future compiler support is planned/ongoing. Please see `TODO`.

### Generating the build system

- Obtain the source code distribution, which should have a toplevel directory called qphix
- Regenerate the build system with commands: 

        cd qphix
        autoreconf

### Configuring the package

The package has some components most of which can be selected with switches to configure:

| Flag | Description |
| --- | --- |
| `--prefix=<install_dir>`  | Install into directory <install_dir> |
| `--enable-clover`         | Enable Clover term related codes (default true) |
| `--enable-proc=PROC`      | Select kernels for processor PROC. Allowed values are SCALAR,AVX,MIC for now	 |
| `--with-qmp<qmp_dir>`     | QMP Library is in <qmp_dir> |
| `--with-qdp=<qdpxx_dir>`  | QDP++ library is in <qdpxx_dir> |
| `--enable-soalen=SOA`     | Set the SOA-length to SOA (SOA=4,8,16 as appropriate)  |
| `--enable-cean`           | Enable C++ Extended Array Notation (Cilk++ notation) in some places (default: disabled) |
| `--enable-mm-malloc`      | Use mm_malloc for allocating aligned arrays if disabled posix_memalign will be used. Currently Xeon Phi will use a mixture of mmap and _mm_malloc irrespective of this value (default: enabled) |

In addition, compiler flags can be passed via variables `CXXFLAGS`, `CFLAGS`
and teh compiler can be selected via vriables `CXX` and `CC`.

When compiling for Xeon Phi in native mode or cross compiling in general it is
useful to switch autoconf into cross compile mode. This can be achieved by
giving a different value for `--host` and `--build` e.g. for Xeon Phi:
`--host=x86_64-linux-gnu --build=none-none-none`

E.g: 

For AVX: 

    configure --prefix=<install_location> \
              --with-qdp=<QDP++ installation>
              --enable-proc=AVX \
              --enable-soalen=8 \
              --enable-clover \
              --enable-openmp \
              --enable-cean \
              --enable-mm-malloc \
              CXXFLAGS="-openmp -g -O2 -finline-functions -fno-alias -std=c++11 -xAVX -vec-report -restrict" \
              CFLAGS="-openmp -g  -O2 -fno-alias -std=c99 -xAVX -vec-report -restrict" \
              CXX="mpiicpc" \
              CC="mpiicc"

or for MIC: 

    configure --prefix=<install_location> \
             --with-qdp=<QDP++ installation>
             --enable-proc=MIC \
             --enable-soalen=8 \
             --enable-clover \ 
             --enable-openmp \
             --enable-cean \
             --enable-mm-malloc \
             CXXFLAGS="-openmp -mmic -vec-report -restrict -mGLOB_default_function_attrs=\"use_gather_scatter_hint=off\" -g -O2 -finline-functions -fno-alias -std=c++0x" \
             CFLAGS="-mmic -vec-report -restrict -mGLOB_default_function_attrs=\"use_gather_scatter_hint=off\" -openmp -g  -O2 -fno-alias -std=c9l9
             CXX="mpiicpc" \
             CC="mpiicc" \
             --host=x86_64-linux-gnu --build=none-none-none

In addition it has been found that GCC and Clang don't like the restrict
keyword in C++ source, but allow defining the `__restrict__` extension. So you
may need to add `-Drestrict=__restrict__` to `CXXFLAGS`.

### Building And Installing

Once configured the package can be installed using `make && make install` the
test programs are build in the `<install_dir>/tests/` directory and are
installed in the `<install_dir>/bin/` directory Since the test programs may
check double precision solves as well it is recommended to build with a double
precision build of QDP++.

### Running the test programs

The test programs need the following typical command line parameters:

| Flag | Description |
| --- | --- |
| `-x Lx -y Ly -z Lz -t Lt` | the dimensions of the lattice |
| `-by By -bz Bz` | (`By` × `Bz`) block size parameters: 4×4 for MIC, 8×8 for Sandy Bridge work well |
| `-pxy Pxy -pxyz Pxyz` | padding factors Pxy Pxyz, Can be set to zero, `Pxy=1` may gain a little for particular lattices (e.g 32³×64) |
| `-c Cores` | Number of cores (all the cores in a system): e.g. 59 for Xeon Phi 5110P, 60 for Xeon Phi 7120, 16 for a dual socket 8 core  per socket Xeon  |
| `-sy Sy -sz Sz` | SMT thread configuration (`Sy=1`, `Sz=4` works well on Knights Corner, `Sy=1`, `Sz=2` works well on Xeon) |
| `-minct  Ct` | Minimum number of blocks in time. This should be the number of sockets ideally. |
| `-compress12` | Employ two row compression (12 number compression) |

For timing routines one can also set: `-i <iters>` Number of timing iters to perform.
For multi-node running one can add: `-geom Px Py Pz Pt` to specify a (`Px Py Pz Pt`) virtual node geometry. 
NB: Cores refers to the number of cores per node.

Some tests may occationally support:

`-prec=PRECISION`: precision to work in `h`=half, `f`=single, `d`=double

E.g. a typical clover dslash test on a dual socket 8 core-per-socket Xeon System would go like:

    ./t_clov_dslash -x 32 -y 32 -z 32 -t 64 -by 8 -bz 8 -pxy 1 -pxyz 0 -c 16 -sy 1 -sz 2 -minct 2 -compress12 -geom 1 1 1 1 -i 500

Whereas on a Xeon Phi 7120 it would be

    ./t_clov_dslash -x 32 -y 32 -z 32 -t 64 -by 4 -bz 4 -pxy 1 -pxyz 0 -c 60 -sy 1 -sz 4 -minct 1 -compress12 -geom 1 1 1 1 -i 500

## Integration with Chroma

### Minimum Chroma Revision

The minimum revision for using this library with Chroma is Git Commit ID:  

    e558743151ff30598b3d0e374b2ccb6fe95a5fbc

from the Chroma distribution on GitHub: https://github.com/JeffersonLab/chroma.git

### Configuring Chroma for building with QPhiX 

QPhiX has been integrated with chroma as a LinOpSysSolver. To use the following
configure options neeed to be passed to Chroma
   
| Flag | Description |
| --- | --- |
| `--with-qphix-solver=<QPhiX Install directory>` | specify the location of the installed QPhiX library |
| `--enable-qphix-solver-arch=PROC` | specify the QPhiX processor architecture (avx,mic) |
| `--enable-qphix-solver-soalen=SOA` | the SOA Length for the problem (4,8,16 etc as appropriate) |
| `--enable-qphix-solver-compress12` | add this option if 12 compression is required |
| `--enable-qphix-solver-inner-type=f` | precision of the inner solver for mixed prec: (`h`=half, `f`=single, `d`=double) |
| `--enable-qphix-solver-inner-soalen=4` | SOA length of inner solver |

### XML Drivers

The XML tag group to use the lin op solver in Propagator applications looks like this for non-mixed precision solvers:

    <InvertParam>
       <invType>QPHIX_CLOVER_INVERTER</invType>
       <SolverType>BICGSTAB</SolverType>
       <MaxIter>1000</MaxIter>		<!-- Maximum iterations before giving up -->
       <RsdTarget>1.0e-7</RsdTarget>		<!-- Desired Target relative residuum -->
       <CloverParams>                         <!-- repeat params in the FermionAction, both Mass and Kappa are OK -->
         <Kappa>0.115</Kappa>
         <clovCoeff>2.0171</clovCoeff>
         <clovCoeffR>2.0171</clovCoeffR>
         <clovCoeffT>0.95</clovCoeffT>
         <AnisoParam>
           <anisoP>true</anisoP>
           <t_dir>3</t_dir>
           <xi_0>2.9</xi_0>
           <nu>0.94</nu>
         </AnisoParam>
       </CloverParams>
       <AntiPeriodicT>false</AntiPeriodicT>    <!-- set to true if BCs are antiperiodic -->
       <Verbose>false</Verbose>
       <NCores>16</NCores>                     <!-- number of cores per node -->
       <ThreadsPerCore>2</ThreadsPerCore>      <!-- number of threads per core -->
       <By>8</By>				  <!-- By and Bz block dimensions -->
       <Bz>8</Bz>
       <Sy>1</Sy>                              <!-- Sy and Sz hyperthread grid dimensions -->
       <Sz>2</Sz>
       <PadXY>1</PadXY>                        <!-- PadXY and PadXYZ padding factors -->
       <PadXYZ>0</PadXYZ>
       <MinCt>2</MinCt>                             <!-- Ct at start -- no of temporal blocks. Set equal to no of sockets -->
       <RsdToleranceFactor>5.0</RsdToleranceFactor> <!-- tolerate slack between iterated and true residua, e,g, due to preconditioning etc -->
       <Tune>true</Tune>                            <!-- tune BLAS routines prior to solve -->
     </InvertParam>

For the Iterative Refinement solver the XML is similar:

    <InvertParam>
      <invType>QPHIX_CLOVER_ITER_REFINE_BICGSTAB_INVERTER</invType>
      <SolverType>BICGSTAB</SolverType>  <!-- Inner solver -->
      <MaxIter>1000</MaxIter>	     <!-- Max outer iters -->
      <RsdTarget>1.0e-7</RsdTarget>      <!-- Target relative residuum for outer solver -->
      <Delta>0.1</Delta>		     <!-- Relative residuum drop factor in inner solve. Ie inner solve reduces residuum this much -->
      <CloverParams>		     <!-- Action params same as in FermionAction -->
        <Kappa>0.115</Kappa>
        <clovCoeff>2.0171</clovCoeff>
        <clovCoeffR>2.0171</clovCoeffR>
        <clovCoeffT>0.95</clovCoeffT>
        <AnisoParam>
          <anisoP>true</anisoP>
          <t_dir>3</t_dir>
          <xi_0>2.9</xi_0>
          <nu>0.94</nu>
        </AnisoParam>
      </CloverParams>
      <AntiPeriodicT>false</AntiPeriodicT> <-- set to true for antiperiodic BCs -->
      <Verbose>false</Verbose>
      <NCores>16</NCores>			<!-- Number of cores per node -->
      <ThreadsPerCore>2</ThreadsPerCore>    <!-- Number of threads per core -->
      <By>8</By>				<!-- Blocking dimensions By and Bz -->
      <Bz>8</Bz>
      <Sy>1</Sy>				<!-- SMT Thread dimensions Sy and Sz -->
      <Sz>2</Sz>
      <PadXY>1</PadXY>			<!-- Padding dimensions PadXY and PadXYZ -->
      <PadXYZ>0</PadXYZ>
      <MinCt>1</MinCt>			<!-- Initial Ct = set to number of sockets -->
      <RsdToleranceFactor>5.0</RsdToleranceFactor>   <!-- Slop factor tolerated between solver and true residua, e.g. due to preconditioning -->
      <Tune>true</Tune>			<!-- set to true to tune BLAS -->
    </InvertParam>

### General pointers about the code

 - this is a C++ code and uses templates
 - The code parameters Floating point type, Vector length, SOA length, and a boolean to indicate 12 compression are typically template parameters e.g.

        Geometry<FT,V,S,compress12>
        Dslash<FT,V,S,compress12> etc...

 -  These parameters typically relate to types or array dimensions and hence need to be set at compile time. This is annoying (especially for the compression) but there it is.

 - The `Geometry<FT,V,S,compress>` class defines the node geometry, blocking, padding etc. It also supplies allocation functions for gauges, spinors, and clover terms.
 - The `Comm<FT,V,S,compress>` class defines multi-node communications.

 - Code generated by QPhiX codegen is copied to the directory: 

        qphix/include/qphix/mic/generated  - for MIC
        qphix/include/qphix/avx/generated  - for AVX

The kernels in the generated directories are included with a combination of:

    qphix/include/qphix/<ARCH>/xxxx_complete_specialization_form.h 

and

    qphix/include/qphix/<ARCH>/xxxx_complete_specialization.h 

where `ARCH=mic` or `ARCH=avx` 
and `xxx` is `dslash_ARCH` or `clov_dslash_ARCH` ie: `dslash_avx clov_dslash_avx` 

The dslash site loops and communications are defined in the functions 

- `DyzPlus` Dslash
- `DyzMinus`: Dslash Dagger 

etc in the files:

    qphix/include/dslash_body.h and qphix/include/clover_dslash_body.h

- Blas linear algebra is implemented by combining site loops which loop over the SOA-s in a vector 
   with Functors to be executed on the vector. The loops are defined in `qphix/include/qphix/site_loops.h`
   and the various functors are defined in `qphix/include/real_functors.h` `qphix/include/complex_functors.h` 

- Solvers are currently implemnented in

    - `invcg.h` -- Conjugate Gradients:   solves `M^\dagger M x = b`
    - `invbicgstab.h` -- BICGStab :       solves `M x = b` 
    - `inv_richardson_multiprec.h` -- Solves `M x = b` using iterative refinement (aka defect correction).

Please see the source code in the tests/ directory on how to set these up:

e.g. in  `tests/testDslashFull.cc` 

The Dslash operators themselves provide two main operations:
      
- `y = D x` (or  `y = A^{-1} D x`)  -- called `Dslash`
- `z = a x - b D y` (or  `z = A x - b D y`) -- so called `achimbdpsi` versions -- `a`, `b` are constants, `A` is the clover term)

From this one can easily construct the Schur preconditioned operators: 
e.g. `y = (A - D A^{-1} D) x`
is constructed as  `z = A^{-1}D x`, `y = Ax - bDz` 

Contact Balint Joo, bjoo@jlab.org
