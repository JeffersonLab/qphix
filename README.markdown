# QPhiX Code Generator

## Licensing Copying and Distribution:

Please see the `LICENSE` file for the License.

## Disclaimers

This is research software and while every effort has been made to ensure it
performs correctly it is not certified or warranted to be free of errors, bugs.
Nor is it certified or warranted that it will not crash, hang or in other way
malfunction. We likewise cannot and will not guarantee any specific performance
on any specific system and use of this software is entirely at your own risk,
as is disclaimed in the `COPYRIGHT` and `LICENSE` files. Further, some targets
are experimental and incomplete. You have been warned.

## Building the Code

Available targets are:

| Target | Description | Example Products |
| --- | --- | --- |
| `mic` | Intel MIC | Intel Xeon Phi Knight's Corner) |
| `avx` | Regular AVX | Intel Xeon Sandy Bridge, Ivy Bridge |
| `avx2` | AVX2 -- untested | Intel Xeon Haswell, Broadwell |
| `avx512` | AVX512 | Intel Xeon Skylake, Kaby Lake; Intel Xeon Phi Knight's Landing |
| `sse`  | For SSE -- untested | Intel Pentium 4, Intel Core, Intel Xeon Nehalem |
| `scalar` | C, scalar code (`SOALEN=1`) | |

Run the code with `make TARGET` this will build the code generator
and generate the code for the target.  The target specific code will be placed
in the subdirectory corresponding to the target (i.e. MIC code goes to `mic/`,
AVX code goes to `avx/` etc).

To build for Intel Sandy Bridge and Intel Ivy Bridge processors, use the
following:

    make avx

The generated files will be dropped in the target directory, i.e.

    ./mic
    ./avx 
 
etc. which will be created in the root directory if they don't yet exist.

In case you forget the `

## How is the target configured?

For each target there is a separate `customMake.<target>` which sets Makefile
variables for the target. In addition some of these variables can be
over-ridden on the command line to `make` (such as `PRECISION`, `SOALEN`, `AVX2`,
`ENABLE_LOW_PRECISION`) as is done in the Makefile

Some env variables which affect the code generation are: 

| Variable | Purpose |
| --- | --- |
| `PREF_L1_SPINOR_IN` | Generate L1 prefetches for input spinors |
| `PREF_L2_SPINOR_IN` | Generate L2 prefetches for input spinors |
| `PREF_L1_SPINOR_OUT` | Generate L1 prefetches for output spinors |
| `PREF_L2_SPINOR_OUT` | Generate L2 prefetches for output spinors |
| `PREF_L1_GAUGE` | Generate L1 prefetches for input gauge fields |
| `PREF_L2_GAUGE` | Generate L2 prefetches for input gauge fields |
| `PREF_L1_CLOVER` | Generate L1 prefetches for input clover terms |
| `PREF_L2_CLOVER` | Generate L2 prefetches for in put clover terms |
| `USE_LDUNPK` | Use Load/Unpack as opposed to gather |
| `USE_PKST` | Use pack-store instead of scatter |
| `USE_SHUFFLES` | Use loads & shuffles to transpose spinor when `SOALEN > 4` |
| `NO_GPREF_L1` | Generate a bunch of normal prefetches instead of gather prefetches for L1  |
| `NO_GPREF_L2` | Generate a bunch of normal prefetches instead of gather prefetches for L2  |
| `ENABLE_STREAMING_STORES` | Enable streaming stores for output spinor |
| `USE_PACKED_GAUGES` | Use 2D xy packing for Gauges -- enabled |
| `USE_PACKED_CLOVER` |  Use 2D xy packing for Clover -- enabled |
| `AVX2` | Generate AVX2 code (experimental)
| `AVX512` | together with `mode=mic`: Generate AVX512 code for Knight's Landing (untested)

## How is the code organized?

There is two parts to the code: 

1. The code generation framework. This works with addresses, registers and
   instructions. Addresses are defined in `address_types.h`. Instructions are
   defined in `instructions.h`.

    NB: Each instruction is an object, belonging to a class. 
    There are two kinds of instructions: those that have memory references and those that do not.
    Each instruction has a `serialize()` method which will print the code (potentially an intrinsic) 
    for that code.

    Then there are simple functions, which usually look like e.g.

    ```cpp
    void loadFVec(InstVector& ivector, const FVec& ret, const Address *s, etc...)  {}
    ```

    These functions typically create an instruction object (in this case a load
    instruction), or several instruction objects, potentially making use of
    input references to vectors and vectors. The created objects are inserted
    into the `InstVector`, which is just a `std::vector` of `Instruction*`s.

    By using functions like `loadFVec`, `gatherFVec`, `scatterFVec` etc, one
    can build linear streams of instructions. 

    The concrete bodies of the instructions, are implemented in the source files:

    | File | Purpose |
    | --- | --- |
    | `inst_scalar.cc` | Scalar C-code (`VECLEN=1`)
    | `inst_dp_vec2.cc` | SSE Double Precision (`VECLEN=2`)
    | `inst_dp_vec4.cc` | AVX Double Precision (`VECLEN=4`)
    | `inst_dp_vec8.cc` | MIC (KNC/AVX512) Double Precision (`VECLEN=8`)
    | `inst_sp_vec4.cc` | SSE Single Precision (`VECLEN=4`)
    | `inst_sp_vec8.cc` | AVX Single Precision (`VECLEN=8`)
    | `inst_sp_vec16.cc` | MIC (KNC/AVX512) Single Precision  (`VECLEN=16`)

2. The generators for the Kernels: 
  
    | File | Purpose |
    | --- | --- |
    | `data_types.h` | defines data_types for Gauge fields, Spinors, and CloverTerm parts along with associated load, store, prefetch, stream etc functions (which use the `Instruction` objects) |
    | `dslash.h` |  defines the code generator functions used to construct a dslash (projections, etc) |
    | `dslash_common.h` | defines the main generator routines for things like projections, mat mults, face packs etc. |
    | `dslash.cc` | implements the main kernels using `dslash_common.h` |
    | `codegen.cc` | is the driver, that just calls `generate_code()` from `dslash.cc` |

The various code-generation options in the `inst_xx_vecX.cc` files, the headers
etc are controlled by compiler `#define`s defined in the `Makefile`s, and
`customMakef.target` files.
