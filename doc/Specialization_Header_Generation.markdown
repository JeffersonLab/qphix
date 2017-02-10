# Specialization Header Generation

If you do not want to read about the history nor the details, you can just skip
to the sections about adding a new architecture or kernel family below.

## Description of the needed code {#desc}

The code generator currently supports seven instruction set architectures (ISA)
(like AVX, AVX2, QPX) as well as four families of kernels (like Dslash, Clover
Dslash, Twisted Mass Dslash). Interface code for every family of kernels and
every architecture has to be written.

For each kernel function (we'll use `dslash_plus_vec` in this example), there
needs to be a general template like the following:

```{.cpp}
template <typename FT, int veclen, int soalen, bool compress12>
inline void dslash_plus_vec(
    const typename Geometry<FT, veclen, soalen, compress12>::FourSpinorBlock
        *xyBase,
    const typename Geometry<FT, veclen, soalen, compress12>::FourSpinorBlock
        *zbBase,
    const typename Geometry<FT, veclen, soalen, compress12>::FourSpinorBlock
        *zfBase,
    const typename Geometry<FT, veclen, soalen, compress12>::FourSpinorBlock
        *tbBase,
    const typename Geometry<FT, veclen, soalen, compress12>::FourSpinorBlock
        *tfBase,
    typename Geometry<FT, veclen, soalen, compress12>::FourSpinorBlock *oBase,
    const typename Geometry<FT, veclen, soalen, compress12>::SU3MatrixBlock
        *gBase,
    const int xbOffs[veclen],
    const int xfOffs[veclen],
    const int ybOffs[veclen],
    const int yfOffs[veclen],
    const int offs[veclen],
    const int gOffs[veclen],
    const int siprefdist1,
    const int siprefdist2,
    const int siprefdist3,
    const int siprefdist4,
    const int gprefdist,
    const int pfyOffs[veclen],
    const typename Geometry<FT, veclen, soalen, compress12>::FourSpinorBlock
        *pfBase2,
    const typename Geometry<FT, veclen, soalen, compress12>::FourSpinorBlock
        *pfBase3,
    const typename Geometry<FT, veclen, soalen, compress12>::FourSpinorBlock
        *pfBase4,
    const unsigned int accumulate[8],
    const FT coeff_s,
    const FT coeff_t_f,
    const FT coeff_t_b) {
    // BASE CASE TEMPLATE. Do nothing for now. Define this in
    // dslash_generated_c.h later
    fprintf(stderr, "Generic veclen and soalen not supported yet.\n");
    abort();
}
```

Then also for each `FPTYPE`, `VEC`, `SOA` and `COMPRESS12`, there needs to be a
template specialization what includes the appropriate kernel

```{.cpp}
template <>
inline void dslash_plus_vec<FPTYPE, VEC, SOA, COMPRESS12>(
    const Geometry<FPTYPE, VEC, SOA, COMPRESS12>::FourSpinorBlock *xyBase,
    const Geometry<FPTYPE, VEC, SOA, COMPRESS12>::FourSpinorBlock *zbBase,
    const Geometry<FPTYPE, VEC, SOA, COMPRESS12>::FourSpinorBlock *zfBase,
    const Geometry<FPTYPE, VEC, SOA, COMPRESS12>::FourSpinorBlock *tbBase,
    const Geometry<FPTYPE, VEC, SOA, COMPRESS12>::FourSpinorBlock *tfBase,
    Geometry<FPTYPE, VEC, SOA, COMPRESS12>::FourSpinorBlock *oBase,
    const Geometry<FPTYPE, VEC, SOA, COMPRESS12>::SU3MatrixBlock *gBase,
    const int xbOffs[VEC],
    const int xfOffs[VEC],
    const int ybOffs[VEC],
    const int yfOffs[VEC],
    const int offs[VEC],
    const int gOffs[VEC],
    const int siprefdist1,
    const int siprefdist2,
    const int siprefdist3,
    const int siprefdist4,
    const int gprefdist,
    const int pfyOffs[VEC],
    const Geometry<FPTYPE, VEC, SOA, COMPRESS12>::FourSpinorBlock *pfBase2,
    const Geometry<FPTYPE, VEC, SOA, COMPRESS12>::FourSpinorBlock *pfBase3,
    const Geometry<FPTYPE, VEC, SOA, COMPRESS12>::FourSpinorBlock *pfBase4,
    const unsigned int accumulate[8],
    const FPTYPE coeff_s,
    const FPTYPE coeff_t_f,
    const FPTYPE coeff_t_b) {
#include INCLUDE_FILE_VAR(qphix/avx512/generated/dslash_plus_body_, FPTYPE, VEC, SOA, COMPRESS_SUFFIX)
}
```

`INCLUDE_FILE_VAR` is a preprocessor macro that builds the filename of the
generated kernel using the `##` string concatenation. 

This code has to be replicated for every single kernel function in a family of
kernels, there is no way around it. The C++ compiler that translates QPhiX will
have to choose the appropriate template specialization depending on the
`./configure` options chosen by the user of QPhiX.

In order to link everything together, the above template specialization is
included with all possible combination of the template parameters into the
QPhiX code, like this (excerpt):

```{.cpp}
#define FPTYPE float
#define VEC 16

#define COMPRESS12 false
#define COMPRESS_SUFFIX _18

#define SOA 4
#include "qphix/avx512/dslash_avx512_complete_specialization_form.h"
#undef SOA

#define SOA 8
#include "qphix/avx512/dslash_avx512_complete_specialization_form.h"
#undef SOA
```

These files are usually hundreds of lines long and contain only few little
information, namely the veclengths and soalength sensible for the ISA at hand.
There has to be such a file for every kernel family and every architecture,
only with minimal changes and a lot of duplication.

## Automatic generation

The specialization generation code has been written by Martin Ueding in order
to generate the “glue code” for twisted mass on AVX2 and later on introduce a
new family of kernels. It uses the [Jinja 2](http://jinja.pocoo.org/) template
engine for Python. This is commonly used to generate HTML but can also generate
C++ code. It is way more powerful than the C++ preprocessor which had been used
before to generate the code. The most important feature are loops. Also it
allows injection of additional variables from a Python script. All blocks with
`{{ ... }}` and `{% ... %}` will be interpreted by the Jinja template library.

### Function declarations and definitions

The large chunks of C++ code that contain the function declarations are now
template files. For the basic Wilson Dslash, this file is
`jinja/dslash_general.h.j2`. The `.j2` file extension shows that it is a
Jinja 2 template, although that particular file does not contain any template
code currently. The file `jinja/dslash_specialization.h.j2` contains the bodies
with the include command for the kernel body. Instead of using the C++
preprocessor `INCLUDE_FILE_VAR` macro, it uses the Jinja macro
`include_generated_kernel` like so:

```{.jinja}
{{ include_generated_kernel(ISA, kernel, "achimbdpsi_plus_body", FPTYPE, VEC, SOA, COMPRESS12) }}
```

The macro emits a `#include` but can build the filename depending on the
current ISA and the kernel, making it a bit more intelligent. Also it can
figure out the suffix for the gauge compression automatically. The definition
of this macro is the following:

```{.jinja}
{% macro include_generated_kernel(isa, kernel, operator, fptype, vec, soa, compress12) %}

{% if compress12 == '' %}
{% set suffix = '' %}
{% else %}
{% set suffix = '_12' if compress12 == 'true' else '_18' %}
{% endif %}

{% set filename = 'qphix/%s/generated/%s_%s_%s_%s_v%d_s%d%s'|format(isa, kernel, operator, fptype, fptype, vec, soa, suffix) %}
#include "{{ filename }}"
{% endmacro %}
```

This allows to generate the template specialization files for every
architecture automatically, therefore consolidating all the code. Instead of
having seven (or more) copies, there is only now now.

### Specialization file

The main advantage is the generation of the explicit instantiations of the
template parameters. These files are the same for each ISA, except that the
name of the kernel family is exchanged. Adding a kernel family had previously
meant to duplicate two files for each ISA and make some trivial changes. Now
there is the following code that emits the needed C++ code:

```{.jinja}
{% for FPTYPE, VEC, soalens in defines %}
{% set kernel = kernel_pattern|format(fptype=FPTYPE) %}
{% for SOA in soalens %}
{% for COMPRESS12 in ['true', 'false'] %}
#define FPTYPE {{ FPTYPE }}
#define VEC {{ VEC }}
#define SOA {{ SOA }}
#define COMPRESS12 {{ COMPRESS12 }}
{% include 'jinja/%s_specialization.h.j2'|format(kernel_base) %}
#undef FPTYPE
#undef VEC
#undef SOA
#undef COMPRESS12
{% endfor %}
{% endfor %}
{% endfor %}
```

What it does is to iterate through all the floating point types (`FPTYPE`), the
associated vector length for that platform. Then it iterates through the
desired SoA length and also over the compression choices. For each it generates
a block with appropriate `#define` statements.

In principle it would be possible to replace all uses of the C++ preprocessor
with Jinja code. However, this would make locating source lines a bit harder.
Also the generated files would be in the order of 20 MB, which is a burden on
the version control system and on text editors. Compilation time seems to be
largely independent of that and needs around 8 minutes on JURECA either way.

### Automake file

Each architecture has its `Makefile.am` which specifies which files are needed
in the compilation. Such a file is generated for each ISA.

### ISA include file

For each family of kernels, there is a `KERNEL_generated.h` file in
`include/qphix/` of the QPhiX repository. These files are also automatically
generated. These files contain blocks that include the appropriate files for
the chosen ISA:

```{.cpp}
#if defined (QPHIX_AVX2_SOURCE)
#warning "including qphix/avx2/dslash_avx2_complete_specialization.h"
#include "qphix/avx2/dslash_avx2_complete_specialization.h"
#endif
```

## Adding a new architecture

The main burden when adding a new architecture is to write the appropriate
intrinsics in the main code generation framework. This is not covered in this
article and it is assumed that this has already been done.

In the file `jinja/isa.js`, which is a simple [JSON](http://json.org/) file,
another block has to bee added. For AVX, the block looks like this:

```{.js}
"avx": {
    "fptypes": {
        "double": {"veclen": 4, "soalens": [2, 4]},
        "float": {"veclen": 8, "soalens": [4, 8]}
    },
    "extra_includes_global": ["immintrin.h"],
    "extra_includes_local": ["qphix/avx/avx_utils.h"]
},
```

It just specifies the parameters used for the different floating point types as
well as auxiliary files that might need to be included.

When you are done, you can run

    python3 generate_files.py NEW_ISA

where `NEW_ISA` is the ISA you have added.

## Adding a new kernel family

When a new family of kernels is implemented, the kernel has to be added to the
list of kernels in `jinja/kernels.js`. This list looks like this:

```{.js}
[
    "clov_%(fptype_underscore)sdslash",
    "dslash",
    "tmf_clov_%(fptype_underscore)sdslash",
    "tmf_dslash"
]
```

The syntax with `%(fptype_underscore)s` must be used if the filenames of the
generated kernels contain the name of the floating point type. This part of the
string will be expanded to `double_` for `double`. This is needed in the clover
operators because the clover term itself has an independent floating point
type.

The two files `KERNEL_general.h.j2` and `KERNEL_specialization.h.j2` need to be
written. The changes there depend on the new implementation.

## Copying files to QPhiX

After the new code has been generated, it needs to be copied to the QPhiX
repository and be checked in there. The following bunches of files need to be
transferred, where `ISA` stands for some ISA:

| Destination in `qphix-codegen` | Target in `qphix` |
| --- | --- |
| `ISA/*` | `include/qphix/ISA/generated/` |
| `jinja/*_generated.h` | `include/qphix/` |
| `jinja/ISA/*` | `include/qphix/ISA/` |

## Closing remark

One says that every problem in computer science can be solved with another
indirection. Except the problem of too many indirections. Hopefully this
additional layer will make it easier to work with this code generator and
extend it.

The layers within the QPhiX Code Generator are:

- C++ Code generator generates C++ code
- Python program generates “glue code” and Automake files

And then in QPhiX the steps are:

- Automake generates the `Makefile.in` that is then used in the compilation
- C++ preprocessor know about the chosen parameters (like `FT` and `VECLEN`)
  and only leaves the parts of the glue code that is relevant.
- Then the generated kernels are included by the preprocessor and everything
  gets compiled.

<!-- vim: set spell tw=79 :-->
