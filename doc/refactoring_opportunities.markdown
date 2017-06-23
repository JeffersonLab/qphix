# Refactoring Opportunities

(Author: Martin Ueding <dev@martin-ueding.de>)

During my work on the non-degenerate twisted mass case, I have discovered
opportunities to refactor the code and make it simpler or easier to maintain. I
do not want to change code which I did not write. Instead, this is a list of
things that might be worth looking at if one wants to refactor the whole
codebase.

## Duplicated plus and minus functions

There are pairs of functions for the plus and minus case. This sign is the
`isign` which states whether the Dslash or its hermitian conjugate is to be
applied. These pairs look like this:

```{.cpp}
void NDTMDyzPlus(int tid,
                 const FourSpinorBlock *psi,
                 FourSpinorBlock *res,
                 const SU3MatrixBlock *u,
                 double mu,
                 double mu_inv,
                 int cb);

void NDTMDyzMinus(int tid,
                  const FourSpinorBlock *psi,
                  FourSpinorBlock *res,
                  const SU3MatrixBlock *u,
                  double mu,
                  double mu_inv,
                  int cb);
```

Looking through both implementations with `diff`, I found that they only differ
on a single line, namely in the kernel they call. In the NDTM case, the caller
looked like this:

```{.cpp}
template <typename FT, int veclen, int soalen, bool compress12>
void TMClovDslash<FT, veclen, soalen, compress12>::dslash(
    FourSpinorBlock *res,
    const FourSpinorBlock *psi,
    const SU3MatrixBlock *u,
    const FullCloverBlock *clov[2],
    int isign,
    int cb)
{
    if (isign == 1) {
        DPsiPlus(u, clov[0], psi, res, cb);
    }

    if (isign == -1) {
        DPsiMinus(u, clov[1], psi, res, cb);
    }
}
```

Now it just passes the `isign` parameter onward and chooses the right part of
the clover term:

```{.cpp}
template <typename FT, int veclen, int soalen, bool compress12>
void NDTMClovDslash<FT, veclen, soalen, compress12>::dslash(
    FourSpinorBlock *res,
    const FourSpinorBlock *psi,
    const SU3MatrixBlock *u,
    const FullCloverBlock *invclov[2],
    int isign,
    int cb)
{
    DPsi(u, invclov[isign == 1 ? 0 : 1], psi, res, cb, isign);
}
```

The functions `DPsiPlus` and `DPsiMinus` have been unified to a single `DPsi`
function where the kernel to use is chosen depending on `isign`:

```{.cpp}
auto kernel =
    (isign == 1 ? ndtm_clov_dslash_plus_vec<FT, veclen, soalen, compress12>
                : ndtm_clov_dslash_minus_vec<FT, veclen, soalen, compress12>);

kernel(xyBase + X, zbBase + X, ..., forw_t_coeff_T, back_t_coeff_T);
```

This way several thousand lines of replicated code can be removed. There should
not really be any change in the performance because the branching has to be
done at some point anyway. The user will call the `dslash` function which will
call all the depending functions.

---

Further, there seems to be a large overlap between the following:

- Plus and Minus versions of the same function
- Dslash and achimbdpsi operators: They use different factors (anisotropy vs.
  beta), but it seems that those are just different names for pretty much the
  same thing. A lot of duplicated code could be saved there.
- There is probably a lot of duplicated code between the non-clover and clover
  operators.
- Since every new operator seems to be a copy of the previous code, a lot of
  common elements (offset calculations, ...) could be refactored into common
  functions.

## Tolerance type traits

The test cases use type traits for the tolerances, like this:

```{.cpp}
template <typename T>
struct tolerance {
    static const Double small; // Always fail
};

template <>
const Double tolerance<half>::small = Double(5.0e-3);

template <>
const Double tolerance<float>::small = Double(1.0e-6);
```

This was duplicated in all the test cases. Now there is a header file
`tolerance_type_traits.h` which just has those type traits and can be included.

## Standardized names for the common template parameters

The four common template parameters, (1) floating point type, (2) vector
length, (3), SoA length, and (4) gauge compression, are named differently for
different classes and functions. Examples are:

```{.cpp}
template <typename FT, int V, int soalen, bool compress>
template <typename FT, int V, int S, bool compress>
template <typename FT, int veclen, int soalen, bool compress>
template <typename FT, int veclen, int soalen, bool compress12>
template <typename T, int V, int S, const bool compressP>
```

The names in the kernel specializations probably need to be different because
there they are set by the preprocessor. In all other code, those names could be
unified to make it easier to do further refactorings.

<!-- vim: set spell tw=79 :-->
