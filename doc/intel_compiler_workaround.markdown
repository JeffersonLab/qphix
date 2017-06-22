# Intel C++ Compiler Workaround                {#intel-cpp-compiler-workaround}

The Intel C++ 17 compiler seems to have a bug in a `const`-conversion. Namely
it does not convert `T **` to `T const *const *`, whereas GCC and LLVM do that.
Also the standard says that this conversion should be done. This article is a
summary of the [Stack Overflow post][SO] where this problem has been
triangulated.

[SO]: http://stackoverflow.com/q/43066405/653152

## Existing One-Flavor Code

In QPhiX, there are BLAS routines that take one or more spinors and work with
them (e.g. scalar product or AXPY). The parameters that are read-only are taken
as “pointer to `const`”, like this:

```{.cpp}
typedef float Spinor[3][4][2][8];

void blas_routine(Spinor *out, const Spinor *in) {
    std::cout << out << " " << in << "\n";
}
```

This has worked just fine, also with the Intel C++ 17 compiler. One full real
example is the `copySpinor` function:

```{.cpp}
template <typename FT, int V, int S, bool compress>
void copySpinor(
    typename Geometry<FT, V, S, compress>::FourSpinorBlock *res,
    const typename Geometry<FT, V, S, compress>::FourSpinorBlock *src,
    const Geometry<FT, V, S, compress> &geom,
    int n_blas_simt)
{
    CopyFunctor<FT, V, S, compress> f(res, src);
    siteLoopNoReduction<FT, V, S, compress, CopyFunctor<FT, V, S, compress>>(
        f, geom, n_blas_simt);
}
```

## Generalization for Multiple Flavors

In the course of generalizing this code to multiple flavors, the existing
functionality was reused. Each function was *overloaded* to take an array of
spinor pointers, like this:

```{.cpp}
template <typename FT, int V, int S, bool compress, int num_flav>
void copySpinor(
    typename Geometry<FT, V, S, compress>::FourSpinorBlock *res[num_flav],
    const typename Geometry<FT, V, S, compress>::FourSpinorBlock *src[num_flav],
    const Geometry<FT, V, S, compress> &geom,
    int n_blas_simt)
{
    for (uint8_t f = 0; f < num_flav; ++f) {
        copySpinor(res[f], src[f], geom, n_blas_simt);
    }
}
```

The only difference in the argument list is the addition of the `[num_flav]` in
the end. So `res` and `src` are now arrays of pointers to spinors. There is one
little caveat: For function parameters, `T[N]` is the same as `T[]` and also
equivalent to `T *`. Although there is this `num_flav`, the function still does
not know the number of elements in the array. This has to be given as an
explicit template argument when calling this function.

## The Forbidden Conversion

The code above does *not* compile. The fundamental problem is that `T **`
cannot be converted to `T const **`. The reason is the following:

Create some array of pointers.

```{.cpp}
T **a;
```

Make the implicit conversion from `T **` to `T const **`. This does not
compile, but let's press on.

```{.cpp}
T const **b = a;
```

Create a pointer to `T const`. The values that `c` points to are really
`const`. The values that `b` point to are mutable, we just have “forgotten”
that for a moment.

```{.cpp}
T const *c;
```

Let the first pointer in `b` point to `c`. This is allowed, since `*b` is of
type `T const *`, just like `c`.

```{.cpp}
b[0] = c;
```

Now try to change something that `a` points to. It should work, there is no
`const` anywhere in `a`. But since we let did `a[0] = c` using the proxy
`b`, `a[0][0]` points to a `T const`!.

```{.cpp}
T some_other_t;
a[0][0] = some_other_t;
```

This is the reason why the above conversion is forbidden.

## The Allowed but Buggy Conversion

The solution to the whole problem is introducing another `const`:

```{.cpp}
T const *const *b = a;
```

Then we are not allowed to do `b[0] = c` any more, nothing bad can happen.

The above code with `T const *const *` works well when one actually has
`num_flav` pointers to `Spinor` and `num_flav` pointers to `Spinor const`.
However, in algorithms like the CG, one has a lot of mutable spinor fields that
sometimes take the role of an input or an output argument. This means that
`Spinor **` is passed into the argument that expects `Spinor const *const *`.

With GCC that is no problem, it happily converts `T **` to `T const *const *`.
Intel C++ 17 does not do that, and this is where some workaround is needed
until that conversion (allowed by the standard) is implemented.

## Template Workaround

There are a few ways to work around this bug:

- Write a `const` wrapper array for each argument before calling the functions.
  This is very tedious work for the caller. The magic should happen in the
  library, not in the client code.

- Remove the whole `const` correctness throughout the code. This of course is a
  bad idea.

- Add a non-`const` overload for each function parameter. There are a few
  functions with three spinor parameters, that would mean eight overloads. Also
  not maintainable.

The least worst solution is to let a template deduce the `const` parameter for
us.

For the `copySpinor` function, it would like the following: `Spinor1` is now a
template argument that is expected be `Spinor` or `Spinor const`. This way
every possible way will be covered:

```{.cpp}
template <typename FT,
          int V,
          int S,
          bool compress,
          int num_flav,
          typename Spinor1>
void copySpinor(
    typename Geometry<FT, V, S, compress>::FourSpinorBlock *res[num_flav],
    Spinor1 *const src[num_flav],
    const Geometry<FT, V, S, compress> &geom,
    int n_blas_simt)
{
    for (uint8_t f = 0; f < num_flav; ++f) {
        copySpinor(res[f], src[f], geom, n_blas_simt);
    }
}
```

This *does* work. The only nuisance is that supplying a type other than
`Spinor` or `Spinor const` will result in incomprehensible template error
messages *within* this function. It will complain about the types that are fed
into the old `copySpinor` functions because that only takes `Spinor const` for
the second argument.

In order to make the error messages more useful to the user, I first used a
static assertion to ensure that those types are the same:

```{.cpp}
static_assert(std::is_same<const Spinor, const S>::value,
              "Template parameter must be `const Spinor` or `Spinor`.");
```

This uses the compile time type traits. `std::is_same<T1, T2>::value` will only
be `true` if the types are the same. Another trick is used is that `const const
S` is just `const S`. Therefore this can directly check whether `S` is `Spinor
const` or `Spinor`.

The error message of this is pretty direct and the user will see the failure
the static assertion.

In the comments of the [Stack Overflow post][SO], it has been suggested to use
`std::enable_if` for this. This uses SFINAE (Substitution Failure Is Not An
Error), which is an emergent property of C++. `enable_if<C, T>` is a `struct`
which has a member `typedef T type`, *iff* the condition `C` is `true`.
Otherwise the `struct` does not have this member. For the return value (which
shall be `void`), I use the `enable_if` together with the `is_same`. If the
user puts in the wrong type for `Spinor1`, the `is_same<>::value` will evaluate
to `false`. Then the `enable_if<false, void>` will *not* have the `typedef void
type` member, therefore `is_same<false, void>::type` will be a *substitution
falue*. But that is not an error, it is just that this overload is taken out
from the overload resolution. Therefore this function does not “exist” any
more. The user will get the error that the function is not viable because it
has been disabled due to the non-matching types of `Spinor1` and what is
required.

This is the final version of the is BLAS routine:

```{.cpp}
template <typename FT,
          int V,
          int S,
          bool compress,
          int num_flav,
          typename Spinor1>
typename std::enable_if<
    std::is_same<const typename Geometry<FT, V, S, compress>::FourSpinorBlock,
                 const Spinor1>::value,
    void>::type
copySpinor(
    typename Geometry<FT, V, S, compress>::FourSpinorBlock *res[num_flav],
    Spinor1 *const src[num_flav],
    const Geometry<FT, V, S, compress> &geom,
    int n_blas_simt)
{
    for (uint8_t f = 0; f < num_flav; ++f) {
        copySpinor(res[f], src[f], geom, n_blas_simt);
    }
}
```

## Virtual Member Functions

Virtual functions and templates do not mix. Therefore this approach does not
work for templates. One has to create an overload with `const` and without
`const` for each argument. These are extra guarded with an `#ifdef
__INTEL_COMPILER`. This way it will be more apparent that they are only needed
to work around the limitation.

It looks as follows in the abstract solver interface:

```{.cpp}
    virtual void
    operator()(Spinor *x[num_flav],
               const Spinor *const rhs[num_flav],
               const double RsdTarget,
               int &niters,
               double &rsd_sq_final,
               unsigned long &site_flops,
               unsigned long &mv_apps,
               int isign,
               bool verboseP) = 0;

#ifdef __INTEL_COMPILER
    virtual void
    operator()(Spinor *x[num_flav],
               Spinor *const rhs[num_flav],
               const double RsdTarget,
               int &niters,
               double &rsd_sq_final,
               unsigned long &site_flops,
               unsigned long &mv_apps,
               int isign,
               bool verboseP) {
        (*this)(x,
                const_cast<Spinor const *const *>(rhs),
                RsdTarget,
                niters,
                rsd_sq_final,
                site_flops,
                mv_apps,
                isign,
                verboseP);
    }
#endif
```

There is an explicit `const_cast` needed, which will add the `const`. Since
there is a second `const`, nothing bad could happen, it just forces the Intel
compiler to do the right thing.

<!-- vim: set spell tw=79 :-->
