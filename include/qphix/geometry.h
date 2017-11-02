#pragma once

#include "qphix/dslash_utils.h"
#include "qphix/print_utils.h"
#include "qphix_codegen/decl_common.h"

#include <cstdlib>
#include <stdexcept>
#include <iostream>
#include <utility>

#if defined(QPHIX_AVX_SOURCE) || defined(QPHIX_AVX2_SOURCE) ||                      \
    defined(QPHIX_AVX512_SOURCE)
#include <immintrin.h>
#endif

namespace QPhiX
{

struct CorePhase {
  int Ct;
  int Cyz;
  int startBlock;
};

#if defined(QPHIX_MIC_SOURCE) || defined(QPHIX_AVX512_SOURCE)
inline float cvtHalf2Float(half val)
{
  float ret;
#if defined(QPHIX_MIC_SOURCE)
  _mm512_mask_packstorelo_ps(
      &ret,
      0x1,
      _mm512_mask_extloadunpacklo_ps(
          _mm512_undefined_ps(), 0x1, &val, _MM_UPCONV_PS_FLOAT16, _MM_HINT_NONE));
#elif defined(QPHIX_AVX512_SOURCE)
  _mm512_mask_storeu_ps(&ret, 0x1, _mm512_cvtph_ps(_mm256_set1_epi16(val)));
#endif // defined(QPHIX_MIC_SOURCE)
  return ret;
}

inline half cvtFloat2Half(float val)
{
  half ret;
#if defined(QPHIX_MIC_SOURCE)
  _mm512_mask_extpackstorelo_ps(
      &ret,
      0x1,
      _mm512_mask_loadunpacklo_ps(_mm512_undefined_ps(), 0x1, &val),
      _MM_DOWNCONV_PS_FLOAT16,
      _MM_HINT_NONE);
#elif defined(QPHIX_AVX512_SOURCE)
  // ret = _mm256_extract_epi16( _mm512_cvt_roundps_ph(_mm512_set1_ps(val),
  // _MM_FROUND_TO_NEAREST_INT), 0);
  unsigned temp;
  _mm512_mask_storeu_epi32(&temp,
                           0x01,
                           _mm512_castsi256_si512(_mm512_cvt_roundps_ph(
                               _mm512_set1_ps(val), _MM_FROUND_TO_NEAREST_INT)));
  ret = (half)temp;
#endif // defined(QPHIX_MIC_SOURCE)
  return ret;
}

// rep: cast 'in' of type T2 into a 'T1' and return it.
template <typename T1, typename T2>
T1 rep(const T2 &in)
{
  if (sizeof(T1) != sizeof(T2)) {

    if (sizeof(T1) == 2) // we are converting float/double to half
      return cvtFloat2Half((float)in);
    else if (sizeof(T2) == 2) // we are converting half to float/double
      return (T1)cvtHalf2Float(in);
    else
      return (T1)in; // we are converting between float and double so just
    // cast is enough
  } else {
    return static_cast<T1>(in); // both T1 and T2 are same
  }
}

#else // defined(QPHIX_MIC_SOURCE) || defined(QPHIX_AVX512_SOURCE)

// rep: cast 'in' of type T2 into a 'T1' and return it.
template <typename T1, typename T2>
T1 rep(const T2 &in)
{
  return (T1)(in);
}

#endif // defined(QPHIX_MIC_SOURCE) || defined(QPHIX_AVX512_SOURCE)

template <typename FT>
inline char const *type_name();

template <>
inline char const *type_name<double>()
{
  return "double";
}

template <>
inline char const *type_name<float>()
{
  return "float";
}

template <>
inline char const *type_name<half>()
{
  return "half";
}

template <typename T, int V, int S, bool compressP>
class Geometry
{
 public:
  typedef T FT;
  int constexpr static veclen = V;
  int constexpr static soalen = S;
  bool constexpr static compress12 = compressP;

  typedef typename Types<T, V, S, compressP>::FourSpinorBlock FourSpinorBlock;
  typedef typename Types<T, V, S, compressP>::TwoSpinorBlock TwoSpinorBlock;
  typedef typename Types<T, V, S, compressP>::SU3MatrixBlock SU3MatrixBlock;
  typedef typename Types<T, V, S, compressP>::CloverBlock CloverBlock;
  typedef typename Types<T, V, S, compressP>::FullCloverBlock FullCloverBlock;

  Geometry(const int latt_size[],
           int By_,
           int Bz_,
           int NCores_,
           int Sy_,
           int Sz_,
           int PadXY_,
           int PadXYZ_,
           int MinCt_,
           bool const verbose = false)
      : Nd(4), By(By_), Bz(Bz_), num_cores(NCores_), Sy(Sy_), Sz(Sz_), PadXY(PadXY_),
        PadXYZ(PadXYZ_), MinCt(MinCt_), nsimt(Sy_ * Sz_),
        num_threads(NCores_ * Sy_ * Sz_)
  {
    Nx_ = latt_size[0];
    Ny_ = latt_size[1];
    Nz_ = latt_size[2];
    Nt_ = latt_size[3];
    Nxh_ = Nx_ / 2;

    // Ensure that blocking is possible.
    if (Ny_ % By_ != 0 || Nz_ % Bz_ != 0) {
      throw std::domain_error("Local lattice size Ny must be divisible by "
                              "blocking length By. Same for Nz and Bz.");
    }

    nvecs_ = Nxh() / S;
    if (Nxh() % S != 0)
      nvecs_++;

    if (V % S != 0) {
      std::cerr << "Error: Geometry constructor: SOALEN=" << S
                << " does not divide V=" << V << std::endl;
      std::abort();
    }
    ngy_ = V / S;

    // Padding constants
    Pxy = (nvecs_ * Ny_ + PadXY);
    Pxyz = (Pxy * Nz_ + PadXYZ);

    // Allos sizes
    spinor_bytes = (Pxyz * Nt_ + 1) * sizeof(FourSpinorBlock);
    gauge_bytes = ((Pxyz * Nt_ * S) / V) * sizeof(SU3MatrixBlock);
    clover_bytes = ((Pxyz * Nt_ * S) / V) * sizeof(CloverBlock);
    full_clover_bytes = ((Pxyz * Nt_ * S) / V) * sizeof(FullCloverBlock);

    // This works out the phase breakdown
    int ly = Ny_ / By;
    int lz = Nz_ / Bz;
    int rem = ly * lz;
    int stblk = 0;
    n_phases = 0;
    int n_cores_per_minct = num_cores / MinCt;
    while (rem > 0) {
      int ctd = n_cores_per_minct / rem;
      int ctu = (n_cores_per_minct + rem - 1) / rem;
      CorePhase &p = getCorePhase(n_phases);
      p.Ct = (ctu <= 4 ? ctu : ctd) * MinCt;
      p.Cyz = num_cores / p.Ct;
      if (p.Cyz > rem)
        p.Cyz = rem;
      p.startBlock = stblk;
      stblk += p.Cyz;
      rem -= p.Cyz;
      //	masterPrintf("Phase %d: Cyz = %d Ct = %d, start = %d\n",
      // n_phases, p.Cyz, p.Ct, p.startBlock);
      n_phases++;
    }

    if (verbose) {
      masterPrintf("Constructed a new Geometry object.\n"
                   "  Template parameters are:\n"
                   "    FT          %s\n"
                   "    veclen      %d\n"
                   "    soalen      %d\n"
                   "    compress12  %s\n",
                   type_name<T>(),
                   veclen,
                   soalen,
                   compress12 ? "true" : "false");
      masterPrintf(
          "  The local sizes for QPhiX data structures in MibiByte (1024**2) are:\n"
          "    Spinor:     %10ld\n"
          "    Gauge:      %10ld\n"
          "    Clover:     %10ld\n"
          "    FullClover: %10ld\n",
          spinor_bytes / 1024 / 1024,
          gauge_bytes / 1024 / 1024,
          clover_bytes / 1024 / 1024,
          full_clover_bytes / 1024 / 1024);

      masterPrintf("  Geometry properties are:\n");
      masterPrintf("    Local Lattice size %d %d %d %d\n", Nx_, Ny_, Nz_, Nt_);
      masterPrintf("    number of vectors  %d\n", nvecs_);
      masterPrintf("    Blocking           %d %d\n", By, Bz);
      masterPrintf("    nGY                %d\n", ngy_);
      masterPrintf("    PadXY PadXYZ       %d %d\n", PadXY, PadXYZ);
      masterPrintf("    MinCt              %d\n", MinCt);
      masterPrintf("    Phases             %d\n", n_phases);
      masterPrintf("  Threading options are:\n");
      masterPrintf("    Cores              %d\n", num_cores);
      masterPrintf("    Sy Sz              %d %d\n", Sy, Sz);
      masterPrintf("    Total threads      %d\n", num_threads);
    }
  }

  ~Geometry() {}

  int Nxh() const { return Nxh_; }
  int Nx() const { return Nx_; }
  int Ny() const { return Ny_; }
  int Nz() const { return Nz_; }
  int Nt() const { return Nt_; }
  int nVecs() const { return nvecs_; }
  int nGY() const { return ngy_; }

  int getBy() const { return By; }
  int getBz() const { return Bz; }
  int getSy() const { return Sy; }
  int getSz() const { return Sz; }
  int getPadXY() const { return PadXY; }
  int getPadXYZ() const { return PadXYZ; }
  int getPxy() const { return Pxy; }
  int getPxyz() const { return Pxyz; }
  int getNSIMT() const { return nsimt; }
  int getNumThreads() const { return num_threads; }
  int getNumCores() const { return num_cores; }
  int getMinCt() const { return MinCt; }
  int getVolCB() const { return Nxh_ * Ny_ * Nz_ * Nt_; }
  int getNumPhases() const { return n_phases; }

  CorePhase &getCorePhase(int i) { return phase[i]; }
  const CorePhase &getCorePhase(int i) const { return phase[i]; }

  /*! \brief Checkerboarded FourSpinor Allocator
  *
  * Allocates a single checkerboard of a Four Spinor.
  * An extra spinor is allocated beyond what is required.
  */
  FourSpinorBlock *allocCBFourSpinor()
  {

    FourSpinorBlock *ret_val = (FourSpinorBlock *)BUFFER_MALLOC(spinor_bytes, 128);
    if (ret_val == (FourSpinorBlock *)0x0) {
      masterPrintf("Failed to allocate FourSpinorBlock\n");
      std::abort();
    }

    // Zero the field.
    // Cast the pointer.
    T *ret_val_ft = (T *)ret_val;

    // change from number of bytes to number of T type elements
    size_t num_ft = spinor_bytes / sizeof(T);

    // Zero it all (including) (especially) the pad regions.
    // FIXME: this is not NUMA friendly necessarily
    for (int i = 0; i < num_ft; i++) {
      ret_val_ft[i] = rep<T, double>(0.0);
    }

    return ret_val + 1;
  }

  void free(FourSpinorBlock *p)
  {
    // The C standard `free` will be a no-operation when the null pointer is passed.
    // This function should exhibit the same behavior.
    if (p == nullptr) {
      return;
    }

    FourSpinorBlock *freeme = p - 1;
    BUFFER_FREE(freeme, spinor_bytes);
  }

  /*! \brief Checkerboard Gauge Field Allocation
  *
  * This function allocates memory for a single checkerboard of
  * a gauge field
  */
  SU3MatrixBlock *allocCBGauge()
  {
    SU3MatrixBlock *ret_val = (SU3MatrixBlock *)BUFFER_MALLOC(gauge_bytes, 128);
    if (ret_val == (SU3MatrixBlock *)0x0) {
      masterPrintf("Failed to allocate SU3MatrixBlock\n");
      std::abort();
    }

    // For AVX we should loop and zero it here....
    // later on.

    // Zero the field.
    // Cast the pointer.
    T *ret_val_ft = (T *)ret_val;

    // change from number of bytes to number of T type elements
    size_t num_ft = gauge_bytes / sizeof(T);

    // Zero it all (including) (especially) the pad regions.
    // FIXME: this is not NUMA friendly necessarily
    for (int i = 0; i < num_ft; i++) {
      ret_val_ft[i] = rep<T, double>(0.0);
    }

    return ret_val;
  }

  void free(SU3MatrixBlock *p) { BUFFER_FREE(p, gauge_bytes); }

  CloverBlock *allocCBClov()
  {
    CloverBlock *ret_val = (CloverBlock *)BUFFER_MALLOC(clover_bytes, 128);
    if (ret_val == (CloverBlock *)0x0) {
      masterPrintf("Failed to allocate CloverBlock\n");
      std::abort();
    }

    // For AVX we should loop and zero it here....
    // later on.

    // Zero the field.
    // Cast the pointer.
    T *ret_val_ft = (T *)ret_val;

    // change from number of bytes to number of T type elements
    size_t num_ft = clover_bytes / sizeof(T);

    // Zero it all (including) (especially) the pad regions.
    // FIXME: this is not NUMA friendly necessarily
    for (int i = 0; i < num_ft; i++) {
      ret_val_ft[i] = rep<T, double>(0.0);
    }

    return ret_val;
  }

  void free(CloverBlock *p) { BUFFER_FREE(p, clover_bytes); }

  FullCloverBlock *allocCBFullClov()
  {
    FullCloverBlock *ret_val =
        (FullCloverBlock *)BUFFER_MALLOC(full_clover_bytes, 128);
    if (ret_val == (FullCloverBlock *)0x0) {
      masterPrintf("Failed to allocate FullCloverBlock\n");
      std::abort();
    }

    // For AVX we should loop and zero it here....
    // later on.

    // Zero the field.
    // Cast the pointer.
    T *ret_val_ft = (T *)ret_val;

    // change from number of bytes to number of T type elements
    size_t num_ft = full_clover_bytes / sizeof(T);

    // Zero it all (including) (especially) the pad regions.
    // FIXME: this is not NUMA friendly necessarily
    for (int i = 0; i < num_ft; i++) {
      ret_val_ft[i] = rep<T, double>(0.0);
    }

    return ret_val;
  }

  void free(FullCloverBlock *p) { BUFFER_FREE(p, full_clover_bytes); }

 private:
  int Nxh_;
  int Nx_;
  int Ny_;
  int Nz_;
  int Nt_;

  const int Nd;
  const int By;
  const int Bz;
  const int Sy;
  const int Sz;
  const int PadXY;
  const int PadXYZ;
  int Pxy;
  int Pxyz;

  /** Minimum number of cores in T dir.

    - MinCt = 1 for single socket/Xeon Phi
    - MinCt = 2 for dual socket
    - MinCt = 4 for quad socket
  */
  int MinCt;

  const int nsimt;
  const int num_threads;
  const int num_cores;

  int nvecs_;
  int ngy_;

  // Dhiraj's new core mapping
  static const int MAX_PHASES = 128;
  CorePhase phase[MAX_PHASES];
  int n_phases;

  size_t gauge_bytes;
  size_t spinor_bytes;
  size_t clover_bytes;
  size_t full_clover_bytes;
};

template <typename FT, int veclen, int soalen, bool compress12>
class FourSpinorHandle
{
 public:
  typedef Geometry<FT, veclen, soalen, compress12> Geom;
  typedef typename Geom::FourSpinorBlock ValueType;

  FourSpinorHandle(Geom &geom) : value(geom.allocCBFourSpinor()), geom(geom) {}

  FourSpinorHandle(FourSpinorHandle<FT, veclen, soalen, compress12> &&other) noexcept
      : value(std::move(other.value)),
        geom(other.geom)
  {
    other.value = nullptr;
  }

  ~FourSpinorHandle() { geom.free(value); }

  ValueType *get() const { return value; }

 private:
  ValueType *value;
  Geom &geom;
};

template <typename FT, int veclen, int soalen, bool compress12>
FourSpinorHandle<FT, veclen, soalen, compress12>
makeFourSpinorHandle(Geometry<FT, veclen, soalen, compress12> &geom)
{
  return {geom};
}

template <typename FT, int veclen, int soalen, bool compress12>
class GaugeHandle
{
 public:
  typedef Geometry<FT, veclen, soalen, compress12> Geom;
  typedef typename Geom::SU3MatrixBlock ValueType;

  GaugeHandle(Geom &geom) : value(geom.allocCBGauge()), geom(geom) {}

  ~GaugeHandle() { geom.free(value); }

  ValueType *get() const { return value; }

 private:
  ValueType *value;
  Geom &geom;
};

template <typename FT, int veclen, int soalen, bool compress12>
GaugeHandle<FT, veclen, soalen, compress12>
makeGaugeHandle(Geometry<FT, veclen, soalen, compress12> &geom)
{
  return {geom};
}

template <typename FT, int veclen, int soalen, bool compress12>
class CloverHandle
{
 public:
  typedef Geometry<FT, veclen, soalen, compress12> Geom;
  typedef typename Geom::CloverBlock ValueType;

  CloverHandle(Geom &geom) : value(geom.allocCBClov()), geom(geom) {}

  ~CloverHandle() { geom.free(value); }

  ValueType *get() const { return value; }

 private:
  ValueType *value;
  Geom &geom;
};

template <typename FT, int veclen, int soalen, bool compress12>
CloverHandle<FT, veclen, soalen, compress12>
makeCloverHandle(Geometry<FT, veclen, soalen, compress12> &geom)
{
  return {geom};
}

template <typename FT, int veclen, int soalen, bool compress12>
class FullCloverHandle
{
 public:
  typedef Geometry<FT, veclen, soalen, compress12> Geom;
  typedef typename Geom::FullCloverBlock ValueType;

  FullCloverHandle(Geom &geom) : value(geom.allocCBFullClov()), geom(geom) {}

  ~FullCloverHandle() { geom.free(value); }

  ValueType *get() const { return value; }

 private:
  ValueType *value;
  Geom &geom;
};

template <typename FT, int veclen, int soalen, bool compress12>
FullCloverHandle<FT, veclen, soalen, compress12>
makeFullCloverHandle(Geometry<FT, veclen, soalen, compress12> &geom)
{
  return {geom};
}

} // Namespace
