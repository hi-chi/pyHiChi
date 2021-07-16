#pragma once

#ifdef _MSC_VER
    #define forceinline __forceinline
#elif defined(__GNUC__)
    #define forceinline inline __attribute__((__always_inline__))
#elif defined(__CLANG__)
    #if __has_attribute(__always_inline__)
        #define forceinline inline __attribute__((__always_inline__))
    #else
        #define forceinline inline
    #endif
#else
    #define forceinline inline
#endif


#ifdef _MSC_VER
    #define PRAGMA(opt) __pragma(opt)
#else
    #define PRAGMA(opt) _Pragma(#opt)
#endif


#if _OPENMP >= 201307
    #define OMP_FOR()  PRAGMA(omp parallel for)
    #define OMP_FOR_COLLAPSE()  PRAGMA(omp parallel for collapse(2))
    #define OMP_FOR_SIMD()  PRAGMA(omp parallel for simd)
    #define OMP_SIMD()  PRAGMA(omp simd)
#else
    #define OMP_FOR()  PRAGMA(omp parallel for)
    #define OMP_FOR_COLLAPSE()  PRAGMA(omp parallel for)
    #define OMP_FOR_SIMD()  PRAGMA(omp parallel for)
    #define OMP_SIMD()  PRAGMA(ivdep)
#endif
