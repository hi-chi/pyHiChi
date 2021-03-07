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


#if _OPENMP >= 201307
    #define OMP_FOR()   __pragma("omp parallel for")
    #define OMP_FOR_COLLAPSE()   __pragma("omp parallel for collapse(2)")
    #define OMP_SIMD()  __pragma("omp simd")
#else
    #define OMP_FOR()   __pragma("omp parallel for")
    #define OMP_FOR_COLLAPSE()   __pragma("omp parallel for")
    #define OMP_SIMD()  __pragma("ivdep")
#endif