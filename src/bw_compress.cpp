// Parallel block compression for the vendored libBigWig writer
// (bwsWriteBlocks() in libBigWig/bwWrite.c).
//
// The OpenMP loop lives here, in C++, so that only C++ code in the package uses
// OpenMP. Writing R Extensions (section 1.2.1.1) notes that using OpenMP from
// more than one of C, C++ and Fortran in one package is not portable, and the
// package is linked by the C++ compiler with $(SHLIB_OPENMP_CXXFLAGS).

#include <zlib.h>

// Compresses src[j] (srcLen[j] bytes) into dst[j] (capacity `bound` bytes) for
// j in [0, n) at zlib level `level`, using up to `nthreads` threads. Stores each
// compressed length in dstLen[j]. Returns 0 on success, 1 if any block failed.
extern "C" int bamscale_compress_blocks(unsigned char **dst, uLongf *dstLen,
                                        uLongf bound, unsigned char **src,
                                        const uLongf *srcLen, long long n,
                                        int level, int nthreads) {
  int err = 0;
#ifdef _OPENMP
#pragma omp parallel for num_threads(nthreads > 1 ? nthreads : 1) schedule(dynamic, 32)
#endif
  for (long long j = 0; j < n; j++) {
    uLongf sz = bound;
    if (compress2(dst[j], &sz, src[j], srcLen[j], level) != Z_OK) {
#ifdef _OPENMP
#pragma omp atomic write
#endif
      err = 1;
    }
    dstLen[j] = sz;
  }
  return err;
}
