/*
 v4_index_transform.h
 Copyright (c) 2026 Terumasa Tadano
 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory
 or http://opensource.org/licenses/mit-license.php for information.
*/
/*
 One index transform of the quartic interaction tensor as a complex GEMM,
 shared by the two V4 builders in scph_v3v4_elements.cpp.
*/
#pragma once

#include <Eigen/Core> // must precede blas_wrapper.h: under EIGEN_USE_BLAS Eigen's blas.h declares zgemm_
#include <complex>
#include <cstddef>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "blas_wrapper.h"

namespace PHON_NS
{
namespace v4_index_transform
{
// The input is a row-major (ns x ncols) buffer whose OUTERMOST mode index is
// contracted with E[out][in] (row-major ns x ns); the result is written as
// (ncols x ns) with the new index innermost, so the flat layout rotates
// [a rest] -> [rest out]:
//     out[col][x] = alpha * sum_a E[x][a] * in[a][col].
// In column-major BLAS terms C(ns x ncols) = E_buf^T * In_buf^T with ldb = ncols
// and ldc = ns. ncols = ns^3 for a whole (k1, k2) slice (k-point builder), ns^2
// for one set (ik_prod, is) of the band builder.
//
// The columns are split into one slab per OpenMP thread, each slab being one
// zgemm on a (ns x ns) x (ns x n_slab) problem, so the product uses the OpenMP
// team even when the BLAS runs single-threaded (a BLAS pinned to one thread
// under MPI, as the CMake note recommends, or a sequential library). MKL and an
// OpenMP-built OpenBLAS run one thread per call inside a parallel region; a
// pthreads OpenBLAS or Accelerate must be pinned (OPENBLAS_NUM_THREADS=1 /
// VECLIB_MAXIMUM_THREADS=1) or they oversubscribe the cores.
inline void transform_index_gemm(const std::complex<double> *evec_row_major, const std::complex<double> *buf_in,
                                 std::complex<double> *buf_out, const std::size_t ns, const std::size_t ncols,
                                 const std::complex<double> alpha_in)
{
    auto slab = [&](const std::size_t c0, const std::size_t c1) {
        if (c1 <= c0) return;
        int m = static_cast<int>(ns);
        int n = static_cast<int>(c1 - c0);
        int k = static_cast<int>(ns);
        int ldb = static_cast<int>(ncols);
        auto alpha = alpha_in;
        auto beta = std::complex<double>(0.0, 0.0);
        zgemm_cpx("T",
                  "T",
                  &m,
                  &n,
                  &k,
                  &alpha,
                  const_cast<std::complex<double> *>(evec_row_major),
                  &m,
                  const_cast<std::complex<double> *>(buf_in + c0),
                  &ldb,
                  &beta,
                  buf_out + c0 * ns,
                  &m);
    };
#ifdef _OPENMP
#pragma omp parallel
    {
        const auto nt = static_cast<std::size_t>(omp_get_num_threads());
        const auto it = static_cast<std::size_t>(omp_get_thread_num());
        slab(ncols * it / nt, ncols * (it + 1) / nt);
    }
#else
    slab(0, ncols);
#endif
}
} // namespace v4_index_transform
} // namespace PHON_NS
