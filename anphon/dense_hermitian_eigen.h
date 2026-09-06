/*
 dense_hermitian_eigen.h

 Copyright (c) 2026 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory
 or http://opensource.org/licenses/mit-license.php for information.
*/

#pragma once

#include <Eigen/Core>
#include <complex>

namespace PHON_NS
{
// Backend seam for dense Hermitian eigenproblems (used by dynamical-matrix
// diagonalization).
//
// v1 backend: LAPACK zheev on the calling rank. Candidate later backends
// behind the same call: ELPA and MAGMA.
//
// mat_in is row-pointer indexed [i][j] and is copied into column-major scratch
// internally. LWORK is fixed at (2n-1)*10 with no workspace query because zheev
// may take a different (blocked vs unblocked) path for different workspace
// sizes and bit-identical results with the historical call are required.
// compute_evec drives JOBZ ('V'/'N'); evec_out (nullable) independently gates
// the eigenvector write-back.
void solve_dense_hermitian(int n, const std::complex<double> *const *mat_in, double *eval_out,
                           std::complex<double> **evec_out, bool compute_evec, char uplo = 'U');

// Divide-and-conquer variant (LAPACK zheevd, workspace queried) on Eigen
// matrices, for the SCPH solver: eigenvalues ascending in eval_out; the
// eigenvectors, when evec_out is given, in its columns (the convention of
// Eigen::SelfAdjointEigenSolver::eigenvectors()). zheevd is several times
// faster than Eigen's tridiagonal QR at n of a few hundred and uses the
// threaded MKL/OpenBLAS kernels; results agree to roundoff (eigenvector phases
// may differ, which the SCPH solver is invariant to).
void solve_dense_hermitian_dc(const Eigen::MatrixXcd &mat_in, Eigen::VectorXd &eval_out, Eigen::MatrixXcd *evec_out);
} // namespace PHON_NS
