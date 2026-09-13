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
// Dense Hermitian eigensolver using LAPACK zheev on the calling rank.
// Copy row-pointer mat_in[i][j] to column-major scratch. Keep
// LWORK = (2n-1)*10 to preserve the historical LAPACK path and results.
// compute_evec selects JOBZ; nullable evec_out controls write-back.
void solve_dense_hermitian(int n, const std::complex<double> *const *mat_in, double *eval_out,
                           std::complex<double> **evec_out, bool compute_evec, char uplo = 'U');

// The same call returning the LAPACK INFO (0 on success) instead of exiting, for
// use inside OpenMP regions: exit() aborts through MPI and must be called from
// the master thread after the region.
int solve_dense_hermitian_info(int n, const std::complex<double> *const *mat_in, double *eval_out,
                               std::complex<double> **evec_out, bool compute_evec, char uplo = 'U');

// LAPACK zheevd for SCPH Eigen matrices, with queried workspace.
// Return ascending eigenvalues and, when requested, eigenvectors in columns.
// Eigenvector phases may differ from Eigen's solver.
void solve_dense_hermitian_dc(const Eigen::MatrixXcd &mat_in, Eigen::VectorXd &eval_out, Eigen::MatrixXcd *evec_out);
} // namespace PHON_NS
