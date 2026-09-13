/*
 dense_symmetric_eigen.h

 Copyright (c) 2026 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory
 or http://opensource.org/licenses/mit-license.php for information.
*/

#pragma once

#include <vector>

namespace PHON_NS
{
// Dense symmetric eigensolver using LAPACK dsyev on the calling rank.
// Overwrite column-major A[n][n] with eigenvectors in columns; return
// ascending eigenvalues in w. num_lowest requests a subset, but this
// backend computes the full spectrum for the caller to truncate.
void solve_dense_symmetric(int n, std::vector<double> &A, std::vector<double> &w, int num_lowest = -1);
} // namespace PHON_NS
