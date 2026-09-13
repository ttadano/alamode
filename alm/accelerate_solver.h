// Apple Accelerate sparse KKT solver. Keep Eigen/AccelerateSupport in its own
// translation unit to avoid conflicting BLAS/LAPACK declarations from
// blas_wrapper.h and lapack_wrapper.h.
#pragma once

#ifdef USE_ACCEL_BACKEND
#include <Eigen/SparseCore>

// Solve the symmetric-indefinite KKT system K x = rhs with Apple Accelerate's sparse LDL^T
// factorization. Returns true on success (factorization and solve both reported Eigen::Success).
bool solve_kkt_accelerate_ldlt(const Eigen::SparseMatrix<double> &K, const Eigen::VectorXd &rhs, Eigen::VectorXd &sol);
#endif
