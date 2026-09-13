//
// Created by Terumasa Tadano on 25/06/04.
//

#pragma once

#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <cstddef> // for size_t
#include "constraint.h"

/**
 * Minimize ||A x - b||_2 using SVD. amat is column-major M x N;
 * bvec has length M and param_out length N. verbosity > 0 prints diagnostics.
 * Returns 0 on success, a nonzero solver status on failure.
 */
auto least_squares_svd(const size_t N, const size_t M, double *amat, const double *bvec, double *param_out,
                       const int verbosity) -> int;

/**
 * Minimize ||A x - b||_2 subject to C x = d using rank reduction and dgglse.
 * amat is column-major M x N; cmat is a P x N row-pointer array.
 * bvec, dvec, and param_out have lengths M, P, and N.
 * verbosity > 0 prints diagnostics. Returns the LAPACK status (0 on success).
 */
auto least_squares_with_constraints_gqr(const size_t N, const size_t M, const size_t P, double *amat,
                                        const double *bvec, double *param_out, const double *const *cmat,
                                        const double *dvec, const int verbosity) -> int;

/**
 * Minimize ||A x - b||_2 subject to C x = d using SVD and the
 * null-space projector I - C^+ C, starting from x0 = C^+ d.
 * verbosity > 0 prints diagnostics. Returns 0 on success, nonzero on failure.
 */
auto least_squares_with_constraints_svd(const size_t N, const size_t M, const size_t P,
                                        double *amat,              // A: (M×N) column-major, can be overwritten
                                        double *bvec,              // b: (M) vector, can be overwritten
                                        double *param_out,         // output x (length N)
                                        const double *const *cmat, // C[i][j] pointer array (not contiguous)
                                        const double *dvec_orig,   // d: (P) vector, will be copied locally
                                        const int verbosity) -> int;


/**
 * Minimize ||A x - b||_2 with the selected sparse solver.
 * sp_mat is M x N, sp_bvec has length M, and x_out has length N.
 * tolerance_iteration and maxnum_iteration control iterative convergence.
 * Returns 0 on success, nonzero on failure.
 */
auto least_squares_eigen_sparse_solver(const Eigen::SparseMatrix<double> &sp_mat, const Eigen::VectorXd &sp_bvec,
                                       Eigen::VectorXd &x_out, const std::string &solver_type,
                                       const double tolerance_iteration, const int maxnum_iteration) -> int;


/**
 * Select independent rows of cmat[P][N] using pivoted QR of C^T.
 * Return column-major C_red[r][N], matching d_red[r], and numerical rank r.
 * verbosity > 0 prints diagnostics. Returns 0 on success or LAPACK INFO.
 */
auto get_independent_rows(const size_t N, const size_t P, const double *const *cmat, const double *dvec,
                          const int verbosity, std::vector<double> &C_red, std::vector<double> &d_red, int &r) -> int;


/**
 * Solve sparse min ||A x - b||_2 subject to C x = d.
 * A is M x N and C is P x N; b and d have lengths M and P.
 * Return x[N] and Lagrange multipliers lambda[P].
 */
auto solveGQRSparse(const Eigen::SparseMatrix<double> &A, const Eigen::VectorXd &b,
                    const Eigen::SparseMatrix<double> &C, const Eigen::VectorXd &d, Eigen::VectorXd &x,
                    Eigen::VectorXd &lambda, const int verbosity = 0, const std::string &solver_type = "",
                    const double tolerance_iteration = 1.0e-8, const int maxnum_iteration = 10000) -> void;

/**
 * Find independent rows of column-major A_data[M][N]. Return rank and
 * zero-based row pivots in pivot order. For tol <= 0, use
 * max(M,N)*|R(0,0)|*epsilon. verbosity > 0 prints diagnostics.
 * Returns 0 on success or LAPACK INFO.
 */
auto find_independent_rows_dense(int M, int N, double *A_data, double tol, int &rank, std::vector<int> &pivots,
                                 const int verbosity = 0) -> int;


// Sentinel for find_independent_rows_dense(): use the existing LAPACK-style auto tolerance
// max(M, N) * |R(0,0)| * eps. This preserves the current QR rank policy.
constexpr double rank_tolerance_auto = -1.0;


/// Extract independent rows of C_sparse (P x N) using pivoted QR of C^T.
///
/// @param C_sparse  Input (P×N), row-major sparse
/// @param dvec      Input length-P
/// @param verbosity >1 prints debug
/// @param tolerance Rank tolerance. Use rank_tolerance_auto to preserve the auto policy.
/// @param C_red     Output (r×N) row-major sparse of independent rows
/// @param d_red     Output length-r of corresponding dvec entries
/// @param r         Output numerical row-rank
/// @returns 0 on success, non-zero LAPACK INFO on failure
auto get_independent_rows_lapack_sparse(const Eigen::SparseMatrix<double> &C_sparse, const Eigen::VectorXd &dvec,
                                        const int verbosity, const double tolerance, Eigen::SparseMatrix<double> &C_red,
                                        Eigen::VectorXd &d_red, int &r) -> int;


auto get_independent_rows_lapack_sparse(const size_t ncols, ConstraintSparseForm &C_sparse, const int verbosity,
                                        const double tolerance, int &r) -> int;
