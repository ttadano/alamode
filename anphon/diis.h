//
// Created by Terumasa Tadano on 2025/11/26.
//

// Header file for performing the generalized Direct Inversion in the Iterative Subspace (DIIS)

#pragma once

#include <Eigen/Core>
#include <deque>
#include <vector>

/// Generalized DIIS optimizer: extrapolate trial vectors to minimize the residual.
class GDIIS
{
public:
    /// Store up to max_history vectors with mixing_beta in (0, 1].
    /// verbosity > 0 enables logging.
    GDIIS(int max_history = 10, double mixing_beta = 0.5, int verbosity = 0);

    ~GDIIS() = default;

    /// Add a trial vector and its residual, e.g. f(x) - x, to the history.
    void push(const Eigen::VectorXd &x_trial, const Eigen::VectorXd &error);

    /// Write the DIIS extrapolation to x_new; return true on success.
    bool extrapolate(Eigen::VectorXd &x_new);

    /// Clear the DIIS history.
    void clear();

    /// Return the number of stored vectors.
    [[nodiscard]] int size() const
    {
        return static_cast<int>(history_x.size());
    }

    /// Return true when at least two vectors are stored.
    [[nodiscard]] bool is_ready() const
    {
        return size() >= 2;
    }

    /// Set the maximum history size.
    void set_max_history(int max_hist)
    {
        max_history_ = max_hist;
    }

    /// Set the mixing parameter beta in (0, 1].
    void set_mixing_beta(double beta)
    {
        mixing_beta_ = beta;
    }

private:
    int max_history_;                          ///< Maximum number of vectors to keep
    double mixing_beta_;                       ///< Mixing parameter for simple mixing fallback
    int verbosity_;                            ///< Verbosity level for logging
    std::deque<Eigen::VectorXd> history_x;     ///< History of trial vectors
    std::deque<Eigen::VectorXd> history_error; ///< History of error vectors

    /// Solve for DIIS coefficients; return true on success.
    bool solve_diis_equations(Eigen::VectorXd &coeffs);
};
