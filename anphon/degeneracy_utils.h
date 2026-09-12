/*
 degeneracy_utils.h

 Copyright (c) 2026 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#pragma once

#include <cmath>
#include <vector>
#include "constants.h"

namespace PHON_NS
{
inline void find_degenerate_groups(const unsigned int ns, const double *eval_at_k, std::vector<int> &degeneracy_out,
                                   const double tol_omega = 1.0e-7)
{
    degeneracy_out.clear();

    auto omega_prev = eval_at_k[0];
    auto ideg = 1;

    for (unsigned int is = 1; is < ns; ++is) {
        const auto omega_now = eval_at_k[is];

        if (std::abs(omega_now - omega_prev) < tol_omega) {
            ++ideg;
        } else {
            degeneracy_out.push_back(ideg);
            ideg = 1;
            omega_prev = omega_now;
        }
    }
    degeneracy_out.push_back(ideg);
}

// Replace per-branch data by its average over each degenerate subspace at
// one k point. eval_at_k holds the ns sorted eigenvalues; data is [ns][width]
// row-major (width = 1 for scalars, 3 for Cartesian vectors, ntemp for
// temperature rows) and is averaged column-wise within each group of
// consecutive eigenvalues closer than tol_omega (~0.01 cm^-1 by default).
inline void average_over_degenerate_modes(const int ns, const double *eval_at_k, const int width, double *data,
                                          const double tol_omega = 1.0e-7)
{
    std::vector<int> degeneracy_at_k;
    find_degenerate_groups(ns, eval_at_k, degeneracy_at_k, tol_omega);

    std::vector<double> data_sum(width);

    int is = 0;
    for (const auto ideg_now: degeneracy_at_k) {
        if (ideg_now > 1) {
            for (int l = 0; l < width; ++l) data_sum[l] = 0.0;
            for (int k = is; k < is + ideg_now; ++k) {
                for (int l = 0; l < width; ++l) {
                    data_sum[l] += data[k * width + l];
                }
            }
            for (int k = is; k < is + ideg_now; ++k) {
                for (int l = 0; l < width; ++l) {
                    data[k * width + l] = data_sum[l] / static_cast<double>(ideg_now);
                }
            }
        }
        is += ideg_now;
    }
}

// ---------------------------------------------------------------------------
// Transport degeneracy blocks.
//
// A block replaces the Wigner pair weight by the band-like (Peierls) limit. That is
// legitimate only where the eigenvectors are numerically ambiguous: there the
// individual velocities are undefined and only the block trace Tr(P V^u P V^v P) is.
//
// TOL_CM is an EMPIRICAL numerical criterion, not a physical one, and no universal
// correctness follows from it:
//   * Branches closer than TOL_CM are TREATED as degenerate. A rough scale for the
//     eigenvalue noise of the diagonalisation is d(omega) ~ ns * eps * omega_max^2 / (2 omega)
//     (lambda = omega^2 with absolute error ~ ns * eps * lambda_max); for ns = 300,
//     omega_max = 1000 cm^-1 that is ~3e-8 cm^-1 at omega = 1 cm^-1 but ~3e-6 cm^-1 at
//     omega = 0.01 cm^-1, above TOL_CM. It is a scale estimate, not a bound, and the eps8
//     frequency guard (~1.1e-3 cm^-1) does not remove every soft mode. So a genuine
//     splitting below TOL_CM may well be numerically resolvable; the criterion merges it
//     anyway.
//   * Merging a pair with true splitting dw and summed HWHM G overestimates its weight by
//     the factor 1 + (dw/G)^2. Conductivity::compute_kappa reports the worst dw/G among
//     merged blocks so that this approximation is visible rather than silent.
//   * Constant lifetime within a block (the precondition for the trace form) is guaranteed
//     by SUBDIVIDING the damping-averaging groups below, not by the size of TOL_CM.
inline double transport_block_tol_cm()
{
    return 1.0e-6;
}

// Block bounds [lo, hi) for every branch at one k point. Damping groups from
// find_degenerate_groups are subdivided; within a group a new block starts whenever a
// branch lies further than tol_cm from the block's FIRST member (anchored, so a chain of
// individually close branches cannot grow into a wide block).
inline void transport_block_bounds(const unsigned int ns, const double *eval_at_k, const double tol_cm,
                                   std::vector<int> &lo_out, std::vector<int> &hi_out)
{
    lo_out.assign(ns, 0);
    hi_out.assign(ns, 0);

    std::vector<int> damping_groups;
    find_degenerate_groups(ns, eval_at_k, damping_groups);

    auto gbegin = 0u;
    for (const auto ndeg: damping_groups) {
        const auto gend = gbegin + static_cast<unsigned int>(ndeg);
        auto is = gbegin;
        while (is < gend) {
            const auto anchor = in_kayser(eval_at_k[is]);
            auto hi = is + 1;
            while (hi < gend && std::abs(in_kayser(eval_at_k[hi]) - anchor) < tol_cm) ++hi;
            for (auto k = is; k < hi; ++k) {
                lo_out[k] = static_cast<int>(is);
                hi_out[k] = static_cast<int>(hi);
            }
            is = hi;
        }
        gbegin = gend;
    }
}
} // namespace PHON_NS
