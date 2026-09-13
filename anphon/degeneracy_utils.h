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

// Transport blocks use the Peierls limit Tr(P V^u P V^v P) for nearly
// degenerate modes. TOL_CM is an empirical threshold and can merge
// numerically resolvable splittings. For splitting dw and summed HWHM G,
// this overestimates the pair weight by 1 + (dw/G)^2; compute_kappa
// reports the worst dw/G. Blocks subdivide the damping-averaging groups
// to keep lifetimes constant; damping values do not determine membership.
inline double transport_block_tol_cm()
{
    return 1.0e-6;
}

// Return block bounds [lo, hi) for each branch. Subdivide the frequency
// groups from find_degenerate_groups (1e-7 Ry) using distance from each
// block's first frequency, preventing chains of close modes from merging.
// Staying within damping-averaging groups ensures constant block lifetimes.
inline void transport_block_bounds(const unsigned int ns, const double *eval_at_k, const double tol_cm,
                                   std::vector<int> &lo_out, std::vector<int> &hi_out)
{
    lo_out.assign(ns, 0);
    hi_out.assign(ns, 0);

    std::vector<int> freq_groups;
    find_degenerate_groups(ns, eval_at_k, freq_groups);

    auto gbegin = 0u;
    for (const auto ndeg: freq_groups) {
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
