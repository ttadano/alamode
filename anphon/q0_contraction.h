/*
 q0_contraction.h

 One multi-threaded sweep over the reciprocal-space quartic coefficients

   v4[ik*nk_dense + jk][a*ns + b][c*ns + d]

 (ik: coarse irreducible k, jk: dense k, row (a,b): modes at k and -k, column
 (c,d): modes at jk and -jk) that produces every contraction with the
 Gamma-point displacement q0 needed by the SCPH/QHA structural optimization:

   v3_renorm[jk][b][c,d] = v3_with_umn[jk][b][c,d]
                           + sum_a v4[g*nk_dense + jk][a,b][c,d] q0[a]      (all dense jk)
   q4_q0[ik][a][b]       = sum_{c,d} v4[ik*nk_dense + jg][a,b][c,d] q0[c] q0[d]  (all coarse ik)

 with g the coarse irreducible index of Gamma and jg its dense index. The v1,
 v2 and v0 renormalizations by q0 are O(ns^2) contractions of q4_q0
 (Relaxation::renormalize_v1_from_q0 etc.), so v4 is read exactly once per
 structure step: nk_dense + nk_irred_coarse - 1 slices, each element loaded
 from DRAM a single time. (Before this kernel the four renormalizations swept
 v4 nk_dense + nk_irred_coarse + 2 times, three of them single-threaded.)

 Distributed layout: v4 may be split over MPI ranks by contiguous ranges of
 units u = ik_prod * ns + a (the ns rows (a, b) of one slice). Options::unit_begin /
 unit_end restrict the sweep to the rows this process owns; the outputs of the
 ranks then add up to the single-process result (v3_with_umn is seeded on one
 rank only, Options::seed_v3_with_umn; q4_q0 is reset on every rank and filled
 for the owned units only). The kernel itself stays MPI-free.

 Layout of the sweep: a thread's owner unit is one (jk, b) pair, i.e. the rows
 (a, b), a over the owned range, of the slice (g, jk). The columns are processed in tiles
 so that the accumulator segment of v3_renorm[jk][b] and the segment of
 w[c*ns+d] = q0[c] q0[d] stay in cache across the a-loop, and the row segment
 is used for both outputs from a single load. Every unit writes its own
 v3_renorm row and (at jk == jg) its own column of q4_q0, so there is no
 reduction across threads and the result is independent of the thread count.

 Header-only; OpenMP only (no MPI, no Eigen), so that test_q0_contraction.cpp
 can exercise it against the naive loops.
*/

#pragma once

#include <algorithm>
#include <complex>
#include <cstddef>
#include <limits>
#include <vector>

namespace PHON_NS::q0_contraction
{

// Column tile: 16384 complex = 256 KB of v4 per row segment; the matching
// accumulator segment (256 KB) and w segment (128 KB) fit in L2. Measured on a
// 64-core cluster node (ns = 256..320): 216-222 GB/s at 256 KB against 170-183 GB/s
// at 64 KB and 191-217 GB/s at 1-4 MB; the M2 Max is insensitive to the tile.
constexpr std::size_t default_tile = 16384;

struct Options
{
    std::size_t tile = default_tile;
    // owned units [unit_begin, unit_end) of u = ik_prod * ns + a; the default is everything
    std::size_t unit_begin = 0;
    std::size_t unit_end = std::numeric_limits<std::size_t>::max();
    // false on the ranks that must not add v3_with_umn (it is added exactly once)
    bool seed_v3_with_umn = true;
};

// v4 may be nullptr (no quartic terms): v3_renorm = v3_with_umn (or zero), q4_q0 = 0.
// The same fast path is taken when every q0 component is exactly zero.
inline void contract_v4_with_q0(const std::size_t ns, const std::size_t nk_dense, const std::size_t nk_irred_coarse,
                                const std::size_t ik_gamma_irred, const std::size_t jk_gamma_dense,
                                const std::complex<double> *const *const *v4, const double *q0,
                                const std::complex<double> *const *const *v3_with_umn,
                                std::complex<double> ***v3_renorm, std::complex<double> ***q4_q0,
                                const Options &opt = Options{})
{
    const std::size_t ns2 = ns * ns;
    const std::size_t tile = opt.tile;
    constexpr std::complex<double> czero(0.0, 0.0);

    const std::size_t nunits = nk_irred_coarse * nk_dense * ns;
    const std::size_t unit_begin = std::min(opt.unit_begin, nunits);
    const std::size_t unit_end = std::min(opt.unit_end, nunits);

    // owned rows a of slice ik_prod: [a0, a1)
    auto owned_a_range = [&](const std::size_t ik_prod, std::size_t &a0, std::size_t &a1) {
        const std::size_t base = ik_prod * ns;
        a0 = unit_begin > base ? std::min(unit_begin - base, ns) : 0;
        a1 = unit_end > base ? std::min(unit_end - base, ns) : 0;
        if (a1 < a0) {
            a1 = a0;
        }
    };

    auto seed_row = [&](const std::size_t jk, const std::size_t b) {
        if (opt.seed_v3_with_umn) {
            std::copy(v3_with_umn[jk][b], v3_with_umn[jk][b] + ns2, v3_renorm[jk][b]);
        } else {
            std::fill(v3_renorm[jk][b], v3_renorm[jk][b] + ns2, czero);
        }
    };

    // q4_q0 is rebuilt from scratch on every call: unowned entries stay zero
    for (std::size_t ik = 0; ik < nk_irred_coarse; ++ik) {
        for (std::size_t a = 0; a < ns; ++a) {
            std::fill(q4_q0[ik][a], q4_q0[ik][a] + ns, czero);
        }
    }

    auto q0_is_zero = true;
    for (std::size_t a = 0; a < ns; ++a) {
        if (q0[a] != 0.0) {
            q0_is_zero = false;
            break;
        }
    }

    if (v4 == nullptr || q0_is_zero) {
        for (std::size_t jk = 0; jk < nk_dense; ++jk) {
            for (std::size_t b = 0; b < ns; ++b) {
                seed_row(jk, b);
            }
        }
        return;
    }

    // w[c*ns + d] = q0[c] q0[d] (real)
    std::vector<double> w(ns2);
    for (std::size_t c = 0; c < ns; ++c) {
        for (std::size_t d = 0; d < ns; ++d) {
            w[c * ns + d] = q0[c] * q0[d];
        }
    }

    // ---- Row-Gamma slices (g, jk), all dense jk: the v3 correction, plus
    //      q4_q0[g] from the slice (g, jg). The inner loops run on the
    //      (re, im) doubles of the rows so that they vectorize as plain
    //      real axpy / dot operations (std::complex<double> is layout-
    //      compatible with double[2]).
    const long nthread_units = static_cast<long>(nk_dense * ns);

#pragma omp parallel
    {
        std::vector<std::complex<double>> q4_local(ns);

        // dynamic: the units of the slice (g, jg) do twice the work of the others
#pragma omp for schedule(dynamic, 1)
        for (long iu = 0; iu < nthread_units; ++iu) {
            const std::size_t jk = static_cast<std::size_t>(iu) / ns;
            const std::size_t b = static_cast<std::size_t>(iu) % ns;
            const std::size_t ik_prod = ik_gamma_irred * nk_dense + jk;
            const auto slice = v4[ik_prod];
            const bool at_gamma = (jk == jk_gamma_dense);
            std::size_t a0, a1;
            owned_a_range(ik_prod, a0, a1);

            std::complex<double> *acc_row = v3_renorm[jk][b];
            seed_row(jk, b);
            if (at_gamma) {
                std::fill(q4_local.begin(), q4_local.end(), czero);
            }

            for (std::size_t c0 = 0; c0 < ns2; c0 += tile) {
                const std::size_t len = std::min(tile, ns2 - c0);
                auto *acc = reinterpret_cast<double *>(acc_row + c0);
                const double *wt = w.data() + c0;

                for (std::size_t a = a0; a < a1; ++a) {
                    const double qa = q0[a];
                    const auto *row = reinterpret_cast<const double *>(slice[a * ns + b] + c0);

                    if (at_gamma) {
                        double dot_re = 0.0, dot_im = 0.0;
                        for (std::size_t c = 0; c < len; ++c) {
                            const double vr = row[2 * c];
                            const double vi = row[2 * c + 1];
                            acc[2 * c] += qa * vr;
                            acc[2 * c + 1] += qa * vi;
                            dot_re += vr * wt[c];
                            dot_im += vi * wt[c];
                        }
                        q4_local[a] += std::complex<double>(dot_re, dot_im);
                    } else {
                        for (std::size_t c = 0; c < 2 * len; ++c) {
                            acc[c] += qa * row[c];
                        }
                    }
                }
            }

            if (at_gamma) {
                for (std::size_t a = a0; a < a1; ++a) {
                    q4_q0[ik_gamma_irred][a][b] = q4_local[a];
                }
            }
        }
    }

    // ---- Column-Gamma slices (ik, jg), ik != g: q4_q0[ik] only, owned rows only.
    for (std::size_t ik = 0; ik < nk_irred_coarse; ++ik) {
        if (ik == ik_gamma_irred) {
            continue;
        }
        const std::size_t ik_prod = ik * nk_dense + jk_gamma_dense;
        const auto slice = v4[ik_prod];
        std::size_t a0, a1;
        owned_a_range(ik_prod, a0, a1);
        const long r_begin = static_cast<long>(a0 * ns);
        const long r_end = static_cast<long>(a1 * ns);

#pragma omp parallel for schedule(static)
        for (long r = r_begin; r < r_end; ++r) {
            const auto *row = reinterpret_cast<const double *>(slice[r]);
            double dot_re = 0.0, dot_im = 0.0;
            for (std::size_t c = 0; c < ns2; ++c) {
                dot_re += row[2 * c] * w[c];
                dot_im += row[2 * c + 1] * w[c];
            }
            q4_q0[ik][static_cast<std::size_t>(r) / ns][static_cast<std::size_t>(r) % ns] =
                std::complex<double>(dot_re, dot_im);
        }
    }
}

} // namespace PHON_NS::q0_contraction
