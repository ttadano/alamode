// Unit test for anphon/v4_index_transform.h: the rotating complex GEMM
// out[col][x] = alpha * sum_a E[x][a] in[a][col] against the naive loops, for
// ncols = ns^2 (band builder) and ns^3 (k-point builder), a non-unit complex
// prefactor and several OpenMP thread counts (uneven slab boundaries, and
// empty slabs when the team is larger than the column count). beta = 0 must
// overwrite whatever the output holds. The team size actually obtained is
// checked, so a build without OpenMP (or a thread limit) reports the missing
// parallel coverage instead of passing silently.
#include "v4_index_transform.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <random>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

int main()
{
    using cplx = std::complex<double>;
    std::mt19937 gen(12345);
    std::uniform_real_distribution<double> uni(-1.0, 1.0);
    int nfail = 0;
#ifdef _OPENMP
    omp_set_dynamic(0);
#endif
    for (const std::size_t ns: {std::size_t(1), std::size_t(5), std::size_t(7)}) {
        for (const std::size_t ncols: {ns * ns, ns * ns * ns}) {
            std::vector<cplx> E(ns * ns), in(ns * ncols), out(ncols * ns), ref(ncols * ns);
            for (auto &x: E) x = cplx(uni(gen), uni(gen));
            for (auto &x: in) x = cplx(uni(gen), uni(gen));
            const cplx alpha(0.37, -0.21);
            for (std::size_t c = 0; c < ncols; ++c) {
                for (std::size_t x = 0; x < ns; ++x) {
                    cplx sum(0.0, 0.0);
                    for (std::size_t a = 0; a < ns; ++a) sum += E[x * ns + a] * in[a * ncols + c];
                    ref[c * ns + x] = alpha * sum;
                }
            }
            for (const int nt: {1, 3, 8}) {
                int nt_actual = 1;
#ifdef _OPENMP
                omp_set_num_threads(nt);
#pragma omp parallel
                {
#pragma omp single
                    nt_actual = omp_get_num_threads();
                }
#endif
                std::fill(out.begin(), out.end(), cplx(7.0, -7.0));
                PHON_NS::v4_index_transform::transform_index_gemm(E.data(), in.data(), out.data(), ns, ncols, alpha);
                double maxdiff = 0.0;
                for (std::size_t i = 0; i < out.size(); ++i) maxdiff = std::max(maxdiff, std::abs(out[i] - ref[i]));
                const bool ok = maxdiff < 1e-13;
                const bool covered = (nt_actual == nt);
                std::printf("ns=%zu ncols=%zu threads=%d (got %d) max|diff|=%.2e %s%s\n", ns, ncols, nt, nt_actual,
                            maxdiff, ok ? "ok" : "FAILED", covered ? "" : " [team size not obtained: slab coverage missing]");
                if (!ok || !covered) ++nfail;
            }
        }
    }
    std::printf(nfail == 0 ? "test_v4_index_transform: all passed\n" : "test_v4_index_transform: %d FAILED\n", nfail);
    return nfail == 0 ? 0 : 1;
}
