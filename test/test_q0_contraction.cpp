/*
 test_q0_contraction.cpp

 Unit test for anphon/q0_contraction.h: the single sweep over the quartic
 coefficients v4 that yields the v3 renormalization by q0 and the quartic
 contraction q4_q0[ik][a][b] = sum_{c,d} v4[ik][a,b][c,d] q0[c] q0[d] feeding
 the v2 / v1 / v0 renormalizations (Relaxation::renormalize_*_from_q0).

 The kernel is compared with naive re-implementations of the four loops it
 replaces (the pre-fusion renormalize_v3/v2/v1/v0_from_q0 bodies) on random
 complex data with distinct values in every slice, for Gamma at irreducible /
 dense index 0 and at nonzero indices, with several column tiles per row,
 after a second call with a different q0 (outputs must not accumulate), for
 the nullptr / all-zero fast paths, for 1 versus several OpenMP threads
 (results must be identical: every owner unit reduces sequentially), and for
 the distributed layout (contiguous unit ranges whose partial results must add
 up to the full sweep, with v3_with_umn seeded on one part only).

 Built by the anphon CMake project as `test_q0_contraction`. Exits 0 on
 success and prints the first failing check otherwise.
*/

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <random>
#include <string>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "ndarray.h"
#include "q0_contraction.h"

using cplx = std::complex<double>;
using PHON_NS::q0_contraction::contract_v4_with_q0;
using PHON_NS::q0_contraction::Options;

namespace
{
int nfail = 0;

void check(const bool ok, const std::string &what)
{
    if (!ok) {
        ++nfail;
        std::printf("FAIL: %s\n", what.c_str());
    }
}

bool close(const cplx &a, const cplx &b, const double tol)
{
    return std::abs(a - b) <= tol * (1.0 + std::abs(b));
}

struct Case
{
    std::size_t ns, nk_dense, nk_irred, g, jg;
    std::size_t tile;
    unsigned seed;
};

// Naive reference: the pre-fusion loops (relaxation.cpp before the sweep).
struct Reference
{
    NDArray<cplx, 3> v3_renorm;   // [nk_dense][ns][ns^2]
    NDArray<cplx, 3> q4;          // [nk_irred][ns][ns]
    std::vector<cplx> v1_quartic; // sum_{b,c,d} v4[g][a,b][c,d] q0[b] q0[c] q0[d]
    cplx v0_quartic;              // sum_{a,b,c,d} v4[g][b,a][c,d] q0[a] q0[b] q0[c] q0[d]
};

Reference compute_reference(const Case &cs, const NDArray<cplx, 3> &v4, const std::vector<double> &q0,
                            const NDArray<cplx, 3> &v3_with_umn)
{
    const auto ns = cs.ns;
    const auto ns2 = ns * ns;
    Reference ref;
    ref.v3_renorm.resize(cs.nk_dense, ns, ns2);
    ref.q4.resize(cs.nk_irred, ns, ns);
    ref.v1_quartic.assign(ns, cplx(0.0, 0.0));
    ref.v0_quartic = cplx(0.0, 0.0);

    // renormalize_v3_from_q0
    for (std::size_t jk = 0; jk < cs.nk_dense; ++jk) {
        for (std::size_t b = 0; b < ns; ++b) {
            for (std::size_t cd = 0; cd < ns2; ++cd) {
                cplx sum = v3_with_umn[jk][b][cd];
                for (std::size_t a = 0; a < ns; ++a) {
                    sum += v4[cs.g * cs.nk_dense + jk][a * ns + b][cd] * q0[a];
                }
                ref.v3_renorm[jk][b][cd] = sum;
            }
        }
    }
    // quartic part of renormalize_v2_from_q0 (before rotation / symmetrization)
    for (std::size_t ik = 0; ik < cs.nk_irred; ++ik) {
        for (std::size_t a = 0; a < ns; ++a) {
            for (std::size_t b = 0; b < ns; ++b) {
                cplx sum(0.0, 0.0);
                for (std::size_t c = 0; c < ns; ++c) {
                    for (std::size_t d = 0; d < ns; ++d) {
                        sum += v4[ik * cs.nk_dense + cs.jg][a * ns + b][c * ns + d] * q0[c] * q0[d];
                    }
                }
                ref.q4[ik][a][b] = sum;
            }
        }
    }
    // quartic part of renormalize_v1_from_q0
    for (std::size_t a = 0; a < ns; ++a) {
        for (std::size_t b = 0; b < ns; ++b) {
            for (std::size_t c = 0; c < ns; ++c) {
                for (std::size_t d = 0; d < ns; ++d) {
                    ref.v1_quartic[a] += v4[cs.g * cs.nk_dense + cs.jg][a * ns + b][c * ns + d] * q0[b] * q0[c] * q0[d];
                }
            }
        }
    }
    // quartic part of renormalize_v0_from_q0 (note the transposed row index)
    for (std::size_t a = 0; a < ns; ++a) {
        for (std::size_t b = 0; b < ns; ++b) {
            for (std::size_t c = 0; c < ns; ++c) {
                for (std::size_t d = 0; d < ns; ++d) {
                    ref.v0_quartic +=
                        v4[cs.g * cs.nk_dense + cs.jg][b * ns + a][c * ns + d] * q0[a] * q0[b] * q0[c] * q0[d];
                }
            }
        }
    }
    return ref;
}

void fill_random(NDArray<cplx, 3> &arr, std::mt19937 &rng)
{
    std::uniform_real_distribution<double> dist(-1.0, 1.0);
    const auto n = arr.size();
    auto *p = &arr[0][0][0];
    for (std::size_t i = 0; i < n; ++i) {
        p[i] = cplx(dist(rng), dist(rng));
    }
}

void run_case(const Case &cs, const std::string &label)
{
    const auto ns = cs.ns;
    const auto ns2 = ns * ns;
    const double tol = 1.0e-13;
    std::mt19937 rng(cs.seed);
    std::uniform_real_distribution<double> dist(-1.0, 1.0);

    NDArray<cplx, 3> v4(cs.nk_irred * cs.nk_dense, ns2, ns2);
    NDArray<cplx, 3> v3_with_umn(cs.nk_dense, ns, ns2);
    fill_random(v4, rng);
    fill_random(v3_with_umn, rng);

    std::vector<double> q0(ns);
    for (auto &q: q0) {
        q = dist(rng);
    }
    q0[ns / 2] = 0.0; // one exactly-zero component must not trigger the fast path

    NDArray<cplx, 3> v3_renorm(cs.nk_dense, ns, ns2);
    NDArray<cplx, 3> q4(cs.nk_irred, ns, ns);

    auto compare = [&](const Reference &ref, const std::string &tag) {
        for (std::size_t jk = 0; jk < cs.nk_dense; ++jk) {
            for (std::size_t b = 0; b < ns; ++b) {
                for (std::size_t cd = 0; cd < ns2; ++cd) {
                    if (!close(v3_renorm[jk][b][cd], ref.v3_renorm[jk][b][cd], tol)) {
                        check(false,
                              label + tag + ": v3_renorm mismatch at jk=" + std::to_string(jk) +
                                  " b=" + std::to_string(b) + " cd=" + std::to_string(cd));
                        return;
                    }
                }
            }
        }
        for (std::size_t ik = 0; ik < cs.nk_irred; ++ik) {
            for (std::size_t a = 0; a < ns; ++a) {
                for (std::size_t b = 0; b < ns; ++b) {
                    if (!close(q4[ik][a][b], ref.q4[ik][a][b], tol)) {
                        check(false,
                              label + tag + ": q4 mismatch at ik=" + std::to_string(ik) + " a=" + std::to_string(a) +
                                  " b=" + std::to_string(b));
                        return;
                    }
                }
            }
        }
        // the derived contractions used by renormalize_v1_from_q0 / renormalize_v0_from_q0
        for (std::size_t a = 0; a < ns; ++a) {
            cplx v1(0.0, 0.0);
            for (std::size_t b = 0; b < ns; ++b) {
                v1 += q4[cs.g][a][b] * q0[b];
            }
            if (!close(v1, ref.v1_quartic[a], tol)) {
                check(false, label + tag + ": v1 quartic mismatch at a=" + std::to_string(a));
                return;
            }
        }
        cplx v0(0.0, 0.0);
        for (std::size_t a = 0; a < ns; ++a) {
            for (std::size_t b = 0; b < ns; ++b) {
                v0 += q4[cs.g][b][a] * q0[a] * q0[b];
            }
        }
        check(close(v0, ref.v0_quartic, tol), label + tag + ": v0 quartic mismatch");
    };

    // 1. against the naive loops
    contract_v4_with_q0(ns,
                        cs.nk_dense,
                        cs.nk_irred,
                        cs.g,
                        cs.jg,
                        v4,
                        q0.data(),
                        v3_with_umn,
                        v3_renorm,
                        q4,
                        Options{cs.tile});
    compare(compute_reference(cs, v4, q0, v3_with_umn), " [first call]");

    // 2. a second call with a different q0 must not accumulate on the previous outputs
    for (auto &q: q0) {
        q = dist(rng);
    }
    contract_v4_with_q0(ns,
                        cs.nk_dense,
                        cs.nk_irred,
                        cs.g,
                        cs.jg,
                        v4,
                        q0.data(),
                        v3_with_umn,
                        v3_renorm,
                        q4,
                        Options{cs.tile});
    compare(compute_reference(cs, v4, q0, v3_with_umn), " [second call]");

#ifdef _OPENMP
    // 3. thread-count independence (bitwise: each owner unit reduces sequentially)
    NDArray<cplx, 3> v3_renorm_1(cs.nk_dense, ns, ns2);
    NDArray<cplx, 3> q4_1(cs.nk_irred, ns, ns);
    const int nthreads_saved = omp_get_max_threads();
    omp_set_num_threads(1);
    contract_v4_with_q0(ns,
                        cs.nk_dense,
                        cs.nk_irred,
                        cs.g,
                        cs.jg,
                        v4,
                        q0.data(),
                        v3_with_umn,
                        v3_renorm_1,
                        q4_1,
                        Options{cs.tile});
    omp_set_num_threads(std::max(nthreads_saved, 4));
    contract_v4_with_q0(ns,
                        cs.nk_dense,
                        cs.nk_irred,
                        cs.g,
                        cs.jg,
                        v4,
                        q0.data(),
                        v3_with_umn,
                        v3_renorm,
                        q4,
                        Options{cs.tile});
    omp_set_num_threads(nthreads_saved);
    // bitwise comparison of the double representations (a numerical == would
    // also accept e.g. differently signed zeros)
    const bool same = std::memcmp(&v3_renorm[0][0][0], &v3_renorm_1[0][0][0], v3_renorm.size() * sizeof(cplx)) == 0 &&
                      std::memcmp(&q4[0][0][0], &q4_1[0][0][0], q4.size() * sizeof(cplx)) == 0;
    check(same, label + ": 1-thread and multi-thread results differ");
#endif

    // 4. fast paths: nullptr v4 and exactly-zero q0 give v3_with_umn and zeros
    auto expect_trivial = [&](const std::string &tag) {
        for (std::size_t jk = 0; jk < cs.nk_dense; ++jk) {
            for (std::size_t b = 0; b < ns; ++b) {
                for (std::size_t cd = 0; cd < ns2; ++cd) {
                    if (v3_renorm[jk][b][cd] != v3_with_umn[jk][b][cd]) {
                        check(false, label + tag + ": v3_renorm != v3_with_umn");
                        return;
                    }
                }
            }
        }
        for (std::size_t i = 0; i < q4.size(); ++i) {
            if ((&q4[0][0][0])[i] != cplx(0.0, 0.0)) {
                check(false, label + tag + ": q4 != 0");
                return;
            }
        }
    };
    contract_v4_with_q0(ns,
                        cs.nk_dense,
                        cs.nk_irred,
                        cs.g,
                        cs.jg,
                        nullptr,
                        q0.data(),
                        v3_with_umn,
                        v3_renorm,
                        q4,
                        Options{cs.tile});
    expect_trivial(" [v4 == nullptr]");
    // poison the outputs, then the zero-q0 path must overwrite them
    fill_random(v3_renorm, rng);
    fill_random(q4, rng);
    std::vector<double> q0_zero(ns, 0.0);
    contract_v4_with_q0(ns,
                        cs.nk_dense,
                        cs.nk_irred,
                        cs.g,
                        cs.jg,
                        v4,
                        q0_zero.data(),
                        v3_with_umn,
                        v3_renorm,
                        q4,
                        Options{cs.tile});
    expect_trivial(" [q0 == 0]");

    // 5. distributed layout: contiguous unit ranges [u0, u1) of u = ik_prod * ns + a, the
    //    outputs of the parts must add up to the full sweep (v3_with_umn seeded on the first
    //    part only, q4 reset on every part)
    for (auto &q: q0) {
        q = dist(rng);
    }
    NDArray<cplx, 3> v3_full(cs.nk_dense, ns, ns2);
    NDArray<cplx, 3> q4_full(cs.nk_irred, ns, ns);
    contract_v4_with_q0(ns,
                        cs.nk_dense,
                        cs.nk_irred,
                        cs.g,
                        cs.jg,
                        v4,
                        q0.data(),
                        v3_with_umn,
                        v3_full,
                        q4_full,
                        Options{cs.tile});
    const std::size_t nunits = cs.nk_irred * cs.nk_dense * ns;
    for (const int nparts: {1, 2, 3, 5}) {
        // boundaries at random positions (parts may be empty or cut through a slice)
        std::vector<std::size_t> bounds{0};
        std::uniform_int_distribution<std::size_t> pick(0, nunits);
        for (int ip = 1; ip < nparts; ++ip) {
            bounds.push_back(pick(rng));
        }
        bounds.push_back(nunits);
        std::sort(bounds.begin(), bounds.end());
        NDArray<cplx, 3> v3_sum(cs.nk_dense, ns, ns2);
        NDArray<cplx, 3> q4_sum(cs.nk_irred, ns, ns);
        std::fill(&v3_sum[0][0][0], &v3_sum[0][0][0] + v3_sum.size(), cplx(0.0, 0.0));
        std::fill(&q4_sum[0][0][0], &q4_sum[0][0][0] + q4_sum.size(), cplx(0.0, 0.0));
        for (int ip = 0; ip < nparts; ++ip) {
            Options opt{cs.tile};
            opt.unit_begin = bounds[ip];
            opt.unit_end = bounds[ip + 1];
            opt.seed_v3_with_umn = (ip == 0);
            fill_random(v3_renorm, rng); // poisoned outputs must be overwritten
            fill_random(q4, rng);
            contract_v4_with_q0(ns,
                                cs.nk_dense,
                                cs.nk_irred,
                                cs.g,
                                cs.jg,
                                v4,
                                q0.data(),
                                v3_with_umn,
                                v3_renorm,
                                q4,
                                opt);
            for (std::size_t i = 0; i < v3_sum.size(); ++i) {
                (&v3_sum[0][0][0])[i] += (&v3_renorm[0][0][0])[i];
            }
            for (std::size_t i = 0; i < q4_sum.size(); ++i) {
                (&q4_sum[0][0][0])[i] += (&q4[0][0][0])[i];
            }
        }
        bool ok = true;
        for (std::size_t i = 0; i < v3_sum.size() && ok; ++i) {
            ok = close((&v3_sum[0][0][0])[i], (&v3_full[0][0][0])[i], tol);
        }
        for (std::size_t i = 0; i < q4_sum.size() && ok; ++i) {
            ok = close((&q4_sum[0][0][0])[i], (&q4_full[0][0][0])[i], tol);
        }
        check(ok, label + " [partition into " + std::to_string(nparts) + " unit ranges]: partial sums != full sweep");
    }
    // deterministic partition: the first part holds no Gamma-family row of the row family
    // (it ends inside the very first slice), the cuts go through the Gamma-Gamma slice and
    // through a column-Gamma slice, one part is a whole slice
    if (cs.nk_irred > 1 && cs.nk_dense > 1) {
        const std::size_t u_gg = (cs.g * cs.nk_dense + cs.jg) * ns + ns / 2; // inside (g, jg)
        const std::size_t ik_other = cs.g == 0 ? 1 : 0;
        const std::size_t u_cg = (ik_other * cs.nk_dense + cs.jg) * ns + ns / 3 + 1; // inside (ik_other, jg)
        std::vector<std::size_t> bounds{0, ns / 2, ns, u_gg, u_cg, nunits}; // [ns/2, ns) and [0,ns/2): no Gamma rows
        std::sort(bounds.begin(), bounds.end());
        NDArray<cplx, 3> v3_sum(cs.nk_dense, ns, ns2);
        NDArray<cplx, 3> q4_sum(cs.nk_irred, ns, ns);
        std::fill(&v3_sum[0][0][0], &v3_sum[0][0][0] + v3_sum.size(), cplx(0.0, 0.0));
        std::fill(&q4_sum[0][0][0], &q4_sum[0][0][0] + q4_sum.size(), cplx(0.0, 0.0));
        for (std::size_t ip = 0; ip + 1 < bounds.size(); ++ip) {
            Options opt{cs.tile};
            opt.unit_begin = bounds[ip];
            opt.unit_end = bounds[ip + 1];
            opt.seed_v3_with_umn = (ip == 0);
            fill_random(v3_renorm, rng);
            fill_random(q4, rng);
            contract_v4_with_q0(ns,
                                cs.nk_dense,
                                cs.nk_irred,
                                cs.g,
                                cs.jg,
                                v4,
                                q0.data(),
                                v3_with_umn,
                                v3_renorm,
                                q4,
                                opt);
#ifdef _OPENMP
            // the restricted sweep must not depend on the thread count (bitwise)
            {
                NDArray<cplx, 3> v3_1(cs.nk_dense, ns, ns2);
                NDArray<cplx, 3> q4_1(cs.nk_irred, ns, ns);
                const int nthreads_saved = omp_get_max_threads();
                omp_set_num_threads(1);
                contract_v4_with_q0(ns,
                                    cs.nk_dense,
                                    cs.nk_irred,
                                    cs.g,
                                    cs.jg,
                                    v4,
                                    q0.data(),
                                    v3_with_umn,
                                    v3_1,
                                    q4_1,
                                    opt);
                omp_set_num_threads(nthreads_saved);
                check(std::memcmp(&v3_1[0][0][0], &v3_renorm[0][0][0], v3_1.size() * sizeof(cplx)) == 0 &&
                          std::memcmp(&q4_1[0][0][0], &q4[0][0][0], q4_1.size() * sizeof(cplx)) == 0,
                      label + " [deterministic partition, part " + std::to_string(ip) + "]: thread-count dependence");
            }
#endif
            for (std::size_t i = 0; i < v3_sum.size(); ++i) {
                (&v3_sum[0][0][0])[i] += (&v3_renorm[0][0][0])[i];
            }
            for (std::size_t i = 0; i < q4_sum.size(); ++i) {
                (&q4_sum[0][0][0])[i] += (&q4[0][0][0])[i];
            }
        }
        bool ok = true;
        for (std::size_t i = 0; i < v3_sum.size() && ok; ++i) {
            ok = close((&v3_sum[0][0][0])[i], (&v3_full[0][0][0])[i], tol);
        }
        for (std::size_t i = 0; i < q4_sum.size() && ok; ++i) {
            ok = close((&q4_sum[0][0][0])[i], (&q4_full[0][0][0])[i], tol);
        }
        check(ok, label + " [deterministic partition through both Gamma families]: partial sums != full sweep");
    }
    // fast paths without seeding must give zeros
    {
        Options opt{cs.tile};
        opt.seed_v3_with_umn = false;
        fill_random(v3_renorm, rng);
        fill_random(q4, rng);
        contract_v4_with_q0(ns,
                            cs.nk_dense,
                            cs.nk_irred,
                            cs.g,
                            cs.jg,
                            nullptr,
                            q0.data(),
                            v3_with_umn,
                            v3_renorm,
                            q4,
                            opt);
        bool zero = true;
        for (std::size_t i = 0; i < v3_renorm.size() && zero; ++i) {
            zero = (&v3_renorm[0][0][0])[i] == cplx(0.0, 0.0);
        }
        for (std::size_t i = 0; i < q4.size() && zero; ++i) {
            zero = (&q4[0][0][0])[i] == cplx(0.0, 0.0);
        }
        check(zero, label + " [v4 == nullptr, no seed]: outputs not zero");
    }
    // a part that owns nothing and does not seed must return zeros
    {
        Options opt{cs.tile};
        opt.unit_begin = nunits;
        opt.unit_end = nunits;
        opt.seed_v3_with_umn = false;
        fill_random(v3_renorm, rng);
        fill_random(q4, rng);
        contract_v4_with_q0(ns, cs.nk_dense, cs.nk_irred, cs.g, cs.jg, v4, q0.data(), v3_with_umn, v3_renorm, q4, opt);
        bool zero = true;
        for (std::size_t i = 0; i < v3_renorm.size() && zero; ++i) {
            zero = (&v3_renorm[0][0][0])[i] == cplx(0.0, 0.0);
        }
        for (std::size_t i = 0; i < q4.size() && zero; ++i) {
            zero = (&q4[0][0][0])[i] == cplx(0.0, 0.0);
        }
        check(zero, label + " [empty unit range without seed]: outputs not zero");
    }
}
} // namespace

int main()
{
    // ns^2 = 49 columns with tile 16: four tiles, the last one partial
    run_case({7, 4, 3, 0, 0, 16, 11u}, "gamma-first");
    // Gamma at irreducible index 1 and dense index 2; ns^2 = 25 with tile 8
    run_case({5, 3, 2, 1, 2, 8, 23u}, "gamma-shifted");
    // Gamma-only (the large-cell configuration): one slice serves both families
    run_case({6, 1, 1, 0, 0, 4096, 37u}, "gamma-only");
    // rows shorter than one tile, several dense k and coarse irreducible k
    run_case({4, 5, 4, 2, 3, 4096, 41u}, "small-rows");

    if (nfail == 0) {
        std::printf("test_q0_contraction: all checks passed\n");
        return EXIT_SUCCESS;
    }
    std::printf("test_q0_contraction: %d check(s) failed\n", nfail);
    return EXIT_FAILURE;
}
