/*
 phonon_thermodynamics.cpp

 Copyright (c) 2014, 2015, 2016 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#include "mpi_common.h"
#include "thermodynamics.h"
#include "anharmonic_core.h"
#include "constants.h"
#include "dynamical.h"
#include "kpoint.h"
#include "mathfunctions.h"
#include "memory.h"
#include "pointers.h"
#include "system.h"
#include <iostream>
#include <complex>

using namespace PHON_NS;

Thermodynamics::Thermodynamics(PHON *phon): Pointers(phon)
{
    T_to_Ryd = k_Boltzmann / Ryd;
    calc_FE_bubble = false;
    FE_bubble = nullptr;
}

Thermodynamics::~Thermodynamics()
{
    if (FE_bubble) {
        memory->deallocate(FE_bubble);
    }
};

void Thermodynamics::setup()
{
    MPI_Bcast(&classical, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD);
}

double Thermodynamics::Cv(const double omega,
                          const double temp_in) const
{
    if (std::abs(temp_in) < eps) return 0.0;

    const auto x = omega / (T_to_Ryd * temp_in);
    return k_Boltzmann * std::pow(x / (2.0 * sinh(0.5 * x)), 2);
}

double Thermodynamics::Cv_classical(const double omega,
                                    const double temp_in) const
{
    if (std::abs(temp_in) < eps) return 0.0;

    return k_Boltzmann;
}

double Thermodynamics::fB(const double omega,
                          const double temp_in) const
{
    if (std::abs(temp_in) < eps || omega < eps8) return 0.0;

    const auto x = omega / (T_to_Ryd * temp_in);
    return 1.0 / (std::exp(x) - 1.0);
}

double Thermodynamics::fC(const double omega,
                          const double temp_in) const
{
    if (std::abs(temp_in) < eps || omega < eps8) return 0.0;

    const auto x = omega / (T_to_Ryd * temp_in);
    return 1.0 / x;
}

double Thermodynamics::Cv_tot(const double temp_in,
                              const unsigned int nk_irred,
                              const unsigned int ns,
                              const std::vector<std::vector<KpointList>> &kp_irred,
                              double *weight_k_irred,
                              double **eval_in) const
{
    int i;
    unsigned int ik, is;
    double omega;
    auto ret = 0.0;

    const auto N = nk_irred * ns;
    int ik_irred;

    if (classical) {
#pragma omp parallel for private(ik_irred, ik, is, omega), reduction(+:ret)
        for (i = 0; i < N; ++i) {
            ik_irred = i / ns;
            ik = kp_irred[ik_irred][0].knum;
            is = i % ns;

            omega = eval_in[ik][is];

            if (omega < eps8) continue;

            ret += Cv_classical(omega, temp_in) * weight_k_irred[ik_irred];
        }
    } else {
#pragma omp parallel for private(ik_irred, ik, is, omega), reduction(+:ret)
        for (i = 0; i < N; ++i) {
            ik_irred = i / ns;
            ik = kp_irred[ik_irred][0].knum;
            is = i % ns;

            omega = eval_in[ik][is];

            if (omega < eps8) continue;

            ret += Cv(omega, temp_in) * weight_k_irred[ik_irred];
        }
    }

    return ret;
}

double Thermodynamics::internal_energy(const double temp_in,
                                       const unsigned int nk_irred,
                                       const unsigned int ns,
                                       const std::vector<std::vector<KpointList>> &kp_irred,
                                       double *weight_k_irred,
                                       double **eval_in) const
{
    int i;
    unsigned int ik, is;
    double omega;
    auto ret = 0.0;

    const auto N = nk_irred * ns;
    int ik_irred;

    if (classical) {
#pragma omp parallel for private(ik_irred, ik, is, omega), reduction(+:ret)
        for (i = 0; i < N; ++i) {
            ik_irred = i / ns;
            ik = kp_irred[ik_irred][0].knum;
            is = i % ns;
            omega = eval_in[ik][is];

            if (omega < eps8) continue;

            ret += T_to_Ryd * temp_in * weight_k_irred[ik_irred];
        }
        ret *= 2.0;

    } else {
#pragma omp parallel for private(ik_irred, ik, is, omega), reduction(+:ret)
        for (i = 0; i < N; ++i) {
            ik_irred = i / ns;
            ik = kp_irred[ik_irred][0].knum;
            is = i % ns;
            omega = eval_in[ik][is];

            if (omega < eps8) continue;

            ret += omega * coth_T(omega, temp_in) * weight_k_irred[ik_irred];
        }

    }
    return ret * 0.5;
}

double Thermodynamics::vibrational_entropy(const double temp_in,
                                           const unsigned int nk_irred,
                                           const unsigned int ns,
                                           const std::vector<std::vector<KpointList>> &kp_irred,
                                           double *weight_k_irred,
                                           double **eval_in) const
{
    int i;
    unsigned int ik, is;
    double omega, x;
    auto ret = 0.0;

    const auto N = nk_irred * ns;
    int ik_irred;

    if (classical) {
#pragma omp parallel for private(ik_irred, ik, is, omega, x), reduction(+:ret)
        for (i = 0; i < N; ++i) {
            ik_irred = i / ns;
            ik = kp_irred[ik_irred][0].knum;
            is = i % ns;
            omega = eval_in[ik][is];

            if (omega < eps8 || std::abs(temp_in) < eps) continue;

            x = omega / (temp_in * T_to_Ryd);
            ret += (std::log(x) - 1.0) * weight_k_irred[ik_irred];
        }

    } else {
#pragma omp parallel for private(ik_irred, ik, is, omega, x), reduction(+:ret)
        for (i = 0; i < N; ++i) {
            ik_irred = i / ns;
            ik = kp_irred[ik_irred][0].knum;
            is = i % ns;
            omega = eval_in[ik][is];

            if (omega < eps8 || std::abs(temp_in) < eps) continue;

            x = omega / (temp_in * T_to_Ryd);
            ret += (std::log(1.0 - std::exp(-x)) - x / (std::exp(x) - 1.0)) * weight_k_irred[ik_irred];
        }
    }
    return -k_Boltzmann * ret;
}

double Thermodynamics::free_energy_QHA(const double temp_in,
                                       const unsigned int nk_irred,
                                       const unsigned int ns,
                                       const std::vector<std::vector<KpointList>> &kp_irred,
                                       double *weight_k_irred,
                                       double **eval_in) const
{
    int i;
    unsigned int ik, is;
    double omega, x;
    auto ret = 0.0;

    const auto N = nk_irred * ns;
    int ik_irred;

    if (classical) {
#pragma omp parallel for private(ik_irred, ik, is, omega, x), reduction(+:ret)
        for (i = 0; i < N; ++i) {
            ik_irred = i / ns;
            ik = kp_irred[ik_irred][0].knum;
            is = i % ns;
            omega = eval_in[ik][is];

            if (omega < eps8) continue;

            if (std::abs(temp_in) > eps) {
                x = omega / (temp_in * T_to_Ryd);
                ret += std::log(x) * weight_k_irred[ik_irred];
            }
        }

        return temp_in * T_to_Ryd * ret;

    }
#pragma omp parallel for private(ik_irred, ik, is, omega, x), reduction(+:ret)
    for (i = 0; i < N; ++i) {
        ik_irred = i / ns;
        ik = kp_irred[ik_irred][0].knum;
        is = i % ns;
        omega = eval_in[ik][is];

        if (omega < eps8) continue;

        if (std::abs(temp_in) < eps) {
            ret += 0.5 * omega * weight_k_irred[ik_irred];
        } else {
            x = omega / (temp_in * T_to_Ryd);
            ret += (0.5 * x + std::log(1.0 - std::exp(-x))) * weight_k_irred[ik_irred];
        }
    }

    if (std::abs(temp_in) < eps) return ret;

    return temp_in * T_to_Ryd * ret;
}

double Thermodynamics::disp2_avg(const double T,
                                 const unsigned int ns1,
                                 const unsigned int ns2) const
{
    int i;
    double ret = 0.0;
    unsigned int ik, is;
    const auto nk = kpoint->nk;
    const auto ns = dynamical->neval;
    double omega;

    int N = nk * ns;

    if (classical) {
#pragma omp parallel for private(ik, is, omega), reduction(+:ret)
        for (i = 0; i < N; ++i) {
            ik = i / ns;
            is = i % ns;
            omega = dynamical->eval_phonon[ik][is];

            // Skip when omega is almost zero. 
            // (neglect divergent contributions from acoustic modes at gamma point)
            if (omega < eps8) continue;

            ret += real(dynamical->evec_phonon[ik][is][ns1]
                    * std::conj(dynamical->evec_phonon[ik][is][ns2]))
                * T * T_to_Ryd / (omega * omega);
        }
    } else {
#pragma omp parallel for private(ik, is, omega), reduction(+:ret)
        for (i = 0; i < N; ++i) {
            ik = i / ns;
            is = i % ns;
            omega = dynamical->eval_phonon[ik][is];

            // Skip when omega is almost zero. 
            // (neglect divergent contributions from acoustic modes at gamma point)
            if (omega < eps8) continue;

            ret += real(dynamical->evec_phonon[ik][is][ns1]
                    * std::conj(dynamical->evec_phonon[ik][is][ns2]))
                * (fB(omega, T) + 0.5) / omega;
        }
    }


    ret *= 1.0 / (static_cast<double>(nk)
        * std::sqrt(system->mass[system->map_p2s[ns1 / 3][0]]
            * system->mass[system->map_p2s[ns2 / 3][0]]));

    return ret;
}

double Thermodynamics::disp_corrfunc(const double T_in,
                                     const unsigned int ncrd1,
                                     const unsigned int ncrd2,
                                     const double cell_shift[3],
                                     const unsigned int nk,
                                     const unsigned int ns,
                                     double **xk_in,
                                     double **eval_in,
                                     std::complex<double> ***evec_in) const
{
    int i;
    int N = nk * ns;
    unsigned int ik, is;
    double omega;
    double phase;
    const std::complex<double> im(0.0, 1.0);
    double ret = 0.0;

    if (classical) {
#pragma omp parallel for private(ik, is, omega, phase), reduction(+:ret)
        for (i = 0; i < N; ++i) {
            ik = i / ns;
            is = i % ns;
            omega = eval_in[ik][is];

            if (omega < eps8) continue;

            phase = 2.0 * pi * (xk_in[ik][0] * cell_shift[0]
                + xk_in[ik][1] * cell_shift[1]
                + xk_in[ik][2] * cell_shift[2]);

            ret += real(std::conj(evec_in[ik][is][ncrd1])
                    * evec_in[ik][is][ncrd2]
                    * std::exp(phase))
                * T_in * T_to_Ryd / (omega * omega);

        }

    } else {
#pragma omp parallel for private(ik, is, omega, phase), reduction(+:ret)
        for (i = 0; i < N; ++i) {
            ik = i / ns;
            is = i % ns;
            omega = eval_in[ik][is];

            if (omega < eps8) continue;

            phase = 2.0 * pi * (xk_in[ik][0] * cell_shift[0]
                + xk_in[ik][1] * cell_shift[1]
                + xk_in[ik][2] * cell_shift[2]);

            ret += real(std::conj(evec_in[ik][is][ncrd1])
                    * evec_in[ik][is][ncrd2]
                    * std::exp(im * phase))
                * (fB(omega, T_in) + 0.5) / omega;
        }
    }

    ret *= 1.0 / (static_cast<double>(nk)
        * std::sqrt(system->mass[system->map_p2s[ncrd1 / 3][0]]
            * system->mass[system->map_p2s[ncrd2 / 3][0]]));

    return ret;
}


double Thermodynamics::coth_T(const double omega,
                              const double T) const
{
    // This function returns coth(hbar*omega/2*kB*T)

    // if T = 0.0 and omega > 0, coth(hbar*omega/(2*kB*T)) = 1.0
    if (T < eps) return 1.0;

    const auto x = omega / (T_to_Ryd * T * 2.0);
    return 1.0 + 2.0 / (std::exp(2.0 * x) - 1.0);
}


void Thermodynamics::compute_free_energy_bubble()
{
    const auto NT = static_cast<unsigned int>((system->Tmax - system->Tmin) / system->dT) + 1;

    if (mympi->my_rank == 0) {
        std::cout << std::endl;
        std::cout << " -----------------------------------------------------------------"
            << std::endl;
        std::cout << " Calculating the vibrational free energy from the Bubble diagram " << std::endl;
    }

    memory->allocate(FE_bubble, NT);

    compute_FE_bubble(dynamical->eval_phonon,
                      dynamical->evec_phonon,
                      FE_bubble);

    if (mympi->my_rank == 0) {
        std::cout << " done!" << std::endl << std::endl;
    }
}


void Thermodynamics::compute_FE_bubble(double **eval,
                                       std::complex<double> ***evec,
                                       double *FE_bubble) const
{
    // This function calculates the free energy of the bubble diagram
    int i;
    double omega_sum[2];
    double nsum[2];
    const auto nk = kpoint->nk;
    const auto nk_reduced = kpoint->kpoint_irred_all.size();
    const auto ns = dynamical->neval;
    unsigned int i0, iT;
    unsigned int arr_cubic[3];
    const int nks0 = nk_reduced * ns;
    const auto NT = static_cast<unsigned int>((system->Tmax - system->Tmin) / system->dT) + 1;
    const auto factor = -1.0 / (static_cast<double>(nk * nk) * 48.0);

    double n0, n1, n2;
    double *FE_local;
    double *FE_tmp;

    memory->allocate(FE_local, NT);
    memory->allocate(FE_tmp, NT);
    std::vector<KsListGroup> triplet;

    std::vector<int> vks_l;
    vks_l.clear();

    for (i0 = 0; i0 < nks0; ++i0) {
        if (i0 % mympi->nprocs == mympi->my_rank) {
            vks_l.push_back(i0);
        }
    }

    unsigned int nk_tmp;

    if (nks0 % mympi->nprocs != 0) {
        nk_tmp = nks0 / mympi->nprocs + 1;
    } else {
        nk_tmp = nks0 / mympi->nprocs;
    }
    if (vks_l.size() < nk_tmp) {
        vks_l.push_back(-1);
    }

    for (iT = 0; iT < NT; ++iT) FE_local[iT] = 0.0;

    for (i0 = 0; i0 < nk_tmp; ++i0) {

        if (vks_l[i0] != -1) {

            const auto ik0 = kpoint->kpoint_irred_all[vks_l[i0] / ns][0].knum;
            const auto is0 = vks_l[i0] % ns;

            kpoint->get_unique_triplet_k(vks_l[i0] / ns,
                                         anharmonic_core->use_triplet_symmetry,
                                         true,
                                         triplet, 1);

            const int npair_uniq = triplet.size();

            arr_cubic[0] = ns * ik0 + is0;

            for (iT = 0; iT < NT; ++iT) FE_tmp[iT] = 0.0;

            for (int ik = 0; ik < npair_uniq; ++ik) {
                int multi = triplet[ik].group.size();

                arr_cubic[0] = ns * ik0 + is0;

                const unsigned int ik1 = triplet[ik].group[0].ks[0];
                const unsigned int ik2 = triplet[ik].group[0].ks[1];

                for (unsigned int is1 = 0; is1 < ns; ++is1) {
                    arr_cubic[1] = ns * ik1 + is1;

                    for (unsigned int is2 = 0; is2 < ns; ++is2) {
                        arr_cubic[2] = ns * ik2 + is2;

                        const auto omega0 = eval[ik0][is0];
                        const auto omega1 = eval[ik1][is1];
                        const auto omega2 = eval[ik2][is2];

                        if (omega0 < eps8 || omega1 < eps8 || omega2 < eps8) continue;

                        omega_sum[0] = 1.0 / (omega0 + omega1 + omega2);
                        omega_sum[1] = 1.0 / (-omega0 + omega1 + omega2);

                        const auto v3_tmp = std::norm(anharmonic_core->V3(arr_cubic))
                            * static_cast<double>(multi);


                        for (iT = 0; iT < NT; ++iT) {
                            const auto temp = system->Tmin + static_cast<double>(iT) * system->dT;

                            if (classical) {
                                n0 = fC(omega0, temp);
                                n1 = fC(omega1, temp);
                                n2 = fC(omega2, temp);

                                nsum[0] = n0 * (n1 + n2) + n1 * n2;
                                nsum[1] = n0 * (n1 + n2) - n1 * n2;
                            } else {
                                n0 = fB(omega0, temp);
                                n1 = fB(omega1, temp);
                                n2 = fB(omega2, temp);

                                nsum[0] = (1.0 + n0) * (1.0 + n1 + n2) + n1 * n2;
                                nsum[1] = n0 * n1 - n1 * n2 + n2 * n0 + n0;
                            }

                            FE_tmp[iT] += v3_tmp * (nsum[0] * omega_sum[0] + 3.0 * nsum[1] * omega_sum[1]);
                        }
                    }
                }
            }
            const auto weight = static_cast<double>(kpoint->kpoint_irred_all[vks_l[i0] / ns].size());
            for (iT = 0; iT < NT; ++iT) FE_local[iT] += FE_tmp[iT] * weight;
        }
    }

    MPI_Allreduce(&FE_local[0], &FE_bubble[0], NT, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    for (iT = 0; iT < NT; ++iT) {
        FE_bubble[iT] *= factor;
    }

    memory->deallocate(FE_local);
    memory->deallocate(FE_tmp);
}


void Thermodynamics::compute_FE_bubble_SCPH(double ***eval_in,
                                            std::complex<double> ****evec_in,
                                            double *FE_bubble)
{
    // This function calculates the free energy from the bubble diagram
    // at the given temperature and lattice dynamics wavefunction
    int i;
    int ik;
    int multi;
    double omega0, omega1, omega2, omega_sum[2];
    double n0, n1, n2, nsum[2];
    unsigned int nk = kpoint->nk;
    unsigned int nk_reduced = kpoint->kpoint_irred_all.size();
    unsigned int ns = dynamical->neval;
    double v3_tmp;
    unsigned int ik0, ik1, ik2, is0, is1, is2, i0, iT;
    unsigned int arr_cubic[3];
    int nks0 = nk_reduced * ns;
    unsigned int NT = static_cast<unsigned int>((system->Tmax - system->Tmin) / system->dT) + 1;
    double temp;
    double factor = -1.0 / (static_cast<double>(nk * nk) * 48.0);

    double *FE_local;
    double *FE_tmp;

    memory->allocate(FE_local, NT);
    memory->allocate(FE_tmp, NT);
    std::vector<KsListGroup> triplet;

    std::vector<int> vks_l;
    vks_l.clear();

    for (i0 = 0; i0 < nks0; ++i0) {
        if (i0 % mympi->nprocs == mympi->my_rank) {
            vks_l.push_back(i0);
        }
    }

    unsigned int nk_tmp;

    if (nks0 % mympi->nprocs != 0) {
        nk_tmp = nks0 / mympi->nprocs + 1;
    } else {
        nk_tmp = nks0 / mympi->nprocs;
    }
    if (vks_l.size() < nk_tmp) {
        vks_l.push_back(-1);
    }

    for (iT = 0; iT < NT; ++iT) FE_local[iT] = 0.0;

    if (mympi->my_rank == 0) {
        std::cout << " Total number of modes per MPI process: " << nk_tmp << std::endl;
    }

    for (i0 = 0; i0 < nk_tmp; ++i0) {

        if (vks_l[i0] != -1) {

            ik0 = kpoint->kpoint_irred_all[vks_l[i0] / ns][0].knum;
            is0 = vks_l[i0] % ns;

            kpoint->get_unique_triplet_k(vks_l[i0] / ns,
                                         anharmonic_core->use_triplet_symmetry,
                                         true,
                                         triplet, 1);

            int npair_uniq = triplet.size();

            arr_cubic[0] = ns * ik0 + is0;

            for (iT = 0; iT < NT; ++iT) FE_tmp[iT] = 0.0;

            for (ik = 0; ik < npair_uniq; ++ik) {
                multi = static_cast<double>(triplet[ik].group.size());

                arr_cubic[0] = ns * ik0 + is0;

                ik1 = triplet[ik].group[0].ks[0];
                ik2 = triplet[ik].group[0].ks[1];


                for (is1 = 0; is1 < ns; ++is1) {
                    arr_cubic[1] = ns * ik1 + is1;

                    for (is2 = 0; is2 < ns; ++is2) {
                        arr_cubic[2] = ns * ik2 + is2;

                        for (iT = 0; iT < NT; ++iT) {

                            temp = system->Tmin + static_cast<double>(iT) * system->dT;

                            omega0 = eval_in[iT][ik0][is0];
                            omega1 = eval_in[iT][ik1][is1];
                            omega2 = eval_in[iT][ik2][is2];

                            if (omega0 < eps8 || omega1 < eps8 || omega2 < eps8) continue;


                            omega_sum[0] = 1.0 / (omega0 + omega1 + omega2);
                            omega_sum[1] = 1.0 / (-omega0 + omega1 + omega2);

                            v3_tmp = std::norm(anharmonic_core->V3(arr_cubic, eval_in[iT], evec_in[iT]))
                                * static_cast<double>(multi);

                            if (classical) {
                                n0 = fC(omega0, temp);
                                n1 = fC(omega1, temp);
                                n2 = fC(omega2, temp);

                                nsum[0] = n0 * (n1 + n2) + n1 * n2;
                                nsum[1] = n0 * (n1 + n2) - n1 * n2;
                            } else {
                                n0 = fB(omega0, temp);
                                n1 = fB(omega1, temp);
                                n2 = fB(omega2, temp);

                                nsum[0] = (1.0 + n0) * (1.0 + n1 + n2) + n1 * n2;
                                nsum[1] = n0 * n1 - n1 * n2 + n2 * n0 + n0;
                            }

                            FE_tmp[iT] += v3_tmp * (nsum[0] * omega_sum[0] + 3.0 * nsum[1] * omega_sum[1]);
                        }
                    }
                }
            }
            double weight = static_cast<double>(kpoint->kpoint_irred_all[vks_l[i0] / ns].size());
            for (iT = 0; iT < NT; ++iT) FE_local[iT] += FE_tmp[iT] * weight;
        }
    }

    MPI_Allreduce(&FE_local[0], &FE_bubble[0], NT, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    for (iT = 0; iT < NT; ++iT) {
        FE_bubble[iT] *= factor;
    }

    memory->deallocate(FE_local);
    memory->deallocate(FE_tmp);
}


double Thermodynamics::FE_scph_correction(unsigned int iT,
                                          double **eval,
                                          std::complex<double> ***evec) const
{
    const auto nk = kpoint->nk;
    const auto ns = dynamical->neval;
    const auto temp = system->Tmin + static_cast<double>(iT) * system->dT;
    const auto N = nk * ns;

    double ret = 0.0;

#pragma omp parallel for reduction(+ : ret)
    for (int i = 0; i < N; ++i) {
        int ik = i / ns;
        int is = i % ns;
        double omega = eval[ik][is];
        if (std::abs(omega) < eps6) continue;

        std::complex<double> tmp_c = std::complex<double>(0.0, 0.0);

        for (int js = 0; js < ns; ++js) {
            double omega2_harm = dynamical->eval_phonon[ik][js];
            if (omega2_harm >= 0.0) {
                omega2_harm = std::pow(omega2_harm, 2);
            } else {
                omega2_harm = -std::pow(omega2_harm, 2);
            }

            for (int ks = 0; ks < ns; ++ks) {
                for (int ls = 0; ls < ns; ++ls) {
                    tmp_c += omega2_harm
                        * dynamical->evec_phonon[ik][js][ks]
                        * std::conj(dynamical->evec_phonon[ik][js][ls])
                        * std::conj(evec[ik][is][ks]) * evec[ik][is][ls];
                }
            }
        }

        if (thermodynamics->classical) {
            ret += (tmp_c.real() - omega * omega) * thermodynamics->fC(omega, temp) / (4.0 * omega);
        } else {
            ret += (tmp_c.real() - omega * omega) * (1.0 + 2.0 * thermodynamics->fB(omega, temp)) / (8.0 * omega);
        }
    }

    return ret / static_cast<double>(nk);
}
