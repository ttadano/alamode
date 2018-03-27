/*
phonon_dos.cpp

Copyright (c) 2014, 2015, 2016 Terumasa Tadano

This file is distributed under the terms of the MIT license.
Please see the file 'LICENCE.txt' in the root directory 
or http://opensource.org/licenses/mit-license.php for information.
*/

#include "mpi_common.h"
#include "phonon_dos.h"
#include "kpoint.h"
#include "constants.h"
#include "error.h"
#include "system.h"
#include "memory.h"
#include "dynamical.h"
#include "write_phonons.h"
#include "integration.h"
#include <algorithm>
#include <vector>
#include <fstream>
#include <iomanip>
#include "parsephon.h"
#include "symmetry_core.h"
#include "thermodynamics.h"

using namespace PHON_NS;

Dos::Dos(PHON *phon): Pointers(phon)
{
    set_default_variables();
}

Dos::~Dos()
{
    deallocate_variables();
}

void Dos::set_default_variables()
{
    flag_dos = false;
    compute_dos = true;
    projected_dos = false;
    two_phonon_dos = false;
    scattering_phase_space = 0;
    energy_dos = nullptr;
    dos_phonon = nullptr;
    pdos_phonon = nullptr;
    dos2_phonon = nullptr;
    sps3_mode = nullptr;
    sps3_with_bose = nullptr;
    kmap_irreducible = nullptr;
}

void Dos::deallocate_variables()
{
    if (energy_dos) {
        memory->deallocate(energy_dos);
    }
    if (dos_phonon) {
        memory->deallocate(dos_phonon);
    }
    if (pdos_phonon) {
        memory->deallocate(pdos_phonon);
    }
    if (dos2_phonon) {
        memory->deallocate(dos2_phonon);
    }
    if (sps3_mode) {
        memory->deallocate(sps3_mode);
    }
    if (sps3_with_bose) {
        memory->deallocate(sps3_mode);
    }
    if (kmap_irreducible) {
        memory->deallocate(kmap_irreducible);
    }
}


void Dos::setup()
{
    // This function must not be called before dynamica->setup_dynamical()

    int i;

    MPI_Bcast(&emin, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&emax, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&delta_e, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&compute_dos, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD);
    MPI_Bcast(&projected_dos, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD);
    MPI_Bcast(&two_phonon_dos, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD);
    MPI_Bcast(&scattering_phase_space, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (kpoint->kpoint_mode == 2) {
        flag_dos = true;
    } else {
        flag_dos = false;
    }

    if (flag_dos && delta_e < eps12)
        error->exit("Dos::setup()", "Too small delta_e");

    if (flag_dos) {
        n_energy = static_cast<int>((emax - emin) / delta_e);
        memory->allocate(energy_dos, n_energy);

        for (i = 0; i < n_energy; ++i) {
            energy_dos[i] = emin + delta_e * static_cast<double>(i);
        }

        if (compute_dos) {
            memory->allocate(dos_phonon, n_energy);
        }

        if (projected_dos) {
            memory->allocate(pdos_phonon, system->natmin, n_energy);
        }

        if (two_phonon_dos) {
            memory->allocate(dos2_phonon, kpoint->nk_irred, n_energy, 4);
        }

        if (scattering_phase_space == 1) {
            memory->allocate(sps3_mode, kpoint->nk_irred, dynamical->neval, 2);
        } else if (scattering_phase_space == 2) {
            double Tmin = system->Tmin;
            double Tmax = system->Tmax;
            double dT = system->dT;
            unsigned int NT = static_cast<unsigned int>((Tmax - Tmin) / dT) + 1;

            memory->allocate(sps3_with_bose, kpoint->nk_irred,
                             dynamical->neval, NT, 2);
        }

        int ***symmetry_tmp;

        memory->allocate(kmap_irreducible, kpoint->nk);
        memory->allocate(symmetry_tmp, symmetry->nsym, 3, 3);

        for (i = 0; i < symmetry->nsym; ++i) {
            for (int j = 0; j < 3; ++j) {
                for (int k = 0; k < 3; ++k) {
                    symmetry_tmp[i][j][k] = symmetry->SymmList[i].rot[j][k];
                }
            }
        }

        kpoint->generate_irreducible_kmap(kmap_irreducible, nk_irreducible, k_irreducible,
                                          kpoint->nkx, kpoint->nky, kpoint->nkz,
                                          kpoint->xk, symmetry->nsym, symmetry_tmp);

        memory->deallocate(symmetry_tmp);
    }
}

void Dos::calc_dos_all()
{
    unsigned int j, k;
    unsigned int nk = kpoint->nk;
    unsigned int neval = dynamical->neval;
    double **eval;

    memory->allocate(eval, neval, nk);

    for (j = 0; j < nk; ++j) {
        for (k = 0; k < neval; ++k) {
            eval[k][j] = writes->in_kayser(dynamical->eval_phonon[j][k]);
        }
    }

    if (compute_dos) {
        calc_dos(nk_irreducible, kmap_irreducible, eval, n_energy, energy_dos,
                 dos_phonon, neval, integration->ismear, kpoint->kpoint_irred_all);
    }

    if (projected_dos) {
        calc_atom_projected_dos(nk, eval, n_energy, energy_dos,
                                pdos_phonon, neval, system->natmin,
                                integration->ismear, dynamical->evec_phonon);
    }
    memory->deallocate(eval);

    if (two_phonon_dos) {
        calc_two_phonon_dos(n_energy, energy_dos, dos2_phonon,
                            integration->ismear, kpoint->kpoint_irred_all);
    }

    if (scattering_phase_space == 1) {
        calc_total_scattering_phase_space(dynamical->eval_phonon, integration->ismear,
                                          kpoint->kpoint_irred_all, sps3_mode, total_sps3);
    } else if (scattering_phase_space == 2) {
        calc_scattering_phase_space_with_Bose(dynamical->eval_phonon, integration->ismear,
                                              kpoint->kpoint_irred_all, sps3_with_bose);
    }
}

void Dos::calc_dos(const unsigned int nk_irreducible,
                   int *map_k,
                   double **eval,
                   const unsigned int n,
                   double *energy,
                   double *ret,
                   const unsigned int neval,
                   const int smearing_method,
                   std::vector<std::vector<KpointList>> &kpinfo)
{
    int i, j, k;
    double *weight;

    if (mympi->my_rank == 0) std::cout << " Calculating phonon DOS ...";
#ifdef _OPENMP
#pragma omp parallel private (weight, k)
#endif
    {
        memory->allocate(weight, nk_irreducible);
#ifdef _OPENMP
#pragma omp for
#endif
        for (i = 0; i < n; ++i) {

            ret[i] = 0.0;

            for (k = 0; k < neval; ++k) {
                if (smearing_method == -1) {
                    integration->calc_weight_tetrahedron(nk_irreducible, map_k,
                                                         weight, eval[k], energy[i]);
                } else {
                    integration->calc_weight_smearing(kpinfo, weight,
                                                      eval[k], energy[i], smearing_method);
                }

                for (j = 0; j < nk_irreducible; ++j) {
                    ret[i] += weight[j];
                }
            }
        }
        memory->deallocate(weight);
    }

    if (mympi->my_rank == 0) std::cout << " done." << std::endl;
}

void Dos::calc_atom_projected_dos(const unsigned int nk,
                                  double **eval,
                                  const unsigned int n,
                                  double *energy,
                                  double **ret,
                                  const unsigned int neval,
                                  const unsigned int natmin,
                                  const int smearing_method,
                                  std::complex<double> ***evec)
{
    // Calculate atom projected phonon-DOS

    int i;
    unsigned int j, k;
    int *kmap_identity;
    double *weight;
    double **proj;

    if (mympi->my_rank == 0)
        std::cout << " PDOS = 1 : Calculating atom-projected phonon DOS ... ";

    memory->allocate(kmap_identity, nk);
    memory->allocate(proj, neval, nk);


    for (i = 0; i < nk; ++i) kmap_identity[i] = i;

    for (unsigned int iat = 0; iat < natmin; ++iat) {

        for (unsigned int imode = 0; imode < neval; ++imode) {
            for (i = 0; i < nk; ++i) {

                proj[imode][i] = 0.0;

                for (unsigned int icrd = 0; icrd < 3; ++icrd) {
                    proj[imode][i] += std::norm(evec[i][imode][3 * iat + icrd]);
                }
            }
        }
#ifdef _OPENMP
#pragma omp parallel private (weight, k, j)
#endif
        {
            memory->allocate(weight, nk);
#ifdef _OPENMP
#pragma omp for
#endif
            for (i = 0; i < n; ++i) {
                ret[iat][i] = 0.0;

                for (k = 0; k < neval; ++k) {
                    if (smearing_method == -1) {
                        integration->calc_weight_tetrahedron(nk, kmap_identity,
                                                             weight, eval[k], energy[i]);
                    } else {
                        integration->calc_weight_smearing(nk, nk, kmap_identity,
                                                          weight, eval[k], energy[i],
                                                          smearing_method);
                    }

                    for (j = 0; j < nk; ++j) {
                        ret[iat][i] += proj[k][j] * weight[j];
                    }
                }
            }

            memory->deallocate(weight);
        }
    }
    memory->deallocate(proj);
    memory->deallocate(kmap_identity);

    if (mympi->my_rank == 0) std::cout << " done!" << std::endl;
}


void Dos::calc_two_phonon_dos(const unsigned int n,
                              double *energy,
                              double ***ret,
                              const int smearing_method,
                              const std::vector<std::vector<KpointList>> &kpinfo)
{
    int i, j;
    int is, js, ik, jk;
    int k;
    int ib;
    int knum;

    unsigned int nk = kpoint->nk;
    unsigned int ns = dynamical->neval;
    unsigned int nk_reduced = kpoint->nk_irred;

    int ns2 = ns * ns;

    int *kmap_identity;

    double **e_tmp;
    double **weight;
    double emax2 = 2.0 * emax;

    double xk_tmp[3];

    int loc;
    int *k_pair;

    if (mympi->my_rank == 0) {
        std::cout << " TDOS = 1 : Calculating two-phonon DOS for all irreducible k points." << std::endl;
        std::cout << "            This may take a while ... ";
    }

    memory->allocate(kmap_identity, nk);
    memory->allocate(e_tmp, 2, nk);
    memory->allocate(weight, n, nk);
    memory->allocate(k_pair, nk);


    for (i = 0; i < nk; ++i) kmap_identity[i] = i;

    for (ik = 0; ik < nk_reduced; ++ik) {

        knum = kpinfo[ik][0].knum;

        for (jk = 0; jk < nk; ++jk) {
            for (i = 0; i < 3; ++i) xk_tmp[i] = kpoint->xk[knum][i] + kpoint->xk[jk][i];
            k_pair[jk] = kpoint->get_knum(xk_tmp[0], xk_tmp[1], xk_tmp[2]);
        }

        for (i = 0; i < n; ++i) {
            for (j = 0; j < 2; ++j) {
                ret[ik][i][j] = 0.0;
            }
        }

        for (ib = 0; ib < ns2; ++ib) {

            is = ib / ns;
            js = ib % ns;
#ifdef _OPENMP
#pragma omp parallel for private(loc)
#endif
            for (jk = 0; jk < nk; ++jk) {
                loc = k_pair[jk];
                e_tmp[0][jk]
                    = writes->in_kayser(dynamical->eval_phonon[jk][is]
                        + dynamical->eval_phonon[loc][js]);
                e_tmp[1][jk]
                    = writes->in_kayser(dynamical->eval_phonon[jk][is]
                        - dynamical->eval_phonon[loc][js]);
            }


            if (smearing_method == -1) {

                for (j = 0; j < 2; ++j) {
#ifdef _OPENMP
#pragma omp parallel for private(k)
#endif
                    for (i = 0; i < n; ++i) {
                        integration->calc_weight_tetrahedron(nk, kmap_identity,
                                                             weight[i], e_tmp[j], energy[i]);
                        for (k = 0; k < nk; ++k) {
                            ret[ik][i][j] += weight[i][k];
                        }
                    }
                }

            } else {

                for (j = 0; j < 2; ++j) {
#ifdef _OPENMP
#pragma omp parallel for private(k)
#endif
                    for (i = 0; i < n; ++i) {
                        integration->calc_weight_smearing(nk, nk, kmap_identity,
                                                          weight[i], e_tmp[j], energy[i],
                                                          smearing_method);
                        for (k = 0; k < nk; ++k) {
                            ret[ik][i][j] += weight[i][k];
                        }
                    }
                }
            }
        }
    }

    memory->deallocate(e_tmp);
    memory->deallocate(weight);
    memory->deallocate(kmap_identity);
    memory->deallocate(k_pair);

    if (mympi->my_rank == 0) {
        std::cout << "done!" << std::endl;
    }
}

void Dos::calc_total_scattering_phase_space(double **omega,
                                            const int smearing_method,
                                            const std::vector<std::vector<KpointList>> &kpinfo,
                                            double ***ret_mode,
                                            double &ret)
{
    int i, j;
    int is;
    int knum;

    unsigned int nk = kpoint->nk;
    unsigned int ns = dynamical->neval;
    int ns2 = ns * ns;
    int ib;

    int *kmap_identity;

    double multi;
    double omega0;
    double sps_tmp1, sps_tmp2;

    if (mympi->my_rank == 0) {
        std::cout << " SPS = 1 : Calculating three-phonon scattering phase space ... ";
    }

    memory->allocate(kmap_identity, nk);

    for (i = 0; i < nk; ++i) kmap_identity[i] = i;

    ret = 0.0;
    double sps_sum1 = 0.0;
    double sps_sum2 = 0.0;

    for (int ik = 0; ik < kpinfo.size(); ++ik) {

        knum = kpinfo[ik][0].knum;
        multi = static_cast<double>(kpinfo[ik].size()) / static_cast<double>(nk);

        for (is = 0; is < ns; ++is) {

            omega0 = writes->in_kayser(omega[knum][is]);

            sps_tmp1 = 0.0;
            sps_tmp2 = 0.0;
#ifdef _OPENMP
#pragma omp parallel
#endif
            {
                double **e_tmp;
                double *weight;
                int js, ks;
                int jk, loc;
                double xk_tmp[3];

                memory->allocate(weight, nk);
                memory->allocate(e_tmp, 2, nk);
#ifdef _OPENMP
#pragma omp for private(i, j), reduction(+: sps_tmp1, sps_tmp2)
#endif
                for (ib = 0; ib < ns2; ++ib) {

                    js = ib / ns;
                    ks = ib % ns;

                    for (jk = 0; jk < nk; ++jk) {

                        for (i = 0; i < 3; ++i) xk_tmp[i] = kpoint->xk[knum][i] + kpoint->xk[jk][i];
                        loc = kpoint->get_knum(xk_tmp[0], xk_tmp[1], xk_tmp[2]);

                        e_tmp[0][jk] = writes->in_kayser(omega[jk][js] + omega[loc][ks]);
                        e_tmp[1][jk] = writes->in_kayser(omega[jk][js] - omega[loc][ks]);

                    }

                    if (smearing_method == -1) {

                        integration->calc_weight_tetrahedron(nk, kmap_identity,
                                                             weight, e_tmp[0], omega0);
                        for (j = 0; j < nk; ++j) sps_tmp1 += weight[j];
                        integration->calc_weight_tetrahedron(nk, kmap_identity,
                                                             weight, e_tmp[1], omega0);
                        for (j = 0; j < nk; ++j) sps_tmp2 += weight[j];


                    } else {

                        integration->calc_weight_smearing(nk, nk, kmap_identity,
                                                          weight, e_tmp[0], omega0,
                                                          smearing_method);
                        for (j = 0; j < nk; ++j) sps_tmp1 += weight[j];
                        integration->calc_weight_smearing(nk, nk, kmap_identity,
                                                          weight, e_tmp[1], omega0,
                                                          smearing_method);
                        for (j = 0; j < nk; ++j) sps_tmp2 += weight[j];

                    }

                }
                memory->deallocate(e_tmp);
                memory->deallocate(weight);
            }
            sps_sum1 += multi * sps_tmp1;
            sps_sum2 += multi * sps_tmp2;

            ret_mode[ik][is][0] = sps_tmp1;
            ret_mode[ik][is][1] = sps_tmp2;
        }
    }

    memory->deallocate(kmap_identity);

    ret = (sps_sum1 + 2.0 * sps_sum2)
        / (3.0 * static_cast<double>(std::pow(ns, 3.0)));

    if (mympi->my_rank == 0) {
        std::cout << "done!" << std::endl;
    }
}

void Dos::calc_dos_scph(double ***eval_anharm,
                        double **dos_scph)
{
    unsigned int j, k;
    unsigned int nk = kpoint->nk;
    unsigned int neval = dynamical->neval;
    double **eval;

    double Tmin = system->Tmin;
    double Tmax = system->Tmax;
    double dT = system->dT;
    unsigned int NT = static_cast<unsigned int>((Tmax - Tmin) / dT) + 1;

    memory->allocate(eval, neval, nk);

    for (unsigned int iT = 0; iT < NT; ++iT) {

        std::cout << " T = " << std::setw(5) << Tmin + static_cast<double>(iT) * dT << std::endl;

        for (j = 0; j < nk; ++j) {
            for (k = 0; k < neval; ++k) {
                eval[k][j] = writes->in_kayser(eval_anharm[iT][j][k]);
            }
        }

        calc_dos(nk_irreducible, kmap_irreducible, eval, n_energy, energy_dos,
                 dos_scph[iT], neval, integration->ismear, kpoint->kpoint_irred_all);
    }
}

void Dos::calc_scattering_phase_space_with_Bose(double **eval,
                                                const int smearing_method,
                                                const std::vector<std::vector<KpointList>> &kp_info,
                                                double ****ret)
{
    unsigned int i, j, k;
    unsigned int knum, snum;
    double xk_tmp[3];
    double **ret_mode;
    double omega0;
    double Tmin = system->Tmin;
    double Tmax = system->Tmax;
    double dT = system->dT;
    double *temperature;
    int N;
    int ik, iT;
    unsigned int nk_irred = kpoint->nk_irred;
    unsigned int nk = kpoint->nk;
    unsigned int ns = dynamical->neval;
    unsigned int k1, k2;
    unsigned int imode;
    unsigned int *k2_arr;
    unsigned int ns2 = ns * ns;
    double **recv_buf;
    double omega_max = emax;
    double omega_min = emin;

    std::vector<int> ks_g, ks_l;
    int iks;

    if (mympi->my_rank == 0) {
        std::cout << " SPS = 2 : Calculating three-phonon scattering phase space" << std::endl;
        std::cout << "           with the Bose distribution function ...";
    }

    N = static_cast<int>((Tmax - Tmin) / dT) + 1;
    memory->allocate(temperature, N);
    for (i = 0; i < N; ++i) temperature[i] = Tmin + static_cast<double>(i) * dT;

    memory->allocate(k2_arr, nk);

    for (i = 0; i < nk_irred; ++i) {
        for (j = 0; j < ns; ++j) {
            for (k = 0; k < N; ++k) {
                ret[i][j][k][0] = 0.0;
                ret[i][j][k][1] = 0.0;
            }
        }
    }

    memory->allocate(ret_mode, N, 2);

    unsigned int nks_total = nk_irred * ns;
    unsigned int nks_each_thread = nks_total / mympi->nprocs;
    unsigned int nrem = nks_total - nks_each_thread * mympi->nprocs;

    if (nrem > 0) {
        memory->allocate(recv_buf, (nks_each_thread + 1) * mympi->nprocs, 2 * N);
    } else {
        memory->allocate(recv_buf, nks_total, 2 * N);
    }

    ks_g.clear();
    for (ik = 0; ik < nk_irred; ++ik) {

        knum = kp_info[ik][0].knum;

        for (imode = 0; imode < ns; ++imode) {

            omega0 = writes->in_kayser(eval[knum][imode]);
            if (omega0 < omega_min || omega0 > omega_max) continue;

            ks_g.push_back(ik * ns + imode);
        }
    }

    ks_l.clear();
    unsigned int count = 0;
    for (auto it = ks_g.begin(); it != ks_g.end(); ++it) {
        if (count % mympi->nprocs == mympi->my_rank) {
            ks_l.push_back(*it);
        }
        ++count;
    }

    unsigned int nks_tmp;
    if (ks_g.size() % mympi->nprocs > 0) {
        nks_tmp = ks_g.size() / mympi->nprocs + 1;
    } else {
        nks_tmp = ks_g.size() / mympi->nprocs;
    }
    if (ks_l.size() < nks_tmp) ks_l.push_back(-1);

    for (i = 0; i < nks_tmp; ++i) {

        iks = ks_l[i];

        if (iks == -1) {

            for (iT = 0; iT < N; ++iT) {
                ret_mode[iT][0] = 0.0;
                ret_mode[iT][1] = 0.0;
            }

        } else {

            knum = kp_info[iks / ns][0].knum;
            snum = iks % ns;

            for (k1 = 0; k1 < nk; ++k1) {
                for (j = 0; j < 3; ++j) xk_tmp[j] = kpoint->xk[knum][j] - kpoint->xk[k1][j];
                k2 = kpoint->get_knum(xk_tmp[0], xk_tmp[1], xk_tmp[2]);
                k2_arr[k1] = k2;
            }

            omega0 = eval[knum][snum];
            calc_scattering_phase_space_with_Bose_mode(nk, ns, N, omega0, eval,
                                                       temperature, k2_arr,
                                                       smearing_method, ret_mode);
        }
        MPI_Gather(&ret_mode[0][0], 2 * N, MPI_DOUBLE, &recv_buf[mympi->nprocs * i][0],
                   2 * N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

    count = 0;
    for (ik = 0; ik < nk_irred; ++ik) {

        knum = kp_info[ik][0].knum;

        for (imode = 0; imode < ns; ++imode) {

            omega0 = writes->in_kayser(eval[knum][imode]);
            if (omega0 < omega_min || omega0 > omega_max) continue;

            for (iT = 0; iT < N; ++iT) {
                ret[ik][imode][iT][0] = recv_buf[count][2 * iT];
                ret[ik][imode][iT][1] = recv_buf[count][2 * iT + 1];
            }
            ++count;
        }
    }

    memory->deallocate(ret_mode);
    memory->deallocate(k2_arr);
    memory->deallocate(recv_buf);
    memory->deallocate(temperature);

    if (mympi->my_rank == 0) {
        std::cout << " done!" << std::endl;
    }
}

void Dos::calc_scattering_phase_space_with_Bose_mode(const unsigned int nk,
                                                     const unsigned int ns,
                                                     const unsigned int N,
                                                     const double omega,
                                                     double **eval,
                                                     double *temperature,
                                                     unsigned int *k_pair,
                                                     const int smearing_method,
                                                     double **ret)
{
    int ib;
    unsigned int i, is, js, k1, k2;
    unsigned int iT;
    unsigned int ns2 = ns * ns;
    double omega1, omega2;
    double temp;
    double ret1, ret2;
    double n1, n2, f1, f2;

    double **energy_tmp;
    double **weight;
    double ***delta_arr;

    int *kmap_identity;

    memory->allocate(delta_arr, nk, ns2, 2);

    memory->allocate(kmap_identity, nk);
    for (i = 0; i < nk; ++i) kmap_identity[i] = i;

    double omega0 = writes->in_kayser(omega);

#ifdef _OPENMP
#pragma omp parallel private(i, is, js, k1, k2, omega1, omega2, energy_tmp, weight)
#endif
    {
        memory->allocate(energy_tmp, 2, nk);
        memory->allocate(weight, 2, nk);
#ifdef _OPENMP
#pragma omp for
#endif
        for (ib = 0; ib < ns2; ++ib) {
            is = ib / ns;
            js = ib % ns;

            for (k1 = 0; k1 < nk; ++k1) {

                k2 = k_pair[k1];

                omega1 = eval[k1][is];
                omega2 = eval[k2][js];

                energy_tmp[0][k1] = writes->in_kayser(omega1 + omega2);
                energy_tmp[1][k1] = writes->in_kayser(omega1 - omega2);
            }

            if (smearing_method == -1) {
                for (i = 0; i < 2; ++i) {
                    integration->calc_weight_tetrahedron(nk, kmap_identity,
                                                         weight[i], energy_tmp[i], omega0);
                }
            } else {
                for (i = 0; i < 2; ++i) {
                    integration->calc_weight_smearing(nk, nk, kmap_identity,
                                                      weight[i], energy_tmp[i], omega0,
                                                      smearing_method);
                }
            }


            for (k1 = 0; k1 < nk; ++k1) {
                delta_arr[k1][ib][0] = weight[0][k1];
                delta_arr[k1][ib][1] = weight[1][k1];
            }
        }

        memory->deallocate(energy_tmp);
        memory->deallocate(weight);
    }


    for (iT = 0; iT < N; ++iT) {
        temp = temperature[iT];
        ret1 = 0.0;
        ret2 = 0.0;
#ifdef _OPENMP
#pragma omp parallel for private(k1, k2, is, js, omega1, omega2, n1, n2, f1, f2), reduction(+:ret1, ret2)
#endif
        for (ib = 0; ib < ns2; ++ib) {

            is = ib / ns;
            js = ib % ns;

            for (k1 = 0; k1 < nk; ++k1) {

                k2 = k_pair[k1];

                omega1 = eval[k1][is];
                omega2 = eval[k2][js];

                if (omega1 < eps12 || omega2 < eps12) continue;

                if (thermodynamics->classical) {
                    f1 = thermodynamics->fC(omega1, temp);
                    f2 = thermodynamics->fC(omega2, temp);
                    n1 = f1 + f2;
                    n2 = f1 - f2;
                } else {
                    f1 = thermodynamics->fB(omega1, temp);
                    f2 = thermodynamics->fB(omega2, temp);
                    n1 = f1 + f2 + 1.0;
                    n2 = f1 - f2;
                }


                ret1 += delta_arr[k1][ib][0] * n1;
                ret2 += -delta_arr[k1][ib][1] * n2;
            }
        }
        ret[iT][0] = ret1;
        ret[iT][1] = ret2;
    }

    memory->deallocate(delta_arr);
    memory->deallocate(kmap_identity);
}
