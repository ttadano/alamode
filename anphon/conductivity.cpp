/*
 conductivity.cpp

 Copyright (c) 2014, 2015, 2016 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#include "mpi_common.h"
#include "conductivity.h"
#include "constants.h"
#include "dynamical.h"
#include "error.h"
#include "integration.h"
#include "isotope.h"
#include "kpoint.h"
#include "mathfunctions.h"
#include "memory.h"
#include "phonon_dos.h"
#include "thermodynamics.h"
#include "phonon_velocity.h"
#include "anharmonic_core.h"
#include "system.h"
#include "write_phonons.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <set>
#include <vector>

using namespace PHON_NS;

Conductivity::Conductivity(PHON *phon): Pointers(phon)
{
    set_default_variables();
}

Conductivity::~Conductivity()
{
    deallocate_variables();
};

void Conductivity::set_default_variables()
{
    calc_kappa_spec = 0;
    ntemp = 0;
    damping3 = nullptr;
    kappa = nullptr;
    kappa_spec = nullptr;
    kappa_nondiagonal = nullptr;
    Temperature = nullptr;
    vel = nullptr;
    velmat = nullptr;
    nondiag_hflux = 1;
}

void Conductivity::deallocate_variables()
{
    if (damping3) {
        memory->deallocate(damping3);
    }
    if (kappa) {
        memory->deallocate(kappa);
    }
    if (kappa_spec) {
        memory->deallocate(kappa_spec);
    }
    if (kappa_nondiagonal) {
        memory->deallocate(kappa_nondiagonal);
    }
    if (Temperature) {
        memory->deallocate(Temperature);
    }
    if (vel) {
        memory->deallocate(vel);
    }
    if (velmat) {
        memory->deallocate(velmat);
    }
}


void Conductivity::setup_kappa()
{
    unsigned int i, j, k;

    nk = kpoint->nk;
    ns = dynamical->neval;

    ntemp = static_cast<unsigned int>((system->Tmax - system->Tmin) / system->dT) + 1;
    memory->allocate(Temperature, ntemp);

    for (i = 0; i < ntemp; ++i) {
        Temperature[i] = system->Tmin + static_cast<double>(i) * system->dT;
    }

    const auto nks_total = kpoint->nk_irred * ns;
    const auto nks_each_thread = nks_total / mympi->nprocs;
    const auto nrem = nks_total - nks_each_thread * mympi->nprocs;

    if (nrem > 0) {
        memory->allocate(damping3, (nks_each_thread + 1) * mympi->nprocs, ntemp);
    } else {
        memory->allocate(damping3, nks_total, ntemp);
    }

    if (mympi->my_rank == 0) {
        memory->allocate(vel, nk, ns, 3);

        for (i = 0; i < nk; ++i) {
            phonon_velocity->phonon_vel_k(kpoint->xk[i], vel[i]);

            // Generate phonon velocity in Cartesian coordinate
            for (j = 0; j < ns; ++j) {
                rotvec(vel[i][j], vel[i][j], system->lavec_p);
                for (k = 0; k < 3; ++k) vel[i][j][k] /= 2.0 * pi;
                for (k = 0; k < 3; ++k) vel[i][j][k] *= Bohr_in_Angstrom * 1.0e-10 / time_ry;
            }
        }

        // nondiag_hflux = 1;

        if (nondiag_hflux) {
            memory->allocate(velmat, nk, ns, ns, 3);
            std::cout << std::scientific;

            for (i = 0; i < nk; ++i) {
     //           std::cout << "k = " << i << std::endl;
                phonon_velocity->velocity_matrix_analytic(kpoint->xk[i],
                                                          fcs_phonon->fc2_ext,
                                                          dynamical->eval_phonon[i],
                                                          dynamical->evec_phonon[i],
                                                          velmat[i]);

                for (j = 0; j < ns; ++j) {
                    for (k = 0; k < ns; ++k) {
                        rotvec(velmat[i][j][k], velmat[i][j][k], system->lavec_p);
                        for (auto mu = 0; mu < 3; ++mu) {
                            velmat[i][j][k][mu] *= Bohr_in_Angstrom * 1.0e-10 / (time_ry * 2.0 * pi);
                        }
                    }
                }

                /*for (auto mu = 0; mu < 3; ++mu) {
                    std::cout << "mu = " << mu << std::endl;

                    std::cout << "Diagonal:\n";

                    for (j = 0; j < ns; ++j) {
                        std::cout << std::setw(20) << vel[i][j][mu] << std::endl;
                    }

                    std::cout << "Full:\n";
                    for (j = 0; j < ns; ++j) {
                        for (k = 0; k < ns; ++k) {
                            std::cout << std::setw(20) << velmat[i][j][k][mu].real() << std::setw(15) << velmat[i][j][k]
                            [
                                mu].imag();
                        }
                        std::cout << std::endl;
                    }
                    std::cout << std::endl;
                }
                std::cout << std::endl;*/
            }
        }


    }

    vks_job.clear();

    for (i = 0; i < kpoint->nk_irred; ++i) {
        for (j = 0; j < ns; ++j) {
            vks_job.insert(i * ns + j);
        }
    }
}

void Conductivity::prepare_restart()
{
    // Write phonon frequency to result file

    int i;
    std::string line_tmp;
    unsigned int nk_tmp, ns_tmp;
    unsigned int multiplicity;
    int nks_done, *arr_done;

    double vel_dummy[3];

    nshift_restart = 0;

    vks_done.clear();

    if (mympi->my_rank == 0) {

        if (!phon->restart_flag) {

            writes->fs_result << "##Phonon Frequency" << std::endl;
            writes->fs_result << "#K-point (irreducible), Branch, Omega (cm^-1)" << std::endl;

            for (i = 0; i < kpoint->nk_irred; ++i) {
                int ik = kpoint->kpoint_irred_all[i][0].knum;
                for (int is = 0; is < dynamical->neval; ++is) {
                    writes->fs_result << std::setw(6) << i + 1 << std::setw(6) << is + 1;
                    writes->fs_result << std::setw(15) << writes->in_kayser(dynamical->eval_phonon[ik][is]) << std::
                        endl;
                }
            }

            writes->fs_result << "##END Phonon Frequency" << std::endl << std::endl;
            writes->fs_result << "##Phonon Relaxation Time" << std::endl;

        } else {

            while (writes->fs_result >> line_tmp) {

                if (line_tmp == "#GAMMA_EACH") {

                    writes->fs_result >> nk_tmp >> ns_tmp;
                    writes->fs_result >> multiplicity;

                    unsigned int nks_tmp = (nk_tmp - 1) * ns + ns_tmp - 1;

                    for (i = 0; i < multiplicity; ++i) {
                        writes->fs_result >> vel_dummy[0] >> vel_dummy[1] >> vel_dummy[2];
                    }

                    for (i = 0; i < ntemp; ++i) {
                        writes->fs_result >> damping3[nks_tmp][i];
                        damping3[nks_tmp][i] *= time_ry / Hz_to_kayser;
                    }
                    vks_done.push_back(nks_tmp);
                }
            }
        }

        writes->fs_result.close();
        writes->fs_result.open(writes->file_result.c_str(), std::ios::app | std::ios::out);
    }

    // Add vks_done list here

    if (mympi->my_rank == 0) {
        nks_done = vks_done.size();
    }
    MPI_Bcast(&nks_done, 1, MPI_INT, 0, MPI_COMM_WORLD);
    nshift_restart = nks_done;

    if (nks_done > 0) {
        memory->allocate(arr_done, nks_done);

        if (mympi->my_rank == 0) {
            for (i = 0; i < nks_done; ++i) {
                arr_done[i] = vks_done[i];
            }
        }
        MPI_Bcast(&arr_done[0], nks_done, MPI_INT, 0, MPI_COMM_WORLD);

        // Remove vks_done elements from vks_job

        for (i = 0; i < nks_done; ++i) {

            std::set<int>::iterator it_set = vks_job.find(arr_done[i]);

            if (it_set == vks_job.end()) {
                std::cout << " rank = " << mympi->my_rank
                    << " arr_done = " << arr_done[i] << std::endl;
                error->exit("prepare_restart", "This cannot happen");
            } else {
                vks_job.erase(it_set);
            }
        }
        memory->deallocate(arr_done);
    }
    vks_done.clear();
}


void Conductivity::calc_anharmonic_imagself()
{
    unsigned int i;
    unsigned int *nks_thread;
    double *damping3_loc;


    // Distribute (k,s) to individual MPI threads

    unsigned int nks_g = vks_job.size();
    vks_l.clear();

    unsigned int icount = 0;

    for (auto it = vks_job.begin(); it != vks_job.end(); ++it) {
        if (icount % mympi->nprocs == mympi->my_rank) {
            vks_l.push_back(*it);
        }
        ++icount;
    }

    if (mympi->my_rank == 0) {
        memory->allocate(nks_thread, mympi->nprocs);
    }

    unsigned int nks_tmp = vks_l.size();
    MPI_Gather(&nks_tmp, 1, MPI_UNSIGNED, &nks_thread[mympi->my_rank],
               1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

    if (mympi->my_rank == 0) {
        std::cout << std::endl;
        std::cout << " Start calculating anharmonic phonon self-energies ... " << std::endl;
        std::cout << " Total Number of phonon modes to be calculated : " << nks_g << std::endl;
        std::cout << " All modes are distributed to MPI threads as the following :" << std::endl;
        for (i = 0; i < mympi->nprocs; ++i) {
            std::cout << " RANK: " << std::setw(5) << i + 1;
            std::cout << std::setw(8) << "MODES: " << std::setw(5) << nks_thread[i] << std::endl;
        }
        std::cout << std::endl << std::flush;

        memory->deallocate(nks_thread);
    }

    unsigned int nk_tmp = nks_g / mympi->nprocs + 1;

    if (nks_g % mympi->nprocs != 0) {
        nk_tmp = nks_g / mympi->nprocs + 1;
    } else {
        nk_tmp = nks_g / mympi->nprocs;
    }

    if (vks_l.size() < nk_tmp) {
        vks_l.push_back(-1);
    }

    memory->allocate(damping3_loc, ntemp);


    for (i = 0; i < nk_tmp; ++i) {

        int iks = vks_l[i];

        if (iks == -1) {

            for (unsigned int j = 0; j < ntemp; ++j) damping3_loc[j] = eps; // do nothing

        } else {

            unsigned int knum = kpoint->kpoint_irred_all[iks / ns][0].knum;
            unsigned int snum = iks % ns;

            double omega = dynamical->eval_phonon[knum][snum];

            if (integration->ismear == 0 || integration->ismear == 1) {
                anharmonic_core->calc_damping_smearing(ntemp,
                                                       Temperature,
                                                       omega,
                                                       iks / ns,
                                                       snum,
                                                       damping3_loc);
            } else if (integration->ismear == -1) {
                anharmonic_core->calc_damping_tetrahedron(ntemp,
                                                          Temperature,
                                                          omega,
                                                          iks / ns,
                                                          snum,
                                                          damping3_loc);
            }
        }

        MPI_Gather(&damping3_loc[0], ntemp, MPI_DOUBLE,
                   damping3[nshift_restart + i * mympi->nprocs], ntemp,
                   MPI_DOUBLE, 0, MPI_COMM_WORLD);

        if (mympi->my_rank == 0) {
            write_result_gamma(i, nshift_restart, vel, damping3);
            std::cout << " MODE " << std::setw(5) << i + 1 << " done." << std::endl << std::flush;
        }
    }

    memory->deallocate(damping3_loc);
}

void Conductivity::write_result_gamma(const unsigned int ik,
                                      const unsigned int nshift,
                                      double ***vel_in,
                                      double **damp_in) const
{
    unsigned int np = mympi->nprocs;
    unsigned int k;

    for (unsigned int j = 0; j < np; ++j) {

        unsigned int iks_g = ik * np + j + nshift;

        if (iks_g >= kpoint->nk_irred * ns) break;

        writes->fs_result << "#GAMMA_EACH" << std::endl;
        writes->fs_result << iks_g / ns + 1 << " " << iks_g % ns + 1 << std::endl;

        unsigned int nk_equiv = kpoint->kpoint_irred_all[iks_g / ns].size();

        writes->fs_result << nk_equiv << std::endl;
        for (k = 0; k < nk_equiv; ++k) {
            unsigned int ktmp = kpoint->kpoint_irred_all[iks_g / ns][k].knum;
            writes->fs_result << std::setw(15) << vel_in[ktmp][iks_g % ns][0];
            writes->fs_result << std::setw(15) << vel_in[ktmp][iks_g % ns][1];
            writes->fs_result << std::setw(15) << vel_in[ktmp][iks_g % ns][2] << std::endl;
        }

        for (k = 0; k < ntemp; ++k) {
            writes->fs_result << std::setw(15)
                << damp_in[iks_g][k] * Hz_to_kayser / time_ry << std::endl;
        }
        writes->fs_result << "#END GAMMA_EACH" << std::endl;
    }
}

void Conductivity::compute_kappa()
{
    unsigned int i;
    unsigned int ik, is;
    unsigned int iks;

    const auto factor_toSI = 1.0e+18 / (std::pow(Bohr_in_Angstrom, 3) * system->volume_p);

    if (mympi->my_rank == 0) {

        std::string file_kl;
        std::ofstream ofs_kl;
        double damp_tmp;

        double **lifetime;
        double **gamma_total;
        double ****kappa_mode;

        memory->allocate(lifetime, kpoint->nk_irred * ns, ntemp);
        memory->allocate(kappa_mode, ntemp, 9, ns, kpoint->nk_irred);
        memory->allocate(gamma_total, kpoint->nk_irred * ns, ntemp);

        average_self_energy_at_degenerate_point(kpoint->nk_irred * ns, ntemp, damping3);

        if (isotope->include_isotope) {
            for (iks = 0; iks < kpoint->nk_irred * ns; ++iks) {
                unsigned int snum = iks % ns;
                if (dynamical->is_imaginary[iks / ns][snum]) {
                    for (i = 0; i < ntemp; ++i) {
                        lifetime[iks][i] = 0.0;
                        gamma_total[iks][i] = 1.0e+100; // very big number
                    }
                } else {
                    for (i = 0; i < ntemp; ++i) {
                        damp_tmp = damping3[iks][i] + isotope->gamma_isotope[iks / ns][snum];
                        gamma_total[iks][i] = damp_tmp;
                        if (damp_tmp > 1.0e-100) {
                            lifetime[iks][i] = 1.0e+12 * time_ry * 0.5 / damp_tmp;
                        } else {
                            lifetime[iks][i] = 0.0;
                        }
                    }
                }
            }
        } else {
            for (iks = 0; iks < kpoint->nk_irred * ns; ++iks) {

                if (dynamical->is_imaginary[iks / ns][iks % ns]) {
                    for (i = 0; i < ntemp; ++i) {
                        lifetime[iks][i] = 0.0;
                        gamma_total[iks][i] = 1.0e+100; // very big number
                    }
                } else {
                    for (i = 0; i < ntemp; ++i) {
                        damp_tmp = damping3[iks][i];
                        gamma_total[iks][i] = damp_tmp;
                        if (damp_tmp > 1.0e-100) {
                            lifetime[iks][i] = 1.0e+12 * time_ry * 0.5 / damp_tmp;
                        } else {
                            lifetime[iks][i] = 0.0;
                        }
                    }
                }
            }
        }

        memory->allocate(kappa, ntemp, 3, 3);

        for (i = 0; i < ntemp; ++i) {
            for (unsigned int j = 0; j < 3; ++j) {
                for (unsigned int k = 0; k < 3; ++k) {

                    if (Temperature[i] < eps) {
                        // Set kappa as zero when T = 0.
                        for (is = 0; is < ns; ++is) {
                            for (ik = 0; ik < kpoint->nk_irred; ++ik) {
                                kappa_mode[i][3 * j + k][is][ik] = 0.0;
                            }
                        }
                    } else {
                        for (is = 0; is < ns; ++is) {
                            for (ik = 0; ik < kpoint->nk_irred; ++ik) {
                                const auto knum = kpoint->kpoint_irred_all[ik][0].knum;
                                const auto omega = dynamical->eval_phonon[knum][is];
                                auto vv_tmp = 0.0;
                                const auto nk_equiv = kpoint->kpoint_irred_all[ik].size();

                                // Accumulate group velocity (diad product) for the reducible k points
                                for (auto ieq = 0; ieq < nk_equiv; ++ieq) {
                                    const auto ktmp = kpoint->kpoint_irred_all[ik][ieq].knum;
                                    vv_tmp += vel[ktmp][is][j] * vel[ktmp][is][k];
                                }

                                if (thermodynamics->classical) {
                                    kappa_mode[i][3 * j + k][is][ik] = thermodynamics->Cv_classical(
                                            omega, Temperature[i])
                                        * vv_tmp * lifetime[ns * ik + is][i];
                                } else {
                                    kappa_mode[i][3 * j + k][is][ik] = thermodynamics->Cv(omega, Temperature[i])
                                        * vv_tmp * lifetime[ns * ik + is][i];
                                }


                                // Convert to SI unit
                                kappa_mode[i][3 * j + k][is][ik] *= factor_toSI;

                            }
                        }
                    }

                    kappa[i][j][k] = 0.0;

                    for (is = 0; is < ns; ++is) {
                        for (ik = 0; ik < kpoint->nk_irred; ++ik) {
                            kappa[i][j][k] += kappa_mode[i][3 * j + k][is][ik];
                        }
                    }

                    kappa[i][j][k] /= static_cast<double>(nk);
                }
            }
        }

        unsigned int js;

        if (nondiag_hflux) {
            memory->allocate(kappa_nondiagonal, ntemp, 3, 3);

            for (i = 0; i < ntemp; ++i) {
                for (unsigned int j = 0; j < 3; ++j) {
                    for (unsigned int k = 0; k < 3; ++k) {
                        kappa_nondiagonal[i][j][k] = std::complex<double>(0.0, 0.0);

                        if (Temperature[i] > eps) {

                            for (is = 0; is < ns; ++is) {
                                for (js = 0; js < ns; ++js) {
                                    if (js == is) continue; // skip the diagonal component

                                    for (ik = 0; ik < kpoint->nk_irred; ++ik) {
                                        const auto knum = kpoint->kpoint_irred_all[ik][0].knum;

                                        const auto omega1 = dynamical->eval_phonon[knum][is];
                                        const auto omega2 = dynamical->eval_phonon[knum][js];

                                        if (omega1 < eps8 || omega2 < eps8) continue;
                                        auto vv_tmp = std::complex<double>(0.0, 0.0);
                                        const auto nk_equiv = kpoint->kpoint_irred_all[ik].size();

                                        // Accumulate group velocity (diad product) for the reducible k points
                                        for (auto ieq = 0; ieq < nk_equiv; ++ieq) {
                                            const auto ktmp = kpoint->kpoint_irred_all[ik][ieq].knum;
                                            vv_tmp += velmat[ktmp][is][js][j] * velmat[ktmp][js][is][k];
                                        }

                                        kappa_nondiagonal[i][j][k]
                                            += 2.0 * (omega1 * omega2) / (omega1 + omega2)
                                            * (thermodynamics->Cv(omega1, Temperature[i]) / omega1
                                                + thermodynamics->Cv(omega2, Temperature[i]) / omega2)
                                            * (gamma_total[ik * ns + is][i] + gamma_total[ik * ns + js][i])
                                            / (4.0 * std::pow(omega1 - omega2, 2.0)
                                                + std::pow(gamma_total[ik * ns + is][i]
                                                           + gamma_total[ik * ns + js][i],
                                                           2.0))
                                            * vv_tmp;
                                    }

                                }
                            }
                        }

                        kappa_nondiagonal[i][j][k] *= factor_toSI * 1.0e+12 * time_ry / static_cast<double>(nk);
                    }
                }
            }

        }


        memory->deallocate(lifetime);
        memory->deallocate(gamma_total);

        if (calc_kappa_spec)
            compute_frequency_resolved_kappa(ntemp, kappa_mode, integration->ismear);

        memory->deallocate(kappa_mode);

        std::cout << "#kappa nondiagonal\n";
        for (i = 0; i < ntemp; ++i) {
            for (unsigned int j = 0; j < 3; ++j) {
                for (unsigned int k = 0; k < 3; ++k) {
                    std::cout << std::setw(20) << kappa_nondiagonal[i][j][k].real() << std::setw(15) <<
                        kappa_nondiagonal[i][j][k].imag();
                }
            }
            std::cout << '\n';
        }
        std::cout << '\n';
    }
}


void Conductivity::average_self_energy_at_degenerate_point(const int n,
                                                           const int m,
                                                           double **damping) const
{
    int j, k, l;
    int nkr = kpoint->nk_irred;

    double *eval_tmp;
    double tol_omega = 1.0e-7; // Approximately equal to 0.01 cm^{-1}

    std::vector<int> degeneracy_at_k;

    memory->allocate(eval_tmp, ns);

    double *damping_sum;

    memory->allocate(damping_sum, m);

    for (int i = 0; i < nkr; ++i) {
        int ik = kpoint->kpoint_irred_all[i][0].knum;

        for (j = 0; j < ns; ++j) eval_tmp[j] = dynamical->eval_phonon[ik][j];

        degeneracy_at_k.clear();

        double omega_prev = eval_tmp[0];
        int ideg = 1;

        for (j = 1; j < ns; ++j) {
            double omega_now = eval_tmp[j];

            if (std::abs(omega_now - omega_prev) < tol_omega) {
                ++ideg;
            } else {
                degeneracy_at_k.push_back(ideg);
                ideg = 1;
                omega_prev = omega_now;
            }
        }
        degeneracy_at_k.push_back(ideg);

        int is = 0;
        for (j = 0; j < degeneracy_at_k.size(); ++j) {
            ideg = degeneracy_at_k[j];

            if (ideg > 1) {

                for (l = 0; l < m; ++l) damping_sum[l] = 0.0;

                for (k = is; k < is + ideg; ++k) {
                    for (l = 0; l < m; ++l) {
                        damping_sum[l] += damping[ns * i + k][l];
                    }
                }

                for (k = is; k < is + ideg; ++k) {
                    for (l = 0; l < m; ++l) {
                        damping[ns * i + k][l] = damping_sum[l] / static_cast<double>(ideg);
                    }
                }
            }

            is += ideg;
        }
    }
    memory->deallocate(damping_sum);
}

void Conductivity::compute_frequency_resolved_kappa(const int ntemp,
                                                    double ****kappa_mode,
                                                    const int smearing_method)
{
    int i, j;
    int *kmap_identity;
    double **eval;

    std::cout << std::endl;
    std::cout << " KAPPA_SPEC = 1 : Calculating thermal conductivity spectra ... ";

    memory->allocate(kappa_spec, dos->n_energy, ntemp, 3);
    memory->allocate(kmap_identity, nk);
    memory->allocate(eval, ns, nk);

    for (i = 0; i < nk; ++i) kmap_identity[i] = i;

    for (i = 0; i < nk; ++i) {
        for (j = 0; j < ns; ++j) {
            eval[j][i] = writes->in_kayser(dynamical->eval_phonon[i][j]);
        }
    }


#ifdef _OPENMP
#pragma omp parallel private (j)
#endif
    {
        int k;
        int knum;
        double *weight;
        memory->allocate(weight, nk);

#ifdef _OPENMP
#pragma omp for
#endif
        for (i = 0; i < dos->n_energy; ++i) {

            for (j = 0; j < ntemp; ++j) {
                for (k = 0; k < 3; ++k) {
                    kappa_spec[i][j][k] = 0.0;
                }
            }

            for (int is = 0; is < ns; ++is) {
                if (smearing_method == -1) {
                    integration->calc_weight_tetrahedron(nk, kmap_identity, weight,
                                                         eval[is], dos->energy_dos[i]);
                } else {
                    integration->calc_weight_smearing(nk, nk, kmap_identity, weight,
                                                      eval[is], dos->energy_dos[i],
                                                      smearing_method);
                }

                for (j = 0; j < ntemp; ++j) {
                    for (k = 0; k < 3; ++k) {
                        for (int ik = 0; ik < kpoint->nk_irred; ++ik) {
                            knum = kpoint->kpoint_irred_all[ik][0].knum;
                            kappa_spec[i][j][k] += kappa_mode[j][3 * k + k][is][ik] * weight[knum];
                        }
                    }
                }
            }
        }
        memory->deallocate(weight);
    }

    memory->deallocate(kmap_identity);
    memory->deallocate(eval);

    std::cout << " done!" << std::endl;
}
