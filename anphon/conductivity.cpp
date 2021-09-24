/*
 conductivity.cpp

 Copyright (c) 2014 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#include "mpi_common.h"
#include "conductivity.h"
#include "constants.h"
#include "dynamical.h"
#include "error.h"
#include "ewald.h"
#include "integration.h"
#include "parsephon.h"
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
#include <vector>

// DONE: check the number of unique pairs 
// DONE: implement the restriction of energy in pair generation with a default cutoff,
// TODO: implement adaptive smearing in integration module.

using namespace PHON_NS;

Conductivity::Conductivity(PHON *phon) : Pointers(phon)
{
    set_default_variables();
}

Conductivity::~Conductivity()
{
    deallocate_variables();
}

void Conductivity::set_default_variables()
{
    calc_kappa_spec = 0;
    ntemp = 0;
    damping3 = nullptr;
    damping4 = nullptr;
    kappa = nullptr;
    kappa_spec = nullptr;
    kappa_coherent = nullptr;
    temperature = nullptr;
    vel = nullptr;
    velmat = nullptr;
    calc_coherent = 0;
    file_coherent_elems = "";
    nshift_restart = 0;
    nshift_restart4 = 0;
    kmesh_4ph = nullptr;
    fph_rta = 0;
    dymat_4ph = nullptr;
    phase_storage_4ph = nullptr;
}

void Conductivity::deallocate_variables()
{
    if (damping3) {
        deallocate(damping3);
    }
    if (damping4) {
        deallocate(damping4);
    }
    if (kappa) {
        deallocate(kappa);
    }
    if (kappa_spec) {
        deallocate(kappa_spec);
    }
    if (kappa_coherent) {
        deallocate(kappa_coherent);
    }
    if (temperature) {
        deallocate(temperature);
    }
    if (vel) {
        deallocate(vel);
    }
    if (velmat) {
        deallocate(velmat);
    }
    if (kmesh_4ph) {
        delete kmesh_4ph;
    }
    if (dymat_4ph) {
        delete dymat_4ph;
    }
    if (phase_storage_4ph) {
        delete phase_storage_4ph;
    }
}

void Conductivity::setup_kappa()
{
    MPI_Bcast(&calc_coherent, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&nk_coarse[0], 3, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

    unsigned int i, j, k;

    nk = dos->kmesh_dos->nk;
    ns = dynamical->neval;

    ntemp = static_cast<unsigned int>((system->Tmax - system->Tmin) / system->dT) + 1;
    allocate(temperature, ntemp);

    for (i = 0; i < ntemp; ++i) {
        temperature[i] = system->Tmin + static_cast<double>(i) * system->dT;
    }

    const auto nks_total = dos->kmesh_dos->nk_irred * ns;
    const auto nks_each_thread = nks_total / mympi->nprocs;
    const auto nrem = nks_total - nks_each_thread * mympi->nprocs;


   // if (fph_rta < 2) {
    if (nrem > 0) {
        allocate(damping3, (nks_each_thread + 1) * mympi->nprocs, ntemp);
    } else {
        allocate(damping3, nks_total, ntemp);
    }

// This if statement was changed as (fph_rta > 0) in Wenhao's implementation
    if (anharmonic_core->quartic_mode) {

        if (mympi->my_rank == 0) {
            std::cout << " QUARTIC = 1: Four-phonon scattering rate will be calculated.\n";

            std::cout << "  KMESH for 4-ph:\n";
            std::cout << "   nk1 : " << std::setw(5) << nk_coarse[0] << '\n';
            std::cout << "   nk2 : " << std::setw(5) << nk_coarse[1] << '\n';
            std::cout << "   nk3 : " << std::setw(5) << nk_coarse[2] << '\n';

        }
        if (nrem > 0) {
            allocate(damping4, (nks_each_thread + 1) * mympi->nprocs, ntemp);
        } else {
            allocate(damping4, nks_total, ntemp);
        }

        double **eval_tmp;
        std::complex<double> ***evec_tmp;

        if (nk_coarse[0] * nk_coarse[1] * nk_coarse[2] > 0) {
            const auto neval = dynamical->neval;

            kmesh_4ph = new KpointMeshUniform(nk_coarse);
            kmesh_4ph->setup(symmetry->SymmList, system->rlavec_p);
            auto nk_4ph = kmesh_4ph->nk;
            dymat_4ph = new DymatEigenValue(true, false, nk_4ph, neval);

            allocate(eval_tmp, nk_4ph, neval);
            allocate(evec_tmp, nk_4ph, neval, neval);

            dynamical->get_eigenvalues_dymat(nk_4ph,
                                             kmesh_4ph->xk,
                                             kmesh_4ph->kvec_na,
                                             fcs_phonon->fc2_ext,
                                             ewald->fc2_without_dipole,
                                             true,
                                             eval_tmp,
                                             evec_tmp);

            if (!dynamical->get_projection_directions().empty()) {
                if (mympi->my_rank == 0) {
                    for (auto ik = 0; ik < nk_4ph; ++ik) {
                        dynamical->project_degenerate_eigenvectors(system->lavec_p,
                                                                   fcs_phonon->fc2_ext,
                                                                   kmesh_4ph->xk[ik],
                                                                   dynamical->get_projection_directions(),
                                                                   evec_tmp[ik]);
                    }
                }

                MPI_Bcast(&evec_tmp[0][0][0],
                          nk_4ph * neval * neval,
                          MPI_CXX_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
            }

            dymat_4ph->set_eigenvals_and_eigenvecs(nk_4ph,
                                                   eval_tmp,
                                                   evec_tmp);
            deallocate(eval_tmp);
            deallocate(evec_tmp);


        } else {
            const auto neval = dynamical->neval;
            kmesh_4ph = new KpointMeshUniform(dos->kmesh_dos->nk_i);
            kmesh_4ph->setup(symmetry->SymmList, system->rlavec_p);
            auto nk_4ph = kmesh_4ph->nk;
            dymat_4ph = new DymatEigenValue(true, false, nk_4ph, neval);

            allocate(eval_tmp, nk_4ph, neval);
            allocate(evec_tmp, nk_4ph, neval, neval);

            dynamical->get_eigenvalues_dymat(nk_4ph,
                                             kmesh_4ph->xk,
                                             kmesh_4ph->kvec_na,
                                             fcs_phonon->fc2_ext,
                                             ewald->fc2_without_dipole,
                                             true,
                                             eval_tmp,
                                             evec_tmp);

            if (!dynamical->get_projection_directions().empty()) {
                if (mympi->my_rank == 0) {
                    for (auto ik = 0; ik < nk_4ph; ++ik) {
                        dynamical->project_degenerate_eigenvectors(system->lavec_p,
                                                                   fcs_phonon->fc2_ext,
                                                                   kmesh_4ph->xk[ik],
                                                                   dynamical->get_projection_directions(),
                                                                   evec_tmp[ik]);
                    }
                }

                MPI_Bcast(&evec_tmp[0][0][0],
                          nk_4ph * neval * neval,
                          MPI_CXX_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
            }

            dymat_4ph->set_eigenvals_and_eigenvecs(nk_4ph,
                                                   eval_tmp,
                                                   evec_tmp);
            deallocate(eval_tmp);
            deallocate(evec_tmp);
        }

        phase_storage_4ph = new PhaseFactorStorage(kmesh_4ph->nk_i);
        phase_storage_4ph->create(true);
    }

    const auto factor = Bohr_in_Angstrom * 1.0e-10 / time_ry;

    if (mympi->my_rank == 0) {
        allocate(vel, nk, ns, 3);
        if (calc_coherent) {
            allocate(velmat, nk, ns, ns, 3);
        }
    } else {
        allocate(vel, 1, 1, 1);
        if (calc_coherent) {
            allocate(velmat, 1, 1, 1, 3);
        }
    }

    phonon_velocity->get_phonon_group_velocity_mesh_mpi(*dos->kmesh_dos,
                                                        system->lavec_p,
                                                        fcs_phonon->fc2_ext,
                                                        vel);
    if (mympi->my_rank == 0) {
        for (i = 0; i < nk; ++i) {
            for (j = 0; j < ns; ++j) {
                for (k = 0; k < 3; ++k) vel[i][j][k] *= factor;
            }
        }
    }

    if (calc_coherent) {
        phonon_velocity->calc_phonon_velmat_mesh(velmat);
        if (calc_coherent == 2) {
            file_coherent_elems = input->job_title + ".kc_elem";
        }
    }

    vks_job.clear();
    vks_job4.clear();

    for (i = 0; i < dos->kmesh_dos->nk_irred; ++i) {
        for (j = 0; j < ns; ++j) {
            vks_job.insert(i * ns + j);
        }
    }

    for (i = 0; i < kmesh_4ph->nk_irred; ++i) {
        for (j = 0; j < ns; ++j) {
            vks_job4.insert(i * ns + j);
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
    vks_done4.clear();

    if (mympi->my_rank == 0) {

        if (!phon->restart_flag) {

            if (fph_rta < 2) {
                writes->fs_result << "##Phonon Frequency" << std::endl;
                writes->fs_result << "#K-point (irreducible), Branch, Omega (cm^-1)" << std::endl;

                for (i = 0; i < dos->kmesh_dos->nk_irred; ++i) {
                    const auto ik = dos->kmesh_dos->kpoint_irred_all[i][0].knum;
                    for (auto is = 0; is < dynamical->neval; ++is) {
                        writes->fs_result << std::setw(6) << i + 1 << std::setw(6) << is + 1;
                        writes->fs_result << std::setw(15) << writes->in_kayser(dos->dymat_dos->get_eigenvalues()[ik][is])
                        << std::
                        endl;
                    }

                    writes->fs_result << "##END Phonon Frequency" << std::endl << std::endl;
                    writes->fs_result << "##Phonon Relaxation Time" << std::endl;
                }
                if (fph_rta > 0) {
                    writes->fs_result4 << "##Phonon Frequency" << std::endl;
                    writes->fs_result4 << "#K-point (irreducible), Branch, Omega (cm^-1)" << std::endl;

                    for (i = 0; i < kmesh_4ph->nk_irred; ++i) {
                        const int ik = kmesh_4ph->kpoint_irred_all[i][0].knum;
                        for (auto is = 0; is < dynamical->neval; ++is) {
                            writes->fs_result4 << std::setw(6) << i + 1 << std::setw(6) << is + 1;
                            writes->fs_result4 << std::setw(15) << writes->in_kayser(dymat_4ph->get_eigenvalues()[ik][is]) << std::
                            endl;
                        }
                    }

                    writes->fs_result4 << "##END Phonon Frequency" << std::endl << std::endl;
                    writes->fs_result4 << "##Phonon Relaxation Time" << std::endl;
                }

            } else {

                if (fph_rta < 2) {
                    while (writes->fs_result >> line_tmp) {

                        if (line_tmp == "#GAMMA_EACH") {

                            writes->fs_result >> nk_tmp >> ns_tmp;
                            writes->fs_result >> multiplicity;

                            const auto nks_tmp = (nk_tmp - 1) * ns + ns_tmp - 1;

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
                if (fph_rta > 0) {
                    while (writes->fs_result4 >> line_tmp) {

                        if (line_tmp == "#GAMMA_EACH") {

                            writes->fs_result4 >> nk_tmp >> ns_tmp;
                            writes->fs_result4 >> multiplicity;

                            const auto nks_tmp = (nk_tmp - 1) * ns + ns_tmp - 1;

                            for (i = 0; i < multiplicity; ++i) {
                                writes->fs_result4 >> vel_dummy[0] >> vel_dummy[1] >> vel_dummy[2];
                            }

                            for (i = 0; i < ntemp; ++i) {
                                writes->fs_result4 >> damping4[nks_tmp][i];
                                damping3[nks_tmp][i] *= time_ry / Hz_to_kayser;
                            }
                            vks_done4.push_back(nks_tmp);
                        }
                    }
                }
            }

            if (fph_rta < 2) {
                writes->fs_result.close();
                writes->fs_result.open(writes->file_result.c_str(), std::ios::app | std::ios::out);
            }
            if (fph_rta > 0) {
                writes->fs_result4.close();
                writes->fs_result4.open(writes->file_result4.c_str(), std::ios::app | std::ios::out);
            }
        }

        // Add vks_done, nshift_restart
        if (fph_rta < 2) {
            if (mympi->my_rank == 0) {
                nks_done = vks_done.size();
            }
            MPI_Bcast(&nks_done, 1, MPI_INT, 0, MPI_COMM_WORLD);
            nshift_restart = nks_done;

            if (nks_done > 0) {
                allocate(arr_done, nks_done);

                if (mympi->my_rank == 0) {
                    for (i = 0; i < nks_done; ++i) {
                        arr_done[i] = vks_done[i];
                    }
                }
                MPI_Bcast(&arr_done[0], nks_done, MPI_INT, 0, MPI_COMM_WORLD);

                // Remove vks_done elements from vks_job

                for (i = 0; i < nks_done; ++i) {

                    const auto it_set = vks_job.find(arr_done[i]);

                    if (it_set == vks_job.end()) {
                        std::cout << " rank = " << mympi->my_rank
                        << " arr_done = " << arr_done[i] << std::endl;
                        exit("prepare_restart", "This cannot happen");
                    } else {
                        vks_job.erase(it_set);
                    }
                }
                deallocate(arr_done);
            }
            vks_done.clear();
        }
        // vks_done4, nshift_restart4
        if (fph_rta > 0) {
            if (mympi->my_rank == 0) {
                nks_done = vks_done4.size();
            }
            MPI_Bcast(&nks_done, 1, MPI_INT, 0, MPI_COMM_WORLD);
            nshift_restart4 = nks_done;

            if (nks_done > 0) {
                allocate(arr_done, nks_done);

                if (mympi->my_rank == 0) {
                    for (i = 0; i < nks_done; ++i) {
                        arr_done[i] = vks_done4[i];
                    }
                }
                MPI_Bcast(&arr_done[0], nks_done, MPI_INT, 0, MPI_COMM_WORLD);

                // Remove vks_done elements from vks_job

                for (i = 0; i < nks_done; ++i) {

                    const auto it_set = vks_job4.find(arr_done[i]);

                    if (it_set == vks_job4.end()) {
                        std::cout << " rank = " << mympi->my_rank
                        << " arr_done = " << arr_done[i] << std::endl;
                        exit("prepare_restart", "This cannot happen");
                    } else {
                        vks_job4.erase(it_set);
                    }
                }
                deallocate(arr_done);
            }
            vks_done.clear();
        }
    }
}

void Conductivity::calc_anharmonic_imagself3()
{
    unsigned int i;
    unsigned int *nks_thread = nullptr;
    double *damping3_loc = nullptr;

    // Distribute (k,s) to individual MPI threads

    const auto nks_g = vks_job.size();
    vks_l.clear();

    unsigned int icount = 0;

    for (const auto &it: vks_job) {
        if (icount % mympi->nprocs == mympi->my_rank) {
            vks_l.push_back(it);
        }
        ++icount;
    }

    if (mympi->my_rank == 0) {
        allocate(nks_thread, mympi->nprocs);
    }

    auto nks_tmp = vks_l.size();
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

        deallocate(nks_thread);
    }

    unsigned int nk_tmp;

    if (nks_g % mympi->nprocs != 0) {
        nk_tmp = nks_g / mympi->nprocs + 1;
    } else {
        nk_tmp = nks_g / mympi->nprocs;
    }

    if (vks_l.size() < nk_tmp) {
        vks_l.push_back(-1);
    }

    allocate(damping3_loc, ntemp);

    for (i = 0; i < nk_tmp; ++i) {

        const auto iks = vks_l[i];

        if (iks == -1) {

            for (unsigned int j = 0; j < ntemp; ++j) damping3_loc[j] = eps; // do nothing

        } else {

            const auto knum = dos->kmesh_dos->kpoint_irred_all[iks / ns][0].knum;
            const auto snum = iks % ns;

            const auto omega = dos->dymat_dos->get_eigenvalues()[knum][snum];

            if (integration->ismear >= 0) {
                anharmonic_core->calc_damping_smearing(ntemp,
                                                       temperature,
                                                       omega,
                                                       iks / ns,
                                                       snum,
                                                       dos->kmesh_dos,
                                                       dos->dymat_dos->get_eigenvalues(),
                                                       dos->dymat_dos->get_eigenvectors(),
                                                       damping3_loc);
            } else if (integration->ismear == -1) {
                anharmonic_core->calc_damping_tetrahedron(ntemp,
                                                          temperature,
                                                          omega,
                                                          iks / ns,
                                                          snum,
                                                          dos->kmesh_dos,
                                                          dos->dymat_dos->get_eigenvalues(),
                                                          dos->dymat_dos->get_eigenvectors(),
                                                          damping3_loc);
            }
        }

        MPI_Gather(&damping3_loc[0], ntemp, MPI_DOUBLE,
                   damping3[nshift_restart + i * mympi->nprocs], ntemp,
                   MPI_DOUBLE, 0, MPI_COMM_WORLD);

        if (mympi->my_rank == 0) {
            write_result_gamma(i, nshift_restart, vel, damping3, 1);
            std::cout << " MODE " << std::setw(5) << i + 1 << " done." << std::endl << std::flush;
        }
    }
    deallocate(damping3_loc);

}

void Conductivity::calc_anharmonic_imagself4()
{
    unsigned int i;
    unsigned int *nks_thread;

    // Distribute (k,s) to individual MPI threads

    const auto nks_g = vks_job4.size();
    vks_l.clear();

    unsigned int icount = 0;

    for (const auto &it: vks_job4) {
        if (icount % mympi->nprocs == mympi->my_rank) {
            vks_l.push_back(it);
        }
        ++icount;
    }

    if (mympi->my_rank == 0) {
        allocate(nks_thread, mympi->nprocs);
    }

    double *damping4_loc;

    auto nks_tmp = vks_l.size();
    MPI_Gather(&nks_tmp, 1, MPI_UNSIGNED, &nks_thread[mympi->my_rank],
               1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

    if (mympi->my_rank == 0) {
        std::cout << std::endl;
        std::cout << " Computing 4-phonon scattering amplitude ... " << std::endl;
        std::cout << " WARNING: This is very very expensive!! Please be patient." << std::endl;
        std::cout << " Total Number of phonon modes to be calculated : " << nks_g << std::endl;
        std::cout << " All modes are distributed to MPI thread/s as the following :" << std::endl;
        for (i = 0; i < mympi->nprocs; ++i) {
            std::cout << " RANK: " << std::setw(5) << i + 1;
            std::cout << std::setw(8) << "MODES: " << std::setw(5) << nks_thread[i] << std::endl;
        }
        std::cout << std::endl << std::flush;

        deallocate(nks_thread);
    }

    unsigned int nk_tmp;

    if (nks_g % mympi->nprocs != 0) {
        nk_tmp = nks_g / mympi->nprocs + 1;
    } else {
        nk_tmp = nks_g / mympi->nprocs;
    }

    if (vks_l.size() < nk_tmp) {
        vks_l.push_back(-1);
    }

    allocate(damping4_loc, ntemp);

    for (i = 0; i < nk_tmp; ++i) {

        const auto iks = vks_l[i];

        if (iks == -1) {

            for (unsigned int j = 0; j < ntemp; ++j) damping4_loc[j] = eps; // do nothing

        } else {

            const auto knum = kmesh_4ph->kpoint_irred_all[iks / ns][0].knum;
            const auto snum = iks % ns;

            const auto omega = dymat_4ph->get_eigenvalues()[knum][snum];

            if (integration->ismear_4ph == 0 || integration->ismear_4ph == 1) {
                anharmonic_core->calc_damping4_smearing(ntemp,
                                                        temperature,
                                                        omega,
                                                        iks / ns,
                                                        snum,
                                                        kmesh_4ph,
                                                        dymat_4ph->get_eigenvalues(),
                                                        dymat_4ph->get_eigenvectors(),
                                                        phase_storage_4ph,
                                                        damping4_loc);

            } else if (integration->ismear_4ph == -1) {
//                anharmonic_core->calc_damping_tetrahedron(ntemp,
//                                                          Temperature,
//                                                          omega,
//                                                          iks / ns,
//                                                          snum,
//                                                          damping4_loc);
            }
        }

        MPI_Gather(&damping4_loc[0], ntemp, MPI_DOUBLE,
                   damping4[nshift_restart4 + i * mympi->nprocs], ntemp,
                   MPI_DOUBLE, 0, MPI_COMM_WORLD);

        if (mympi->my_rank == 0) {
            write_result_gamma(i, nshift_restart4, vel, damping4, -1);

            //for (auto j = 0; j < ntemp; ++j) {
            //    std::cout << std::scientific << std::setprecision(15)
            //    << damping4[nshift_restart4 + i * mympi->nprocs][j] << std::endl;
            //}
            std::cout << " MODE " << std::setw(5) << i + 1 << " done." << std::endl << std::flush;
        }

    }
    deallocate(damping4_loc);
}

void Conductivity::calc_anharmonic_imagself()
{

    if (fph_rta < 2) {
        calc_anharmonic_imagself3();
    }
    
    if (fph_rta > 0) {
        calc_anharmonic_imagself4();
    }
}

void Conductivity::write_result_gamma(const unsigned int ik,
                                      const unsigned int nshift,
                                      double ***vel_in,
                                      double **damp_in,
                                      int mode) const
{
    const unsigned int np = mympi->nprocs;
    unsigned int k;

    if (mode == 1) {
        // damping 3
        for (unsigned int j = 0; j < np; ++j) {

            const auto iks_g = ik * np + j + nshift;

            if (iks_g >= dos->kmesh_dos->nk_irred * ns) break;

            writes->fs_result << "#GAMMA_EACH" << std::endl;
            writes->fs_result << iks_g / ns + 1 << " " << iks_g % ns + 1 << std::endl;

            const auto nk_equiv = dos->kmesh_dos->kpoint_irred_all[iks_g / ns].size();

            writes->fs_result << nk_equiv << std::endl;
            for (k = 0; k < nk_equiv; ++k) {
                const auto ktmp = dos->kmesh_dos->kpoint_irred_all[iks_g / ns][k].knum;
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
    } else if (mode == -1 ) {
        // damping 4
        for (unsigned int j = 0; j < np; ++j) {

            const auto iks_g = ik * np + j + nshift;

            if (iks_g >= kmesh_4ph->nk_irred * ns) break;

            writes->fs_result4 << "#GAMMA_EACH" << std::endl;
            writes->fs_result4 << iks_g / ns + 1 << " " << iks_g % ns + 1 << std::endl;

            const auto nk_equiv = kmesh_4ph->kpoint_irred_all[iks_g / ns].size();

            writes->fs_result4 << nk_equiv << std::endl;
            for (k = 0; k < nk_equiv; ++k) {
                const auto ktmp = kmesh_4ph->kpoint_irred_all[iks_g / ns][k].knum;
                writes->fs_result4 << std::setw(15) << vel_in[ktmp][iks_g % ns][0];
                writes->fs_result4 << std::setw(15) << vel_in[ktmp][iks_g % ns][1];
                writes->fs_result4 << std::setw(15) << vel_in[ktmp][iks_g % ns][2] << std::endl;
            }

            for (k = 0; k < ntemp; ++k) {
                writes->fs_result4 << std::setw(15)
                                << damp_in[iks_g][k] * Hz_to_kayser / time_ry << std::endl;
            }
            writes->fs_result4 << "#END GAMMA_EACH" << std::endl;
        }
    }

}

void Conductivity::write_result_gamma(const unsigned int ik,
                                      const unsigned int nshift,
                                      double ***vel_in,
                                      double **damp_in) const
{
    const unsigned int np = mympi->nprocs;
    unsigned int k;

    for (unsigned int j = 0; j < np; ++j) {

        const auto iks_g = ik * np + j + nshift;

        if (iks_g >= dos->kmesh_dos->nk_irred * ns) break;

        writes->fs_result << "#GAMMA_EACH" << std::endl;
        writes->fs_result << iks_g / ns + 1 << " " << iks_g % ns + 1 << std::endl;

        const auto nk_equiv = dos->kmesh_dos->kpoint_irred_all[iks_g / ns].size();

        writes->fs_result << nk_equiv << std::endl;
        for (k = 0; k < nk_equiv; ++k) {
            const auto ktmp = dos->kmesh_dos->kpoint_irred_all[iks_g / ns][k].knum;
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
    unsigned int iks;

    if (mympi->my_rank == 0) {

        std::string file_kl;
        std::ofstream ofs_kl;
        double damp_tmp;

        double **lifetime;
        double **gamma_total;

        allocate(lifetime, dos->kmesh_dos->nk_irred * ns, ntemp);
        allocate(gamma_total, dos->kmesh_dos->nk_irred * ns, ntemp);

        for (iks = 0; iks < dos->kmesh_dos->nk_irred * ns; ++iks) {
            for (i = 0; i < ntemp; ++i) {
                gamma_total[iks][i] = damping3[iks][i];
                if (fph_rta == 1) {
                    gamma_total[iks][i] += damping4[iks][i];
                       }
            }
        }
         average_self_energy_at_degenerate_point(dos->kmesh_dos->nk_irred * ns,
                                                ntemp,
                                                dos->kmesh_dos,
                                                dos->dymat_dos->get_eigenvalues(),
                                                damping3);

               if (isotope->include_isotope) {
            for (iks = 0; iks < dos->kmesh_dos->nk_irred * ns; ++iks) {
                const auto snum = iks % ns;
                if (dynamical->is_imaginary[iks / ns][snum]) {
                    for (i = 0; i < ntemp; ++i) {
                        lifetime[iks][i] = 0.0;
                        gamma_total[iks][i] = 1.0e+100; // very big number
                    }
                } else {
                    for (i = 0; i < ntemp; ++i) {
                        //damp_tmp = damping3[iks][i] + isotope->gamma_isotope[iks / ns][snum];
                        gamma_total[iks][i] += isotope->gamma_isotope[iks / ns][snum];
                        damp_tmp = gamma_total[iks][i];
                        if (damp_tmp > 1.0e-100) {
                            lifetime[iks][i] = 1.0e+12 * time_ry * 0.5 / damp_tmp;
                        } else {
                            lifetime[iks][i] = 0.0;
                        }
                    }
                }
            }
        } else {
            for (iks = 0; iks < dos->kmesh_dos->nk_irred * ns; ++iks) {

                if (dynamical->is_imaginary[iks / ns][iks % ns]) {
                    for (i = 0; i < ntemp; ++i) {
                        lifetime[iks][i] = 0.0;
                        gamma_total[iks][i] = 1.0e+100; // very big number
                    }
                } else {
                    for (i = 0; i < ntemp; ++i) {
                        //damp_tmp = damping3[iks][i];
                        damp_tmp = gamma_total[iks][i];
                        if (damp_tmp > 1.0e-100) {
                            lifetime[iks][i] = 1.0e+12 * time_ry * 0.5 / damp_tmp;
                        } else {
                            lifetime[iks][i] = 0.0;
                            gamma_total[iks][i] = 1.0e+100;
                        }
                    }
                }
            }
        }

        allocate(kappa, ntemp, 3, 3);

        if (calc_kappa_spec) {
            allocate(kappa_spec, dos->n_energy, ntemp, 3);
        }

        compute_kappa_intraband(dos->kmesh_dos,
                                dos->dymat_dos->get_eigenvalues(),
                                lifetime,
                                kappa,
                                kappa_spec);
        deallocate(lifetime);

        if (calc_coherent) {
            allocate(kappa_coherent, ntemp, 3, 3);
            compute_kappa_coherent(dos->kmesh_dos,
                                   dos->dymat_dos->get_eigenvalues(),
                                   gamma_total,
                                   kappa_coherent);
        }

        deallocate(gamma_total);
    }
}

void Conductivity::average_self_energy_at_degenerate_point(const int n,
                                                           const int m,
                                                           const KpointMeshUniform *kmesh_in,
                                                           const double *const *eval_in,
                                                           double **damping) const
{
    int j, k, l;
    const auto nkr = kmesh_in->nk_irred;

    double *eval_tmp;
    const auto tol_omega = 1.0e-7; // Approximately equal to 0.01 cm^{-1}

    std::vector<int> degeneracy_at_k;

    allocate(eval_tmp, ns);

    double *damping_sum;

    allocate(damping_sum, m);

    for (auto i = 0; i < nkr; ++i) {
        const auto ik = kmesh_in->kpoint_irred_all[i][0].knum;

        for (j = 0; j < ns; ++j) eval_tmp[j] = eval_in[ik][j];

        degeneracy_at_k.clear();

        auto omega_prev = eval_tmp[0];
        auto ideg = 1;

        for (j = 1; j < ns; ++j) {
            const auto omega_now = eval_tmp[j];

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
    deallocate(damping_sum);
}

void Conductivity::compute_kappa_intraband(const KpointMeshUniform *kmesh_in,
                                           const double *const *eval_in,
                                           const double *const *lifetime,
                                           double ***kappa_intra,
                                           double ***kappa_spec_out) const
{
    int i, is, ik;
    double ****kappa_mode;
    const auto factor_toSI = 1.0e+18 / (std::pow(Bohr_in_Angstrom, 3) * system->volume_p);

    const auto nk_irred = kmesh_in->nk_irred;
    allocate(kappa_mode, ntemp, 9, ns, nk_irred);

    for (i = 0; i < ntemp; ++i) {
        for (unsigned int j = 0; j < 3; ++j) {
            for (unsigned int k = 0; k < 3; ++k) {

                if (temperature[i] < eps) {
                    // Set kappa as zero when T = 0.
                    for (is = 0; is < ns; ++is) {
                        for (ik = 0; ik < nk_irred; ++ik) {
                            kappa_mode[i][3 * j + k][is][ik] = 0.0;
                        }
                    }
                } else {
                    for (is = 0; is < ns; ++is) {
                        for (ik = 0; ik < nk_irred; ++ik) {
                            const auto knum = kmesh_in->kpoint_irred_all[ik][0].knum;
                            const auto omega = eval_in[knum][is];
                            //std::cout << omega << std::endl;
                            auto vv_tmp = 0.0;
                            const auto nk_equiv = kmesh_in->kpoint_irred_all[ik].size();

                            // Accumulate group velocity (diad product) for the reducible k points
                            for (auto ieq = 0; ieq < nk_equiv; ++ieq) {
                                const auto ktmp = kmesh_in->kpoint_irred_all[ik][ieq].knum;
                                vv_tmp += vel[ktmp][is][j] * vel[ktmp][is][k];
                            }

                            if (thermodynamics->classical) {
                                kappa_mode[i][3 * j + k][is][ik]
                                      = thermodynamics->Cv_classical(omega, temperature[i])
                                      * vv_tmp * lifetime[ns * ik + is][i];
                            } else {
                                kappa_mode[i][3 * j + k][is][ik]
                                      = thermodynamics->Cv(omega, temperature[i])
                                      * vv_tmp * lifetime[ns * ik + is][i];
                            }

                            // Convert to SI unit
                            kappa_mode[i][3 * j + k][is][ik] *= factor_toSI;

                        }
                    }
                }

                kappa_intra[i][j][k] = 0.0;

                for (is = 0; is < ns; ++is) {
                    for (ik = 0; ik < nk_irred; ++ik) {
                        kappa_intra[i][j][k] += kappa_mode[i][3 * j + k][is][ik];
                    }
                }

                kappa_intra[i][j][k] /= static_cast<double>(nk);
            }
        }
    }

    if (calc_kappa_spec) {
        //allocate(kappa_spec_out, dos->n_energy, ntemp, 3);
        compute_frequency_resolved_kappa(ntemp,
                                         integration->ismear,
                                         dos->kmesh_dos,
                                         dos->dymat_dos->get_eigenvalues(),
                                         kappa_mode,
                                         kappa_spec_out);
    }

    deallocate(kappa_mode);
}

void Conductivity::compute_kappa_coherent(const KpointMeshUniform *kmesh_in,
                                          const double *const *eval_in,
                                          const double *const *gamma_total,
                                          double ***kappa_coherent_out) const
{
    // Compute the coherent part of thermal conductivity
    // based on the Michelle's paper.
    int ib;
    const auto factor_toSI = 1.0e+18 / (std::pow(Bohr_in_Angstrom, 3) * system->volume_p);
    const auto common_factor = factor_toSI * 1.0e+12 * time_ry / static_cast<double>(nk);
    const auto common_factor_output = factor_toSI * 1.0e+12 * time_ry;
    const int ns2 = ns * ns;
    const auto czero = std::complex<double>(0.0, 0.0);
    std::vector<std::complex<double>> kappa_tmp(ns2, czero);
    std::complex<double> **kappa_save = nullptr;

    const auto nk_irred = kmesh_in->nk_irred;

    std::ofstream ofs;
    if (calc_coherent == 2) {
        ofs.open(file_coherent_elems.c_str(), std::ios::out);
        if (!ofs) exit("compute_kappa_coherent", "cannot open file_kc");
        ofs << "# Temperature [K], 1st and 2nd xyz components, ibranch, jbranch, ik_irred, "
               "omega1 [cm^-1], omega2 [cm^-1], kappa_elems real, kappa_elems imag" << std::endl;
        allocate(kappa_save, ns2, nk_irred);
    }

    for (auto i = 0; i < ntemp; ++i) {
        for (unsigned int j = 0; j < 3; ++j) {
            for (unsigned int k = 0; k < 3; ++k) {

                kappa_coherent_out[i][j][k] = 0.0;

                if (temperature[i] > eps) {
#pragma omp parallel for
                    for (ib = 0; ib < ns2; ++ib) {
                        kappa_tmp[ib] = czero;
                        const int is = ib / ns;
                        const int js = ib % ns;

                        if (js == is) continue; // skip the diagonal component

                        for (auto ik = 0; ik < nk_irred; ++ik) {
                            const auto knum = kmesh_in->kpoint_irred_all[ik][0].knum;
                            const auto omega1 = eval_in[knum][is];
                            const auto omega2 = eval_in[knum][js];

                            if (omega1 < eps8 || omega2 < eps8) continue;
                            auto vv_tmp = czero;
                            const auto nk_equiv = kmesh_in->kpoint_irred_all[ik].size();

                            // Accumulate group velocity (diad product) for the reducible k points
                            for (auto ieq = 0; ieq < nk_equiv; ++ieq) {
                                const auto ktmp = kmesh_in->kpoint_irred_all[ik][ieq].knum;
                                vv_tmp += velmat[ktmp][is][js][j] * velmat[ktmp][js][is][k];
                            }
                            auto kcelem_tmp = 2.0 * (omega1 * omega2) / (omega1 + omega2)
                                  * (thermodynamics->Cv(omega1, temperature[i]) / omega1
                                        + thermodynamics->Cv(omega2, temperature[i]) / omega2)
                                  * 2.0 * (gamma_total[ik * ns + is][i] + gamma_total[ik * ns + js][i])
                                  / (4.0 * std::pow(omega1 - omega2, 2.0)
                                        + 4.0 * std::pow(gamma_total[ik * ns + is][i]
                                                               + gamma_total[ik * ns + js][i], 2.0))
                                  * vv_tmp;
                            kappa_tmp[ib] += kcelem_tmp;

                            if (calc_coherent == 2 && j == k) {
                                kappa_save[ib][ik] = kcelem_tmp * common_factor_output;
                            }
                        }
                    } // end OpenMP parallelization over ib

                    for (ib = 0; ib < ns2; ++ib) {
                        if (std::abs(kappa_tmp[ib].imag()) > eps10) {
                            warn("compute_kappa_coherent",
                                 "The kappa_coherent_out has imaginary component.");
                        }
                        kappa_coherent_out[i][j][k] += kappa_tmp[ib].real();
                    }

                    if (calc_coherent == 2 && j == k) {
                        for (ib = 0; ib < ns2; ++ib) {

                            const int is = ib / ns;
                            const int js = ib % ns;

                            for (auto ik = 0; ik < nk_irred; ++ik) {
                                if (is == js) kappa_save[ib][ik] = czero;

                                ofs << std::setw(5) << temperature[i];
                                ofs << std::setw(3) << j + 1 << std::setw(3) << k + 1;
                                ofs << std::setw(4) << is + 1;
                                ofs << std::setw(4) << js + 1;
                                ofs << std::setw(6) << ik + 1;
                                const auto knum = kmesh_in->kpoint_irred_all[ik][0].knum;
                                const auto omega1 = eval_in[knum][is];
                                const auto omega2 = eval_in[knum][js];
                                ofs << std::setw(15) << writes->in_kayser(omega1);
                                ofs << std::setw(15) << writes->in_kayser(omega2);
                                ofs << std::setw(15) << kappa_save[ib][ik].real();
                                ofs << std::setw(15) << kappa_save[ib][ik].imag();
                                ofs << '\n';
                            }
                        }
                        ofs << '\n';
                    }
                }
                kappa_coherent_out[i][j][k] *= common_factor;
            }
        }
    }

    if (calc_coherent == 2) {
        ofs.close();
        deallocate(kappa_save);
    }
}

void Conductivity::compute_frequency_resolved_kappa(const int ntemp,
                                                    const int smearing_method,
                                                    const KpointMeshUniform *kmesh_in,
                                                    const double *const *eval_in,
                                                    const double *const *const *const *kappa_mode,
                                                    double ***kappa_spec_out) const
{
    int i, j;
    unsigned int *kmap_identity;
    double **eval;

    std::cout << std::endl;
    std::cout << " KAPPA_SPEC = 1 : Calculating thermal conductivity spectra ... ";

    allocate(kmap_identity, nk);
    allocate(eval, ns, nk);

    for (i = 0; i < nk; ++i) kmap_identity[i] = i;

    for (i = 0; i < nk; ++i) {
        for (j = 0; j < ns; ++j) {
            eval[j][i] = writes->in_kayser(eval_in[i][j]);
        }
    }

#ifdef _OPENMP
#pragma omp parallel private (j)
#endif
    {
        int k;
        int knum;
        double *weight;
        allocate(weight, nk);

#ifdef _OPENMP
#pragma omp for
#endif
        for (i = 0; i < dos->n_energy; ++i) {

            for (j = 0; j < ntemp; ++j) {
                for (k = 0; k < 3; ++k) {
                    kappa_spec_out[i][j][k] = 0.0;
                }
            }

            for (int is = 0; is < ns; ++is) {
                if (smearing_method == -1) {
                    integration->calc_weight_tetrahedron(nk, kmap_identity,
                                                         eval[is], dos->energy_dos[i],
                                                         dos->tetra_nodes_dos->get_ntetra(),
                                                         dos->tetra_nodes_dos->get_tetras(),
                                                         weight);
                } else {
                    integration->calc_weight_smearing(nk, nk, kmap_identity,
                                                      eval[is], dos->energy_dos[i],
                                                      smearing_method, weight);
                }

                for (j = 0; j < ntemp; ++j) {
                    for (k = 0; k < 3; ++k) {
                        for (int ik = 0; ik < kmesh_in->nk_irred; ++ik) {
                            knum = kmesh_in->kpoint_irred_all[ik][0].knum;
                            kappa_spec_out[i][j][k] += kappa_mode[j][3 * k + k][is][ik] * weight[knum];
                        }
                    }
                }
            }
        }
        deallocate(weight);
    }

    deallocate(kmap_identity);
    deallocate(eval);

    std::cout << " done!" << std::endl;
}

void Conductivity::set_kmesh_coarse(const unsigned int *nk_in)
{
    for (auto i = 0; i < 3; ++i) nk_coarse[i] = nk_in[i];
}

KpointMeshUniform *Conductivity::get_kmesh_coarse() const
{
    return kmesh_4ph;
}
