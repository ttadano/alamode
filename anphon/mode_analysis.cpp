/*
 mode_analysis.cpp

Copyright (c) 2018 Terumasa Tadano

This file is distributed under the terms of the MIT license.
Please see the file 'LICENCE.txt' in the root directory
or http://opensource.org/licenses/mit-license.php for information.
*/

#include "mpi_common.h"
#include "mode_analysis.h"
#include "anharmonic_core.h"
#include "dynamical.h"
#include "error.h"
#include "kpoint.h"
#include "mathfunctions.h"
#include "memory.h"
#include "integration.h"
#include "parsephon.h"
#include "phonon_dos.h"
#include "selfenergy.h"
#include "symmetry_core.h"
#include "system.h"
#include "thermodynamics.h"
#include "write_phonons.h"
#include <iostream>
#include <iomanip>
#include <fstream>


using namespace PHON_NS;

ModeAnalysis::ModeAnalysis(PHON *phon) : Pointers(phon)
{
    set_default_variables();
}

ModeAnalysis::~ModeAnalysis()
{
    deallocate_variables();
}

void ModeAnalysis::set_default_variables()
{
    ks_analyze_mode = false;
    calc_imagpart = true;
    calc_realpart = false;
    calc_fstate_omega = false;
    calc_fstate_k = false;
    print_V3 = false;
    spectral_func = false;
}

void ModeAnalysis::deallocate_variables() {}


void ModeAnalysis::setup_mode_analysis()
{
    // Judge if ks_analyze_mode should be turned on or not.

    unsigned int i;

    if (mympi->my_rank == 0) {
        if (!ks_input.empty()) {
            ks_analyze_mode = true;

            std::cout << std::endl;
            std::cout << " KS_INPUT-tag is given : Analysis on the specified phonon modes" << std::endl;
            std::cout << " will be performed instead of thermal conductivity calculation." << std::endl;
            std::cout << std::endl;

            std::ifstream ifs_ks;
            ifs_ks.open(ks_input.c_str(), std::ios::in);
            if (!ifs_ks)
                error->exit("setup_mode_analysis",
                            "Cannot open file KS_INPUT");

            unsigned int nlist;
            double ktmp[3];
            unsigned int snum_tmp;

            ifs_ks >> nlist;

            if (nlist <= 0)
                error->exit("setup_mode_analysis",
                            "First line in KS_INPUT files should be a positive integer.");

            if (calc_fstate_k) {
                kslist_fstate_k.clear();

                for (i = 0; i < nlist; ++i) {
                    ifs_ks >> ktmp[0] >> ktmp[1] >> ktmp[2] >> snum_tmp;

                    if (snum_tmp <= 0 || snum_tmp > dynamical->neval) {
                        error->exit("setup_mode_analysis",
                                    "Mode index out of range.");
                    }

                    kslist_fstate_k.emplace_back(ktmp, snum_tmp - 1);
                }
                std::cout << " The number of entries = "
                    << kslist_fstate_k.size() << std::endl;
            } else {
                kslist.clear();
                for (i = 0; i < nlist; ++i) {
                    ifs_ks >> ktmp[0] >> ktmp[1] >> ktmp[2] >> snum_tmp;
                    int knum_tmp = kpoint->get_knum(ktmp[0], ktmp[1], ktmp[2]);

                    if (knum_tmp == -1)
                        error->exit("setup_mode_analysis",
                                    "Given kpoint does not exist in given k-point grid.");
                    if (snum_tmp <= 0 || snum_tmp > dynamical->neval) {
                        error->exit("setup_mode_analysis", "Mode index out of range.");
                    }
                    kslist.push_back(knum_tmp * dynamical->neval + snum_tmp - 1);
                }
                std::cout << " The number of entries = " << kslist.size() << std::endl;
            }

            ifs_ks.close();
        } else {
            ks_analyze_mode = false;
        }
    }

    MPI_Bcast(&ks_analyze_mode, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD);
    MPI_Bcast(&calc_realpart, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD);
    MPI_Bcast(&calc_fstate_k, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD);
    MPI_Bcast(&calc_fstate_omega, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD);

    unsigned int nlist;

    if (kpoint->kpoint_mode == 3) {

        if (!ks_analyze_mode) {
            error->exit("setup_mode_analysis", "KPMODE = 3 must be used with FSTATE_K = 1");
        }
        int j;
        double **vec_tmp;
        unsigned int *mode_tmp;


        // Broadcast kslist_fstate_k

        nlist = kslist_fstate_k.size();
        MPI_Bcast(&nlist, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

        memory->allocate(vec_tmp, nlist, 3);
        memory->allocate(mode_tmp, nlist);

        if (mympi->my_rank == 0) {
            for (i = 0; i < nlist; ++i) {
                for (j = 0; j < 3; ++j) {
                    vec_tmp[i][j] = kslist_fstate_k[i].xk[j];
                }
                mode_tmp[i] = kslist_fstate_k[i].nmode;
            }
        }

        MPI_Bcast(&vec_tmp[0][0], 3 * nlist, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&mode_tmp[0], nlist, MPI_INT, 0, MPI_COMM_WORLD);

        if (mympi->my_rank > 0) {
            kslist_fstate_k.clear();

            for (i = 0; i < nlist; ++i) {
                kslist_fstate_k.emplace_back(vec_tmp[i], mode_tmp[i]);
            }
        }

        memory->deallocate(vec_tmp);
        memory->deallocate(mode_tmp);

    } else {

        unsigned int *kslist_arr;
        nlist = kslist.size();

        // Broadcast kslist

        MPI_Bcast(&nlist, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
        memory->allocate(kslist_arr, nlist);

        if (mympi->my_rank == 0) {
            for (i = 0; i < nlist; ++i) kslist_arr[i] = kslist[i];
        }
        MPI_Bcast(&kslist_arr[0], nlist, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

        if (mympi->my_rank > 0) {
            kslist.clear();
            for (i = 0; i < nlist; ++i) kslist.push_back(kslist_arr[i]);
        }
        memory->deallocate(kslist_arr);
    }


    if (ks_analyze_mode) {
        if (kpoint->kpoint_mode == 2 && anharmonic_core->use_triplet_symmetry) {
            anharmonic_core->use_triplet_symmetry = false;
            if (mympi->my_rank == 0) {
                std::cout << std::endl;
                std::cout << " TRISYM was automatically set to 0." << std::endl;
                std::cout << std::endl;
            }
        }

        if (anharmonic_core->quartic_mode > 0) {
            // This is for quartic vertexes.

            if (mympi->my_rank == 0) {
                std::cout << " QUARTIC = 1 : Frequency shift due to the loop diagram associated with" << std::endl;
                std::cout << "               quartic anharmonicity will be calculated." << std::endl;
                std::cout << "               Please check the accuracy of the quartic IFCs " << std::endl;
                std::cout << "               before doing serious calculations." << std::endl;
                std::cout << std::endl;
            }
        }

        if (calc_realpart && integration->ismear != 0) {
            error->exit("setup_mode_analysis",
                        "Sorry. REALPART = 1 can be used only with ISMEAR = 0");
        }

        if (spectral_func && integration->ismear != -1) {
            error->exit("setup_mode_analysis",
                        "Sorry. SELF_W = 1 can be used only with the tetrahedron method (ISMEAR = -1).");
        }

        if (calc_fstate_k && kpoint->kpoint_mode != 3) {
            error->exit("setup_mode_analysis",
                        "KPMODE should be 3 when FSTATE_K = 1.");
        }
        if (!calc_fstate_k && kpoint->kpoint_mode == 3) {
            error->exit("setup_mode_analysis",
                        "KPMODE = 3 works only when FSTATE_K = 1");
        }

        if (calc_fstate_k && (calc_fstate_omega || print_V3 || spectral_func || calc_realpart)) {
            error->warn("setup_mode_analysis",
                        "FSTATE_K = 1 shouldn't be set with the followings: PRINTV3=1, REALPART=1, FSTATE_W=1, SELF_W=1");
        }

        dynamical->modify_eigenvectors();
    }
}

void ModeAnalysis::run_mode_analysis()
{
    double Tmax = system->Tmax;
    double Tmin = system->Tmin;
    double dT = system->dT;
    double *T_arr;

    unsigned int NT = static_cast<unsigned int>((Tmax - Tmin) / dT) + 1;
    memory->allocate(T_arr, NT);
    for (unsigned int i = 0; i < NT; ++i) T_arr[i] = Tmin + static_cast<double>(i) * dT;

    double epsilon = integration->epsilon;

    if (calc_fstate_k) {

        // Momentum-resolved final state amplitude
        print_momentum_resolved_final_state(NT, T_arr, epsilon);

    } else {

        print_selfenergy(NT, T_arr);

        if (print_V3) print_V3_elements();
        //        if (print_V3) print_Phi3_elements();

        if (calc_fstate_omega) print_frequency_resolved_final_state(NT, T_arr);

        if (spectral_func) print_spectral_function(NT, T_arr);

    }

    memory->deallocate(T_arr);
}


void ModeAnalysis::print_selfenergy(const int NT,
                                    double *T_arr)
{
    int ns = dynamical->neval;
    int knum, snum;
    int i, j;
    double omega;
    double *damping_a;
    double omega_shift;

    std::ofstream ofs_linewidth, ofs_shift;
    std::string file_linewidth, file_shift;

    std::complex<double> *self_tadpole;
    std::complex<double> *self_a, *self_b, *self_c, *self_d, *self_e;
    std::complex<double> *self_f, *self_g, *self_h, *self_i, *self_j;

    damping_a = nullptr;
    self_tadpole = nullptr;
    self_a = nullptr;
    self_b = nullptr;
    self_c = nullptr;
    self_d = nullptr;
    self_e = nullptr;
    self_f = nullptr;
    self_g = nullptr;
    self_h = nullptr;
    self_i = nullptr;
    self_j = nullptr;

    if (mympi->my_rank == 0) {
        std::cout << std::endl;
        std::cout << " Calculate the line width (FWHM) of phonons" << std::endl;
        std::cout << " due to 3-phonon interactions for given "
            << kslist.size() << " modes." << std::endl;

        if (calc_realpart) {
            if (anharmonic_core->quartic_mode == 1) {
                std::cout << " REALPART = 1 and " << std::endl;
                std::cout << " QUARTIC  = 1     : Additionally, frequency shift of phonons due to 3-phonon" << std::
                    endl;
                std::cout << "                    and 4-phonon interactions will be calculated." << std::endl;
            } else {
                std::cout << " REALPART = 1 : Additionally, frequency shift of phonons due to 3-phonon" << std::
                    endl;
                std::cout << "                interactions will be calculated." << std::endl;
            }
        }

        if (anharmonic_core->quartic_mode == 2) {
            std::cout << std::endl;
            std::cout << " QUARTIC = 2 : Additionally, phonon line width due to 4-phonon" << std::endl;
            std::cout << "               interactions will be calculated." << std::endl;
            std::cout << " WARNING     : This is very very expensive." << std::endl;
        }
    }

    memory->allocate(damping_a, NT);
    memory->allocate(self_a, NT);
    memory->allocate(self_b, NT);
    memory->allocate(self_tadpole, NT);

    if (anharmonic_core->quartic_mode == 2) {
        memory->allocate(self_c, NT);
        memory->allocate(self_d, NT);
        memory->allocate(self_e, NT);
        memory->allocate(self_f, NT);
        memory->allocate(self_g, NT);
        memory->allocate(self_h, NT);
        memory->allocate(self_i, NT);
        memory->allocate(self_j, NT);
    }

    for (i = 0; i < kslist.size(); ++i) {
        knum = kslist[i] / ns;
        snum = kslist[i] % ns;

        omega = dynamical->eval_phonon[knum][snum];

        if (mympi->my_rank == 0) {
            std::cout << std::endl;
            std::cout << " Number : " << std::setw(5) << i + 1 << std::endl;
            std::cout << "  Phonon at k = (";
            for (j = 0; j < 3; ++j) {
                std::cout << std::setw(10) << std::fixed << kpoint->xk[knum][j];
                if (j < 2) std::cout << ",";
            }
            std::cout << ")" << std::endl;
            std::cout << "  Mode index = " << std::setw(5) << snum + 1 << std::endl;
            std::cout << "  Frequency (cm^-1) : "
                << std::setw(15) << writes->in_kayser(omega) << std::endl;
        }

        int ik_irred = kpoint->kmap_to_irreducible[knum];

        if (integration->ismear == -1) {
            anharmonic_core->calc_damping_tetrahedron(NT, T_arr, omega, ik_irred, snum, damping_a);
        } else {
            selfenergy->selfenergy_a(NT, T_arr, omega, knum, snum, self_a);
            for (j = 0; j < NT; ++j) damping_a[j] = self_a[j].imag();
        }
        if (anharmonic_core->quartic_mode == 2) {
            selfenergy->selfenergy_c(NT, T_arr, omega, knum, snum, self_c);
            //   selfenergy->selfenergy_d(NT, T_arr, omega, knum, snum, self_d);
            //   selfenergy->selfenergy_e(NT, T_arr, omega, knum, snum, self_e);
            //   selfenergy->selfenergy_f(NT, T_arr, omega, knum, snum, self_f);
            //                 selfenergy->selfenergy_g(NT, T_arr, omega, knum, snum, self_g);
            //                 selfenergy->selfenergy_h(NT, T_arr, omega, knum, snum, self_h);
            //                 selfenergy->selfenergy_i(NT, T_arr, omega, knum, snum, self_i);
            //                 selfenergy->selfenergy_j(NT, T_arr, omega, knum, snum, self_j);
        }

        if (mympi->my_rank == 0) {
            file_linewidth = input->job_title + ".Gamma." + std::to_string(i + 1);
            ofs_linewidth.open(file_linewidth.c_str(), std::ios::out);
            if (!ofs_linewidth)
                error->exit("print_selfenergy",
                            "Cannot open file file_linewidth");

            ofs_linewidth << "# xk = ";

            for (j = 0; j < 3; ++j) {
                ofs_linewidth << std::setw(15) << kpoint->xk[knum][j];
            }
            ofs_linewidth << std::endl;
            ofs_linewidth << "# mode = " << snum + 1 << std::endl;
            ofs_linewidth << "# Frequency = " << writes->in_kayser(omega) << std::endl;
            ofs_linewidth << "## Temperature dependence of 2*Gamma (FWHM) for the given mode" << std::endl;
            ofs_linewidth << "## T[K], 2*Gamma3 (cm^-1) (bubble)";
            if (anharmonic_core->quartic_mode == 2) ofs_linewidth << ", 2*Gamma4(cm^-1) <-- specific diagram only";
            ofs_linewidth << std::endl;

            for (j = 0; j < NT; ++j) {
                ofs_linewidth << std::setw(10) << T_arr[j]
                    << std::setw(15) << writes->in_kayser(2.0 * damping_a[j]);

                if (anharmonic_core->quartic_mode == 2) {
                    //							ofs_mode_tau << std::setw(15) << writes->in_kayser(damp4[j]);
                    ofs_linewidth << std::setw(15) << writes->in_kayser(2.0 * self_c[j].imag());
                    ofs_linewidth << std::setw(15) << writes->in_kayser(2.0 * self_d[j].imag());
                    ofs_linewidth << std::setw(15) << writes->in_kayser(2.0 * self_e[j].imag());
                    ofs_linewidth << std::setw(15) << writes->in_kayser(2.0 * self_f[j].imag());
                    ofs_linewidth << std::setw(15) << writes->in_kayser(2.0 * self_g[j].imag());
                    ofs_linewidth << std::setw(15) << writes->in_kayser(2.0 * self_h[j].imag());
                    ofs_linewidth << std::setw(15) << writes->in_kayser(2.0 * self_i[j].imag());
                    ofs_linewidth << std::setw(15) << writes->in_kayser(2.0 * self_j[j].imag());
                }

                ofs_linewidth << std::endl;
            }
            ofs_linewidth.close();
            std::cout << "  Phonon line-width is printed in " << file_linewidth << std::endl;
        }


        if (calc_realpart) {
            selfenergy->selfenergy_tadpole(NT, T_arr, omega, knum, snum, self_tadpole);
            //                selfenergy->selfenergy_a(NT, T_arr, omega, knum, snum, self_a);

            if (anharmonic_core->quartic_mode == 1) {
                selfenergy->selfenergy_b(NT, T_arr, omega, knum, snum, self_b);
            }

            if (mympi->my_rank == 0) {
                file_shift = input->job_title + ".Shift." + std::to_string(i + 1);
                ofs_shift.open(file_shift.c_str(), std::ios::out);
                if (!ofs_shift)
                    error->exit("print_selfenergy",
                                "Cannot open file file_shift");

                ofs_shift << "# xk = ";

                for (j = 0; j < 3; ++j) {
                    ofs_shift << std::setw(15) << kpoint->xk[knum][j];
                }
                ofs_shift << std::endl;
                ofs_shift << "# mode = " << snum + 1 << std::endl;
                ofs_shift << "# Frequency = " << writes->in_kayser(omega) << std::endl;
                ofs_shift << "## T[K], Shift3 (cm^-1) (tadpole), Shift3 (cm^-1) (bubble)";
                if (anharmonic_core->quartic_mode == 1) ofs_shift << ", Shift4 (cm^-1) (loop)";
                ofs_shift << ", Shifted frequency (cm^-1)";
                ofs_shift << std::endl;


                for (j = 0; j < NT; ++j) {
                    ofs_shift << std::setw(10) << T_arr[j];
                    ofs_shift << std::setw(15) << writes->in_kayser(-self_tadpole[j].real());
                    ofs_shift << std::setw(15) << writes->in_kayser(-self_a[j].real());

                    omega_shift = omega - self_tadpole[j].real() - self_a[j].real();

                    if (anharmonic_core->quartic_mode == 1) {
                        ofs_shift << std::setw(15) << writes->in_kayser(-self_b[j].real());
                        omega_shift -= self_b[j].real();
                    }
                    ofs_shift << std::setw(15) << writes->in_kayser(omega_shift);
                    ofs_shift << std::endl;
                }

                ofs_shift.close();
                std::cout << "  Phonon frequency shift is printed in " << file_shift << std::endl;
            }
        }
    }

    if (damping_a) {
        memory->deallocate(damping_a);
    }
    if (self_tadpole) {
        memory->deallocate(self_tadpole);
    }
    if (self_a) {
        memory->deallocate(self_a);
    }
    if (self_b) {
        memory->deallocate(self_b);
    }
    if (self_c) {
        memory->deallocate(self_c);
    }
    if (self_d) {
        memory->deallocate(self_d);
    }
    if (self_e) {
        memory->deallocate(self_e);
    }
    if (self_f) {
        memory->deallocate(self_f);
    }
    if (self_g) {
        memory->deallocate(self_g);
    }
    if (self_h) {
        memory->deallocate(self_h);
    }
    if (self_i) {
        memory->deallocate(self_i);
    }
    if (self_j) {
        memory->deallocate(self_j);
    }
}


void ModeAnalysis::print_frequency_resolved_final_state(const unsigned int NT,
                                                        double *T_arr)
{
    int i, j;
    unsigned int knum, snum;
    double omega, omega0;
    double ***gamma_final;
    double *freq_array;
    int ienergy;
    std::ofstream ofs_omega;
    std::string file_omega;
    int ns = dynamical->neval;

    memory->allocate(gamma_final, NT, dos->n_energy, 2);
    memory->allocate(freq_array, dos->n_energy);

    for (i = 0; i < dos->n_energy; ++i) {
        freq_array[i] = dos->energy_dos[i] * time_ry / Hz_to_kayser;
    }

    if (mympi->my_rank == 0) {
        std::cout << std::endl;
        std::cout << " FSTATE_W = 1 : Calculate the frequency-resolved final state amplitude" << std::endl;
        std::cout << "                due to 3-phonon interactions." << std::endl;
    }

    for (i = 0; i < kslist.size(); ++i) {
        knum = kslist[i] / ns;
        snum = kslist[i] % ns;

        omega0 = dynamical->eval_phonon[knum][snum];

        if (mympi->my_rank == 0) {
            std::cout << std::endl;
            std::cout << " Number : " << std::setw(5) << i + 1 << std::endl;
            std::cout << "  Phonon at k = (";
            for (j = 0; j < 3; ++j) {
                std::cout << std::setw(10) << std::fixed << kpoint->xk[knum][j];
                if (j < 2) std::cout << ",";
            }
            std::cout << ")" << std::endl;
            std::cout << "  Mode index = " << std::setw(5) << snum + 1 << std::endl;
            std::cout << "  Frequency (cm^-1) : " << std::setw(15)
                << writes->in_kayser(omega0) << std::endl;
        }

        if (integration->ismear == -1) {
            calc_frequency_resolved_final_state_tetrahedron(NT,
                                                            T_arr,
                                                            omega0,
                                                            dos->n_energy,
                                                            freq_array,
                                                            kpoint->kmap_to_irreducible[knum],
                                                            snum,
                                                            gamma_final);
        } else {
            calc_frequency_resolved_final_state(NT,
                                                T_arr,
                                                omega0,
                                                dos->n_energy,
                                                freq_array,
                                                kpoint->kmap_to_irreducible[knum],
                                                snum,
                                                gamma_final);
        }


        if (mympi->my_rank == 0) {
            file_omega = input->job_title + ".fw." + std::to_string(i + 1);
            ofs_omega.open(file_omega.c_str(), std::ios::out);
            if (!ofs_omega)
                error->exit("print_frequency_resolved_final_state",
                            "Cannot open file file_omega");

            ofs_omega << "# xk = ";

            for (j = 0; j < 3; ++j) {
                ofs_omega << std::setw(15) << kpoint->xk[knum][j];
            }
            ofs_omega << std::endl;
            ofs_omega << "# mode = " << snum << std::endl;
            ofs_omega << "# Frequency = " << writes->in_kayser(omega0) << std::endl;

            ofs_omega << "## Frequency-resolved final state amplitude for given modes" << std::endl;
            ofs_omega << "## Gamma[omega][temperature] (absorption, emission)";
            ofs_omega << std::endl;

            ofs_omega << "## ";
            for (j = 0; j < NT; ++j) {
                ofs_omega << std::setw(10) << T_arr[j];
            }
            ofs_omega << std::endl;
            for (ienergy = 0; ienergy < dos->n_energy; ++ienergy) {
                omega = dos->energy_dos[ienergy];

                ofs_omega << std::setw(10) << omega;
                for (j = 0; j < NT; ++j) {
                    ofs_omega << std::setw(15) << gamma_final[j][ienergy][1];
                    ofs_omega << std::setw(15) << gamma_final[j][ienergy][0];
                }

                ofs_omega << std::endl;
            }
            ofs_omega.close();
            std::cout << "  Frequency-resolved final state amplitude is printed in " << file_omega << std::endl;
        }
    }

    memory->deallocate(freq_array);
    memory->deallocate(gamma_final);
}


void ModeAnalysis::calc_frequency_resolved_final_state(const unsigned int N,
                                                       double *T,
                                                       const double omega0,
                                                       const unsigned int M,
                                                       const double *omega,
                                                       const unsigned int ik_in,
                                                       const unsigned int snum,
                                                       double ***ret)
{
    int i, j;

    double multi;
    int knum, knum_minus;
    int k1, k2;
    int is, js;
    unsigned int arr[3];
    double omega_inner[2];
    double v3_tmp;
    double T_tmp;
    double n1, n2;
    double f1, f2;
    double prod_tmp[2];
    double ***ret_mpi;
    int nk = kpoint->nk;
    int ns = dynamical->neval;

    double epsilon = integration->epsilon;

    std::vector<KsListGroup> triplet;

    kpoint->get_unique_triplet_k(ik_in,
                                 anharmonic_core->use_triplet_symmetry,
                                 false,
                                 triplet);
    memory->allocate(ret_mpi, N, M, 2);

    for (i = 0; i < N; ++i) {
        for (j = 0; j < M; ++j) {
            ret_mpi[i][j][0] = 0.0;
            ret_mpi[i][j][1] = 0.0;
        }
    }

    for (int ik = mympi->my_rank; ik < triplet.size(); ik += mympi->nprocs) {
        multi = static_cast<double>(triplet[ik].group.size());
        knum = kpoint->kpoint_irred_all[ik_in][0].knum;
        knum_minus = kpoint->knum_minus[knum];

        arr[0] = ns * knum_minus + snum;
        k1 = triplet[ik].group[0].ks[0];
        k2 = triplet[ik].group[0].ks[1];

        for (is = 0; is < ns; ++is) {
            for (js = 0; js < ns; ++js) {
                arr[1] = ns * k1 + is;
                arr[2] = ns * k2 + js;

                omega_inner[0] = dynamical->eval_phonon[k1][is];
                omega_inner[1] = dynamical->eval_phonon[k2][js];

                v3_tmp = std::norm(anharmonic_core->V3(arr));

                for (i = 0; i < N; ++i) {
                    T_tmp = T[i];

                    if (thermodynamics->classical) {
                        f1 = thermodynamics->fC(omega_inner[0], T_tmp);
                        f2 = thermodynamics->fC(omega_inner[1], T_tmp);
                        n1 = f1 + f2;
                        n2 = f1 - f2;
                    } else {
                        f1 = thermodynamics->fB(omega_inner[0], T_tmp);
                        f2 = thermodynamics->fB(omega_inner[1], T_tmp);
                        n1 = f1 + f2 + 1.0;
                        n2 = f1 - f2;
                    }

                    if (integration->ismear == 0) {
                        prod_tmp[0] = n1
                            * (delta_lorentz(omega0 - omega_inner[0] - omega_inner[1], epsilon)
                                - delta_lorentz(omega0 + omega_inner[0] + omega_inner[1], epsilon));
                        prod_tmp[1] = n2
                            * (delta_lorentz(omega0 + omega_inner[0] - omega_inner[1], epsilon)
                                - delta_lorentz(omega0 - omega_inner[0] + omega_inner[1], epsilon));

                        for (j = 0; j < M; ++j) {
                            ret_mpi[i][j][0] += v3_tmp * multi
                                * delta_lorentz(omega[j] - omega_inner[0], epsilon)
                                * prod_tmp[0];
                            ret_mpi[i][j][1] += v3_tmp * multi
                                * delta_lorentz(omega[j] - omega_inner[0], epsilon)
                                * prod_tmp[1];
                        }
                    } else if (integration->ismear == 1) {
                        prod_tmp[0] = n1
                            * (delta_gauss(omega0 - omega_inner[0] - omega_inner[1], epsilon)
                                - delta_gauss(omega0 + omega_inner[0] + omega_inner[1], epsilon));
                        prod_tmp[1] = n2
                            * (delta_gauss(omega0 + omega_inner[0] - omega_inner[1], epsilon)
                                - delta_gauss(omega0 - omega_inner[0] + omega_inner[1], epsilon));

                        for (j = 0; j < M; ++j) {
                            ret_mpi[i][j][0] += v3_tmp * multi
                                * delta_gauss(omega[j] - omega_inner[0], epsilon)
                                * prod_tmp[0];
                            ret_mpi[i][j][1] += v3_tmp * multi
                                * delta_gauss(omega[j] - omega_inner[0], epsilon)
                                * prod_tmp[1];
                        }
                    }
                }
            }
        }
    }
    for (i = 0; i < N; ++i) {
        for (j = 0; j < M; ++j) {
            ret_mpi[i][j][0] *= pi * std::pow(0.5, 4) / static_cast<double>(nk);
            ret_mpi[i][j][1] *= pi * std::pow(0.5, 4) / static_cast<double>(nk);
        }
    }

    MPI_Reduce(&ret_mpi[0][0][0], &ret[0][0][0], 2 * N * M, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    memory->deallocate(ret_mpi);
    triplet.clear();
}


void ModeAnalysis::calc_frequency_resolved_final_state_tetrahedron(const unsigned int N,
                                                                   double *T,
                                                                   const double omega0,
                                                                   const unsigned int M,
                                                                   const double *omega,
                                                                   const unsigned int ik_in,
                                                                   const unsigned int snum,
                                                                   double ***ret)
{
    int i, j;
    int ik, ib;
    unsigned int jk;

    unsigned int is, js;
    unsigned int k1, k2;
    unsigned int arr[3];

    double omega_inner[2];
    double n1, n2;
    double f1, f2;
    double xk_tmp[3];
    double v3_tmp;
    int nk = kpoint->nk;
    int ns = dynamical->neval;
    int ns2 = ns * ns;

    int *kmap_identity;
    double **energy_tmp;
    double **weight_tetra;
    double **v3_arr;
    double ***delta_arr;
    double prod_tmp[2];

    double epsilon = integration->epsilon;

    std::vector<KsListGroup> triplet;

    kpoint->get_unique_triplet_k(ik_in,
                                 anharmonic_core->use_triplet_symmetry,
                                 false,
                                 triplet);

    for (i = 0; i < N; ++i) {
        for (j = 0; j < M; ++j) {
            ret[i][j][0] = 0.0;
            ret[i][j][1] = 0.0;
        }
    }

    unsigned int npair_uniq = triplet.size();

    memory->allocate(v3_arr, npair_uniq, ns2);
    memory->allocate(delta_arr, npair_uniq, ns2, 2);

    int knum = kpoint->kpoint_irred_all[ik_in][0].knum;
    int knum_minus = kpoint->knum_minus[knum];

    memory->allocate(kmap_identity, nk);

    for (i = 0; i < nk; ++i) kmap_identity[i] = i;


#ifdef _OPENMP
#pragma omp parallel private(is, js, k1, k2, xk_tmp, energy_tmp, i, weight_tetra, ik, jk, arr)
#endif
    {
        memory->allocate(energy_tmp, 3, nk);
        memory->allocate(weight_tetra, 3, nk);

#ifdef _OPENMP
#pragma omp for
#endif
        for (ib = 0; ib < ns2; ++ib) {
            is = ib / ns;
            js = ib % ns;

            for (k1 = 0; k1 < nk; ++k1) {
                // Prepare two-phonon frequency for the tetrahedron method

                for (i = 0; i < 3; ++i) xk_tmp[i] = kpoint->xk[knum][i] - kpoint->xk[k1][i];

                k2 = kpoint->get_knum(xk_tmp[0], xk_tmp[1], xk_tmp[2]);

                energy_tmp[0][k1] = dynamical->eval_phonon[k1][is] + dynamical->eval_phonon[k2][js];
                energy_tmp[1][k1] = dynamical->eval_phonon[k1][is] - dynamical->eval_phonon[k2][js];
                energy_tmp[2][k1] = -energy_tmp[1][k1];
            }

            for (i = 0; i < 3; ++i) {
                integration->calc_weight_tetrahedron(nk,
                                                     kmap_identity,
                                                     weight_tetra[i],
                                                     energy_tmp[i],
                                                     omega0);
            }

            // Loop for irreducible k points
            for (ik = 0; ik < npair_uniq; ++ik) {
                delta_arr[ik][ib][0] = 0.0;
                delta_arr[ik][ib][1] = 0.0;

                for (i = 0; i < triplet[ik].group.size(); ++i) {
                    jk = triplet[ik].group[i].ks[0];
                    delta_arr[ik][ib][0] += weight_tetra[0][jk];
                    delta_arr[ik][ib][1] += weight_tetra[1][jk] - weight_tetra[2][jk];
                }

                // Calculate the matrix element V3 only when the weight is nonzero.
                if (delta_arr[ik][ib][0] > 0.0 || std::abs(delta_arr[ik][ib][1]) > 0.0) {
                    k1 = triplet[ik].group[0].ks[0];
                    k2 = triplet[ik].group[0].ks[1];

                    arr[0] = ns * knum_minus + snum;
                    arr[1] = ns * k1 + is;
                    arr[2] = ns * k2 + js;

                    v3_arr[ik][ib] = std::norm(anharmonic_core->V3(arr));
                } else {
                    v3_arr[ik][ib] = 0.0;
                }
            }
        }

        memory->deallocate(energy_tmp);
        memory->deallocate(weight_tetra);
    }


    for (ik = 0; ik < npair_uniq; ++ik) {
        for (is = 0; is < ns; ++is) {
            for (js = 0; js < ns; ++js) {
                v3_tmp = v3_arr[ik][ns * is + js];

                if (v3_tmp > eps) {
                    k1 = triplet[ik].group[0].ks[0];
                    k2 = triplet[ik].group[0].ks[1];
                    omega_inner[0] = dynamical->eval_phonon[k1][is];
                    omega_inner[1] = dynamical->eval_phonon[k2][js];

#ifdef _OPENMP
#pragma omp parallel for private(f1, f2, n1, n2, prod_tmp, j)
#endif
                    for (i = 0; i < N; ++i) {
                        if (thermodynamics->classical) {
                            f1 = thermodynamics->fC(omega_inner[0], T[i]);
                            f2 = thermodynamics->fC(omega_inner[1], T[i]);

                            n1 = f1 + f2;
                            n2 = f1 - f2;
                        } else {
                            f1 = thermodynamics->fB(omega_inner[0], T[i]);
                            f2 = thermodynamics->fB(omega_inner[1], T[i]);

                            n1 = f1 + f2 + 1.0;
                            n2 = f1 - f2;
                        }

                        prod_tmp[0] = v3_tmp * n1 * delta_arr[ik][ns * is + js][0];
                        prod_tmp[1] = -v3_tmp * n2 * delta_arr[ik][ns * is + js][1];

                        for (j = 0; j < M; ++j) {
                            ret[i][j][0] += prod_tmp[0] * delta_gauss(omega[j] - omega_inner[0], epsilon);
                            ret[i][j][1] += prod_tmp[1] * delta_gauss(omega[j] - omega_inner[0], epsilon);
                        }
                    }
                }
            }
        }
    }

    for (i = 0; i < N; ++i) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (j = 0; j < M; ++j) {
            ret[i][j][0] *= pi * std::pow(0.5, 4);
            ret[i][j][1] *= pi * std::pow(0.5, 4);
        }
    }

    memory->deallocate(v3_arr);
    memory->deallocate(delta_arr);
    memory->deallocate(kmap_identity);

    triplet.clear();
}


void ModeAnalysis::print_momentum_resolved_final_state(const unsigned int NT,
                                                       double *T_arr,
                                                       double epsilon)
{
    int i, j, k, l, m;
    int iT;
    int is, js;
    int nklist;
    int mode;
    double xk1[3], xk2[3], xk3[3];
    double kvec[3];
    double f1, f2, n1, n2;
    double norm, T_tmp;
    double V3norm;
    double **eval, **eval2;
    double **gamma_k;
    int ns = dynamical->neval;

    std::complex<double> ***evec;

    std::ofstream ofs_mode_tau;
    std::string file_mode_tau;

    if (mympi->my_rank == 0) {
        std::cout << std::endl;
        if (integration->ismear == -1) {
            std::cout << " ISMEAR = -1: Tetrahedron method will be used." << std::endl;
            std::cout << " Sorry. Currently, ISMEAR = -1 cannot be used with FSTATE_K = 1.";
            error->exit("calc_momentum_resolved_final_state", "exit.");
        } else if (integration->ismear == 0) {
            std::cout << " ISMEAR = 0: Lorentzian broadening with epsilon = "
                << std::fixed << std::setprecision(2) << epsilon << " (cm^-1)" << std::endl;
        } else if (integration->ismear == 1) {
            std::cout << " ISMEAR = 1: Gaussian broadening with epsilon = "
                << std::fixed << std::setprecision(2) << epsilon << " (cm^-1)" << std::endl;
        } else {
            error->exit("print_momentum_resolved_final_state", "Invalid ISMEAR");
        }

        std::cout << std::endl;
        std::cout << " FSTATE_K = 1 : Calculate the momentum-resolved final state amplitude" << std::endl;
        std::cout << "                due to 3-phonon interactions for given "
            << kslist_fstate_k.size() << " entries." << std::endl;
        std::cout << std::endl;
    }

    double **xk_plane, **xk_plane2;
    double **kvec_plane;
    int nk1_plane, nk2_plane, nk_plane;
    double xk_vec1[3], xk_vec2[3];
    double div1, div2;
    double *eval_tmp;
    double omega_sum[3];
    double frac;
    int knum_triangle[3];
    std::vector<std::vector<double>> ***kplist_conserved;
    std::vector<KpointListWithCoordinate> ***kplist_for_target_mode;
    std::vector<double> xk_vec;
    double xk_norm[3], xk_tmp[3];
    double norm1, norm2, dprod, norm_ref;
    double theta, theta_ref;

    memory->allocate(kplist_conserved, ns, ns, 2);
    memory->allocate(kplist_for_target_mode, ns, ns, kslist_fstate_k.size());

    theta_ref = 0.0;

    // Loop over k point planes

    for (i = 0; i < kpoint->kp_plane_geometry.size(); ++i) {
        nk1_plane = kpoint->kp_plane_geometry[i].npoints[0];
        nk2_plane = kpoint->kp_plane_geometry[i].npoints[1];

        div1 = 1.0 / static_cast<double>(nk1_plane - 1);
        div2 = 1.0 / static_cast<double>(nk2_plane - 1);

        nk_plane = nk1_plane * nk2_plane;

        for (j = 0; j < 3; ++j) {
            xk_vec1[j] = kpoint->kp_plane_geometry[i].xk_edges[0][j]
                - kpoint->kp_plane_geometry[i].xk_origin[j];
            xk_vec2[j] = kpoint->kp_plane_geometry[i].xk_edges[1][j]
                - kpoint->kp_plane_geometry[i].xk_origin[j];
        }


        for (j = 0; j < 3; ++j) {
            xk_norm[j] = xk_vec1[j];
        }

        rotvec(xk_norm, xk_norm, system->rlavec_p, 'T');
        norm_ref = std::sqrt(xk_norm[0] * xk_norm[0]
            + xk_norm[1] * xk_norm[1]
            + xk_norm[2] * xk_norm[2]);

        memory->allocate(xk_plane, nk_plane, 3);
        memory->allocate(xk_plane2, nk_plane, 3);
        memory->allocate(kvec_plane, nk_plane, 3);
        memory->allocate(eval, nk_plane, ns);
        memory->allocate(eval2, nk_plane, ns);
        memory->allocate(eval_tmp, ns);
        memory->allocate(evec, 1, 1, 1);

        // Constructing xk's for the plane
        m = 0;
        for (j = 0; j < nk1_plane; ++j) {
            for (k = 0; k < nk2_plane; ++k) {
                for (l = 0; l < 3; ++l) {
                    xk_plane[m][l] = kpoint->kp_plane_geometry[i].xk_origin[l]
                        + xk_vec1[l] * static_cast<double>(j) * div1
                        + xk_vec2[l] * static_cast<double>(k) * div2;
                }
                ++m;
            }
        }

        // Get frequencies of each k point

        for (j = 0; j < nk_plane; ++j) {
            for (k = 0; k < 3; ++k) kvec_plane[j][k] = dynamical->fold(xk_plane[j][k]);
            rotvec(kvec_plane[j], kvec_plane[j], system->rlavec_p, 'T');
            norm = std::sqrt(kvec_plane[j][0] * kvec_plane[j][0]
                + kvec_plane[j][1] * kvec_plane[j][1]
                + kvec_plane[j][2] * kvec_plane[j][2]);

            if (norm > eps) {
                for (k = 0; k < 3; ++k) kvec_plane[j][k] /= norm;
            }
        }

        for (j = 0; j < nk_plane; ++j) {
            dynamical->eval_k(xk_plane[j],
                              kvec_plane[j],
                              fcs_phonon->fc2_ext,
                              eval[j],
                              evec[0],
                              false);
            for (k = 0; k < ns; ++k) {
                eval[j][k] = dynamical->freq(eval[j][k]);
            }
        }

        // Loop over k points to analyze the final state amplitude

        for (j = 0; j < kslist_fstate_k.size(); ++j) {
            for (k = 0; k < 3; ++k) xk1[k] = kslist_fstate_k[j].xk[k];
            mode = kslist_fstate_k[j].nmode;

            for (k = 0; k < 3; ++k) kvec[k] = kslist_fstate_k[j].xk[k];
            rotvec(kvec, kvec, system->rlavec_p, 'T');
            norm = std::sqrt(kvec[0] * kvec[0]
                + kvec[1] * kvec[1]
                + kvec[2] * kvec[2]);

            if (norm > eps) {
                for (k = 0; k < 3; ++k) kvec[k] /= norm;
            }

            // for (k = 0; k < 3; ++k) xk1[k] = dynamical->fold(xk1[k]);

            dynamical->eval_k(xk1,
                              kvec,
                              fcs_phonon->fc2_ext,
                              eval_tmp,
                              evec[0],
                              false);
            for (k = 0; k < ns; ++k) eval_tmp[k] = dynamical->freq(eval_tmp[k]);

            // Calculate xk's for the third index that satisfy the momentum conservation

            for (k = 0; k < nk_plane; ++k) {
                for (l = 0; l < 3; ++l) {
                    xk_plane2[k][l] = dynamical->fold(xk1[l] - xk_plane[k][l]);
                }
            }

            // Get frequencies of each k point

            for (k = 0; k < nk_plane; ++k) {
                for (l = 0; l < 3; ++l) kvec_plane[k][l] = xk_plane2[k][l];
                rotvec(kvec_plane[k], kvec_plane[k], system->rlavec_p, 'T');
                norm = std::sqrt(kvec_plane[k][0] * kvec_plane[k][0]
                    + kvec_plane[k][1] * kvec_plane[k][1]
                    + kvec_plane[k][2] * kvec_plane[k][2]);

                if (norm > eps) {
                    for (l = 0; l < 3; ++l) kvec_plane[k][l] /= norm;
                }
            }

            for (k = 0; k < nk_plane; ++k) {
                dynamical->eval_k(xk_plane2[k],
                                  kvec_plane[k],
                                  fcs_phonon->fc2_ext,
                                  eval2[k],
                                  evec[0],
                                  false);
                for (l = 0; l < ns; ++l) {
                    eval2[k][l] = dynamical->freq(eval2[k][l]);
                }
            }

            // Find a list of k points which satisfy the energy conservation

            for (const auto &it : kpoint->kp_planes_tri[i]) {
                // K point indexes for each triangle
                for (k = 0; k < 3; ++k) knum_triangle[k] = it.knum[k];

                //for (k=0; k<3; ++k) std::cout << " " << knum_triangle[k];
                //std::cout << std::endl;


                //for (k = 0; k < 3; ++k) {
                //    std::cout << "xk1 = " << xk1[0] << " " << xk1[1] << " " << xk1[2] << std::endl;
                //    std::cout << "xk2 = " << xk_plane[knum_triangle[k]][0] << " " << xk_plane[knum_triangle[k]][1] << " " << xk_plane[knum_triangle[k]][2] << std::endl;
                //    std::cout << "xk3 = " << xk_plane2[knum_triangle[k]][0] << " " << xk_plane2[knum_triangle[k]][1] << " " << xk_plane2[knum_triangle[k]][2] << std::endl;
                //    std::cout << std::endl;
                //}

                for (is = 0; is < ns; ++is) {
                    for (js = 0; js < ns; ++js) {
                        // The case of delta(w1 - w2 - w3) 

                        //     std::cout << "is = " << is << " js = " << js << std::endl;
                        for (k = 0; k < 3; ++k) {
                            omega_sum[k] = eval_tmp[mode]
                                - eval[knum_triangle[k]][is]
                                - eval2[knum_triangle[k]][js];
                            //std::cout << " " << writes->in_kayser(eval_tmp[mode]) 
                            //<< " " << writes->in_kayser(eval[knum_triangle[k]][is])
                            //<< " " << writes->in_kayser(eval2[knum_triangle[k]][js]) << std::endl;
                        }
                        if ((omega_sum[0] > 0.0 && omega_sum[1] > 0.0 && omega_sum[2] > 0.0) ||
                            (omega_sum[0] < 0.0 && omega_sum[1] < 0.0 && omega_sum[2] < 0.0))
                            continue;

                        if (omega_sum[0] * omega_sum[1] < 0.0) {
                            xk_vec.clear();

                            frac = -omega_sum[0] / (omega_sum[1] - omega_sum[0]);

                            for (k = 0; k < 3; ++k) {
                                xk_vec.push_back((1.0 - frac) * xk_plane[knum_triangle[0]][k]
                                    + frac * xk_plane[knum_triangle[1]][k]);
                            }
                            kplist_conserved[is][js][0].push_back(xk_vec);
                        }

                        if (omega_sum[0] * omega_sum[2] < 0.0) {
                            xk_vec.clear();

                            frac = -omega_sum[0] / (omega_sum[2] - omega_sum[0]);

                            for (k = 0; k < 3; ++k) {
                                xk_vec.push_back((1.0 - frac) * xk_plane[knum_triangle[0]][k]
                                    + frac * xk_plane[knum_triangle[2]][k]);
                            }
                            kplist_conserved[is][js][0].push_back(xk_vec);
                        }

                        if (omega_sum[1] * omega_sum[2] < 0.0) {
                            xk_vec.clear();

                            frac = -omega_sum[1] / (omega_sum[2] - omega_sum[1]);

                            for (k = 0; k < 3; ++k) {
                                xk_vec.push_back((1.0 - frac) * xk_plane[knum_triangle[1]][k]
                                    + frac * xk_plane[knum_triangle[2]][k]);
                            }
                            kplist_conserved[is][js][0].push_back(xk_vec);
                        }

                        // The case of delta(w1 - w2 + w3)

                        for (k = 0; k < 3; ++k) {
                            omega_sum[k] = eval_tmp[mode]
                                - eval[knum_triangle[k]][is]
                                + eval2[knum_triangle[k]][js];
                        }
                        if ((omega_sum[0] > 0.0 && omega_sum[1] > 0.0 && omega_sum[2] > 0.0) ||
                            (omega_sum[0] < 0.0 && omega_sum[1] < 0.0 && omega_sum[2] < 0.0))
                            continue;

                        if (omega_sum[0] * omega_sum[1] < 0.0) {
                            xk_vec.clear();

                            frac = -omega_sum[0] / (omega_sum[1] - omega_sum[0]);

                            for (k = 0; k < 3; ++k) {
                                xk_vec.push_back((1.0 - frac) * xk_plane[knum_triangle[0]][k]
                                    + frac * xk_plane[knum_triangle[1]][k]);
                            }
                            kplist_conserved[is][js][1].push_back(xk_vec);
                        }

                        if (omega_sum[0] * omega_sum[2] < 0.0) {
                            xk_vec.clear();

                            frac = -omega_sum[0] / (omega_sum[2] - omega_sum[0]);

                            for (k = 0; k < 3; ++k) {
                                xk_vec.push_back((1.0 - frac) * xk_plane[knum_triangle[0]][k]
                                    + frac * xk_plane[knum_triangle[2]][k]);
                            }
                            kplist_conserved[is][js][1].push_back(xk_vec);
                        }

                        if (omega_sum[1] * omega_sum[2] < 0.0) {
                            xk_vec.clear();

                            frac = -omega_sum[1] / (omega_sum[2] - omega_sum[1]);

                            for (k = 0; k < 3; ++k) {
                                xk_vec.push_back((1.0 - frac) * xk_plane[knum_triangle[1]][k]
                                    + frac * xk_plane[knum_triangle[2]][k]);
                            }
                            kplist_conserved[is][js][1].push_back(xk_vec);
                        }
                    }
                }
            }

            for (is = 0; is < ns; ++is) {
                for (js = 0; js < ns; ++js) {
                    for (auto it2 = kplist_conserved[is][js][0].begin();
                         it2 != kplist_conserved[is][js][0].end(); ++it2) {
                        for (k = 0; k < 3; ++k) {
                            xk_tmp[k] = (*it2)[k];
                        }
                        rotvec(xk_tmp, xk_tmp, system->rlavec_p, 'T');

                        norm1 = 0.0;
                        dprod = 0.0;
                        for (k = 0; k < 3; ++k) {
                            norm1 += xk_tmp[k] * xk_tmp[k];
                            dprod += xk_tmp[k] * xk_norm[k];
                        }
                        theta = std::acos(dprod / (norm_ref * std::sqrt(norm1)));

                        kplist_for_target_mode[is][js][j].emplace_back(*it2,
                                                                       std::cos(theta + theta_ref) * std::sqrt(norm1),
                                                                       std::sin(theta + theta_ref) * std::sqrt(norm1),
                                                                       i, 0);
                    }

                    for (auto it2 = kplist_conserved[is][js][1].begin();
                         it2 != kplist_conserved[is][js][1].end(); ++it2) {
                        for (k = 0; k < 3; ++k) {
                            xk_tmp[k] = (*it2)[k];
                        }
                        rotvec(xk_tmp, xk_tmp, system->rlavec_p, 'T');

                        norm1 = 0.0;
                        dprod = 0.0;
                        for (k = 0; k < 3; ++k) {
                            norm1 += xk_tmp[k] * xk_tmp[k];
                            dprod += xk_tmp[k] * xk_norm[k];
                        }
                        theta = std::acos(dprod / (norm_ref * std::sqrt(norm1)));

                        kplist_for_target_mode[is][js][j].emplace_back(*it2,
                                                                       std::cos(theta + theta_ref) * std::sqrt(norm1),
                                                                       std::sin(theta + theta_ref) * std::sqrt(norm1),
                                                                       i, 1);
                    }

                    kplist_conserved[is][js][0].clear();
                    kplist_conserved[is][js][1].clear();
                }
            }
        }

        memory->deallocate(xk_plane);
        memory->deallocate(xk_plane2);
        memory->deallocate(kvec_plane);
        memory->deallocate(eval);
        memory->deallocate(eval2);
        memory->deallocate(eval_tmp);
        memory->deallocate(evec);

        rotvec(xk_vec1, xk_vec1, system->rlavec_p, 'T');
        rotvec(xk_vec2, xk_vec2, system->rlavec_p, 'T');

        norm1 = 0.0;
        norm2 = 0.0;
        dprod = 0.0;
        for (j = 0; j < 3; ++j) {
            norm1 += xk_vec1[j] * xk_vec1[j];
            norm2 += xk_vec2[j] * xk_vec2[j];
            dprod += xk_vec1[j] * xk_vec2[j];
        }
        theta = std::acos(dprod / std::sqrt(norm1 * norm2));

        theta_ref += theta;
    }

    memory->deallocate(kplist_conserved);


    std::vector<std::vector<double>> **final_state_xy;
    std::vector<double> triplet_xyG;
    std::vector<int> small_group_k;
    double pos_x, pos_y;
    int selection_type;

    int isym;


    double srot[3][3];
    double xk_sym[3];
    double srot_inv[3][3], srot_inv_t[3][3];
    double ***symop_k;
    double diff;

    memory->allocate(symop_k, symmetry->nsym, 3, 3);
    memory->allocate(final_state_xy, kslist_fstate_k.size(), NT);

    memory->allocate(eval, 3, ns);
    memory->allocate(evec, 3, ns, ns);

    for (i = 0; i < kslist_fstate_k.size(); ++i) {
        for (j = 0; j < 3; ++j) xk1[j] = -kslist_fstate_k[i].xk[j];
        mode = kslist_fstate_k[i].nmode;

        for (j = 0; j < 3; ++j) kvec[j] = dynamical->fold(xk1[j]);
        rotvec(kvec, kvec, system->rlavec_p, 'T');
        norm = std::sqrt(kvec[0] * kvec[0] + kvec[1] * kvec[1] + kvec[2] * kvec[2]);

        if (norm > eps) for (j = 0; j < 3; ++j) kvec[j] /= norm;
        for (j = 0; j < 3; ++j) xk1[j] = dynamical->fold(xk1[j]);

        dynamical->eval_k(xk1, kvec, fcs_phonon->fc2_ext, eval[0], evec[0], true);
        for (j = 0; j < ns; ++j) eval[0][j] = dynamical->freq(eval[0][j]);

        small_group_k.clear();

        for (isym = 0; isym < symmetry->nsym; ++isym) {
            for (j = 0; j < 3; ++j) {
                for (k = 0; k < 3; ++k) {
                    srot[j][k] = static_cast<double>(symmetry->SymmList[isym].rot[j][k]);
                }
            }

            invmat3(srot_inv, srot);
            transpose3(srot_inv_t, srot_inv);

            for (j = 0; j < 3; ++j) {
                for (k = 0; k < 3; ++k) {
                    symop_k[isym][j][k] = srot_inv_t[j][k];
                }
            }

            rotvec(xk_sym, xk1, symop_k[isym]);

            for (j = 0; j < 3; ++j) {
                xk_sym[j] = xk_sym[j] - nint(xk_sym[j]);
                xk_sym[j] = std::fmod(xk_sym[j], 1.0);
            }

            diff = 0.0;
            for (j = 0; j < 3; ++j) {
                diff += std::pow(std::fmod(xk_sym[j] - xk1[j], 1.0), 2);
            }
            if (std::sqrt(diff) < eps8) {
                small_group_k.push_back(isym);
            }
        }

        if (mympi->my_rank == 0) {
            std::cout << " Number : " << std::setw(5) << i + 1 << std::endl;
            std::cout << "  Phonon at k = (";
            for (j = 0; j < 3; ++j) {
                std::cout << std::setw(10) << std::fixed << kslist_fstate_k[i].xk[j];
                if (j < 2) std::cout << ",";
            }
            std::cout << ")" << std::endl;
            std::cout << "  Mode index = " << std::setw(5) << mode + 1 << std::endl;
            std::cout << "  Frequency (cm^-1) : " << std::setw(15)
                << writes->in_kayser(eval[0][mode]) << std::endl;

            int count_kp = 0;

            for (is = 0; is < ns; ++is) {
                for (js = 0; js < ns; ++js) {
                    count_kp += kplist_for_target_mode[is][js][i].size();
                }
            }
            std::cout << "  Number of k points satisfying the selection rule : "
                << count_kp << std::endl;
            std::cout << "  Number of symmetry operations at k point : "
                << small_group_k.size() << std::endl << std::endl;
        }

        for (is = 0; is < ns; ++is) {
            for (js = 0; js < ns; ++js) {
                nklist = kplist_for_target_mode[is][js][i].size();

                if (nklist == 0) continue;

                memory->allocate(gamma_k, nklist, NT);

                for (k = 0; k < nklist; ++k) {
                    for (l = 0; l < NT; ++l) {
                        gamma_k[k][l] = 0.0;
                    }
                }

                for (k = 0; k < nklist; ++k) {
                    for (l = 0; l < 3; ++l)
                        xk2[l] = dynamical->fold(kplist_for_target_mode[is][js][i][k].xk[l]);

                    for (isym = 0; isym < small_group_k.size(); ++isym) {
                        rotvec(xk_sym, xk2, symop_k[small_group_k[isym]]);

                        for (l = 0; l < 3; ++l) xk3[l] = dynamical->fold(-xk1[l] - xk_sym[l]);

                        for (l = 0; l < 3; ++l) kvec[l] = xk_sym[l];
                        rotvec(kvec, kvec, system->rlavec_p, 'T');
                        norm = std::sqrt(kvec[0] * kvec[0] + kvec[1] * kvec[1] + kvec[2] * kvec[2]);

                        if (norm > eps) for (l = 0; l < 3; ++l) kvec[l] /= norm;

                        dynamical->eval_k(xk_sym, kvec, fcs_phonon->fc2_ext, eval[1], evec[1], true);

                        for (l = 0; l < 3; ++l) kvec[l] = xk3[l];
                        rotvec(kvec, kvec, system->rlavec_p, 'T');
                        norm = std::sqrt(kvec[0] * kvec[0] + kvec[1] * kvec[1] + kvec[2] * kvec[2]);

                        if (norm > eps) for (l = 0; l < 3; ++l) kvec[l] /= norm;

                        dynamical->eval_k(xk3, kvec, fcs_phonon->fc2_ext, eval[2], evec[2], true);

                        for (l = 0; l < ns; ++l) {
                            eval[1][l] = dynamical->freq(eval[1][l]);
                            eval[2][l] = dynamical->freq(eval[2][l]);
                        }

                        V3norm = std::norm(anharmonic_core->V3_mode(mode, xk_sym, xk3, is, js, eval, evec));

                        for (iT = 0; iT < NT; ++iT) {
                            T_tmp = T_arr[iT];

                            if (thermodynamics->classical) {
                                f1 = thermodynamics->fC(eval[1][is], T_tmp);
                                f2 = thermodynamics->fC(eval[2][js], T_tmp);
                                n1 = f1 + f2;
                                n2 = f1 - f2;
                            } else {
                                f1 = thermodynamics->fB(eval[1][is], T_tmp);
                                f2 = thermodynamics->fB(eval[2][js], T_tmp);
                                n1 = f1 + f2 + 1.0;
                                n2 = f1 - f2;
                            }


                            if (selection_type == 0) {
                                gamma_k[k][iT] += V3norm * n1;
                            } else if (selection_type == 1) {
                                gamma_k[k][iT] += V3norm * n2;
                            }
                        }
                    }

                    for (iT = 0; iT < NT; ++iT)
                        gamma_k[k][iT] *= pi * std::pow(0.5, 4) / static_cast<double>(small_group_k.size());

                    pos_x = kplist_for_target_mode[is][js][i][k].x;
                    pos_y = kplist_for_target_mode[is][js][i][k].y;
                    selection_type = kplist_for_target_mode[is][js][i][k].selection_type;

                    for (iT = 0; iT < NT; ++iT) {
                        triplet_xyG.clear();
                        triplet_xyG.push_back(pos_x);
                        triplet_xyG.push_back(pos_y);
                        triplet_xyG.push_back(gamma_k[k][iT]);
                        final_state_xy[i][iT].push_back(triplet_xyG);
                    }
                }
                memory->deallocate(gamma_k);
            }
        }

        if (mympi->my_rank == 0) {
            file_mode_tau = input->job_title + ".fk." + std::to_string(i + 1);
            ofs_mode_tau.open(file_mode_tau.c_str(), std::ios::out);
            if (!ofs_mode_tau)
                error->exit("compute_mode_tau",
                            "Cannot open file file_mode_tau");

            ofs_mode_tau << "## Momentum-resolved final state amplitude" << std::endl;

            ofs_mode_tau << "# " << "Gamma at ";
            for (l = 0; l < 3; ++l) ofs_mode_tau << std::setw(10) << -xk1[l];
            ofs_mode_tau << " , mode = " << mode + 1 << std::endl;
            ofs_mode_tau << " # Temperature [K], coordinate in FBZ, final state amplitude" << std::endl;
            for (iT = 0; iT < NT; ++iT) {
                for (k = 0; k < final_state_xy[i][iT].size(); ++k) {
                    ofs_mode_tau << std::setw(10) << T_arr[iT];
                    ofs_mode_tau << std::setw(15) << final_state_xy[i][iT][k][0];
                    ofs_mode_tau << std::setw(15) << final_state_xy[i][iT][k][1];
                    ofs_mode_tau << std::setw(15) << final_state_xy[i][iT][k][2] << std::endl;
                }
            }

            ofs_mode_tau.close();
            std::cout << "  The result is saved in " << file_mode_tau << std::endl;
            std::cout << std::endl;
        }
    }

    memory->deallocate(kplist_for_target_mode);
    memory->deallocate(final_state_xy);
    memory->deallocate(symop_k);
    memory->deallocate(eval);
    memory->deallocate(evec);
}

void ModeAnalysis::print_V3_elements()
{
    int i, j;
    int knum;
    int snum;
    int ns = dynamical->neval;
    double omega;
    double **v3norm;
    std::string file_V3;
    std::ofstream ofs_V3;

    int ik_irred, multi;
    unsigned int nk_size;
    unsigned int ib, is, js, k1, k2;
    std::vector<KsListGroup> triplet;

    for (i = 0; i < kslist.size(); ++i) {
        knum = kslist[i] / ns;
        snum = kslist[i] % ns;

        omega = dynamical->eval_phonon[knum][snum];

        if (mympi->my_rank == 0) {
            std::cout << std::endl;
            std::cout << " Number : " << std::setw(5) << i + 1 << std::endl;
            std::cout << "  Phonon at k = (";
            for (j = 0; j < 3; ++j) {
                std::cout << std::setw(10) << std::fixed << kpoint->xk[knum][j];
                if (j < 2) std::cout << ",";
            }
            std::cout << ")" << std::endl;
            std::cout << "  Mode index = " << std::setw(5) << snum + 1 << std::endl;
            std::cout << "  Frequency (cm^-1) : "
                << std::setw(15) << writes->in_kayser(omega) << std::endl;
        }

        ik_irred = kpoint->kmap_to_irreducible[knum];

        kpoint->get_unique_triplet_k(ik_irred,
                                     anharmonic_core->use_triplet_symmetry,
                                     true,
                                     triplet);
        nk_size = triplet.size();

        memory->allocate(v3norm, nk_size, ns * ns);

        calc_V3norm2(ik_irred, snum, v3norm);

        if (mympi->my_rank == 0) {
            file_V3 = input->job_title + ".V3." + std::to_string(i + 1);
            ofs_V3.open(file_V3.c_str(), std::ios::out);
            if (!ofs_V3)
                error->exit("run_mode_analysis",
                            "Cannot open file file_V3");

            ofs_V3 << "# xk = ";

            for (j = 0; j < 3; ++j) {
                ofs_V3 << std::setw(15) << kpoint->xk[knum][j];
            }
            ofs_V3 << std::endl;
            ofs_V3 << "# mode = " << snum + 1 << std::endl;
            ofs_V3 << "# Frequency = " << writes->in_kayser(omega) << std::endl;
            ofs_V3 << "## Matrix elements |V3|^2 for given mode" << std::endl;
            ofs_V3 << "## q', j', omega(q'j') (cm^-1), q'', j'', ";
            ofs_V3 << "omega(q''j'') (cm^-1), |V3(-qj,q'j',q''j'')|^2 (cm^-2), multiplicity" << std::endl;

            for (j = 0; j < nk_size; ++j) {
                multi = static_cast<double>(triplet[j].group.size());
                k1 = triplet[j].group[0].ks[0];
                k2 = triplet[j].group[0].ks[1];

                ib = 0;

                for (is = 0; is < ns; ++is) {
                    for (js = 0; js < ns; ++js) {
                        ofs_V3 << std::setw(5) << k1 + 1 << std::setw(5) << is + 1;
                        ofs_V3 << std::setw(15)
                            << writes->in_kayser(dynamical->eval_phonon[k1][is]);
                        ofs_V3 << std::setw(5) << k2 + 1 << std::setw(5) << js + 1;
                        ofs_V3 << std::setw(15)
                            << writes->in_kayser(dynamical->eval_phonon[k2][js]);
                        ofs_V3 << std::setw(15) << v3norm[j][ib];
                        ofs_V3 << std::setw(5) << multi;
                        ofs_V3 << std::endl;

                        ++ib;
                    }
                    ofs_V3 << std::endl;
                }
            }

            ofs_V3.close();
        }
        memory->deallocate(v3norm);
    }
}


void ModeAnalysis::calc_V3norm2(const unsigned int ik_in,
                                const unsigned int snum,
                                double **ret)
{
    int ib;
    unsigned int ik;
    unsigned int is, js;
    unsigned int k1, k2;
    unsigned int arr[3];
    unsigned int knum, knum_minus;
    int ns = dynamical->neval;

    int ns2 = ns * ns;

    double factor = std::pow(0.5, 3) * std::pow(Hz_to_kayser / time_ry, 2);
    std::vector<KsListGroup> triplet;

    knum = kpoint->kpoint_irred_all[ik_in][0].knum;
    knum_minus = kpoint->knum_minus[knum];

    kpoint->get_unique_triplet_k(ik_in,
                                 anharmonic_core->use_triplet_symmetry,
                                 true,
                                 triplet);
#ifdef _OPENMP
#pragma omp parallel for private(is, js, ik, k1, k2, arr)
#endif
    for (ib = 0; ib < ns2; ++ib) {
        is = ib / ns;
        js = ib % ns;

        for (ik = 0; ik < triplet.size(); ++ik) {
            k1 = triplet[ik].group[0].ks[0];
            k2 = triplet[ik].group[0].ks[1];

            arr[0] = ns * knum_minus + snum;
            arr[1] = ns * k1 + is;
            arr[2] = ns * k2 + js;

            ret[ik][ib] = std::norm(anharmonic_core->V3(arr)) * factor;
        }
    }
}


void ModeAnalysis::print_Phi3_elements()
{
    int i, j;
    int knum;
    int snum;
    int ns = dynamical->neval;
    double omega;
    std::complex<double> **phi3;
    std::string file_V3;
    std::ofstream ofs_V3;

    int ik_irred, multi;
    unsigned int nk_size;
    unsigned int ib, is, js, k1, k2;
    std::vector<KsListGroup> triplet;

    for (i = 0; i < kslist.size(); ++i) {
        knum = kslist[i] / ns;
        snum = kslist[i] % ns;

        omega = dynamical->eval_phonon[knum][snum];

        if (mympi->my_rank == 0) {
            std::cout << std::endl;
            std::cout << " Number : " << std::setw(5) << i + 1 << std::endl;
            std::cout << "  Phonon at k = (";
            for (j = 0; j < 3; ++j) {
                std::cout << std::setw(10) << std::fixed << kpoint->xk[knum][j];
                if (j < 2) std::cout << ",";
            }
            std::cout << ")" << std::endl;
            std::cout << "  Mode index = " << std::setw(5) << snum + 1 << std::endl;
            std::cout << "  Frequency (cm^-1) : "
                << std::setw(15) << writes->in_kayser(omega) << std::endl;
        }

        ik_irred = kpoint->kmap_to_irreducible[knum];

        kpoint->get_unique_triplet_k(ik_irred,
                                     anharmonic_core->use_triplet_symmetry,
                                     true,
                                     triplet, 1);
        nk_size = triplet.size();

        memory->allocate(phi3, nk_size, ns * ns);

        calc_Phi3(knum, snum, triplet, phi3);

        if (mympi->my_rank == 0) {
            file_V3 = input->job_title + ".Phi3." + std::to_string(i + 1);
            ofs_V3.open(file_V3.c_str(), std::ios::out);
            if (!ofs_V3)
                error->exit("print_phi3_element",
                            "Cannot open file file_V3");

            ofs_V3 << "# xk = ";

            for (j = 0; j < 3; ++j) {
                ofs_V3 << std::setw(15) << kpoint->xk[knum][j];
            }
            ofs_V3 << std::endl;
            ofs_V3 << "# mode = " << snum + 1 << std::endl;
            ofs_V3 << "# Frequency = " << writes->in_kayser(omega) << std::endl;
            ofs_V3 << "## Matrix elements Phi3 for given mode" << std::endl;
            ofs_V3 << "## q', j', omega(q'j') (cm^-1), q'', j'', omega(q''j'') (cm^-1), ";
            ofs_V3 << "Phi3(qj,q'j',q''j'') (Ry/(u^{1/2}Bohr)^{3}), multiplicity" << std::endl;

            for (j = 0; j < nk_size; ++j) {
                multi = static_cast<double>(triplet[j].group.size());
                k1 = triplet[j].group[0].ks[0];
                k2 = triplet[j].group[0].ks[1];

                ib = 0;

                for (is = 0; is < ns; ++is) {
                    for (js = 0; js < ns; ++js) {
                        ofs_V3 << std::setw(5) << k1 + 1 << std::setw(5) << is + 1;
                        ofs_V3 << std::setw(15)
                            << writes->in_kayser(dynamical->eval_phonon[k1][is]);
                        ofs_V3 << std::setw(5) << k2 + 1 << std::setw(5) << js + 1;
                        ofs_V3 << std::setw(15)
                            << writes->in_kayser(dynamical->eval_phonon[k2][js]);
                        ofs_V3 << std::setw(15) << phi3[j][ib].real();
                        ofs_V3 << std::setw(15) << phi3[j][ib].imag();
                        ofs_V3 << std::setw(5) << multi;
                        ofs_V3 << std::endl;

                        ++ib;
                    }
                    ofs_V3 << std::endl;
                }
            }

            ofs_V3.close();
        }
        memory->deallocate(phi3);
    }
}


void ModeAnalysis::calc_Phi3(const unsigned int knum,
                             const unsigned int snum,
                             const std::vector<KsListGroup> &triplet,
                             std::complex<double> **ret)
{
    int ib;
    unsigned int ik;
    unsigned int is, js;
    unsigned int k1, k2;
    unsigned int arr[3];
    int ns = dynamical->neval;
    int ns2 = ns * ns;
    double omega[3];

    double factor = std::pow(amu_ry, 1.5);

#ifdef _OPENMP
#pragma omp parallel for private(is, js, ik, k1, k2, arr, omega)
#endif
    for (ib = 0; ib < ns2; ++ib) {
        is = ib / ns;
        js = ib % ns;

        for (ik = 0; ik < triplet.size(); ++ik) {
            k1 = triplet[ik].group[0].ks[0];
            k2 = triplet[ik].group[0].ks[1];

            arr[0] = ns * knum + snum;
            arr[1] = ns * k1 + is;
            arr[2] = ns * k2 + js;

            omega[0] = dynamical->eval_phonon[knum][snum];
            omega[1] = dynamical->eval_phonon[k1][is];
            omega[2] = dynamical->eval_phonon[k2][js];

            ret[ik][ib] = anharmonic_core->V3(arr) * factor
                * std::sqrt(omega[0] * omega[1] * omega[2]);
        }
    }
}


void ModeAnalysis::print_spectral_function(const int NT,
                                           double *T_arr)
{
    int ns = dynamical->neval;
    int i, j;
    int knum, snum;
    int ik_irred, iomega, iT;
    double **self3_imag, **self3_real;
    std::string file_self;
    std::ofstream ofs_self;
    double *omega_array;
    double Omega_min = dos->emin;
    double Omega_max = dos->emax;
    double delta_omega = dos->delta_e;
    double T_now, omega2[2];
    double omega;

    int nomega = static_cast<unsigned int>((Omega_max - Omega_min) / delta_omega) + 1;

    memory->allocate(omega_array, nomega);
    memory->allocate(self3_imag, NT, nomega);
    memory->allocate(self3_real, NT, nomega);

    for (i = 0; i < nomega; ++i) {
        omega_array[i] = Omega_min + delta_omega * static_cast<double>(i);
        omega_array[i] *= time_ry / Hz_to_kayser;
    }

    for (i = 0; i < kslist.size(); ++i) {
        knum = kslist[i] / ns;
        snum = kslist[i] % ns;
        ik_irred = kpoint->kmap_to_irreducible[knum];

        if (mympi->my_rank == 0) {
            std::cout << std::endl;
            std::cout << " SELF_W = 1: Calculate bubble selfenergy with frequency dependency" << std::endl;
            std::cout << " for given " << kslist.size() << " modes." << std::endl;
            std::cout << std::endl;
            std::cout << " Number : " << std::setw(5) << i + 1 << std::endl;
            std::cout << "  Phonon at k = (";
            for (j = 0; j < 3; ++j) {
                std::cout << std::setw(10) << std::fixed << kpoint->xk[knum][j];
                if (j < 2) std::cout << ",";
            }
            std::cout << ")" << std::endl;
            std::cout << "  Mode index = " << std::setw(5) << snum + 1 << std::endl;

            file_self = input->job_title + ".Self." + std::to_string(i + 1);
            ofs_self.open(file_self.c_str(), std::ios::out);
            if (!ofs_self) error->exit("run_mode_analysis", "Cannot open file file_shift");

            ofs_self << "# xk = ";

            for (j = 0; j < 3; ++j) {
                ofs_self << std::setw(15) << kpoint->xk[knum][j];
            }
            ofs_self << std::endl;
            ofs_self << "# mode = " << snum + 1 << std::endl;
            ofs_self << "## T[K], Freq (cm^-1), omega (cm^-1), Self.real (cm^-1), Self.imag (cm^-1)";
            ofs_self << std::endl;
        }

        for (iT = 0; iT < NT; ++iT) {
            T_now = T_arr[iT];
            omega = dynamical->eval_phonon[knum][snum];

            if (mympi->my_rank == 0) {
                std::cout << "  Temperature (K) : " << std::setw(15) << T_now << std::endl;
                std::cout << "  Frequency (cm^-1) : " << std::setw(15) << writes->in_kayser(omega) << std::endl;
            }

            anharmonic_core->calc_self3omega_tetrahedron(T_now,
                                                         dynamical->eval_phonon,
                                                         dynamical->evec_phonon,
                                                         ik_irred,
                                                         snum,
                                                         nomega,
                                                         omega_array,
                                                         self3_imag[iT]);

            // Calculate real part of the self-energy by Kramers-Kronig relation
            for (iomega = 0; iomega < nomega; ++iomega) {
                double self_tmp = 0.0;
                omega2[0] = omega_array[iomega] * omega_array[iomega];
                for (int jomega = 0; jomega < nomega; ++jomega) {
                    if (jomega == iomega) continue;
                    omega2[1] = omega_array[jomega] * omega_array[jomega];
                    self_tmp += omega_array[jomega] * self3_imag[iT][jomega] / (omega2[1] - omega2[0]);
                }
                self3_real[iT][iomega] = 2.0 * delta_omega * time_ry * self_tmp / (pi * Hz_to_kayser);
            }

            if (mympi->my_rank == 0) {
                for (iomega = 0; iomega < nomega; ++iomega) {
                    ofs_self << std::setw(10) << T_now << std::setw(15) << writes->in_kayser(omega);
                    ofs_self << std::setw(10) << writes->in_kayser(omega_array[iomega])
                        << std::setw(15) << writes->in_kayser(self3_real[iT][iomega])
                        << std::setw(15) << writes->in_kayser(self3_imag[iT][iomega]) << std::endl;
                }
                ofs_self << std::endl;
            }
        }
        if (mympi->my_rank == 0) ofs_self.close();
    }

    memory->deallocate(omega_array);
    memory->deallocate(self3_imag);
    memory->deallocate(self3_real);
}
