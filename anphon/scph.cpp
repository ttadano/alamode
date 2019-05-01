/*
 scph.cpp

 Copyright (c) 2015 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#include "mpi_common.h"
#include "scph.h"
#include "dynamical.h"
#include "kpoint.h"
#include "anharmonic_core.h"
#include "ewald.h"
#include "memory.h"
#include "thermodynamics.h"
#include "write_phonons.h"
#include "constants.h"
#include "system.h"
#include "error.h"
#include "mathfunctions.h"
#include "parsephon.h"
#include "phonon_dos.h"
#include "symmetry_core.h"
#include <iostream>
#include <iomanip>
#include <complex>
#include <algorithm>
#include <fftw3.h>
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <boost/lexical_cast.hpp>
#include "timer.h"
#include <cmath>
#include <cstdlib>
#include <vector>

#if defined(WIN32) || defined(_WIN32)
#pragma comment(lib, "libfftw3-3.lib")
#pragma comment(lib, "libfftw3f-3.lib")
#pragma comment(lib, "libfftw3l-3.lib")
#endif

using namespace PHON_NS;

Scph::Scph(PHON *phon): Pointers(phon)
{
    set_default_variables();
}

Scph::~Scph()
{
    deallocate_variables();
}


void Scph::set_default_variables()
{
    im = std::complex<double>(0.0, 1.0);
    restart_scph = false;
    warmstart_scph = false;
    lower_temp = true;
    tolerance_scph = 1.0e-10;
    mixalpha = 0.1;
    maxiter = 100;
    print_self_consistent_fc2 = false;
    selfenergy_offdiagonal = true;
    relax_coordinate = false;

    xk_scph = nullptr;
    kvec_na_scph = nullptr;
    xk_interpolate = nullptr;
    vec_for_v3 = nullptr;
    vec_for_v4 = nullptr;
    invmass_for_v3 = nullptr;
    invmass_for_v4 = nullptr;
    evec_index3 = nullptr;
    evec_index4 = nullptr;
    kmap_interpolate_to_scph = nullptr;
    evec_harmonic = nullptr;
    fcs_group = nullptr;
    fcs_group2 = nullptr;
    omega2_harmonic = nullptr;
    mat_transform_sym = nullptr;
    small_group_at_k = nullptr;
    symop_minus_at_k = nullptr;
    kpoint_map_symmetry = nullptr;
    exp_phase = nullptr;
    exp_phase3 = nullptr;
    mindist_list_scph = nullptr;
}


void Scph::deallocate_variables()
{
    if (xk_scph) {
        memory->deallocate(xk_scph);
    }
    if (kvec_na_scph) {
        memory->deallocate(kvec_na_scph);
    }
    if (xk_interpolate) {
        memory->deallocate(xk_interpolate);
    }
    if (kmap_interpolate_to_scph) {
        memory->deallocate(kmap_interpolate_to_scph);
    }
    if (mindist_list_scph) {
        memory->deallocate(mindist_list_scph);
    }
    if (evec_harmonic) {
        memory->deallocate(evec_harmonic);
    }
    if (omega2_harmonic) {
        memory->deallocate(omega2_harmonic);
    }
    if (vec_for_v3) {
        memory->deallocate(vec_for_v3);
    }
    if (vec_for_v4) {
        memory->deallocate(vec_for_v4);
    }
    if (invmass_for_v3) {
        memory->deallocate(invmass_for_v3);
    }
    if (invmass_for_v4) {
        memory->deallocate(invmass_for_v4);
    }
    if (evec_index3) {
        memory->deallocate(evec_index3);
    }
    if (evec_index4) {
        memory->deallocate(evec_index4);
    }
    if (fcs_group) {
        memory->deallocate(fcs_group);
    }
    if (fcs_group2) {
        memory->deallocate(fcs_group2);
    }
    if (exp_phase) {
        memory->deallocate(exp_phase);
    }
    if (exp_phase3) {
        memory->deallocate(exp_phase3);
    }
    if (mat_transform_sym) {
        memory->deallocate(mat_transform_sym);
    }
    if (small_group_at_k) {
        memory->deallocate(small_group_at_k);
    }
    if (symop_minus_at_k) {
        memory->deallocate(symop_minus_at_k);
    }
    if (kpoint_map_symmetry) {
        memory->deallocate(kpoint_map_symmetry);
    }
}

void Scph::setup_scph()
{
    relax_coordinate = false;
    MPI_Bcast(&relax_coordinate, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD);

    setup_kmesh();
    setup_eigvecs();
    setup_transform_ifc();
    setup_pp_interaction();
    setup_transform_symmetry();
}

void Scph::exec_scph()
{
    const auto nk_ref = kpoint->nk;
    const auto ns = dynamical->neval;
    const auto Tmin = system->Tmin;
    const auto Tmax = system->Tmax;
    const auto dT = system->dT;

    std::complex<double> ****delta_dymat_scph = nullptr;

    const auto NT = static_cast<unsigned int>((Tmax - Tmin) / dT) + 1;

    MPI_Bcast(&restart_scph, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD);
    MPI_Bcast(&selfenergy_offdiagonal, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD);
    MPI_Bcast(&ialgo, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

    memory->allocate(delta_dymat_scph, NT, ns, ns, nk_interpolate);

    if (restart_scph) {

        // Read anharmonic correction to the dynamical matrix from the existing file
        load_scph_dymat_from_file(delta_dymat_scph);

    } else {

        if (dynamical->nonanalytic == 3) {
            error->exit("exec_scph",
                        "Sorry, NONANALYTIC=3 can't be used for the main loop of the SCPH calculation.");
        }
        // Solve the SCPH equation and obtain the correction to the dynamical matrix
        exec_scph_main(delta_dymat_scph);

        if (mympi->my_rank == 0) {
            store_scph_dymat_to_file(delta_dymat_scph);
            write_anharmonic_correction_fc2(delta_dymat_scph, NT);
        }
    }

    if (kpoint->kpoint_mode == 2) {
        if (thermodynamics->calc_FE_bubble) {
            compute_free_energy_bubble_SCPH(kmesh_interpolate,
                                            delta_dymat_scph);
        }
    }

    postprocess(delta_dymat_scph);

    memory->deallocate(delta_dymat_scph);
}

void Scph::postprocess(std::complex<double> ****delta_dymat_scph)
{
    double ***eval_anharm = nullptr;
    const auto nk_ref = kpoint->nk;
    const auto ns = dynamical->neval;
    const auto Tmin = system->Tmin;
    const auto Tmax = system->Tmax;
    const auto dT = system->dT;
    const auto NT = static_cast<unsigned int>((Tmax - Tmin) / dT) + 1;

    if (mympi->my_rank == 0) {

        std::cout << '\n';
        std::cout << " Running postprocess of SCPH (calculation of free energy, MSD, DOS)" << std::endl;
        std::cout << " The number of temperature points: " << std::setw(4) << NT << std::endl;
        std::cout << "   ";

        std::complex<double> ***evec_tmp = nullptr;

        memory->allocate(eval_anharm, NT, nk_ref, ns);
        memory->allocate(evec_tmp, nk_ref, ns, ns);

        double **dos_scph = nullptr;
        double ***pdos_scph = nullptr;
        double *heat_capacity = nullptr;
        double *FE_QHA = nullptr;
        double *dFE_scph = nullptr;
        double **msd_scph = nullptr;
        double ***ucorr_scph = nullptr;

        if (kpoint->kpoint_mode == 2) {
            if (dos->compute_dos) {
                memory->allocate(dos_scph, NT, dos->n_energy);

                if (dos->projected_dos) {
                    memory->allocate(pdos_scph, NT, ns, dos->n_energy);
                }
            }
            memory->allocate(heat_capacity, NT);
            memory->allocate(FE_QHA, NT);
            memory->allocate(dFE_scph, NT);

            if (writes->print_msd) {
                memory->allocate(msd_scph, NT, ns);
            }
            if (writes->print_ucorr) {
                memory->allocate(ucorr_scph, NT, ns, ns);
            }
        }

        for (auto iT = 0; iT < NT; ++iT) {
            auto T = Tmin + dT * static_cast<double>(iT);

            exec_interpolation(kmesh_interpolate,
                               delta_dymat_scph[iT],
                               nk_ref,
                               kpoint->xk,
                               kpoint->kvec_na,
                               eval_anharm[iT],
                               evec_tmp);

            if (kpoint->kpoint_mode == 2) {

                if (dos->compute_dos) {
                    dos->calc_dos_from_given_frequency(eval_anharm[iT],
                                                       dos_scph[iT]);
                }

                heat_capacity[iT] = thermodynamics->Cv_tot(T,
                                                           kpoint->nk_irred,
                                                           ns,
                                                           kpoint->kpoint_irred_all,
                                                           &kpoint->weight_k[0],
                                                           eval_anharm[iT]);

                FE_QHA[iT] = thermodynamics->free_energy_QHA(T,
                                                             kpoint->nk_irred,
                                                             ns,
                                                             kpoint->kpoint_irred_all,
                                                             &kpoint->weight_k[0],
                                                             eval_anharm[iT]);

                dFE_scph[iT] = thermodynamics->FE_scph_correction(iT,
                                                                  eval_anharm[iT],
                                                                  evec_tmp);

                if (writes->print_msd) {
                    double shift[3]{0.0, 0.0, 0.0};

                    for (auto is = 0; is < ns; ++is) {
                        msd_scph[iT][is] = thermodynamics->disp_corrfunc(T, is, is,
                                                                         shift, kpoint->nk,
                                                                         ns,
                                                                         kpoint->xk,
                                                                         eval_anharm[iT],
                                                                         evec_tmp);
                    }
                }

                if (writes->print_ucorr) {
                    double shift[3];
                    for (auto i = 0; i < 3; ++i) shift[i] = static_cast<double>(writes->shift_ucorr[i]);

                    for (auto is = 0; is < ns; ++is) {
                        for (auto js = 0; js < ns; ++js) {
                            ucorr_scph[iT][is][js] = thermodynamics->disp_corrfunc(T, is, js,
                                                                                   shift, kpoint->nk,
                                                                                   ns,
                                                                                   kpoint->xk,
                                                                                   eval_anharm[iT],
                                                                                   evec_tmp);
                        }
                    }
                }
            }

            std::cout << '.' << std::flush;
            if (iT % 25 == 24) {
                std::cout << std::endl;
                std::cout << std::setw(3);
            }
        }
        std::cout << "\n\n";

        if (kpoint->kpoint_mode == 0) {
            writes->write_scph_energy(eval_anharm);
        } else if (kpoint->kpoint_mode == 1) {
            writes->write_scph_bands(eval_anharm);
        } else if (kpoint->kpoint_mode == 2) {
            if (dos->compute_dos) {
                writes->write_scph_dos(dos_scph);
            }
            writes->write_scph_thermodynamics(heat_capacity, FE_QHA, dFE_scph);
            if (writes->print_msd) {
                writes->write_scph_msd(msd_scph);
            }
            if (writes->print_ucorr) {
                writes->write_scph_ucorr(ucorr_scph);
            }
        }
        memory->deallocate(eval_anharm);
        memory->deallocate(evec_tmp);

        if (dos_scph) memory->deallocate(dos_scph);
        if (pdos_scph) memory->deallocate(pdos_scph);
        if (heat_capacity) memory->deallocate(heat_capacity);
        if (FE_QHA) memory->deallocate(FE_QHA);
        if (dFE_scph) memory->deallocate(dFE_scph);

        //std::cout << "\n done.\n" << std::flush;
    }
}


void Scph::load_scph_dymat_from_file(std::complex<double> ****dymat_out)
{
    const auto ns = dynamical->neval;
    const auto Tmin = system->Tmin;
    const auto Tmax = system->Tmax;
    const auto dT = system->dT;
    auto NT = static_cast<unsigned int>((Tmax - Tmin) / dT) + 1;
    std::vector<double> Temp_array(NT);

    for (int i = 0; i < NT; ++i) {
        Temp_array[i] = Tmin + dT * static_cast<double>(i);
    }

    if (mympi->my_rank == 0) {

        auto consider_offdiagonal = selfenergy_offdiagonal;
        double temp;
        std::ifstream ifs_dymat;
        auto file_dymat = input->job_title + ".scph_dymat";
        bool consider_offdiag_tmp;
        unsigned int nk_interpolate_ref[3];
        unsigned int nk_scph_tmp[3];
        double Tmin_tmp, Tmax_tmp, dT_tmp;
        double dymat_real, dymat_imag;
        std::string str_dummy;
        int nonanalytic_tmp;

        std::cout << " RESTART_SCPH is true." << std::endl;
        std::cout << " Dynamical matrix is read from file ...";

        ifs_dymat.open(file_dymat.c_str(), std::ios::in);

        if (!ifs_dymat) {
            error->exit("load_scph_dymat_from_file",
                        "Cannot open scph_dymat file");
        }

        // Read computational settings from file and check the consistency.
        ifs_dymat >> nk_interpolate_ref[0] >> nk_interpolate_ref[1] >> nk_interpolate_ref[2];
        ifs_dymat >> nk_scph_tmp[0] >> nk_scph_tmp[1] >> nk_scph_tmp[2];
        ifs_dymat >> Tmin_tmp >> Tmax_tmp >> dT_tmp;
        ifs_dymat >> nonanalytic_tmp >> consider_offdiag_tmp;

        if (nk_interpolate_ref[0] != kmesh_interpolate[0] ||
            nk_interpolate_ref[1] != kmesh_interpolate[1] ||
            nk_interpolate_ref[2] != kmesh_interpolate[2]) {
            error->exit("load_scph_dymat_from_file",
                        "The number of KMESH_INTERPOLATE is not consistent");
        }
        if (nk_scph_tmp[0] != kmesh_scph[0] ||
            nk_scph_tmp[1] != kmesh_scph[1] ||
            nk_scph_tmp[2] != kmesh_scph[2]) {
            error->exit("load_scph_dymat_from_file",
                        "The number of KMESH_SCPH is not consistent");
        }
        if (nonanalytic_tmp != dynamical->nonanalytic) {
            error->warn("load_scph_dymat_from_file",
                        "The NONANALYTIC tag is not consistent");
        }
        if (consider_offdiag_tmp != consider_offdiagonal) {
            error->exit("load_scph_dymat_from_file",
                        "The SELF_OFFDIAG tag is not consistent");
        }

        // Check if the precalculated data for the given temperature range exists
        int NT_ref = static_cast<unsigned int>((Tmax_tmp - Tmin_tmp) / dT_tmp) + 1;
        std::vector<double> Temp_array_ref(NT_ref);
        for (int i = 0; i < NT_ref; ++i) {
            Temp_array_ref[i] = Tmin_tmp + dT_tmp * static_cast<double>(i);
        }
        std::vector<int> flag_load(NT_ref);
        for (int i = 0; i < NT_ref; ++i) {
            flag_load[i] = 0;
            for (int j = 0; j < NT; ++j) {
                if (std::abs(Temp_array_ref[i] - Temp_array[j]) < eps6) {
                    flag_load[i] = 1;
                    break;
                }
            }
        }
        int icount = 0;
        for (int iT = 0; iT < NT_ref; ++iT) {
            ifs_dymat >> str_dummy >> temp;
            for (int is = 0; is < ns; ++is) {
                for (int js = 0; js < ns; ++js) {
                    for (int ik = 0; ik < nk_interpolate; ++ik) {
                        ifs_dymat >> dymat_real >> dymat_imag;
                        if (flag_load[iT]) {
                            dymat_out[icount][is][js][ik]
                                = std::complex<double>(dymat_real, dymat_imag);
                        }
                    }
                }
            }
            if (flag_load[iT]) icount += 1;
        }

        ifs_dymat.close();

        if (icount != NT) {
            error->exit("load_scph_dymat_from_file",
                        "The temperature information is not consistent");
        }
        std::cout << " done." << std::endl;
    }
    // Broadcast to all MPI threads
    mpi_bcast_complex(dymat_out, NT, nk_interpolate, ns);
}

void Scph::store_scph_dymat_to_file(std::complex<double> ****dymat_in)
{
    int i;
    const auto ns = dynamical->neval;
    const auto Tmin = system->Tmin;
    const auto Tmax = system->Tmax;
    const auto dT = system->dT;
    std::ofstream ofs_dymat;
    auto file_dymat = input->job_title + ".scph_dymat";

    const auto NT = static_cast<unsigned int>((Tmax - Tmin) / dT) + 1;

    ofs_dymat.open(file_dymat.c_str(), std::ios::out);

    if (!ofs_dymat) {
        error->exit("store_scph_dymat_to_file",
                    "Cannot open scph_dymat file");
    }
    for (i = 0; i < 3; ++i) {
        ofs_dymat << std::setw(5) << kmesh_interpolate[i];
    }
    ofs_dymat << std::endl;
    for (i = 0; i < 3; ++i) {
        ofs_dymat << std::setw(5) << kmesh_scph[i];
    }
    ofs_dymat << std::endl;
    ofs_dymat << std::setw(10) << Tmin;
    ofs_dymat << std::setw(10) << Tmax;
    ofs_dymat << std::setw(10) << dT << std::endl;
    ofs_dymat << std::setw(5) << dynamical->nonanalytic;
    ofs_dymat << std::setw(5) << selfenergy_offdiagonal << std::endl;

    for (auto iT = 0; iT < NT; ++iT) {
        const auto temp = Tmin + static_cast<double>(iT) * dT;
        ofs_dymat << "# " << temp << std::endl;
        for (auto is = 0; is < ns; ++is) {
            for (auto js = 0; js < ns; ++js) {
                for (auto ik = 0; ik < nk_interpolate; ++ik) {
                    ofs_dymat << std::setprecision(15)
                        << std::setw(25) << dymat_in[iT][is][js][ik].real();
                    ofs_dymat << std::setprecision(15)
                        << std::setw(25) << dymat_in[iT][is][js][ik].imag();
                    ofs_dymat << std::endl;
                }
            }
        }
    }
    ofs_dymat.close();
    std::cout << "  " << std::setw(input->job_title.length() + 12) << std::left << file_dymat;
    std::cout << " : Anharmonic dynamical matrix (restart file)" << std::endl;
}

void Scph::exec_scph_main(std::complex<double> ****dymat_anharm)
{
    int iT, ik, is;
    const auto nk = nk_scph;
    const auto ns = dynamical->neval;
    const auto nk_reduced_scph = kp_irred_scph.size();
    const auto nk_irred_interpolate = kp_irred_interpolate.size();
    const auto Tmin = system->Tmin;
    const auto Tmax = system->Tmax;
    const auto dT = system->dT;
    double ***omega2_anharm;
    std::complex<double> ***evec_anharm_tmp;
    std::complex<double> ***v3_array_all;
    std::complex<double> ***v4_array_all;

    std::vector<double> vec_temp;

    const auto NT = static_cast<unsigned int>((Tmax - Tmin) / dT) + 1;

    // Find the degeneracy at each irreducible k points.

    std::vector<int> *degeneracy_at_k;
    memory->allocate(degeneracy_at_k, nk_reduced_scph);
    find_degeneracy(degeneracy_at_k,
                    nk_reduced_scph,
                    kp_irred_scph,
                    omega2_harmonic);

    // Compute matrix element of 4-phonon interaction

    memory->allocate(omega2_anharm, NT, nk, ns);
    memory->allocate(evec_anharm_tmp, nk, ns, ns);
    memory->allocate(v4_array_all, nk_irred_interpolate * nk_scph,
                     ns * ns, ns * ns);

    // Calculate v4 array. 
    // This operation is the most expensive part of the calculation.
    if (selfenergy_offdiagonal & (ialgo == 1)) {
        compute_V4_elements_mpi_over_band(v4_array_all,
                                          evec_harmonic,
                                          selfenergy_offdiagonal);
    } else {
        compute_V4_elements_mpi_over_kpoint(v4_array_all,
                                            evec_harmonic,
                                            selfenergy_offdiagonal,
                                            relax_coordinate);
    }

    if (relax_coordinate) {
        memory->allocate(v3_array_all, nk, ns, ns * ns);
        compute_V3_elements_mpi_over_kpoint(v3_array_all,
                                            evec_harmonic,
                                            selfenergy_offdiagonal);
    }

    if (mympi->my_rank == 0) {

        std::complex<double> ***cmat_convert;
        memory->allocate(cmat_convert, nk, ns, ns);

        vec_temp.clear();

        if (lower_temp) {
            for (iT = NT - 1; iT >= 0; --iT) {
                vec_temp.push_back(Tmin + static_cast<double>(iT) * dT);
            }
        } else {
            for (iT = 0; iT < NT; ++iT) {
                vec_temp.push_back(Tmin + static_cast<double>(iT) * dT);
            }
        }

        auto converged_prev = false;

        for (int i = 0; i < vec_temp.size(); ++i) {
            const auto temp = vec_temp[i];

            iT = static_cast<unsigned int>((temp - Tmin) / dT);

            // Initialize phonon eigenvectors with harmonic values

            for (ik = 0; ik < nk; ++ik) {
                for (is = 0; is < ns; ++is) {
                    for (int js = 0; js < ns; ++js) {
                        evec_anharm_tmp[ik][is][js] = evec_harmonic[ik][is][js];
                    }
                }
            }
            if (converged_prev) {
                if (lower_temp) {
                    for (ik = 0; ik < nk; ++ik) {
                        for (is = 0; is < ns; ++is) {
                            omega2_anharm[iT][ik][is] = omega2_anharm[iT + 1][ik][is];
                        }
                    }
                } else {
                    for (ik = 0; ik < nk; ++ik) {
                        for (is = 0; is < ns; ++is) {
                            omega2_anharm[iT][ik][is] = omega2_anharm[iT - 1][ik][is];
                        }
                    }
                }
            }

            compute_anharmonic_frequency(v4_array_all,
                                         omega2_anharm[iT],
                                         evec_anharm_tmp,
                                         temp,
                                         degeneracy_at_k,
                                         converged_prev,
                                         cmat_convert,
                                         selfenergy_offdiagonal);

            calc_new_dymat_with_evec(dymat_anharm[iT],
                                     omega2_anharm[iT],
                                     evec_anharm_tmp);

            if (!warmstart_scph) converged_prev = false;
        }

        memory->deallocate(cmat_convert);

    }

    mpi_bcast_complex(dymat_anharm, NT, nk_interpolate, ns);

    memory->deallocate(omega2_anharm);
    memory->deallocate(v4_array_all);
    memory->deallocate(evec_anharm_tmp);
    memory->deallocate(degeneracy_at_k);
}


void Scph::compute_V3_elements_mpi_over_kpoint(std::complex<double> ***v3_out,
                                               std::complex<double> ***evec_in,
                                               const bool self_offdiag)
{
    // Calculate the matrix elements of quartic terms in reciprocal space.
    // This is the most expensive part of the SCPH calculation.

    auto ns = dynamical->neval;
    auto ns2 = ns * ns;
    auto ns3 = ns * ns * ns;
    unsigned int is, js, ks;
    double phase3[3];
    int loc3[3];
    unsigned int **ind;
    unsigned int i, j;
    std::complex<double> ret;
    long int ii;

    const auto inv2pi = 1.0 / (2.0 * pi);
    const auto dnk_represent = static_cast<double>(nk_represent);
    const auto factor = std::pow(0.5, 2) / static_cast<double>(nk_scph);
    static auto complex_zero = std::complex<double>(0.0, 0.0);
    std::complex<double> *v3_array_at_kpair;
    std::complex<double> ***v3_mpi;

    //  nk2_prod = nk_reduced_interpolate * nk_scph;

    if (mympi->my_rank == 0) {
        if (self_offdiag) {
            std::cout << " SELF_OFFDIAG = 1: Calculating all components of v3_array ... ";
        } else {
            std::cout << " SELF_OFFDIAG = 0: Calculating diagonal components of v3_array ... ";
        }
    }

    memory->allocate(v3_array_at_kpair, ngroup);
    memory->allocate(ind, ngroup, 3);
    memory->allocate(v3_mpi, nk_scph, ns, ns2);

    for (unsigned int ik = mympi->my_rank; ik < nk_scph; ik += mympi->nprocs) {

        for (is = 0; is < ngroup; ++is) v3_array_at_kpair[is] = complex_zero;

        unsigned int ielem = 0;

        for (i = 0; i < ngroup; ++i) {

            std::complex<double> sum_tmp = std::complex<double>(0.0, 0.0);
            for (j = 0; j < 3; ++j) ind[i][j] = evec_index3[ielem][j];

            if (tune_type == 0) {
                for (j = 0; j < fcs_group[i].size(); ++j) {
                    const auto phase = xk_scph[ik][0] * (vec_for_v3[ielem][0][0] - vec_for_v3[ielem][1][0])
                        + xk_scph[ik][1] * (vec_for_v3[ielem][0][1] - vec_for_v3[ielem][1][1])
                        + xk_scph[ik][2] * (vec_for_v3[ielem][0][2] - vec_for_v3[ielem][1][2]);

                    const int loc = nint(phase * dnk_represent * inv2pi) % nk_represent + nk_represent - 1;

                    sum_tmp += fcs_group[i][j] * invmass_for_v3[ielem] * exp_phase[loc];

                    ++ielem;
                }
            } else if (tune_type == 1) {
                for (j = 0; j < fcs_group2[i].size(); ++j) {

                    for (ii = 0; ii < 3; ++ii) {
                        phase3[ii] = xk_scph[ik][ii] * (vec_for_v3[ielem][0][ii] - vec_for_v3[ielem][1][ii]);
                        loc3[ii] = nint(phase3[ii] * dnk[ii] * inv2pi) % nk_grid[ii] + nk_grid[ii] - 1;
                    }

                    sum_tmp += fcs_group[i][j] * invmass_for_v3[ielem]
                        * exp_phase3[loc3[0]][loc3[1]][loc3[2]];

                    ++ielem;
                }
            }
            v3_array_at_kpair[i] = sum_tmp;
        }

#pragma omp parallel for private(is)
        for (ii = 0; ii < ns; ++ii) {
            for (is = 0; is < ns2; ++is) {
                v3_mpi[ik][ii][is] = complex_zero;
                v3_out[ik][ii][is] = complex_zero;
            }
        }

        if (self_offdiag) {

            // All matrix elements will be calculated when considering the off-diagonal
            // elements of phonon selfenergy.

#pragma omp parallel for private(is, js, ks, ret, i)
            for (ii = 0; ii < ns3; ++ii) {
                is = ii / ns2;
                js = (ii - ns2 * is) / ns;
                ks = ii % ns;

                ret = std::complex<double>(0.0, 0.0);

                for (i = 0; i < ngroup; ++i) {

                    ret += v3_array_at_kpair[i]
                        * evec_in[0][is][ind[i][0]]
                        * evec_in[ik][js][ind[i][1]]
                        * std::conj(evec_in[ik][ks][ind[i][2]]);
                }

                v3_mpi[ik][is][ns * js + ks] = factor * ret;
            }

        } else {

            // Only diagonal elements will be computed when neglecting the polarization mixing.

            if (ik == 0) {
#pragma omp parallel for private(is, js, ks, ret, i)
                for (ii = 0; ii < ns3; ++ii) {
                    is = ii / ns2;
                    js = (ii - ns2 * is) / ns;
                    ks = ii % ns;

                    ret = std::complex<double>(0.0, 0.0);

                    for (i = 0; i < ngroup; ++i) {

                        ret += v3_array_at_kpair[i]
                            * evec_in[0][is][ind[i][0]]
                            * evec_in[ik][js][ind[i][1]]
                            * std::conj(evec_in[ik][ks][ind[i][2]]);
                    }

                    v3_mpi[ik][is][ns * js + ks] = factor * ret;
                }
            } else {

#pragma omp parallel for private(is, js, ret, i)
                for (ii = 0; ii < ns2; ++ii) {
                    is = ii / ns;
                    js = ii % ns;

                    ret = std::complex<double>(0.0, 0.0);

                    for (i = 0; i < ngroup; ++i) {

                        ret += v3_array_at_kpair[i]
                            * evec_in[0][is][ind[i][0]]
                            * evec_in[ik][js][ind[i][1]]
                            * std::conj(evec_in[ik][js][ind[i][2]]);
                    }

                    v3_mpi[ik][is][(ns + 1) * js] = factor * ret;
                }
            }
        }
    }

    memory->deallocate(v3_array_at_kpair);
    memory->deallocate(ind);
#ifdef MPI_CXX_DOUBLE_COMPLEX
    MPI_Allreduce(&v3_mpi[0][0][0], &v3_out[0][0][0], nk_scph * ns3,
        MPI_CXX_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
#else
    MPI_Allreduce(&v3_mpi[0][0][0], &v3_out[0][0][0], nk_scph * ns3,
                  MPI_COMPLEX16, MPI_SUM, MPI_COMM_WORLD);
#endif

    memory->deallocate(v3_mpi);

    zerofill_elements_acoustic_at_gamma(omega2_harmonic, v3_out, 3);

    if (mympi->my_rank == 0) {
        std::cout << " done !" << std::endl;
        timer->print_elapsed();
    }
}

void Scph::compute_V4_elements_mpi_over_kpoint(std::complex<double> ***v4_out,
                                               std::complex<double> ***evec_in,
                                               const bool self_offdiag,
                                               const bool relax)
{
    // Calculate the matrix elements of quartic terms in reciprocal space.
    // This is the most expensive part of the SCPH calculation.

    const auto nk_reduced_interpolate = kp_irred_interpolate.size();
    const auto ns = dynamical->neval;
    const auto ns2 = ns * ns;
    const auto ns3 = ns * ns * ns;
    const auto ns4 = ns * ns * ns * ns;
    unsigned int is, js, ks, ls;
    double phase3[3];
    int loc3[3];
    unsigned int **ind;
    unsigned int i, j;
    std::complex<double> ret;
    long int ii;

    const auto inv2pi = 1.0 / (2.0 * pi);
    const auto dnk_represent = static_cast<double>(nk_represent);
    const auto factor = std::pow(0.5, 2) / static_cast<double>(nk_scph);
    static auto complex_zero = std::complex<double>(0.0, 0.0);
    std::complex<double> *v4_array_at_kpair;
    std::complex<double> ***v4_mpi;

    const auto nk2_prod = nk_reduced_interpolate * nk_scph;

    if (mympi->my_rank == 0) {
        if (self_offdiag) {
            std::cout << " SELF_OFFDIAG = 1: Calculating all components of v4_array ... ";
        } else {
            std::cout << " SELF_OFFDIAG = 0: Calculating diagonal components of v4_array ... ";
        }
    }

    memory->allocate(v4_array_at_kpair, ngroup2);
    memory->allocate(ind, ngroup2, 4);
    memory->allocate(v4_mpi, nk2_prod, ns2, ns2);

    for (int ik_prod = mympi->my_rank; ik_prod < nk2_prod; ik_prod += mympi->nprocs) {
        const auto ik = ik_prod / nk_scph;
        const auto jk = ik_prod % nk_scph;

        const unsigned int knum = kmap_interpolate_to_scph[kp_irred_interpolate[ik][0].knum];

        for (is = 0; is < ngroup2; ++is) v4_array_at_kpair[is] = complex_zero;

        unsigned int ielem = 0;

        for (i = 0; i < ngroup2; ++i) {

            auto sum_tmp = std::complex<double>(0.0, 0.0);
            for (j = 0; j < 4; ++j) ind[i][j] = evec_index4[ielem][j];

            if (tune_type == 0) {
                for (j = 0; j < fcs_group2[i].size(); ++j) {
                    const auto phase = - xk_scph[knum][0] * vec_for_v4[ielem][0][0]
                        - xk_scph[knum][1] * vec_for_v4[ielem][0][1]
                        - xk_scph[knum][2] * vec_for_v4[ielem][0][2]
                        + xk_scph[jk][0] * (vec_for_v4[ielem][1][0] - vec_for_v4[ielem][2][0])
                        + xk_scph[jk][1] * (vec_for_v4[ielem][1][1] - vec_for_v4[ielem][2][1])
                        + xk_scph[jk][2] * (vec_for_v4[ielem][1][2] - vec_for_v4[ielem][2][2]);

                    const auto loc = nint(phase * dnk_represent * inv2pi) % nk_represent + nk_represent - 1;

                    sum_tmp += fcs_group2[i][j] * invmass_for_v4[ielem] * exp_phase[loc];

                    ++ielem;
                }
            } else if (tune_type == 1) {
                for (j = 0; j < fcs_group2[i].size(); ++j) {

                    for (ii = 0; ii < 3; ++ii) {
                        phase3[ii] = - xk_scph[knum][ii] * vec_for_v4[ielem][0][ii]
                            + xk_scph[jk][ii] * (vec_for_v4[ielem][1][ii] - vec_for_v4[ielem][2][ii]);
                        loc3[ii] = nint(phase3[ii] * dnk[ii] * inv2pi) % nk_grid[ii] + nk_grid[ii] - 1;
                    }

                    sum_tmp += fcs_group2[i][j] * invmass_for_v4[ielem]
                        * exp_phase3[loc3[0]][loc3[1]][loc3[2]];

                    ++ielem;
                }
            }
            v4_array_at_kpair[i] = sum_tmp;
        }

#pragma omp parallel for private(is)
        for (ii = 0; ii < ns2; ++ii) {
            for (is = 0; is < ns2; ++is) {
                v4_mpi[ik_prod][ii][is] = complex_zero;
                v4_out[ik_prod][ii][is] = complex_zero;
            }
        }

        if (self_offdiag) {

            // All matrix elements will be calculated when considering the off-diagonal
            // elements of phonon selfenergy.

#pragma omp parallel for private(is, js, ks, ls, ret, i)
            for (ii = 0; ii < ns4; ++ii) {
                is = ii / ns3;
                js = (ii - ns3 * is) / ns2;
                ks = (ii - ns3 * is - ns2 * js) / ns;
                ls = ii % ns;

                if (is < js) continue;

                ret = std::complex<double>(0.0, 0.0);

                for (i = 0; i < ngroup2; ++i) {

                    ret += v4_array_at_kpair[i]
                        * evec_in[knum][is][ind[i][0]]
                        * std::conj(evec_in[knum][js][ind[i][1]])
                        * evec_in[jk][ks][ind[i][2]]
                        * std::conj(evec_in[jk][ls][ind[i][3]]);
                }

                v4_mpi[ik_prod][ns * is + js][ns * ks + ls] = factor * ret;
            }

        } else {

            // Only diagonal elements will be computed when neglecting the polarization mixing.

            if (relax && (knum == 0 || jk == 0)) {

#pragma omp parallel for private(is, js, ks, ls, ret, i)
                for (ii = 0; ii < ns4; ++ii) {
                    is = ii / ns3;
                    js = (ii - ns3 * is) / ns2;
                    ks = (ii - ns3 * is - ns2 * js) / ns;
                    ls = ii % ns;

                    if (is < js) continue;

                    ret = std::complex<double>(0.0, 0.0);

                    for (i = 0; i < ngroup2; ++i) {

                        ret += v4_array_at_kpair[i]
                            * evec_in[knum][is][ind[i][0]]
                            * std::conj(evec_in[knum][js][ind[i][1]])
                            * evec_in[jk][ks][ind[i][2]]
                            * std::conj(evec_in[jk][ls][ind[i][3]]);
                    }

                    v4_mpi[ik_prod][ns * is + js][ns * ks + ls] = factor * ret;
                }

            } else {

#pragma omp parallel for private(is, js, ret, i)
                for (ii = 0; ii < ns2; ++ii) {
                    is = ii / ns;
                    js = ii % ns;

                    ret = std::complex<double>(0.0, 0.0);

                    for (i = 0; i < ngroup2; ++i) {

                        ret += v4_array_at_kpair[i]
                            * evec_in[knum][is][ind[i][0]]
                            * std::conj(evec_in[knum][is][ind[i][1]])
                            * evec_in[jk][js][ind[i][2]]
                            * std::conj(evec_in[jk][js][ind[i][3]]);
                    }

                    v4_mpi[ik_prod][(ns + 1) * is][(ns + 1) * js] = factor * ret;
                }
            }
        }
    }

    memory->deallocate(v4_array_at_kpair);
    memory->deallocate(ind);
#ifdef MPI_CXX_DOUBLE_COMPLEX
    MPI_Allreduce(&v4_mpi[0][0][0], &v4_out[0][0][0], nk2_prod*ns4, 
        MPI_CXX_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
#else
    MPI_Allreduce(&v4_mpi[0][0][0], &v4_out[0][0][0], nk2_prod * ns4,
                  MPI_COMPLEX16, MPI_SUM, MPI_COMM_WORLD);
#endif

    memory->deallocate(v4_mpi);


    zerofill_elements_acoustic_at_gamma(omega2_harmonic, v4_out, 4);


    if (mympi->my_rank == 0) {
        std::cout << " done !" << std::endl;
        timer->print_elapsed();
    }
}

void Scph::compute_V4_elements_mpi_over_band(std::complex<double> ***v4_out,
                                             std::complex<double> ***evec_in,
                                             const bool self_offdiag)
{
    // Calculate the matrix elements of quartic terms in reciprocal space.
    // This is the most expensive part of the SCPH calculation.

    int ik_prod;
    auto nk_reduced_interpolate = kp_irred_interpolate.size();
    auto ns = dynamical->neval;
    auto ns2 = ns * ns;
    auto ns4 = ns * ns * ns * ns;
    int is, js;
    unsigned int knum;
    double phase3[3];
    int loc3[3];
    unsigned int **ind;
    unsigned int i, j;
    long int *nset_mpi;

    auto inv2pi = 1.0 / (2.0 * pi);
    auto dnk_represent = static_cast<double>(nk_represent);
    auto factor = std::pow(0.5, 2) / static_cast<double>(nk_scph);
    static auto complex_zero = std::complex<double>(0.0, 0.0);
    std::complex<double> *v4_array_at_kpair;
    std::complex<double> ***v4_mpi;

    std::vector<int> ik_vec, jk_vec, is_vec, js_vec;

    auto nk2_prod = nk_reduced_interpolate * nk_scph;

    if (mympi->my_rank == 0) {
        if (self_offdiag) {
            std::cout << " IALGO = 1 : Use different algorithm efficient when ns >> nk" << std::endl;
            std::cout << " SELF_OFFDIAG = 1: Calculating all components of v4_array ... ";
        } else {
            error->exit("compute_V4_elements_mpi_over_kpoint",
                        "This function can be used only when SELF_OFFDIAG = 1");
        }
    }

    memory->allocate(nset_mpi, mympi->nprocs);

    long int nset_tot = nk2_prod * ((ns2 - ns) / 2 + ns);
    long int nset_each = nset_tot / mympi->nprocs;
    long int nres = nset_tot - nset_each * mympi->nprocs;

    for (i = 0; i < mympi->nprocs; ++i) {
        nset_mpi[i] = nset_each;
        if (nres > i) {
            nset_mpi[i] += 1;
        }
    }

    MPI_Bcast(&nset_mpi[0], mympi->nprocs, MPI_LONG, 0, MPI_COMM_WORLD);
    long int nstart = 0;
    for (i = 0; i < mympi->my_rank; ++i) {
        nstart += nset_mpi[i];
    }
    long int nend = nstart + nset_mpi[mympi->my_rank];
    nset_each = nset_mpi[mympi->my_rank];
    memory->deallocate(nset_mpi);

    ik_vec.clear();
    jk_vec.clear();
    is_vec.clear();
    js_vec.clear();

    long int icount = 0;
    for (ik_prod = 0; ik_prod < nk2_prod; ++ik_prod) {
        for (is = 0; is < ns; ++is) {
            for (js = 0; js < ns; ++js) {
                if (is < js) continue;

                if (icount >= nstart && icount < nend) {
                    ik_vec.push_back(ik_prod / nk_scph);
                    jk_vec.push_back(ik_prod % nk_scph);
                    is_vec.push_back(is);
                    js_vec.push_back(js);
                }
                ++icount;
            }
        }
    }

    memory->allocate(v4_array_at_kpair, ngroup2);
    memory->allocate(ind, ngroup2, 4);
    memory->allocate(v4_mpi, nk2_prod, ns2, ns2);

    for (ik_prod = 0; ik_prod < nk2_prod; ++ik_prod) {
#pragma omp parallel for private (js)
        for (is = 0; is < ns2; ++is) {
            for (js = 0; js < ns2; ++js) {
                v4_mpi[ik_prod][is][js] = complex_zero;
                v4_out[ik_prod][is][js] = complex_zero;
            }
        }
    }

    int ik_old = -1;
    int jk_old = -1;

    if (mympi->my_rank == 0) {
        std::cout << " Total number of sets to compute : " << nset_each << std::endl;
    }

    for (long int ii = 0; ii < nset_each; ++ii) {

        auto ik_now = ik_vec[ii];
        auto jk_now = jk_vec[ii];
        auto is_now = is_vec[ii];
        auto js_now = js_vec[ii];

        if (!(ik_now == ik_old && jk_now == jk_old)) {

            // Update v4_array_at_kpair and ind

            knum = kmap_interpolate_to_scph[kp_irred_interpolate[ik_now][0].knum];

            for (is = 0; is < ngroup2; ++is) v4_array_at_kpair[is] = complex_zero;

            unsigned int ielem = 0;

            for (i = 0; i < ngroup2; ++i) {

                std::complex<double> sum_tmp = std::complex<double>(0.0, 0.0);
                for (j = 0; j < 4; ++j) ind[i][j] = evec_index4[ielem][j];

                if (tune_type == 0) {
                    for (j = 0; j < fcs_group2[i].size(); ++j) {
                        auto phase = - xk_scph[knum][0] * vec_for_v4[ielem][0][0]
                            - xk_scph[knum][1] * vec_for_v4[ielem][0][1]
                            - xk_scph[knum][2] * vec_for_v4[ielem][0][2]
                            + xk_scph[jk_now][0] * (vec_for_v4[ielem][1][0] - vec_for_v4[ielem][2][0])
                            + xk_scph[jk_now][1] * (vec_for_v4[ielem][1][1] - vec_for_v4[ielem][2][1])
                            + xk_scph[jk_now][2] * (vec_for_v4[ielem][1][2] - vec_for_v4[ielem][2][2]);

                        auto loc = nint(phase * dnk_represent * inv2pi) % nk_represent + nk_represent - 1;

                        sum_tmp += fcs_group2[i][j] * invmass_for_v4[ielem] * exp_phase[loc];

                        ++ielem;
                    }
                } else if (tune_type == 1) {
                    for (j = 0; j < fcs_group2[i].size(); ++j) {

                        for (unsigned int k = 0; k < 3; ++k) {
                            phase3[k] = - xk_scph[knum][k] * vec_for_v4[ielem][0][k]
                                + xk_scph[jk_now][k] * (vec_for_v4[ielem][1][k] - vec_for_v4[ielem][2][k]);
                            loc3[k] = nint(phase3[k] * dnk[k] * inv2pi) % nk_grid[k] + nk_grid[k] - 1;
                        }

                        sum_tmp += fcs_group2[i][j] * invmass_for_v4[ielem]
                            * exp_phase3[loc3[0]][loc3[1]][loc3[2]];

                        ++ielem;
                    }
                }
                v4_array_at_kpair[i] = sum_tmp;
            }
            ik_old = ik_now;
            jk_old = jk_now;
        }

        ik_prod = ik_now * nk_scph + jk_now;
        int is_prod = ns * is_now + js_now;

#pragma omp parallel for private (i)
        for (js = 0; js < ns2; ++js) {

            unsigned int ks = js / ns;
            unsigned int ls = js % ns;

            auto ret = std::complex<double>(0.0, 0.0);

            for (i = 0; i < ngroup2; ++i) {

                ret += v4_array_at_kpair[i]
                    * evec_in[knum][is_now][ind[i][0]]
                    * std::conj(evec_in[knum][js_now][ind[i][1]])
                    * evec_in[jk_now][ks][ind[i][2]]
                    * std::conj(evec_in[jk_now][ls][ind[i][3]]);
            }

            v4_mpi[ik_prod][is_prod][js] = factor * ret;
        }

        if (mympi->my_rank == 0) {
            std::cout << " SET " << ii + 1 << " done. " << std::endl;
        }

    } // loop over nk2_prod*ns2

    memory->deallocate(v4_array_at_kpair);
    memory->deallocate(ind);
#ifdef MPI_CXX_DOUBLE_COMPLEX
    MPI_Allreduce(&v4_mpi[0][0][0], &v4_out[0][0][0], nk2_prod*ns4, 
        MPI_CXX_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
#else
    MPI_Allreduce(&v4_mpi[0][0][0], &v4_out[0][0][0], nk2_prod * ns4,
                  MPI_COMPLEX16, MPI_SUM, MPI_COMM_WORLD);
#endif

    memory->deallocate(v4_mpi);

    zerofill_elements_acoustic_at_gamma(omega2_harmonic, v4_out, 4);

    if (mympi->my_rank == 0) {
        std::cout << " done !" << std::endl;
        timer->print_elapsed();
    }
}

void Scph::zerofill_elements_acoustic_at_gamma(double **omega2,
                                               std::complex<double> ***v_elems,
                                               const int fc_order) const
{
    // Set V3 or V4 elements involving acoustic modes at Gamma point
    // exactly zero.

    int jk;
    int is, js, ks, ls;
    const auto ns = dynamical->neval;
    bool *is_acoustic;
    memory->allocate(is_acoustic, ns);
    int nacoustic;
    auto threshould = 1.0e-24;
    const auto nk_reduced_interpolate = kp_irred_interpolate.size();
    static auto complex_zero = std::complex<double>(0.0, 0.0);


    if (fc_order < 3 && fc_order > 4) {
        error->exit("zerofill_elements_acoustic_at_gamma",
                    "The fc_order must be either 3 or 4.");
    }

    do {
        nacoustic = 0;
        for (is = 0; is < ns; ++is) {
            if (std::abs(omega2[0][is]) < threshould) {
                is_acoustic[is] = true;
                ++nacoustic;
            } else {
                is_acoustic[is] = false;
            }
        }
        if (nacoustic > 3) {
            error->exit("zerofill_elements_acoustic_at_gamma",
                        "Could not assign acoustic modes at Gamma.");
        }
        threshould *= 2.0;
    } while (nacoustic < 3);


    if (fc_order == 3) {

        // Set V3 to zeros so as to avoid mixing with gamma acoustic modes
        // jk = 0;
        for (is = 0; is < ns; ++is) {
            for (ks = 0; ks < ns; ++ks) {
                for (ls = 0; ls < ns; ++ls) {
                    if (is_acoustic[ks] || is_acoustic[ls]) {
                        v_elems[0][is][ns * ks + ls] = complex_zero;
                    }
                }
            }
        }

        // ik = 0;
        for (jk = 0; jk < nk_scph; ++jk) {
            for (is = 0; is < ns; ++is) {
                if (is_acoustic[is]) {
                    for (ks = 0; ks < ns; ++ks) {
                        for (ls = 0; ls < ns; ++ls) {
                            v_elems[jk][is][ns * ks + ls] = complex_zero;
                        }
                    }
                }
            }
        }

    } else if (fc_order == 4) {
        // Set V4 to zeros so as to avoid mixing with gamma acoustic modes
        // jk = 0;
        for (int ik = 0; ik < nk_reduced_interpolate; ++ik) {
            for (is = 0; is < ns; ++is) {
                for (js = 0; js < ns; ++js) {
                    for (ks = 0; ks < ns; ++ks) {
                        for (ls = 0; ls < ns; ++ls) {
                            if (is_acoustic[ks] || is_acoustic[ls]) {
                                v_elems[nk_scph * ik][ns * is + js][ns * ks + ls] = complex_zero;
                            }
                        }
                    }
                }
            }
        }
        // ik = 0;
        for (jk = 0; jk < nk_scph; ++jk) {
            for (is = 0; is < ns; ++is) {
                for (js = 0; js < ns; ++js) {
                    if (is_acoustic[is] || is_acoustic[js]) {
                        for (ks = 0; ks < ns; ++ks) {
                            for (ls = 0; ls < ns; ++ls) {
                                v_elems[jk][ns * is + js][ns * ks + ls] = complex_zero;
                            }
                        }
                    }
                }
            }
        }

    }

    memory->deallocate(is_acoustic);
}


void Scph::setup_kmesh()
{
    unsigned int ik;
    unsigned int i;
    double xtmp[3];

    // Setup k points for SCPH equation
    MPI_Bcast(&kmesh_scph[0], 3, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    MPI_Bcast(&kmesh_interpolate[0], 3, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

    // Set up k points for Fourier interpolation
    nk_scph = kmesh_scph[0] * kmesh_scph[1] * kmesh_scph[2];
    nk_interpolate = kmesh_interpolate[0] * kmesh_interpolate[1] * kmesh_interpolate[2];

    memory->allocate(xk_scph, nk_scph, 3);
    memory->allocate(xk_interpolate, nk_interpolate, 3);
    memory->allocate(kvec_na_scph, nk_scph, 3);

    kpoint->gen_kmesh(true,
                      kmesh_scph,
                      xk_scph,
                      kp_irred_scph);

    kpoint->gen_kmesh(true,
                      kmesh_interpolate,
                      xk_interpolate,
                      kp_irred_interpolate);


    for (ik = 0; ik < nk_scph; ++ik) {
        for (i = 0; i < 3; ++i) {
            kvec_na_scph[ik][i] = 0.0;
        }
    }

    for (ik = 0; ik < nk_scph; ++ik) {

        for (i = 0; i < 3; ++i) xtmp[i] = xk_scph[ik][i];
        rotvec(xtmp, xtmp, system->rlavec_p, 'T');
        const auto norm = xtmp[0] * xtmp[0] + xtmp[1] * xtmp[1] + xtmp[2] * xtmp[2];

        if (norm > eps) {
            for (i = 0; i < 3; ++i) kvec_na_scph[ik][i] = xk_scph[ik][i] / std::sqrt(norm);
        }
    }

    if (mympi->my_rank == 0) {
        std::cout << " Setting up the SCPH calculations ..." << std::endl << std::endl;
        std::cout << "  Gamma-centered uniform grid with the following mesh density:" << std::endl;
        std::cout << "  nk1:" << std::setw(5) << kmesh_scph[0] << std::endl;
        std::cout << "  nk2:" << std::setw(5) << kmesh_scph[1] << std::endl;
        std::cout << "  nk3:" << std::setw(5) << kmesh_scph[2] << std::endl;
        std::cout << std::endl;
        std::cout << "  Number of k points : " << nk_scph << std::endl;
        std::cout << "  Number of irreducible k points : " << kp_irred_scph.size() << std::endl;
        std::cout << std::endl;
        std::cout << "  Fourier interpolation from reciprocal to real space" << std::endl;
        std::cout << "  will be performed with the following mesh density:" << std::endl;
        std::cout << "  nk1:" << std::setw(5) << kmesh_interpolate[0] << std::endl;
        std::cout << "  nk2:" << std::setw(5) << kmesh_interpolate[1] << std::endl;
        std::cout << "  nk3:" << std::setw(5) << kmesh_interpolate[2] << std::endl;
        std::cout << std::endl;
        std::cout << "  Number of k points : " << nk_interpolate << std::endl;
        std::cout << "  Number of irreducible k points : "
            << kp_irred_interpolate.size() << std::endl;
    }

    memory->allocate(kmap_interpolate_to_scph, nk_interpolate);

    for (ik = 0; ik < nk_interpolate; ++ik) {
        for (i = 0; i < 3; ++i) xtmp[i] = xk_interpolate[ik][i];

        const auto loc = kpoint->get_knum(xtmp, kmesh_scph);

        if (loc == -1)
            error->exit("setup_kmesh",
                        "KMESH_INTERPOLATE should be a integral multiple of KMESH_SCPH");
        kmap_interpolate_to_scph[ik] = loc;
    }
}


void Scph::setup_transform_symmetry()
{
    // Construct small_group_at_k, symop_minus_at_k, and
    // mat_transport_sym.

    unsigned int ik;
    unsigned int is, js;
    unsigned int icrd, jcrd;
    double x1[3], x2[3], k[3], k_minus[3], Sk[3], xtmp[3];
    double S_cart[3][3], S_frac[3][3], S_frac_inv[3][3];
    double S_recip[3][3];
    std::complex<double> im(0.0, 1.0);
    std::complex<double> **gamma_tmp;
    bool *flag;

    const auto natmin = system->natmin;
    const auto ns = dynamical->neval;
    const auto nk_irred_interpolate = kp_irred_interpolate.size();

    memory->allocate(gamma_tmp, ns, ns);
    memory->allocate(mat_transform_sym, nk_irred_interpolate,
                     symmetry->nsym, ns, ns);
    memory->allocate(small_group_at_k, nk_irred_interpolate);
    memory->allocate(symop_minus_at_k, nk_irred_interpolate);
    memory->allocate(kpoint_map_symmetry, nk_interpolate);
    memory->allocate(flag, nk_interpolate);

    for (ik = 0; ik < nk_interpolate; ++ik) {
        flag[ik] = false;
    }

    for (ik = 0; ik < kp_irred_interpolate.size(); ++ik) {

        small_group_at_k[ik].clear();
        symop_minus_at_k[ik].clear();

        const auto knum = kp_irred_interpolate[ik][0].knum;
        for (icrd = 0; icrd < 3; ++icrd) {
            k[icrd] = xk_interpolate[knum][icrd];
            k_minus[icrd] = -k[icrd];
        }
        const auto knum_minus = kpoint->get_knum(k_minus, kmesh_interpolate);

        unsigned int isym = 0;

        for (const auto &it : symmetry->SymmListWithMap) {

            for (icrd = 0; icrd < 3; ++icrd) {
                for (jcrd = 0; jcrd < 3; ++jcrd) {
                    S_cart[icrd][jcrd] = it.rot[3 * icrd + jcrd];
                    S_frac[icrd][jcrd] = it.rot_real[3 * icrd + jcrd];
                    S_recip[icrd][jcrd] = it.rot_reciprocal[3 * icrd + jcrd];
                }
            }

            invmat3(S_frac_inv, S_frac);
            rotvec(Sk, k, S_recip);

            for (auto i = 0; i < 3; ++i) Sk[i] = Sk[i] - nint(Sk[i]);

            const auto knum_sym = kpoint->get_knum(Sk, kmesh_interpolate);
            if (knum_sym == -1)
                error->exit("setup_transform_symmetry",
                            "kpoint not found");

            if (knum_sym == knum) small_group_at_k[ik].push_back(isym);
            if (knum_sym == knum_minus) symop_minus_at_k[ik].push_back(isym);

            if (!flag[knum_sym]) {
                kpoint_map_symmetry[knum_sym].symmetry_op = isym;
                kpoint_map_symmetry[knum_sym].knum_irred_orig = ik;
                kpoint_map_symmetry[knum_sym].knum_orig = knum;
                flag[knum_sym] = true;
            }

            for (is = 0; is < ns; ++is) {
                for (js = 0; js < ns; ++js) {
                    gamma_tmp[is][js] = std::complex<double>(0.0, 0.0);
                }
            }

            for (unsigned int jat = 0; jat < natmin; ++jat) {
                const auto iat = it.mapping[jat];

                // Fractional coordinates of x1 and x2
                for (icrd = 0; icrd < 3; ++icrd) {
                    x1[icrd] = system->xr_p[system->map_p2s[iat][0]][icrd];
                    x2[icrd] = system->xr_p[system->map_p2s[jat][0]][icrd];
                }

                rotvec(xtmp, x1, S_frac_inv);
                for (icrd = 0; icrd < 3; ++icrd) {
                    xtmp[icrd] = xtmp[icrd] - x2[icrd];
                }

                auto phase = 2.0 * pi * (k[0] * xtmp[0] + k[1] * xtmp[1] + k[2] * xtmp[2]);

                for (icrd = 0; icrd < 3; ++icrd) {
                    for (jcrd = 0; jcrd < 3; ++jcrd) {
                        gamma_tmp[3 * iat + icrd][3 * jat + jcrd]
                            = S_cart[icrd][jcrd] * std::exp(im * phase);
                    }
                }
            }

            for (is = 0; is < ns; ++is) {
                for (js = 0; js < ns; ++js) {
                    mat_transform_sym[ik][isym][is][js] = gamma_tmp[is][js];
                }
            }

            ++isym;
        }
    }

    memory->deallocate(gamma_tmp);
    memory->deallocate(flag);
}


void Scph::symmetrize_dynamical_matrix(const unsigned int ik,
                                       Eigen::MatrixXcd &dymat) const
{
    // Symmetrize the dynamical matrix of given index ik.
    using namespace Eigen;
    unsigned int i, isym;
    unsigned int is, js;
    const auto ns = dynamical->neval;
    MatrixXcd dymat_sym = MatrixXcd::Zero(ns, ns);
    MatrixXcd dymat_tmp(ns, ns), gamma(ns, ns);

    const auto nsym_small = small_group_at_k[ik].size();
    const auto nsym_minus = symop_minus_at_k[ik].size();

    for (i = 0; i < nsym_minus; ++i) {
        isym = symop_minus_at_k[ik][i];

        for (is = 0; is < ns; ++is) {
            for (js = 0; js < ns; ++js) {
                gamma(is, js) = mat_transform_sym[ik][isym][is][js];
            }
        }

        dymat_tmp = gamma * dymat * gamma.transpose().conjugate();
        dymat_sym += dymat_tmp.conjugate();
    }

    for (i = 0; i < nsym_small; ++i) {
        isym = small_group_at_k[ik][i];

        for (is = 0; is < ns; ++is) {
            for (js = 0; js < ns; ++js) {
                gamma(is, js) = mat_transform_sym[ik][isym][is][js];
            }
        }

        dymat_tmp = gamma * dymat * gamma.transpose().conjugate();
        dymat_sym += dymat_tmp;
    }


    dymat = dymat_sym / static_cast<double>(nsym_small + nsym_minus);
}

void Scph::replicate_dymat_for_all_kpoints(std::complex<double> ***dymat_inout) const
{
    using namespace Eigen;
    unsigned int i;
    unsigned int is, js;
    const auto ns = dynamical->neval;
    MatrixXcd dymat_tmp(ns, ns), gamma(ns, ns), dymat(ns, ns);

    std::complex<double> ***dymat_all;

    memory->allocate(dymat_all, ns, ns, nk_interpolate);

    for (i = 0; i < nk_interpolate; ++i) {

        const auto ik_irred = kpoint_map_symmetry[i].knum_irred_orig;
        const auto ik_orig = kpoint_map_symmetry[i].knum_orig;
        const auto isym = kpoint_map_symmetry[i].symmetry_op;

        for (is = 0; is < ns; ++is) {
            for (js = 0; js < ns; ++js) {
                gamma(is, js) = mat_transform_sym[ik_irred][isym][is][js];
                dymat(is, js) = dymat_inout[is][js][ik_orig];
            }
        }
        dymat_tmp = gamma * dymat * gamma.transpose().conjugate();

        for (is = 0; is < ns; ++is) {
            for (js = 0; js < ns; ++js) {
                dymat_all[is][js][i] = dymat_tmp(is, js);
            }
        }
    }

    for (is = 0; is < ns; ++is) {
        for (js = 0; js < ns; ++js) {
            for (i = 0; i < nk_interpolate; ++i) {
                dymat_inout[is][js][i] = dymat_all[is][js][i];
            }
        }
    }
    memory->deallocate(dymat_all);
}

void Scph::setup_eigvecs()
{
    const auto ns = dynamical->neval;

    if (mympi->my_rank == 0) {
        std::cout << std::endl
            << " Diagonalizing dynamical matrices for all k points ... ";
    }

    memory->allocate(evec_harmonic, nk_scph, ns, ns);
    memory->allocate(omega2_harmonic, nk_scph, ns);

    // Calculate phonon eigenvalues and eigenvectors for all k-points for scph

#pragma omp parallel for
    for (int ik = 0; ik < nk_scph; ++ik) {

        dynamical->eval_k(xk_scph[ik], kvec_na_scph[ik],
                          fcs_phonon->fc2_ext, omega2_harmonic[ik],
                          evec_harmonic[ik], true);
    }

    if (mympi->my_rank == 0) {
        std::cout << "done !" << std::endl;
    }
}

void Scph::setup_pp_interaction()
{
    // Prepare information for calculating ph-ph interaction coefficients.

    unsigned int i, j;
    double *invsqrt_mass_p;

    if (mympi->my_rank == 0) {
        if (relax_coordinate) {
            std::cout << " Preparing for calculating V3 & V4  ...";
        } else {
            std::cout << " Preparing for calculating V4  ...";
        }
    }

    if (anharmonic_core->quartic_mode != 1) {
        error->exit("setup_pp_interaction",
                    "quartic_mode should be 1 for SCPH");
    }

    memory->allocate(invsqrt_mass_p, system->natmin);

    for (i = 0; i < system->natmin; ++i) {
        invsqrt_mass_p[i] = std::sqrt(1.0 / system->mass[system->map_p2s[i][0]]);
    }

    // Setup for V3 if relax_coordinate = True.
    if (relax_coordinate) {

        std::sort(fcs_phonon->force_constant_with_cell[1].begin(),
                  fcs_phonon->force_constant_with_cell[1].end());

        anharmonic_core->prepare_group_of_force_constants(fcs_phonon->force_constant_with_cell[1], 3,
                                                          ngroup, fcs_group);

        memory->allocate(vec_for_v3, fcs_phonon->force_constant_with_cell[1].size(), 2, 3);
        memory->allocate(invmass_for_v3, fcs_phonon->force_constant_with_cell[1].size());
        memory->allocate(evec_index3, fcs_phonon->force_constant_with_cell[1].size(), 3);


        j = 0;
        for (const auto &it : fcs_phonon->force_constant_with_cell[1]) {
            invmass_for_v3[j]
                = invsqrt_mass_p[it.pairs[0].index / 3]
                * invsqrt_mass_p[it.pairs[1].index / 3]
                * invsqrt_mass_p[it.pairs[2].index / 3];

            ++j;
        }
        anharmonic_core->prepare_relative_vector(fcs_phonon->force_constant_with_cell[1],
                                                 3, vec_for_v3);

        for (i = 0; i < fcs_phonon->force_constant_with_cell[1].size(); ++i) {
            for (j = 0; j < 3; ++j) {
                evec_index3[i][j]
                    = fcs_phonon->force_constant_with_cell[1][i].pairs[j].index;
            }
        }

    }

    std::sort(fcs_phonon->force_constant_with_cell[2].begin(),
              fcs_phonon->force_constant_with_cell[2].end());

    anharmonic_core->prepare_group_of_force_constants(fcs_phonon->force_constant_with_cell[2], 4,
                                                      ngroup2, fcs_group2);

    memory->allocate(vec_for_v4, fcs_phonon->force_constant_with_cell[2].size(), 3, 3);
    memory->allocate(invmass_for_v4, fcs_phonon->force_constant_with_cell[2].size());
    memory->allocate(evec_index4, fcs_phonon->force_constant_with_cell[2].size(), 4);


    j = 0;
    for (const auto &it : fcs_phonon->force_constant_with_cell[2]) {
        invmass_for_v4[j]
            = invsqrt_mass_p[it.pairs[0].index / 3]
            * invsqrt_mass_p[it.pairs[1].index / 3]
            * invsqrt_mass_p[it.pairs[2].index / 3]
            * invsqrt_mass_p[it.pairs[3].index / 3];

        ++j;
    }
    anharmonic_core->prepare_relative_vector(fcs_phonon->force_constant_with_cell[2],
                                             4, vec_for_v4);

    for (i = 0; i < fcs_phonon->force_constant_with_cell[2].size(); ++i) {
        for (j = 0; j < 4; ++j) {
            evec_index4[i][j]
                = fcs_phonon->force_constant_with_cell[2][i].pairs[j].index;
        }
    }


    memory->deallocate(invsqrt_mass_p);

    nk_grid[0] = kmesh_scph[0];
    nk_grid[1] = kmesh_scph[1];
    nk_grid[2] = kmesh_scph[2];

    for (i = 0; i < 3; ++i) dnk[i] = static_cast<double>(nk_grid[i]);

    if (nk_grid[0] == nk_grid[1] && nk_grid[1] == nk_grid[2]) {
        nk_represent = nk_grid[0];
        tune_type = 0;

    } else if (nk_grid[0] == nk_grid[1] && nk_grid[2] == 1) {
        nk_represent = nk_grid[0];
        tune_type = 0;

    } else if (nk_grid[1] == nk_grid[2] && nk_grid[0] == 1) {
        nk_represent = nk_grid[1];
        tune_type = 0;

    } else if (nk_grid[2] == nk_grid[0] && nk_grid[1] == 1) {
        nk_represent = nk_grid[2];
        tune_type = 0;

    } else if (nk_grid[0] == 1 && nk_grid[1] == 1) {
        nk_represent = nk_grid[2];
        tune_type = 0;

    } else if (nk_grid[1] == 1 && nk_grid[2] == 1) {
        nk_represent = nk_grid[0];
        tune_type = 0;

    } else if (nk_grid[2] == 1 && nk_grid[0] == 1) {
        nk_represent = nk_grid[1];
        tune_type = 0;

    } else {
        tune_type = 1;
    }

    int ii;

    if (tune_type == 0) {

        memory->allocate(exp_phase, 2 * nk_represent - 1);
        for (ii = 0; ii < 2 * nk_represent - 1; ++ii) {
            const auto phase = 2.0 * pi * static_cast<double>(ii - nk_represent + 1)
                / static_cast<double>(nk_represent);
            exp_phase[ii] = std::exp(im * phase);
        }

    } else if (tune_type == 1) {
        double phase[3];

        tune_type = 1;
        memory->allocate(exp_phase3, 2 * nk_grid[0] - 1, 2 * nk_grid[1] - 1, 2 * nk_grid[2] - 1);

        for (ii = 0; ii < 2 * nk_grid[0] - 1; ++ii) {
            phase[0] = 2.0 * pi * static_cast<double>(ii - nk_grid[0] + 1) / dnk[0];
            for (int jj = 0; jj < 2 * nk_grid[1] - 1; ++jj) {
                phase[1] = 2.0 * pi * static_cast<double>(jj - nk_grid[1] + 1) / dnk[1];
                for (int kk = 0; kk < 2 * nk_grid[2] - 1; ++kk) {
                    phase[2] = 2.0 * pi * static_cast<double>(kk - nk_grid[2] + 1) / dnk[2];
                    exp_phase3[ii][jj][kk] = std::exp(im * (phase[0] + phase[1] + phase[2]));
                }
            }
        }
    }

    if (mympi->my_rank == 0) {
        std::cout << " done!" << std::endl;
    }
}

void Scph::setup_transform_ifc()
{
    // Compute mindist_list_scph necessary to calculate dynamical matrix
    // from the real-space force constants

    int i, j;
    int ix, iy, iz;
    const auto nk = nk_interpolate;
    const auto nat = system->natmin;
    unsigned int iat;

    int **shift_cell, **shift_cell_super;
    double **xf_p;
    double ****x_all;

    const int nkx = kmesh_interpolate[0];
    const int nky = kmesh_interpolate[1];
    const int nkz = kmesh_interpolate[2];

    const auto ncell = nk;
    const auto ncell_s = 27;

    memory->allocate(shift_cell, ncell, 3);
    memory->allocate(shift_cell_super, ncell_s, 3);
    memory->allocate(xf_p, nat, 3);
    memory->allocate(x_all, ncell_s, ncell, nat, 3);

    unsigned int icell = 0;
    for (ix = 0; ix < nkx; ++ix) {
        for (iy = 0; iy < nky; ++iy) {
            for (iz = 0; iz < nkz; ++iz) {

                shift_cell[icell][0] = ix;
                shift_cell[icell][1] = iy;
                shift_cell[icell][2] = iz;

                ++icell;
            }
        }
    }

    for (i = 0; i < 3; ++i) shift_cell_super[0][i] = 0;
    icell = 1;
    for (ix = -1; ix <= 1; ++ix) {
        for (iy = -1; iy <= 1; ++iy) {
            for (iz = -1; iz <= 1; ++iz) {
                if (ix == 0 && iy == 0 && iz == 0) continue;

                shift_cell_super[icell][0] = ix;
                shift_cell_super[icell][1] = iy;
                shift_cell_super[icell][2] = iz;

                ++icell;
            }
        }
    }

    for (i = 0; i < nat; ++i) {
        rotvec(xf_p[i], system->xr_s[system->map_p2s[i][0]], system->lavec_s);
        rotvec(xf_p[i], xf_p[i], system->rlavec_p);
        for (j = 0; j < 3; ++j) xf_p[i][j] /= 2.0 * pi;
    }

    for (i = 0; i < ncell_s; ++i) {
        for (j = 0; j < ncell; ++j) {
            for (iat = 0; iat < nat; ++iat) {
                x_all[i][j][iat][0] = xf_p[iat][0] + static_cast<double>(shift_cell[j][0])
                    + static_cast<double>(nkx * shift_cell_super[i][0]);
                x_all[i][j][iat][1] = xf_p[iat][1] + static_cast<double>(shift_cell[j][1])
                    + static_cast<double>(nky * shift_cell_super[i][1]);
                x_all[i][j][iat][2] = xf_p[iat][2] + static_cast<double>(shift_cell[j][2])
                    + static_cast<double>(nkz * shift_cell_super[i][2]);

                rotvec(x_all[i][j][iat], x_all[i][j][iat], system->lavec_p);
            }
        }
    }

    double dist;
    std::vector<DistList> dist_tmp;
    ShiftCell shift_tmp{};
    std::vector<int> vec_tmp;

    memory->allocate(mindist_list_scph, nat, nat, ncell);

    for (iat = 0; iat < nat; ++iat) {
        for (unsigned int jat = 0; jat < nat; ++jat) {
            for (icell = 0; icell < ncell; ++icell) {

                dist_tmp.clear();
                for (i = 0; i < ncell_s; ++i) {
                    dist = distance(x_all[0][0][iat], x_all[i][icell][jat]);
                    dist_tmp.emplace_back(i, dist);
                }
                std::sort(dist_tmp.begin(), dist_tmp.end());

                const auto dist_min = dist_tmp[0].dist;
                mindist_list_scph[iat][jat][icell].dist = dist_min;

                for (i = 0; i < ncell_s; ++i) {
                    dist = dist_tmp[i].dist;

                    if (std::abs(dist_min - dist) < eps8) {

                        shift_tmp.sx = shift_cell[icell][0]
                            + nkx * shift_cell_super[dist_tmp[i].cell_s][0];
                        shift_tmp.sy = shift_cell[icell][1]
                            + nky * shift_cell_super[dist_tmp[i].cell_s][1];
                        shift_tmp.sz = shift_cell[icell][2]
                            + nkz * shift_cell_super[dist_tmp[i].cell_s][2];

                        mindist_list_scph[iat][jat][icell].shift.push_back(shift_tmp);
                    }
                }

            }
        }
    }

    memory->deallocate(shift_cell);
    memory->deallocate(shift_cell_super);
    memory->deallocate(xf_p);
    memory->deallocate(x_all);
}


void Scph::exec_interpolation(const unsigned int kmesh_orig[3],
                              std::complex<double> ***dymat_r,
                              const unsigned int nk_dense,
                              double **xk_dense,
                              double **kvec_dense,
                              double **eval_out,
                              std::complex<double> ***evec_out)
{
    unsigned int i, j, is;
    const auto ns = dynamical->neval;
    const auto nk1 = kmesh_orig[0];
    const auto nk2 = kmesh_orig[1];
    const auto nk3 = kmesh_orig[2];

    double *eval_real;
    std::complex<double> **mat_tmp;
    std::complex<double> **mat_harmonic, **mat_harmonic_na;
    std::vector<double> eval_vec(ns);

    memory->allocate(mat_tmp, ns, ns);
    memory->allocate(eval_real, ns);
    memory->allocate(mat_harmonic, ns, ns);

    if (dynamical->nonanalytic) {
        memory->allocate(mat_harmonic_na, ns, ns);
    }

    for (int ik = 0; ik < nk_dense; ++ik) {

        if (dynamical->nonanalytic == 3) {
            dynamical->calc_analytic_k(xk_dense[ik],
                                       ewald->fc2_without_dipole,
                                       mat_harmonic);
        } else {
            dynamical->calc_analytic_k(xk_dense[ik],
                                       fcs_phonon->fc2_ext,
                                       mat_harmonic);
        }

        r2q(xk_dense[ik], nk1, nk2, nk3, ns, dymat_r, mat_tmp);

        for (i = 0; i < ns; ++i) {
            for (j = 0; j < ns; ++j) {
                mat_tmp[i][j] += mat_harmonic[i][j];
            }
        }

        if (dynamical->nonanalytic) {

            if (dynamical->nonanalytic == 1) {
                dynamical->calc_nonanalytic_k(xk_dense[ik],
                                              kvec_dense[ik],
                                              mat_harmonic_na);
            } else if (dynamical->nonanalytic == 2) {
                dynamical->calc_nonanalytic_k2(xk_dense[ik],
                                               kvec_dense[ik],
                                               mat_harmonic_na);

            } else if (dynamical->nonanalytic == 3) {
                ewald->add_longrange_matrix(xk_dense[ik],
                                            kvec_dense[ik],
                                            mat_harmonic_na);
            }

            for (i = 0; i < ns; ++i) {
                for (j = 0; j < ns; ++j) {
                    mat_tmp[i][j] += mat_harmonic_na[i][j];
                }
            }
        }

        diagonalize_interpolated_matrix(mat_tmp, eval_real, evec_out[ik], true);

        for (is = 0; is < ns; ++is) {
            const auto eval_tmp = eval_real[is];

            if (eval_tmp < 0.0) {
                eval_vec[is] = -std::sqrt(-eval_tmp);
            } else {
                eval_vec[is] = std::sqrt(eval_tmp);
            }
        }

        for (is = 0; is < ns; ++is) eval_out[ik][is] = eval_vec[is];

    }

    memory->deallocate(eval_real);
    memory->deallocate(mat_tmp);
    memory->deallocate(mat_harmonic);

    if (dynamical->nonanalytic) {
        memory->deallocate(mat_harmonic_na);
    }
}


void Scph::r2q(const double *xk_in,
               const unsigned int nx,
               const unsigned int ny,
               const unsigned int nz,
               const unsigned int ns,
               std::complex<double> ***dymat_r_in,
               std::complex<double> **dymat_k_out) const
{
    std::complex<double> im(0.0, 1.0);

    const auto ncell = nx * ny * nz;

    for (unsigned int i = 0; i < ns; ++i) {

        const auto iat = i / 3;

        for (unsigned int j = 0; j < ns; ++j) {

            const auto jat = j / 3;

            dymat_k_out[i][j] = std::complex<double>(0.0, 0.0);

            for (unsigned int icell = 0; icell < ncell; ++icell) {

                auto exp_phase = std::complex<double>(0.0, 0.0);

                // This operation is necessary for the Hermiticity of the dynamical matrix.
                for (const auto &it : mindist_list_scph[iat][jat][icell].shift) {

                    auto phase = 2.0 * pi
                        * (static_cast<double>(it.sx) * xk_in[0]
                            + static_cast<double>(it.sy) * xk_in[1]
                            + static_cast<double>(it.sz) * xk_in[2]);

                    exp_phase += std::exp(im * phase);
                }
                exp_phase /= static_cast<double>(mindist_list_scph[iat][jat][icell].shift.size());

                dymat_k_out[i][j] += dymat_r_in[i][j][icell] * exp_phase;
            }

        }
    }
}


void Scph::diagonalize_interpolated_matrix(std::complex<double> **mat_in,
                                           double *eval_out,
                                           std::complex<double> **evec_out,
                                           const bool require_evec) const
{
    unsigned int i, j;
    char JOBZ;
    int INFO;
    double *RWORK;
    std::complex<double> *amat;
    std::complex<double> *WORK;

    int ns = dynamical->neval;

    int LWORK = (2 * ns - 1) * 10;
    memory->allocate(RWORK, 3 * ns - 2);
    memory->allocate(WORK, LWORK);


    if (require_evec) {
        JOBZ = 'V';
    } else {
        JOBZ = 'N';
    }

    char UPLO = 'U';

    memory->allocate(amat, ns * ns);

    unsigned int k = 0;
    for (j = 0; j < ns; ++j) {
        for (i = 0; i < ns; ++i) {
            amat[k++] = mat_in[i][j];
        }
    }

    zheev_(&JOBZ, &UPLO, &ns, amat, &ns, eval_out, WORK, &LWORK, RWORK, &INFO);

    k = 0;

    if (require_evec) {
        // Here we transpose the matrix evec_out so that 
        // evec_out[i] becomes phonon eigenvector of i-th mode.
        for (j = 0; j < ns; ++j) {
            for (i = 0; i < ns; ++i) {
                evec_out[j][i] = amat[k++];
            }
        }
    }

    memory->deallocate(amat);
    memory->deallocate(WORK);
    memory->deallocate(RWORK);
}

void Scph::find_degeneracy(std::vector<int> *degeneracy_out,
                           const unsigned int nk_irred,
                           const std::vector<std::vector<KpointList>> &kp_info,
                           double **eval) const
{
    // eval is omega^2 in atomic unit

    const auto ns = dynamical->neval;
    const auto tol_omega = 1.0e-10;

    for (unsigned int ik = 0; ik < nk_irred; ++ik) {
        const auto knum = kp_info[ik][0].knum;

        degeneracy_out[ik].clear();

        auto omega_prev = eval[knum][0];
        auto ideg = 1;

        for (unsigned int is = 1; is < ns; ++is) {
            const auto omega_now = eval[knum][is];

            if (std::abs(omega_now - omega_prev) < tol_omega) {
                ++ideg;
            } else {
                degeneracy_out[ik].push_back(ideg);
                ideg = 1;
                omega_prev = omega_now;
            }

        }
        degeneracy_out[ik].push_back(ideg);
    }
}

void Scph::calc_new_dymat_with_evec(std::complex<double> ***dymat_out,
                                    double **omega2_in,
                                    std::complex<double> ***evec_in)
{
    std::complex<double> *polarization_matrix, *mat_tmp;
    std::complex<double> *eigval_matrix, *dmat;
    std::complex<double> *beta;
    std::complex<double> ***dymat_q, **dymat_harmonic;
    std::complex<double> im(0.0, 1.0);

    unsigned int ik, is, js;
    int ns = dynamical->neval;

    const unsigned int ns2 = ns * ns;

    auto alpha = std::complex<double>(1.0, 0.0);

    char TRANSA[] = "N";
    char TRANSB[] = "C";

    memory->allocate(polarization_matrix, ns2);
    memory->allocate(mat_tmp, ns2);
    memory->allocate(eigval_matrix, ns2);
    memory->allocate(beta, ns);
    memory->allocate(dmat, ns2);
    memory->allocate(dymat_q, ns, ns, nk_interpolate);
    memory->allocate(dymat_harmonic, ns, ns);

    for (is = 0; is < ns; ++is) beta[is] = std::complex<double>(0.0, 0.0);

    for (ik = 0; ik < nk_interpolate; ++ik) {

        const auto knum = kmap_interpolate_to_scph[ik];

        // create eigval matrix

        for (is = 0; is < ns2; ++is) eigval_matrix[is] = std::complex<double>(0.0, 0.0);

        unsigned int m = 0;
        for (is = 0; is < ns; ++is) {
            for (js = 0; js < ns; ++js) {
                if (is == js) {
                    eigval_matrix[m] = omega2_in[knum][is];
                }
                ++m;
            }
        }

        // create polarization matrix

        m = 0;

        for (is = 0; is < ns; ++is) {
            for (js = 0; js < ns; ++js) {
                polarization_matrix[m++] = evec_in[knum][is][js];
            }
        }

        zgemm_(TRANSA, TRANSB, &ns, &ns, &ns, &alpha,
               eigval_matrix, &ns, polarization_matrix, &ns, beta, mat_tmp, &ns);
        zgemm_(TRANSA, TRANSA, &ns, &ns, &ns, &alpha,
               polarization_matrix, &ns, mat_tmp, &ns, beta, dmat, &ns);

        m = 0;

        for (js = 0; js < ns; ++js) {
            for (is = 0; is < ns; ++is) {
                dymat_q[is][js][ik] = dmat[m];
                ++m;
            }
        }


        // Subtract harmonic contribution
        dynamical->calc_analytic_k(xk_interpolate[ik],
                                   fcs_phonon->fc2_ext,
                                   dymat_harmonic);

        for (is = 0; is < ns; ++is) {
            for (js = 0; js < ns; ++js) {
                dymat_q[is][js][ik] -= dymat_harmonic[is][js];
            }
        }
    }

    memory->deallocate(polarization_matrix);
    memory->deallocate(mat_tmp);
    memory->deallocate(eigval_matrix);
    memory->deallocate(beta);
    memory->deallocate(dmat);
    memory->deallocate(dymat_harmonic);

    const auto nk1 = kmesh_interpolate[0];
    const auto nk2 = kmesh_interpolate[1];
    const auto nk3 = kmesh_interpolate[2];

    std::vector<std::vector<double>> xk_dup;

    int icell = 0;

    for (int ix = 0; ix < nk1; ++ix) {
        for (int iy = 0; iy < nk2; ++iy) {
            for (int iz = 0; iz < nk3; ++iz) {

                for (is = 0; is < ns; ++is) {
                    for (js = 0; js < ns; ++js) {
                        dymat_out[is][js][icell] = std::complex<double>(0.0, 0.0);
                    }
                }

                for (ik = 0; ik < nk_interpolate; ++ik) {

                    duplicate_xk_boundary(xk_interpolate[ik], xk_dup);

                    auto cexp_phase = std::complex<double>(0.0, 0.0);

                    for (const auto &i : xk_dup) {

                        auto phase = 2.0 * pi * (i[0] * static_cast<double>(ix)
                            + i[1] * static_cast<double>(iy)
                            + i[2] * static_cast<double>(iz));
                        cexp_phase += std::exp(-im * phase);

                    }
                    cexp_phase /= static_cast<double>(xk_dup.size());

                    for (is = 0; is < ns; ++is) {
                        for (js = 0; js < ns; ++js) {
                            dymat_out[is][js][icell] += dymat_q[is][js][ik] * cexp_phase;
                        }
                    }


                }
                for (is = 0; is < ns; ++is) {
                    for (js = 0; js < ns; ++js) {
                        dymat_out[is][js][icell] /= static_cast<double>(nk_interpolate);
                    }
                }

                ++icell;
            }
        }
    }


    memory->deallocate(dymat_q);
}


void Scph::compute_anharmonic_frequency(std::complex<double> ***v4_array_all,
                                        double **omega2_out,
                                        std::complex<double> ***evec_anharm_scph,
                                        const double temp,
                                        std::vector<int> *degeneracy_info,
                                        bool &flag_converged,
                                        std::complex<double> ***cmat_convert,
                                        const bool offdiag)
{
    // This is the main function of the SCPH equation.
    // The detailed algorithm can be found in PRB 92, 054301 (2015).
    // Eigen3 library is used for the compact notation of matrix-matrix products.

    using namespace Eigen;

    int ik, jk;
    unsigned int i;
    unsigned int is, js, ks;
    unsigned int kk;
    const auto nk = nk_scph;
    const auto ns = dynamical->neval;
    unsigned int knum, knum_interpolate;
    const auto nk_irred_interpolate = kp_irred_interpolate.size();
    const auto nk1 = kmesh_interpolate[0];
    const auto nk2 = kmesh_interpolate[1];
    const auto nk3 = kmesh_interpolate[2];
    int iloop;

    MatrixXd omega_now(nk, ns), omega_old(nk, ns);
    MatrixXd omega2_HA(nk, ns);
    MatrixXcd mat_tmp(ns, ns), evec_tmp(ns, ns);

    VectorXd eval_tmp(ns);
    MatrixXcd Dymat(ns, ns);
    MatrixXcd Fmat(ns, ns);
    MatrixXcd Qmat = MatrixXcd::Zero(ns, ns);
    MatrixXcd Cmat(ns, ns), Dmat(ns, ns);

    double diff;
    double conv_tol = tolerance_scph;
    double alpha = mixalpha;

    double **eval_interpolate;
    double re_tmp, im_tmp;
    bool has_negative;

    std::complex<double> ctmp;
    std::complex<double> ***mat_omega2_harmonic;
    std::complex<double> ***evec_initial;
    std::complex<double> ***dmat_convert;
    std::complex<double> ***dmat_convert_old;
    std::complex<double> ***evec_new;
    std::complex<double> ***dymat_new, ***dymat_harmonic;
    std::complex<double> ***dymat_q;
    std::complex<double> ***Fmat0;

    const auto complex_one = std::complex<double>(1.0, 0.0);
    const auto complex_zero = std::complex<double>(0.0, 0.0);

    SelfAdjointEigenSolver<MatrixXcd> saes;

    memory->allocate(mat_omega2_harmonic, nk_interpolate, ns, ns);
    memory->allocate(eval_interpolate, nk, ns);
    memory->allocate(evec_initial, nk, ns, ns);
    memory->allocate(evec_new, nk, ns, ns);
    memory->allocate(dmat_convert, nk, ns, ns);
    memory->allocate(dmat_convert_old, nk, ns, ns);
    memory->allocate(dymat_new, ns, ns, nk_interpolate);
    memory->allocate(dymat_q, ns, ns, nk_interpolate);
    memory->allocate(dymat_harmonic, nk_interpolate, ns, ns);
    memory->allocate(Fmat0, nk_irred_interpolate, ns, ns);

    const auto T_in = temp;

    std::cout << " Temperature = " << T_in << " K" << std::endl;

    // Set initial values

    for (ik = 0; ik < nk; ++ik) {
        for (is = 0; is < ns; ++is) {

            if (flag_converged) {
                if (omega2_out[ik][is] < 0.0 && std::abs(omega2_out[ik][is]) > 1.0e-16) {
                    std::cout << "Warning : Large negative frequency detected" << std::endl;
                }

                if (omega2_out[ik][is] < 0.0) {
                    omega_now(ik, is) = std::sqrt(-omega2_out[ik][is]);
                } else {
                    omega_now(ik, is) = std::sqrt(omega2_out[ik][is]);
                }
            } else {
                if (omega2_harmonic[ik][is] < 0.0) {
                    omega_now(ik, is) = std::sqrt(-omega2_harmonic[ik][is]);
                } else {
                    omega_now(ik, is) = std::sqrt(omega2_harmonic[ik][is]);
                }
            }

            omega2_HA(ik, is) = omega2_harmonic[ik][is];

            for (js = 0; js < ns; ++js) {
                evec_initial[ik][is][js] = evec_harmonic[ik][is][js];

                if (!flag_converged) {
                    // Initialize Cmat with identity matrix
                    if (is == js) {
                        cmat_convert[ik][is][js] = complex_one;
                    } else {
                        cmat_convert[ik][is][js] = complex_zero;
                    }
                }

            }
        }
    }

    // Set initial harmonic dymat and eigenvalues

    for (ik = 0; ik < nk_interpolate; ++ik) {
        knum = kmap_interpolate_to_scph[ik];

        for (is = 0; is < ns; ++is) {
            for (js = 0; js < ns; ++js) {
                mat_omega2_harmonic[ik][is][js] = complex_zero;
            }
            mat_omega2_harmonic[ik][is][is] = std::complex<double>(omega2_HA(knum, is), 0.0);
            eval_interpolate[knum][is] = omega2_HA(knum, is);
        }

        dynamical->calc_analytic_k(xk_interpolate[ik], fcs_phonon->fc2_ext,
                                   dymat_harmonic[ik]);
    }

    for (ik = 0; ik < nk_irred_interpolate; ++ik) {

        knum_interpolate = kp_irred_interpolate[ik][0].knum;

        // Fmat harmonic
        for (is = 0; is < ns; ++is) {
            for (js = 0; js < ns; ++js) {
                if (is == js) {
                    Fmat0[ik][is][js] = mat_omega2_harmonic[knum_interpolate][is][is];
                } else {
                    Fmat0[ik][is][js] = complex_zero;
                }
            }
        }
    }

    int icount = 0;

    // Main loop
    for (iloop = 0; iloop < maxiter; ++iloop) {

        for (ik = 0; ik < nk; ++ik) {
            for (is = 0; is < ns; ++is) {
                auto omega1 = omega_now(ik, is);
                if (std::abs(omega1) < eps8) {
                    Qmat(is, is) = complex_zero;
                } else {
                    // Note that the missing factor 2 in the denominator of Qmat is 
                    // already considered in the v4_array_all.
                    if (thermodynamics->classical) {
                        Qmat(is, is) = std::complex<double>(2.0 * T_in * thermodynamics->T_to_Ryd / (omega1 * omega1),
                                                            0.0);
                    } else {
                        auto n1 = thermodynamics->fB(omega1, T_in);
                        Qmat(is, is) = std::complex<double>((2.0 * n1 + 1.0) / omega1, 0.0);
                    }
                }
            }

            for (is = 0; is < ns; ++is) {
                for (js = 0; js < ns; ++js) {
                    Cmat(is, js) = cmat_convert[ik][is][js];
                }
            }

            Dmat = Cmat * Qmat * Cmat.adjoint();

            for (is = 0; is < ns; ++is) {
                for (js = 0; js < ns; ++js) {
                    dmat_convert[ik][is][js] = Dmat(is, js);
                }
            }
        }

        // Mixing dmat
        if (iloop > 0) {
#pragma omp parallel for private(is, js)
            for (ik = 0; ik < nk; ++ik) {
                for (is = 0; is < ns; ++is) {
                    for (js = 0; js < ns; ++js) {
                        dmat_convert[ik][is][js] = alpha * dmat_convert[ik][is][js]
                            + (1.0 - alpha) * dmat_convert_old[ik][is][js];
                    }
                }
            }
        }

        for (ik = 0; ik < nk_irred_interpolate; ++ik) {

            knum_interpolate = kp_irred_interpolate[ik][0].knum;
            knum = kmap_interpolate_to_scph[knum_interpolate];

            // Fmat harmonic

            for (is = 0; is < ns; ++is) {
                for (js = 0; js < ns; ++js) {
                    Fmat(is, js) = Fmat0[ik][is][js];
                }
            }

            // Anharmonic correction to Fmat

            if (!offdiag) {
                for (is = 0; is < ns; ++is) {
                    i = (ns + 1) * is;

                    re_tmp = 0.0;
                    im_tmp = 0.0;

#pragma omp parallel for private(jk,kk,ks), reduction(+:re_tmp, im_tmp)
                    for (jk = 0; jk < nk; ++jk) {

                        kk = nk * ik + jk;

                        for (ks = 0; ks < ns; ++ks) {
                            ctmp = v4_array_all[kk][i][(ns + 1) * ks]
                                * dmat_convert[jk][ks][ks];
                            re_tmp += ctmp.real();
                            im_tmp += ctmp.imag();
                        }
                    }
                    Fmat(is, is) += std::complex<double>(re_tmp, im_tmp);
                }
            } else {

                // Anharmonic correction to Fmat

                for (is = 0; is < ns; ++is) {
                    for (js = 0; js <= is; ++js) {

                        i = ns * is + js;

                        re_tmp = 0.0;
                        im_tmp = 0.0;

#pragma omp parallel for private(jk,kk,ks), reduction(+:re_tmp, im_tmp)
                        for (jk = 0; jk < nk; ++jk) {

                            kk = nk * ik + jk;

                            for (ks = 0; ks < ns; ++ks) {
                                for (unsigned int ls = 0; ls < ns; ++ls) {
                                    ctmp = v4_array_all[kk][i][ns * ks + ls]
                                        * dmat_convert[jk][ks][ls];
                                    re_tmp += ctmp.real();
                                    im_tmp += ctmp.imag();
                                }
                            }
                        }
                        Fmat(is, js) += std::complex<double>(re_tmp, im_tmp);
                    }
                }

            }

            saes.compute(Fmat);
            eval_tmp = saes.eigenvalues();

            for (is = 0; is < ns; ++is) {

                double omega2_tmp = eval_tmp(is);

                if (omega2_tmp < 0.0 && std::abs(omega2_tmp) > 1.0e-16) {

                    std::cout << " Detect imaginary : ";
                    std::cout << "  knum = " << knum + 1 << " is = " << is + 1 << std::endl;
                    for (int j = 0; j < 3; ++j) {
                        std::cout << "  xk = " << std::setw(15) << xk_scph[knum][j];
                    }
                    std::cout << std::endl;

                    if (v4_array_all[nk * ik + knum][(ns + 1) * is][(ns + 1) * is].real() > 0.0) {
                        std::cout << "  onsite V4 is positive" << std::endl;
                        std::cout << std::endl;
                        if (flag_converged) {
                            ++icount;
                            eval_tmp(is) = omega2_out[knum][is] * std::pow(0.99, icount);
                        } else {
                            ++icount;
                            eval_tmp(is) = -eval_tmp(is) * std::pow(0.99, icount);
                        }
                    } else {
                        std::cout << "  onsite V4 is negative" << std::endl;
                        std::cout << std::endl;
                        eval_tmp(is) = std::abs(omega2_tmp);
                    }
                }
            }

            for (is = 0; is < ns; ++is) {
                for (js = 0; js < ns; ++js) {
                    evec_tmp(is, js) = evec_initial[knum][is][js];

                    if (is == js) {
                        Dymat(is, js) = std::complex<double>(eval_tmp(is), 0.0);
                    } else {
                        Dymat(is, js) = complex_zero;
                    }
                }
            }

            // New eigenvector matrix E_{new}= E_{old} * C
            mat_tmp = evec_tmp.transpose() * saes.eigenvectors();
            Dymat = mat_tmp * Dymat * mat_tmp.adjoint();

#ifdef _DEBUG
            Dymat_sym = Dymat;
            symmetrize_dynamical_matrix(ik, Dymat_sym);
            std::complex<double> **dymat_exact;
            memory->allocate(dymat_exact, ns, ns);
            std::cout << "ik = " << ik + 1 << std::endl;
            std::cout << "Dymat" << std::endl;
            std::cout << Dymat << std::endl;
            std::cout << "Dymat_sym" << std::endl;
            std::cout << Dymat_sym << std::endl;
            dynamical->calc_analytic_k(xk_interpolate[knum_interpolate], fcs_phonon->fc2_ext, dymat_exact);
            for (is = 0; is < ns; ++is) {
                for (js = 0; js < ns; ++js) {
                    Dymat_sym(is,js) = dymat_exact[is][js];
                }
            }
            std::cout << "Dymat_exact" << std::endl;
            std::cout << Dymat_sym << std::endl;
            memory->deallocate(dymat_exact);

#endif
            symmetrize_dynamical_matrix(ik, Dymat);

            for (is = 0; is < ns; ++is) {
                for (js = 0; js < ns; ++js) {
                    dymat_q[is][js][knum_interpolate] = Dymat(is, js);
                }
            }
        } // close loop ik

        replicate_dymat_for_all_kpoints(dymat_q);

#ifdef _DEBUG
        for (ik = 0; ik < nk_interpolate; ++ik) {

            knum = kmap_interpolate_to_scph[ik];
                
            for (is = 0; is < ns; ++is) {
                for (js = 0; js < ns; ++js) {
                    Dymat(is,js) = dymat_q[is][js][ik];
                }
            }

            saes.compute(Dymat);
            eval_tmp = saes.eigenvalues();

            for (is = 0; is < ns; ++is) {
                eval_orig(is) = omega2_harmonic(knum,is);
            }

            std::cout << " ik = " << std::setw(4) << ik + 1 << " : ";
            for (i = 0; i < 3; ++i)  std::cout << std::setw(15) << xk_scph[knum][i];
            std::cout << std::endl;

            for (is = 0; is < ns; ++is) {
                std::cout << std::setw(15) << eval_tmp(is);
                std::cout << std::setw(15) << eval_orig(is);
                std::cout << std::setw(15) << eval_tmp(is) - eval_orig(is) << std::endl;
            }
            
        }
#endif

        // Subtract harmonic contribution to the dynamical matrix
        for (ik = 0; ik < nk_interpolate; ++ik) {
            for (is = 0; is < ns; ++is) {
                for (js = 0; js < ns; ++js) {
                    dymat_q[is][js][ik] -= dymat_harmonic[ik][is][js];
                }
            }
        }

        for (is = 0; is < ns; ++is) {
            for (js = 0; js < ns; ++js) {
                fftw_plan plan = fftw_plan_dft_3d(nk1, nk2, nk3,
                                                  reinterpret_cast<fftw_complex*>(dymat_q[is][js]),
                                                  reinterpret_cast<fftw_complex*>(dymat_new[is][js]),
                                                  FFTW_FORWARD, FFTW_ESTIMATE);
                fftw_execute(plan);
                fftw_destroy_plan(plan);

                for (ik = 0; ik < nk_interpolate; ++ik)
                    dymat_new[is][js][ik] /= static_cast<double>(nk_interpolate);
            }
        }

        exec_interpolation(kmesh_interpolate,
                           dymat_new,
                           nk,
                           xk_scph,
                           kvec_na_scph,
                           eval_interpolate, evec_new);

        for (ik = 0; ik < nk; ++ik) {

            for (is = 0; is < ns; ++is) {
                for (js = 0; js < ns; ++js) {
                    mat_tmp(is, js) = evec_initial[ik][js][is];
                    evec_tmp(is, js) = evec_new[ik][js][is];
                }
            }

            Cmat = mat_tmp.adjoint() * evec_tmp;

            for (is = 0; is < ns; ++is) {
                omega_now(ik, is) = eval_interpolate[ik][is];
                for (js = 0; js < ns; ++js) {
                    cmat_convert[ik][is][js] = Cmat(is, js);
                }
            }
        }

        if (iloop == 0) {
            std::cout << "  SCPH ITER " << std::setw(5) << iloop + 1 << " :  DIFF = N/A" << std::endl;

            for (ik = 0; ik < nk; ++ik) {
                for (is = 0; is < ns; ++is) {
                    omega_old(ik, is) = omega_now(ik, is);
                }
            }

        } else {

            std::cout << "  SCPH ITER " << std::setw(5) << iloop + 1 << " : ";
            diff = 0.0;

            for (ik = 0; ik < nk_interpolate; ++ik) {
                knum = kmap_interpolate_to_scph[ik];
                for (is = 0; is < ns; ++is) {
                    diff += std::pow(omega_now(knum, is) - omega_old(knum, is), 2.0);
                }
            }
            diff /= static_cast<double>(nk_interpolate * ns);
            std::cout << " DIFF = " << std::setw(15) << std::sqrt(diff) << std::endl;
            for (ik = 0; ik < nk; ++ik) {
                for (is = 0; is < ns; ++is) {
                    omega_old(ik, is) = omega_now(ik, is);
                }
            }
            if (std::sqrt(diff) < conv_tol) {
                has_negative = false;

                for (ik = 0; ik < nk_interpolate; ++ik) {
                    knum = kmap_interpolate_to_scph[ik];
                    for (is = 0; is < ns; ++is) {
                        if (omega_now(knum, is) < 0.0 && std::abs(omega_now(knum, is)) > eps8) {
                            has_negative = true;
                            break;
                        }
                    }
                }
                if (!has_negative) {
                    std::cout << "  DIFF < SCPH_TOL : break SCPH loop" << std::endl;
                    break;
                }
                std::cout << "  DIFF < SCPH_TOL but a negative frequency is detected." << std::endl;
            }
        }

        for (ik = 0; ik < nk; ++ik) {
            for (is = 0; is < ns; ++is) {
                for (js = 0; js < ns; ++js) {
                    dmat_convert_old[ik][is][js] = dmat_convert[ik][is][js];
                }
            }
        }
    } // end loop iteration

    if (std::sqrt(diff) < conv_tol) {
        std::cout << " Temp = " << T_in;
        std::cout << " : convergence achieved in " << std::setw(5)
            << iloop + 1 << " iterations." << std::endl;
        flag_converged = true;
    } else {
        std::cout << "Temp = " << T_in;
        std::cout << " : not converged." << std::endl;
        flag_converged = false;
    }

    for (ik = 0; ik < nk; ++ik) {
        for (is = 0; is < ns; ++is) {
            if (eval_interpolate[ik][is] < 0.0) {
                if (std::abs(eval_interpolate[ik][is]) <= eps10) {
                    omega2_out[ik][is] = 0.0;
                } else {
                    omega2_out[ik][is] = -std::pow(eval_interpolate[ik][is], 2.0);
                }
            } else {
                omega2_out[ik][is] = std::pow(eval_interpolate[ik][is], 2.0);
            }
            for (js = 0; js < ns; ++js) {
                evec_anharm_scph[ik][is][js] = evec_new[ik][is][js];
            }
        }
    }

    std::cout << "New eigenvalues" << std::endl;
    for (ik = 0; ik < nk_interpolate; ++ik) {
        knum = kmap_interpolate_to_scph[ik];
        for (is = 0; is < ns; ++is) {
            std::cout << " ik_interpolate = " << std::setw(5) << ik + 1;
            std::cout << " is = " << std::setw(5) << is + 1;
            std::cout << " omega2 = " << std::setw(15) << omega2_out[knum][is] << std::endl;
        }
        std::cout << std::endl;
    }

    memory->deallocate(mat_omega2_harmonic);
    memory->deallocate(eval_interpolate);
    memory->deallocate(evec_initial);
    memory->deallocate(dmat_convert);
    memory->deallocate(dmat_convert_old);
    memory->deallocate(evec_new);
    memory->deallocate(dymat_new);
    memory->deallocate(dymat_q);
    memory->deallocate(dymat_harmonic);
    memory->deallocate(Fmat0);
}


void Scph::compute_free_energy_bubble_SCPH(const unsigned int kmesh[3],
                                           std::complex<double> ****delta_dymat_scph)
{
    const auto NT = static_cast<unsigned int>((system->Tmax - system->Tmin) / system->dT) + 1;
    const auto nk_ref = kpoint->nk;
    const auto ns = dynamical->neval;
    double ***eval;
    std::complex<double> ****evec;

    if (mympi->my_rank == 0) {
        std::cout << std::endl;
        std::cout << " -----------------------------------------------------------------"
            << std::endl;
        std::cout << " Calculating the vibrational free energy from the Bubble diagram " << std::endl;
        std::cout << " on top of the SCPH calculation." << std::endl;
        std::cout << '\n';
        std::cout << " This calculation requires allocation of additional memory:" << std::endl;

        size_t nsize = nk_ref * ns * ns * NT * sizeof(std::complex<double>)
            + nk_ref * ns * NT * sizeof(double);

        const auto nsize_dble = static_cast<double>(nsize) / 100000000.0;
        std::cout << "  Estimated memory usage per MPI process: " << std::setw(10)
            << std::fixed << std::setprecision(4) << nsize_dble << " GByte." << std::endl;

        std::cout << "  To avoid possible faults associated with insufficient memory,\n"
            "  please reduce the number of MPI processes per node and/or\n"
            "  the number of temperagure grids.\n\n";
    }

    memory->allocate(thermodynamics->FE_bubble, NT);
    memory->allocate(eval, NT, nk_ref, ns);
    memory->allocate(evec, NT, nk_ref, ns, ns); // This requires lots of RAM

    for (auto iT = 0; iT < NT; ++iT) {
        const auto temp = system->Tmin + system->dT * float(iT);

        exec_interpolation(kmesh,
                           delta_dymat_scph[iT],
                           nk_ref,
                           kpoint->xk,
                           kpoint->kvec_na,
                           eval[iT],
                           evec[iT]);
    }

    thermodynamics->compute_FE_bubble_SCPH(eval, evec, thermodynamics->FE_bubble);

    memory->deallocate(eval);
    memory->deallocate(evec);

    if (mympi->my_rank == 0) {
        std::cout << " done!" << std::endl << std::endl;
    }
}


double Scph::distance(double *x1,
                      double *x2) const
{
    auto dist = std::pow(x1[0] - x2[0], 2)
        + std::pow(x1[1] - x2[1], 2)
        + std::pow(x1[2] - x2[2], 2);
    dist = std::sqrt(dist);

    return dist;
}

void Scph::duplicate_xk_boundary(double *xk_in,
                                 std::vector<std::vector<double>> &vec_xk) const
{
    int i;
    int n[3];
    double sign[3];
    std::vector<double> vec_tmp;

    vec_xk.clear();

    for (i = 0; i < 3; ++i) {
        if (std::abs(std::abs(xk_in[i]) - 0.5) < eps) {
            n[i] = 2;
        } else {
            n[i] = 1;
        }
    }

    for (i = 0; i < n[0]; ++i) {
        sign[0] = 1.0 - 2.0 * static_cast<double>(i);
        for (int j = 0; j < n[1]; ++j) {
            sign[1] = 1.0 - 2.0 * static_cast<double>(j);
            for (int k = 0; k < n[2]; ++k) {
                sign[2] = 1.0 - 2.0 * static_cast<double>(k);

                vec_tmp.clear();
                for (int l = 0; l < 3; ++l) {
                    vec_tmp.push_back(sign[l] * xk_in[l]);
                }
                vec_xk.push_back(vec_tmp);

            }
        }
    }
}

void Scph::write_anharmonic_correction_fc2(std::complex<double> ****delta_dymat,
                                           const unsigned int NT)
{
    unsigned int i, j;
    const auto Tmin = system->Tmin;
    const auto dT = system->dT;
    double ***delta_fc2;
    double **xtmp;
    const auto ns = dynamical->neval;
    unsigned int is, js, icell;
    unsigned int iat, jat;

    std::ofstream ofs_fc2;
    auto file_fc2 = input->job_title + ".scph_dfc2";

    ofs_fc2.open(file_fc2.c_str(), std::ios::out);
    if (!ofs_fc2)
        error->exit("write_anharmonic_correction_fc2",
                    "Cannot open file_fc2");

    const auto ncell = kmesh_interpolate[0] * kmesh_interpolate[1] * kmesh_interpolate[2];

    memory->allocate(delta_fc2, ns, ns, ncell);

    memory->allocate(xtmp, system->natmin, 3);

    ofs_fc2.precision(10);

    for (i = 0; i < system->natmin; ++i) {
        rotvec(xtmp[i], system->xr_s[system->map_p2s[i][0]], system->lavec_s);
        rotvec(xtmp[i], xtmp[i], system->rlavec_p);
        for (j = 0; j < 3; ++j) xtmp[i][j] /= 2.0 * pi;
    }

    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            ofs_fc2 << std::setw(20) << system->lavec_p[j][i];
        }
        ofs_fc2 << std::endl;
    }
    ofs_fc2 << std::setw(5) << system->natmin << std::setw(5) << system->nkd << std::endl;
    for (i = 0; i < system->nkd; ++i) {
        ofs_fc2 << std::setw(5) << system->symbol_kd[i];
    }
    ofs_fc2 << std::endl;

    for (i = 0; i < system->natmin; ++i) {
        for (j = 0; j < 3; ++j) {
            ofs_fc2 << std::setw(15) << xtmp[i][j];
        }
        ofs_fc2 << std::setw(5) << system->kd[system->map_p2s[i][0]] + 1 << std::endl;
    }

    memory->deallocate(xtmp);

    for (unsigned int iT = 0; iT < NT; ++iT) {
        const auto temp = Tmin + dT * static_cast<double>(iT);

        ofs_fc2 << "# Temp = " << temp << std::endl;

        for (is = 0; is < ns; ++is) {
            iat = is / 3;

            for (js = 0; js < ns; ++js) {
                jat = js / 3;

                for (icell = 0; icell < ncell; ++icell) {
                    delta_fc2[is][js][icell]
                        = delta_dymat[iT][is][js][icell].real()
                        * std::sqrt(system->mass[system->map_p2s[iat][0]]
                            * system->mass[system->map_p2s[jat][0]]);
                }

            }
        }


        for (icell = 0; icell < ncell; ++icell) {

            for (is = 0; is < ns; ++is) {
                iat = is / 3;
                const auto icrd = is % 3;

                for (js = 0; js < ns; ++js) {
                    jat = js / 3;
                    const auto jcrd = js % 3;

                    const auto nmulti = mindist_list_scph[iat][jat][icell].shift.size();

                    for (auto it = mindist_list_scph[iat][jat][icell].shift.cbegin();
                         it != mindist_list_scph[iat][jat][icell].shift.cend(); ++it) {

                        ofs_fc2 << std::setw(4) << (*it).sx;
                        ofs_fc2 << std::setw(4) << (*it).sy;
                        ofs_fc2 << std::setw(4) << (*it).sz;
                        ofs_fc2 << std::setw(5) << iat << std::setw(3) << icrd;
                        ofs_fc2 << std::setw(4) << jat << std::setw(3) << jcrd;
                        ofs_fc2 << std::setprecision(15) << std::setw(25)
                            << delta_fc2[is][js][icell] / static_cast<double>(nmulti) << std::endl;

                    }

                }
            }
        }

        ofs_fc2 << std::endl;
    }

    memory->deallocate(delta_fc2);

    ofs_fc2.close();
    std::cout << "  " << std::setw(input->job_title.length() + 12) << std::left << file_fc2;
    std::cout << " : Anharmonic corrections to the second-order IFCs" << std::endl;
}

void Scph::mpi_bcast_complex(std::complex<double> ****data,
                             const int NT,
                             const int nk,
                             const int ns) const
{
#ifdef MPI_COMPLEX16
    MPI_Bcast(&data[0][0][0][0], NT * nk * ns * ns, MPI_COMPLEX16, 0, MPI_COMM_WORLD);
#elif defined MPI_DOUBLE_COMPLEX
    MPI_Bcast(&data[0][0][0][0], NT * nk * ns * ns, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
#else
    unsigned int iT, ik, is, js;
    double ***data_real, ***data_imag;

    memory->allocate(data_real, ns, ns, nk);
    memory->allocate(data_imag, ns, ns, nk);

    for (iT = 0; iT < NT; ++iT) {

        if (mympi->my_rank == 0) {
            for (is = 0; is < ns; ++is) {
                for (js = 0; js < ns; ++js) {
                    for (ik = 0; ik < nk; ++ik) {

                        data_real[is][js][ik] = data[iT][is][js][ik].real();
                        data_imag[is][js][ik] = data[iT][is][js][ik].imag();
                    }
                }
            }
        }

        MPI_Bcast(&data_real[0][0][0], nk * ns * ns, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&data_imag[0][0][0], nk * ns * ns, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        if (mympi->my_rank > 0) {
            for (iT = 0; iT < NT; ++iT) {
                for (is = 0; is < ns; ++is) {
                    for (js = 0; js < ns; ++js) {
                        for (ik = 0; ik < nk; ++ik) {
                            data[iT][is][js][ik]
                                = std::complex<double>(data_real[is][js][ik],
                                                       data_imag[is][js][ik]);
                        }
                    }
                }
            }
        }
    }

    memory->deallocate(data_real);
    memory->deallocate(data_imag);
#endif
}
