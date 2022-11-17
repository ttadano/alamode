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
#include "dielec.h"
#include "ewald.h"
#include "memory.h"
#include "thermodynamics.h"
#include "write_phonons.h"
#include "constants.h"
#include "system.h"
#include "error.h"
#include "mathfunctions.h"
#include "integration.h"
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

Scph::Scph(PHON *phon) : Pointers(phon)
{
    set_default_variables();
}

Scph::~Scph()
{
    deallocate_variables();
}

void Scph::set_default_variables()
{
    restart_scph = false;
    warmstart_scph = false;
    lower_temp = true;
    tolerance_scph = 1.0e-10;
    mixalpha = 0.1;
    maxiter = 100;
    print_self_consistent_fc2 = false;
    selfenergy_offdiagonal = true;
    relax_coordinate = false;

    kmap_interpolate_to_scph = nullptr;
    evec_harmonic = nullptr;
    omega2_harmonic = nullptr;
    mat_transform_sym = nullptr;
    symop_minus_at_k = nullptr;
    kpoint_map_symmetry = nullptr;
    mindist_list_scph = nullptr;

    bubble = 0;
    phi3_reciprocal = nullptr;
    phi4_reciprocal = nullptr;
    compute_Cv_anharmonic = 0;

    kmesh_coarse = nullptr;
    kmesh_dense = nullptr;

    phase_factor_scph = nullptr;
}

void Scph::deallocate_variables()
{
    if (kmap_interpolate_to_scph) {
        deallocate(kmap_interpolate_to_scph);
    }
    if (mindist_list_scph) {
        deallocate(mindist_list_scph);
    }
    if (evec_harmonic) {
        deallocate(evec_harmonic);
    }
    if (omega2_harmonic) {
        deallocate(omega2_harmonic);
    }
    if (mat_transform_sym) {
        deallocate(mat_transform_sym);
    }
    if (symop_minus_at_k) {
        deallocate(symop_minus_at_k);
    }
    if (kpoint_map_symmetry) {
        deallocate(kpoint_map_symmetry);
    }
    if (phi3_reciprocal) {
        deallocate(phi3_reciprocal);
    }
    if (phi4_reciprocal) {
        deallocate(phi4_reciprocal);
    }

    if (kmesh_coarse) delete kmesh_coarse;
    if (kmesh_dense) delete kmesh_dense;

    if (phase_factor_scph) delete phase_factor_scph;
}

void Scph::setup_scph()
{
    relax_coordinate = false;
    MPI_Bcast(&relax_coordinate, 1, MPI_CXX_BOOL, 0, MPI_COMM_WORLD);
    MPI_Bcast(&bubble, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

    setup_kmesh();
    setup_eigvecs();
    setup_transform_ifc();
    setup_pp_interaction();
    setup_transform_symmetry();
}

void Scph::exec_scph()
{
    const auto ns = dynamical->neval;
    const auto Tmin = system->Tmin;
    const auto Tmax = system->Tmax;
    const auto dT = system->dT;

    std::complex<double> ****delta_dymat_scph = nullptr;
    std::complex<double> ****delta_dymat_scph_plus_bubble = nullptr;

    const auto NT = static_cast<unsigned int>((Tmax - Tmin) / dT) + 1;

    MPI_Bcast(&restart_scph, 1, MPI_CXX_BOOL, 0, MPI_COMM_WORLD);
    MPI_Bcast(&selfenergy_offdiagonal, 1, MPI_CXX_BOOL, 0, MPI_COMM_WORLD);
    MPI_Bcast(&ialgo, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

    allocate(delta_dymat_scph, NT, ns, ns, kmesh_coarse->nk);

    if (restart_scph) {

        // Read anharmonic correction to the dynamical matrix from the existing file
        load_scph_dymat_from_file(delta_dymat_scph);

    } else {

        if (dynamical->nonanalytic == 3) {
            exit("exec_scph",
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

    if (bubble) {
        allocate(delta_dymat_scph_plus_bubble, NT, ns, ns, kmesh_coarse->nk);
        bubble_correction(delta_dymat_scph,
                          delta_dymat_scph_plus_bubble);
        if (mympi->my_rank == 0) {
            write_anharmonic_correction_fc2(delta_dymat_scph_plus_bubble, NT, bubble);
        }
    }

    postprocess(delta_dymat_scph,
                delta_dymat_scph_plus_bubble);

    deallocate(delta_dymat_scph);
    if (delta_dymat_scph_plus_bubble) deallocate(delta_dymat_scph_plus_bubble);

}

void Scph::postprocess(std::complex<double> ****delta_dymat_scph,
                       std::complex<double> ****delta_dymat_scph_plus_bubble)
{
    double ***eval_anharm = nullptr;
    const auto ns = dynamical->neval;
    const auto Tmin = system->Tmin;
    const auto Tmax = system->Tmax;
    const auto dT = system->dT;
    const auto NT = static_cast<unsigned int>((Tmax - Tmin) / dT) + 1;

    unsigned int nomega_dielec;

    if (mympi->my_rank == 0) {

        std::cout << '\n';
        std::cout << " Running postprocess of SCPH (calculation of free energy, MSD, DOS)" << std::endl;
        std::cout << " The number of temperature points: " << std::setw(4) << NT << std::endl;
        std::cout << "   ";

        std::complex<double> ***evec_tmp = nullptr;
        double **eval_gam = nullptr;
        std::complex<double> ***evec_gam = nullptr;
        double **xk_gam = nullptr;

        double **dos_scph = nullptr;
        double ***pdos_scph = nullptr;
        double *heat_capacity = nullptr;
        double *heat_capacity_correction = nullptr;
        double *FE_QHA = nullptr;
        double *dFE_scph = nullptr;
        double **msd_scph = nullptr;
        double ***ucorr_scph = nullptr;
        double ****dielec_scph = nullptr;
        double *omega_grid = nullptr;
        double **domega_dt = nullptr;

        if (dos->kmesh_dos) {
            allocate(eval_anharm, NT, dos->kmesh_dos->nk, ns);
            allocate(evec_tmp, dos->kmesh_dos->nk, ns, ns);

            if (dos->compute_dos) {
                allocate(dos_scph, NT, dos->n_energy);

                if (dos->projected_dos) {
                    allocate(pdos_scph, NT, ns, dos->n_energy);
                }
            }
            allocate(heat_capacity, NT);
            allocate(FE_QHA, NT);
            allocate(dFE_scph, NT);

            if (writes->getPrintMSD()) {
                allocate(msd_scph, NT, ns);
            }
            if (writes->getPrintUcorr()) {
                allocate(ucorr_scph, NT, ns, ns);
            }
            if (compute_Cv_anharmonic) {
                allocate(heat_capacity_correction, NT);
                allocate(domega_dt, dos->kmesh_dos->nk, ns);
                if (compute_Cv_anharmonic == 1) {
                    // Use central difference to evaluate temperature derivative of
                    // anharmonic frequencies
                    heat_capacity_correction[0] = 0.0;
                    heat_capacity_correction[NT - 1] = 0.0;
                }
            }

            for (auto iT = 0; iT < NT; ++iT) {
                auto temperature = Tmin + dT * static_cast<double>(iT);

                exec_interpolation(kmesh_interpolate,
                                   delta_dymat_scph[iT],
                                   dos->kmesh_dos->nk,
                                   dos->kmesh_dos->xk,
                                   dos->kmesh_dos->kvec_na,
                                   eval_anharm[iT],
                                   evec_tmp);

                if (dos->compute_dos) {
                    dos->calc_dos_from_given_frequency(dos->kmesh_dos,
                                                       eval_anharm[iT],
                                                       dos->tetra_nodes_dos->get_ntetra(),
                                                       dos->tetra_nodes_dos->get_tetras(),
                                                       dos_scph[iT]);
                }

                heat_capacity[iT] = thermodynamics->Cv_tot(temperature,
                                                           dos->kmesh_dos->nk_irred,
                                                           ns,
                                                           dos->kmesh_dos->kpoint_irred_all,
                                                           &dos->kmesh_dos->weight_k[0],
                                                           eval_anharm[iT]);

                FE_QHA[iT] = thermodynamics->free_energy_QHA(temperature,
                                                             dos->kmesh_dos->nk_irred,
                                                             ns,
                                                             dos->kmesh_dos->kpoint_irred_all,
                                                             &dos->kmesh_dos->weight_k[0],
                                                             eval_anharm[iT]);

                dFE_scph[iT] = thermodynamics->FE_scph_correction(iT,
                                                                  eval_anharm[iT],
                                                                  evec_tmp);

                if (writes->getPrintMSD()) {
                    double shift[3]{0.0, 0.0, 0.0};

                    for (auto is = 0; is < ns; ++is) {
                        msd_scph[iT][is] = thermodynamics->disp_corrfunc(temperature,
                                                                         is, is,
                                                                         shift,
                                                                         dos->kmesh_dos->nk,
                                                                         ns,
                                                                         dos->kmesh_dos->xk,
                                                                         eval_anharm[iT],
                                                                         evec_tmp);
                    }
                }

                if (writes->getPrintUcorr()) {
                    double shift[3];
                    for (auto i = 0; i < 3; ++i) shift[i] = static_cast<double>(writes->getShiftUcorr()[i]);

                    for (auto is = 0; is < ns; ++is) {
                        for (auto js = 0; js < ns; ++js) {
                            ucorr_scph[iT][is][js] = thermodynamics->disp_corrfunc(temperature,
                                                                                   is, js,
                                                                                   shift,
                                                                                   dos->kmesh_dos->nk,
                                                                                   ns,
                                                                                   dos->kmesh_dos->xk,
                                                                                   eval_anharm[iT],
                                                                                   evec_tmp);
                        }
                    }
                }

                if (compute_Cv_anharmonic == 1) {

                    if (iT >= 1 and iT <= NT - 2) {
                        get_derivative_central_diff(dT, dos->kmesh_dos->nk,
                                                    eval_anharm[iT - 1],
                                                    eval_anharm[iT + 1],
                                                    domega_dt);

                        heat_capacity_correction[iT] = thermodynamics->Cv_anharm_correction(temperature,
                                                                                            dos->kmesh_dos->nk_irred,
                                                                                            ns,
                                                                                            dos->kmesh_dos->kpoint_irred_all,
                                                                                            &dos->kmesh_dos->weight_k[0],
                                                                                            eval_anharm[iT],
                                                                                            domega_dt);
                    }

                }

                std::cout << '.' << std::flush;
                if (iT % 25 == 24) {
                    std::cout << std::endl;
                    std::cout << std::setw(3);
                }
            }
            std::cout << "\n\n";

            if (dos->compute_dos) {
                writes->write_scph_dos(dos_scph);
            }
            writes->write_scph_thermodynamics(heat_capacity,
                                              heat_capacity_correction,
                                              FE_QHA,
                                              dFE_scph);
            if (writes->getPrintMSD()) {
                writes->write_scph_msd(msd_scph);
            }
            if (writes->getPrintUcorr()) {
                writes->write_scph_ucorr(ucorr_scph);
            }

            // If delta_dymat_scph_plus_bubble != nullptr, run postprocess again with
            // delta_dymat_scph_plus_bubble.
            if (bubble > 0) {
                std::cout << std::endl;
                std::cout << "   ";

                for (auto iT = 0; iT < NT; ++iT) {
                    auto temperature = Tmin + dT * static_cast<double>(iT);

                    exec_interpolation(kmesh_interpolate,
                                       delta_dymat_scph_plus_bubble[iT],
                                       dos->kmesh_dos->nk,
                                       dos->kmesh_dos->xk,
                                       dos->kmesh_dos->kvec_na,
                                       eval_anharm[iT],
                                       evec_tmp);

                    if (dos->compute_dos) {
                        dos->calc_dos_from_given_frequency(dos->kmesh_dos,
                                                           eval_anharm[iT],
                                                           dos->tetra_nodes_dos->get_ntetra(),
                                                           dos->tetra_nodes_dos->get_tetras(),
                                                           dos_scph[iT]);
                    }

                    heat_capacity[iT] = thermodynamics->Cv_tot(temperature,
                                                               dos->kmesh_dos->nk_irred,
                                                               ns,
                                                               dos->kmesh_dos->kpoint_irred_all,
                                                               &dos->kmesh_dos->weight_k[0],
                                                               eval_anharm[iT]);

                    if (writes->getPrintMSD()) {
                        double shift[3]{0.0, 0.0, 0.0};

                        for (auto is = 0; is < ns; ++is) {
                            msd_scph[iT][is] = thermodynamics->disp_corrfunc(temperature,
                                                                             is, is,
                                                                             shift,
                                                                             dos->kmesh_dos->nk,
                                                                             ns,
                                                                             dos->kmesh_dos->xk,
                                                                             eval_anharm[iT],
                                                                             evec_tmp);
                        }
                    }

                    if (writes->getPrintUcorr()) {
                        double shift[3];
                        for (auto i = 0; i < 3; ++i) shift[i] = static_cast<double>(writes->getShiftUcorr()[i]);

                        for (auto is = 0; is < ns; ++is) {
                            for (auto js = 0; js < ns; ++js) {
                                ucorr_scph[iT][is][js] = thermodynamics->disp_corrfunc(temperature,
                                                                                       is, js,
                                                                                       shift,
                                                                                       dos->kmesh_dos->nk,
                                                                                       ns,
                                                                                       dos->kmesh_dos->xk,
                                                                                       eval_anharm[iT],
                                                                                       evec_tmp);
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

                if (dos->compute_dos) {
                    writes->write_scph_dos(dos_scph, bubble);
                }
                if (writes->getPrintMSD()) {
                    writes->write_scph_msd(msd_scph, bubble);
                }
                if (writes->getPrintUcorr()) {
                    writes->write_scph_ucorr(ucorr_scph, bubble);
                }

            }
            deallocate(eval_anharm);
            eval_anharm = nullptr;
            deallocate(evec_tmp);
            evec_tmp = nullptr;
        }

        if (kpoint->kpoint_general) {
            allocate(eval_anharm, NT, kpoint->kpoint_general->nk, ns);
            allocate(evec_tmp, kpoint->kpoint_general->nk, ns, ns);

            for (auto iT = 0; iT < NT; ++iT) {
                exec_interpolation(kmesh_interpolate,
                                   delta_dymat_scph[iT],
                                   kpoint->kpoint_general->nk,
                                   kpoint->kpoint_general->xk,
                                   kpoint->kpoint_general->kvec_na,
                                   eval_anharm[iT],
                                   evec_tmp);
            }

            writes->write_scph_energy(kpoint->kpoint_general->nk,
                                      eval_anharm);

            if (bubble > 0) {
                for (auto iT = 0; iT < NT; ++iT) {
                    exec_interpolation(kmesh_interpolate,
                                       delta_dymat_scph_plus_bubble[iT],
                                       kpoint->kpoint_general->nk,
                                       kpoint->kpoint_general->xk,
                                       kpoint->kpoint_general->kvec_na,
                                       eval_anharm[iT],
                                       evec_tmp);
                }
                writes->write_scph_energy(kpoint->kpoint_general->nk,
                                          eval_anharm,
                                          bubble);
            }
            deallocate(eval_anharm);
            deallocate(evec_tmp);
            eval_anharm = nullptr;
            evec_tmp = nullptr;
        }

        if (kpoint->kpoint_bs) {
            allocate(eval_anharm, NT, kpoint->kpoint_bs->nk, ns);
            allocate(evec_tmp, kpoint->kpoint_bs->nk, ns, ns);

            for (auto iT = 0; iT < NT; ++iT) {
                exec_interpolation(kmesh_interpolate,
                                   delta_dymat_scph[iT],
                                   kpoint->kpoint_bs->nk,
                                   kpoint->kpoint_bs->xk,
                                   kpoint->kpoint_bs->kvec_na,
                                   eval_anharm[iT],
                                   evec_tmp);
            }

            writes->write_scph_bands(kpoint->kpoint_bs->nk,
                                     kpoint->kpoint_bs->kaxis,
                                     eval_anharm);

            if (bubble > 0) {
                for (auto iT = 0; iT < NT; ++iT) {
                    exec_interpolation(kmesh_interpolate,
                                       delta_dymat_scph_plus_bubble[iT],
                                       kpoint->kpoint_bs->nk,
                                       kpoint->kpoint_bs->xk,
                                       kpoint->kpoint_bs->kvec_na,
                                       eval_anharm[iT],
                                       evec_tmp);
                }
                writes->write_scph_bands(kpoint->kpoint_bs->nk,
                                         kpoint->kpoint_bs->kaxis,
                                         eval_anharm,
                                         bubble);
            }
            deallocate(eval_anharm);
            deallocate(evec_tmp);
            eval_anharm = nullptr;
            evec_tmp = nullptr;
        }

        if (dielec->calc_dielectric_constant) {
            omega_grid = dielec->get_omega_grid(nomega_dielec);
            allocate(dielec_scph, NT, nomega_dielec, 3, 3);
            allocate(eval_gam, 1, ns);
            allocate(evec_gam, 1, ns, ns);
            allocate(xk_gam, 1, 3);
            for (auto i = 0; i < 3; ++i) xk_gam[0][i] = 0.0;

            for (auto iT = 0; iT < NT; ++iT) {
                exec_interpolation(kmesh_interpolate,
                                   delta_dymat_scph[iT],
                                   1,
                                   xk_gam,
                                   xk_gam,
                                   eval_gam,
                                   evec_gam);

                for (auto is = 0; is < ns; ++is) {
                    if (eval_gam[0][is] < 0.0) {
                        eval_gam[0][is] = -std::pow(eval_gam[0][is], 2.0);
                    } else {
                        eval_gam[0][is] = std::pow(eval_gam[0][is], 2.0);
                    }
                }
                dielec->compute_dielectric_function(nomega_dielec,
                                                    omega_grid,
                                                    eval_gam[0],
                                                    evec_gam[0],
                                                    dielec_scph[iT]);
            }
            writes->write_scph_dielec(dielec_scph);
        }

        if (eval_anharm) deallocate(eval_anharm);
        if (evec_tmp) deallocate(evec_tmp);

        if (dos_scph) deallocate(dos_scph);
        if (pdos_scph) deallocate(pdos_scph);
        if (heat_capacity) deallocate(heat_capacity);
        if (heat_capacity_correction) deallocate(heat_capacity_correction);
        if (FE_QHA) deallocate(FE_QHA);
        if (dFE_scph) deallocate(dFE_scph);
        if (dielec_scph) deallocate(dielec_scph);

        if (eval_gam) deallocate(eval_gam);
        if (evec_gam) deallocate(evec_gam);
        if (xk_gam) deallocate(xk_gam);
    }
}

void Scph::load_scph_dymat_from_file(std::complex<double> ****dymat_out)
{
    const auto ns = dynamical->neval;
    const auto Tmin = system->Tmin;
    const auto Tmax = system->Tmax;
    const auto dT = system->dT;
    const auto NT = static_cast<unsigned int>((Tmax - Tmin) / dT) + 1;
    std::vector<double> Temp_array(NT);

    for (int i = 0; i < NT; ++i) {
        Temp_array[i] = Tmin + dT * static_cast<double>(i);
    }

    if (mympi->my_rank == 0) {

        const auto consider_offdiagonal = selfenergy_offdiagonal;
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
            exit("load_scph_dymat_from_file",
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
            exit("load_scph_dymat_from_file",
                 "The number of KMESH_INTERPOLATE is not consistent");
        }
        if (nk_scph_tmp[0] != kmesh_scph[0] ||
            nk_scph_tmp[1] != kmesh_scph[1] ||
            nk_scph_tmp[2] != kmesh_scph[2]) {
            exit("load_scph_dymat_from_file",
                 "The number of KMESH_SCPH is not consistent");
        }
        if (nonanalytic_tmp != dynamical->nonanalytic) {
            warn("load_scph_dymat_from_file",
                 "The NONANALYTIC tag is not consistent");
        }
        if (consider_offdiag_tmp != consider_offdiagonal) {
            exit("load_scph_dymat_from_file",
                 "The SELF_OFFDIAG tag is not consistent");
        }

        // Check if the precalculated data for the given temperature range exists
        const auto NT_ref = static_cast<unsigned int>((Tmax_tmp - Tmin_tmp) / dT_tmp) + 1;
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
                    for (int ik = 0; ik < kmesh_coarse->nk; ++ik) {
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
            exit("load_scph_dymat_from_file",
                 "The temperature information is not consistent");
        }
        std::cout << " done." << std::endl;
    }
    // Broadcast to all MPI threads
    mpi_bcast_complex(dymat_out, NT, kmesh_coarse->nk, ns);
}

void Scph::store_scph_dymat_to_file(const std::complex<double> *const *const *const *dymat_in)
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
        exit("store_scph_dymat_to_file",
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
                for (auto ik = 0; ik < kmesh_coarse->nk; ++ik) {
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
    int ik, is;
    const auto nk = kmesh_dense->nk;
    const auto ns = dynamical->neval;
    const auto nk_irred_interpolate = kmesh_coarse->nk_irred;
    const auto Tmin = system->Tmin;
    const auto Tmax = system->Tmax;
    const auto dT = system->dT;
    double ***omega2_anharm;
    std::complex<double> ***evec_anharm_tmp;
    std::complex<double> ***v3_array_all;
    std::complex<double> ***v4_array_all;

    std::vector<double> vec_temp;

    const auto NT = static_cast<unsigned int>((Tmax - Tmin) / dT) + 1;

    // Compute matrix element of 4-phonon interaction

    allocate(omega2_anharm, NT, nk, ns);
    allocate(evec_anharm_tmp, nk, ns, ns);
    allocate(v4_array_all, nk_irred_interpolate * nk,
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
        allocate(v3_array_all, nk, ns, ns * ns);
        compute_V3_elements_mpi_over_kpoint(v3_array_all,
                                            evec_harmonic,
                                            selfenergy_offdiagonal);
    }

    if (mympi->my_rank == 0) {

        std::complex<double> ***cmat_convert;
        allocate(cmat_convert, nk, ns, ns);

        vec_temp.clear();

        if (lower_temp) {
            for (int i = NT - 1; i >= 0; --i) {
                vec_temp.push_back(Tmin + static_cast<double>(i) * dT);
            }
        } else {
            for (int i = 0; i < NT; ++i) {
                vec_temp.push_back(Tmin + static_cast<double>(i) * dT);
            }
        }

        auto converged_prev = false;

        for (double temp: vec_temp) {
            auto iT = static_cast<unsigned int>((temp - Tmin) / dT);

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
                                         converged_prev,
                                         cmat_convert,
                                         selfenergy_offdiagonal,
                                         writes->getVerbosity());

            calc_new_dymat_with_evec(dymat_anharm[iT],
                                     omega2_anharm[iT],
                                     evec_anharm_tmp);

            if (!warmstart_scph) converged_prev = false;
        }

        deallocate(cmat_convert);

    }

    mpi_bcast_complex(dymat_anharm, NT, kmesh_coarse->nk, ns);

    deallocate(omega2_anharm);
    deallocate(v4_array_all);
    deallocate(evec_anharm_tmp);
}

void Scph::compute_V3_elements_mpi_over_kpoint(std::complex<double> ***v3_out,
                                               const std::complex<double> *const *const *evec_in,
                                               const bool self_offdiag)
{
    // Calculate the matrix elements of quartic terms in reciprocal space.
    // This is the most expensive part of the SCPH calculation.

    auto ns = dynamical->neval;
    auto ns2 = ns * ns;
    auto ns3 = ns * ns * ns;
    unsigned int is, js, ks;
    unsigned int **ind;
    unsigned int i, j;
    std::complex<double> ret;
    long int ii;

    const auto nk_scph = kmesh_dense->nk;
    const auto ngroup_v3 = anharmonic_core->get_ngroup_fcs(3);
    const auto factor = std::pow(0.5, 2) / static_cast<double>(nk_scph);
    static auto complex_zero = std::complex<double>(0.0, 0.0);
    std::complex<double> *v3_array_at_kpair;
    std::complex<double> ***v3_mpi;

    if (mympi->my_rank == 0) {
        if (self_offdiag) {
            std::cout << " SELF_OFFDIAG = 1: Calculating all components of v3_array ... ";
        } else {
            std::cout << " SELF_OFFDIAG = 0: Calculating diagonal components of v3_array ... ";
        }
    }

    allocate(v3_array_at_kpair, ngroup_v3);
    allocate(ind, ngroup_v3, 3);
    allocate(v3_mpi, nk_scph, ns, ns2);

    for (unsigned int ik = mympi->my_rank; ik < nk_scph; ik += mympi->nprocs) {

        anharmonic_core->calc_phi3_reciprocal(kmesh_dense->xk[ik],
                                              kmesh_dense->xk[kmesh_dense->kindex_minus_xk[ik]],
                                              phase_factor_scph,
                                              phi3_reciprocal);

#ifdef _OPENMP
#pragma omp parallel for private(j)
#endif
        for (ii = 0; ii < ngroup_v3; ++ii) {
            v3_array_at_kpair[ii] = phi3_reciprocal[ii] * anharmonic_core->get_invmass_factor(3)[ii];
            for (j = 0; j < 3; ++j) ind[ii][j] = anharmonic_core->get_evec_index(3)[ii][j];
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
            // elements of the phonon self-energy (i.e., when considering polarization mixing).

#pragma omp parallel for private(is, js, ks, ret, i)
            for (ii = 0; ii < ns3; ++ii) {
                is = ii / ns2;
                js = (ii - ns2 * is) / ns;
                ks = ii % ns;

                ret = std::complex<double>(0.0, 0.0);

                for (i = 0; i < ngroup_v3; ++i) {

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

                    for (i = 0; i < ngroup_v3; ++i) {

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

                    for (i = 0; i < ngroup_v3; ++i) {

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

    deallocate(v3_array_at_kpair);
    deallocate(ind);
#ifdef MPI_CXX_DOUBLE_COMPLEX
    MPI_Allreduce(&v3_mpi[0][0][0], &v3_out[0][0][0],
                  static_cast<int>(nk_scph) * ns3,
                  MPI_CXX_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
#else
    MPI_Allreduce(&v3_mpi[0][0][0], &v3_out[0][0][0], static_cast<int>(nk_scph) * ns3,
                  MPI_COMPLEX16, MPI_SUM, MPI_COMM_WORLD);
#endif

    deallocate(v3_mpi);

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

    const size_t nk_reduced_interpolate = kmesh_coarse->nk_irred;
    const size_t ns = dynamical->neval;
    const size_t ns2 = ns * ns;
    const size_t ns3 = ns * ns * ns;
    const size_t ns4 = ns * ns * ns * ns;
    size_t is, js, ks, ls;
    unsigned int **ind;
    unsigned int i, j;
    std::complex<double> ret;
    long int ii;

    const auto nk_scph = kmesh_dense->nk;
    const auto ngroup_v4 = anharmonic_core->get_ngroup_fcs(4);
    const auto factor = std::pow(0.5, 2) / static_cast<double>(nk_scph);
    static auto complex_zero = std::complex<double>(0.0, 0.0);
    std::complex<double> *v4_array_at_kpair;
    std::complex<double> ***v4_mpi;

    const size_t nk2_prod = nk_reduced_interpolate * nk_scph;

    if (mympi->my_rank == 0) {
        if (self_offdiag) {
            std::cout << " SELF_OFFDIAG = 1: Calculating all components of v4_array ... ";
        } else {
            std::cout << " SELF_OFFDIAG = 0: Calculating diagonal components of v4_array ... ";
        }
    }

    allocate(v4_array_at_kpair, ngroup_v4);
    allocate(ind, ngroup_v4, 4);
    allocate(v4_mpi, nk2_prod, ns2, ns2);

    for (size_t ik_prod = mympi->my_rank; ik_prod < nk2_prod; ik_prod += mympi->nprocs) {
        const auto ik = ik_prod / nk_scph;
        const auto jk = ik_prod % nk_scph;

        const unsigned int knum = kmap_interpolate_to_scph[kmesh_coarse->kpoint_irred_all[ik][0].knum];

        anharmonic_core->calc_phi4_reciprocal(kmesh_dense->xk[knum],
                                              kmesh_dense->xk[jk],
                                              kmesh_dense->xk[kmesh_dense->kindex_minus_xk[jk]],
                                              phase_factor_scph,
                                              phi4_reciprocal);

#ifdef _OPENMP
#pragma omp parallel for private(j)
#endif
        for (ii = 0; ii < ngroup_v4; ++ii) {
            v4_array_at_kpair[ii] = phi4_reciprocal[ii] * anharmonic_core->get_invmass_factor(4)[ii];
            for (j = 0; j < 4; ++j) ind[ii][j] = anharmonic_core->get_evec_index(4)[ii][j];
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
            // elements of the phonon self-energy (loop diagram).

#pragma omp parallel for private(is, js, ks, ls, ret, i)
            for (ii = 0; ii < ns4; ++ii) {
                is = ii / ns3;
                js = (ii - ns3 * is) / ns2;
                ks = (ii - ns3 * is - ns2 * js) / ns;
                ls = ii % ns;

                if (is < js) continue;

                ret = std::complex<double>(0.0, 0.0);

                for (i = 0; i < ngroup_v4; ++i) {

                    ret += v4_array_at_kpair[i]
                           * std::conj(evec_in[knum][is][ind[i][0]])
                           * evec_in[knum][js][ind[i][1]]
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

                    for (i = 0; i < ngroup_v4; ++i) {

                        ret += v4_array_at_kpair[i]
                               * std::conj(evec_in[knum][is][ind[i][0]])
                               * evec_in[knum][js][ind[i][1]]
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

                    for (i = 0; i < ngroup_v4; ++i) {

                        ret += v4_array_at_kpair[i]
                               * std::conj(evec_in[knum][is][ind[i][0]])
                               * evec_in[knum][is][ind[i][1]]
                               * evec_in[jk][js][ind[i][2]]
                               * std::conj(evec_in[jk][js][ind[i][3]]);
                    }

                    v4_mpi[ik_prod][(ns + 1) * is][(ns + 1) * js] = factor * ret;
                }
            }
        }
    }

    deallocate(v4_array_at_kpair);
    deallocate(ind);

// Now, communicate the calculated data.
// When the data count is larger than 2^31-1, split it.

    long maxsize = 1;
    maxsize = (maxsize << 31) - 1;

    const size_t count = nk2_prod * ns4;
    const size_t count_sub = ns4;

    if (count <= maxsize) {
#ifdef MPI_CXX_DOUBLE_COMPLEX
        MPI_Allreduce(&v4_mpi[0][0][0], &v4_out[0][0][0], count,
                      MPI_CXX_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
#else
        MPI_Allreduce(&v4_mpi[0][0][0], &v4_out[0][0][0], count,
                      MPI_COMPLEX16, MPI_SUM, MPI_COMM_WORLD);
#endif
    } else if (count_sub <= maxsize) {
        for (size_t ik_prod = 0; ik_prod < nk2_prod; ++ik_prod) {
#ifdef MPI_CXX_DOUBLE_COMPLEX
            MPI_Allreduce(&v4_mpi[ik_prod][0][0], &v4_out[ik_prod][0][0],
                          count_sub,
                          MPI_CXX_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
#else
            MPI_Allreduce(&v4_mpi[ik_prod][0][0], &v4_out[ik_prod][0][0],
                          count_sub,
                          MPI_COMPLEX16, MPI_SUM, MPI_COMM_WORLD);
#endif
        }
    } else {
        for (size_t ik_prod = 0; ik_prod < nk2_prod; ++ik_prod) {
            for (is = 0; is < ns2; ++is) {
#ifdef MPI_CXX_DOUBLE_COMPLEX
                MPI_Allreduce(&v4_mpi[ik_prod][is][0], &v4_out[ik_prod][is][0],
                              ns2,
                              MPI_CXX_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
#else
                MPI_Allreduce(&v4_mpi[ik_prod][is][0], &v4_out[ik_prod][is][0],
                              ns2,
                              MPI_COMPLEX16, MPI_SUM, MPI_COMM_WORLD);
#endif
            }
        }
    }

    deallocate(v4_mpi);

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

    size_t ik_prod;
    const size_t nk_reduced_interpolate = kmesh_coarse->nk_irred;
    const size_t ns = dynamical->neval;
    const size_t ns2 = ns * ns;
    const size_t ns4 = ns * ns * ns * ns;
    int is, js;
    unsigned int knum;
    unsigned int **ind;
    unsigned int i, j;
    long int *nset_mpi;

    const auto nk_scph = kmesh_dense->nk;
    const auto ngroup_v4 = anharmonic_core->get_ngroup_fcs(4);
    auto factor = std::pow(0.5, 2) / static_cast<double>(nk_scph);
    static auto complex_zero = std::complex<double>(0.0, 0.0);
    std::complex<double> *v4_array_at_kpair;
    std::complex<double> ***v4_mpi;

    std::vector<int> ik_vec, jk_vec, is_vec, js_vec;

    auto nk2_prod = nk_reduced_interpolate * nk_scph;

    if (mympi->my_rank == 0) {
        if (self_offdiag) {
            std::cout << " IALGO = 1 : Use different algorithm efficient when nbands >> nk_3ph\n";
            std::cout << " SELF_OFFDIAG = 1: Calculating all components of v4_array ... \n";
        } else {
            exit("compute_V4_elements_mpi_over_kpoint",
                 "This function can be used only when SELF_OFFDIAG = 1");
        }
    }

    allocate(nset_mpi, mympi->nprocs);

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
    deallocate(nset_mpi);

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

    allocate(v4_array_at_kpair, ngroup_v4);
    allocate(ind, ngroup_v4, 4);
    allocate(v4_mpi, nk2_prod, ns2, ns2);

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

            knum = kmap_interpolate_to_scph[kmesh_coarse->kpoint_irred_all[ik_now][0].knum];

            anharmonic_core->calc_phi4_reciprocal(kmesh_dense->xk[knum],
                                                  kmesh_dense->xk[jk_now],
                                                  kmesh_dense->xk[kmesh_dense->kindex_minus_xk[jk_now]],
                                                  phase_factor_scph,
                                                  phi4_reciprocal);

#ifdef _OPENMP
#pragma omp parallel for private(j)
#endif
            for (ii = 0; ii < ngroup_v4; ++ii) {
                v4_array_at_kpair[ii] = phi4_reciprocal[ii] * anharmonic_core->get_invmass_factor(4)[ii];
                for (j = 0; j < 4; ++j) ind[ii][j] = anharmonic_core->get_evec_index(4)[ii][j];
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

            for (i = 0; i < ngroup_v4; ++i) {

                ret += v4_array_at_kpair[i]
                       * std::conj(evec_in[knum][is_now][ind[i][0]])
                       * evec_in[knum][js_now][ind[i][1]]
                       * evec_in[jk_now][ks][ind[i][2]]
                       * std::conj(evec_in[jk_now][ls][ind[i][3]]);
            }

            v4_mpi[ik_prod][is_prod][js] = factor * ret;
        }

        if (mympi->my_rank == 0) {
            std::cout << " SET " << ii + 1 << " done. " << std::endl;
        }

    } // loop over nk2_prod*ns2

    deallocate(v4_array_at_kpair);
    deallocate(ind);

// Now, communicate the calculated data.
// When the data count is larger than 2^31-1, split it.

    long maxsize = 1;
    maxsize = (maxsize << 31) - 1;

    const size_t count = nk2_prod * ns4;
    const size_t count_sub = ns4;

    if (mympi->my_rank == 0) {
        std::cout << "Communicating v4_array over MPI ...";
    }
    if (count <= maxsize) {
#ifdef MPI_CXX_DOUBLE_COMPLEX
        MPI_Allreduce(&v4_mpi[0][0][0], &v4_out[0][0][0], count,
                      MPI_CXX_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
#else
        MPI_Allreduce(&v4_mpi[0][0][0], &v4_out[0][0][0], count,
                      MPI_COMPLEX16, MPI_SUM, MPI_COMM_WORLD);
#endif
    } else if (count_sub <= maxsize) {
        for (size_t ik_prod = 0; ik_prod < nk2_prod; ++ik_prod) {
#ifdef MPI_CXX_DOUBLE_COMPLEX
            MPI_Allreduce(&v4_mpi[ik_prod][0][0], &v4_out[ik_prod][0][0],
                          count_sub,
                          MPI_CXX_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
#else
            MPI_Allreduce(&v4_mpi[ik_prod][0][0], &v4_out[ik_prod][0][0],
                          count_sub,
                          MPI_COMPLEX16, MPI_SUM, MPI_COMM_WORLD);
#endif
        }
    } else {
        for (size_t ik_prod = 0; ik_prod < nk2_prod; ++ik_prod) {
            for (is = 0; is < ns2; ++is) {
#ifdef MPI_CXX_DOUBLE_COMPLEX
                MPI_Allreduce(&v4_mpi[ik_prod][is][0], &v4_out[ik_prod][is][0],
                              ns2,
                              MPI_CXX_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
#else
                MPI_Allreduce(&v4_mpi[ik_prod][is][0], &v4_out[ik_prod][is][0],
                              ns2,
                              MPI_COMPLEX16, MPI_SUM, MPI_COMM_WORLD);
#endif
            }
        }
    }
    if (mympi->my_rank == 0) {
        std::cout << "done.\n";
    }

    deallocate(v4_mpi);

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
    allocate(is_acoustic, ns);
    int nacoustic;
    auto threshould = 1.0e-24;
    const auto nk_reduced_interpolate = kmesh_coarse->nk_irred;
    static auto complex_zero = std::complex<double>(0.0, 0.0);

    if (!(fc_order == 3 || fc_order == 4)) {
        exit("zerofill_elements_acoustic_at_gamma",
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
            exit("zerofill_elements_acoustic_at_gamma",
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
        for (jk = 0; jk < kmesh_dense->nk; ++jk) {
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
                                v_elems[kmesh_dense->nk * ik][ns * is + js][ns * ks + ls] = complex_zero;
                            }
                        }
                    }
                }
            }
        }
        // ik = 0;
        for (jk = 0; jk < kmesh_dense->nk; ++jk) {
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

    deallocate(is_acoustic);
}

void Scph::setup_kmesh()
{
    unsigned int ik;
    unsigned int i;
    double xtmp[3];

    // Setup k points for SCPH equation
    MPI_Bcast(&kmesh_scph[0], 3, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    MPI_Bcast(&kmesh_interpolate[0], 3, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

    kmesh_coarse = new KpointMeshUniform(kmesh_interpolate);
    kmesh_dense = new KpointMeshUniform(kmesh_scph);
    kmesh_coarse->setup(symmetry->SymmList, system->rlavec_p);
    kmesh_dense->setup(symmetry->SymmList, system->rlavec_p);

    if (mympi->my_rank == 0) {
//        if (verbosity > 0) {
        std::cout << " Setting up the SCPH calculations ..." << std::endl << std::endl;
        std::cout << "  Gamma-centered uniform grid with the following mesh density:" << std::endl;
        std::cout << "  nk1:" << std::setw(5) << kmesh_scph[0] << std::endl;
        std::cout << "  nk2:" << std::setw(5) << kmesh_scph[1] << std::endl;
        std::cout << "  nk3:" << std::setw(5) << kmesh_scph[2] << std::endl;
        std::cout << std::endl;
        std::cout << "  Number of k points : " << kmesh_dense->nk << std::endl;
        std::cout << "  Number of irreducible k points : " << kmesh_dense->nk_irred << std::endl;
        std::cout << std::endl;
        std::cout << "  Fourier interpolation from reciprocal to real space" << std::endl;
        std::cout << "  will be performed with the following mesh density:" << std::endl;
        std::cout << "  nk1:" << std::setw(5) << kmesh_interpolate[0] << std::endl;
        std::cout << "  nk2:" << std::setw(5) << kmesh_interpolate[1] << std::endl;
        std::cout << "  nk3:" << std::setw(5) << kmesh_interpolate[2] << std::endl;
        std::cout << std::endl;
        std::cout << "  Number of k points : " << kmesh_coarse->nk << std::endl;
        std::cout << "  Number of irreducible k points : "
                  << kmesh_coarse->nk_irred << std::endl;
//        }
    }

    allocate(kmap_interpolate_to_scph, kmesh_coarse->nk);

    for (ik = 0; ik < kmesh_coarse->nk; ++ik) {
        for (i = 0; i < 3; ++i) xtmp[i] = kmesh_coarse->xk[ik][i];

        const auto loc = kmesh_dense->get_knum(xtmp);

        if (loc == -1)
            exit("setup_kmesh",
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
    const auto nk_irred_interpolate = kmesh_coarse->nk_irred;

    allocate(gamma_tmp, ns, ns);
    allocate(mat_transform_sym, nk_irred_interpolate,
             symmetry->nsym, ns, ns);
    //allocate(small_group_at_k, nk_irred_interpolate);
    allocate(symop_minus_at_k, nk_irred_interpolate);
    allocate(kpoint_map_symmetry, kmesh_coarse->nk);
    allocate(flag, kmesh_coarse->nk);

    for (ik = 0; ik < kmesh_coarse->nk; ++ik) {
        flag[ik] = false;
    }

    for (ik = 0; ik < nk_irred_interpolate; ++ik) {

        symop_minus_at_k[ik].clear();

        const auto knum = kmesh_coarse->kpoint_irred_all[ik][0].knum;
        for (icrd = 0; icrd < 3; ++icrd) {
            k[icrd] = kmesh_coarse->xk[knum][icrd];
            k_minus[icrd] = -k[icrd];
        }
        const auto knum_minus = kmesh_coarse->kindex_minus_xk[knum];

        unsigned int isym = 0;

        for (const auto &it: symmetry->SymmListWithMap) {

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

            const auto knum_sym = kmesh_coarse->get_knum(Sk);
            if (knum_sym == -1)
                exit("setup_transform_symmetry",
                     "kpoint not found");

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

    deallocate(gamma_tmp);
    deallocate(flag);
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

    const auto nsym_small = kmesh_coarse->small_group_of_k[ik].size();
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
        isym = kmesh_coarse->small_group_of_k[ik][i];

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

    allocate(dymat_all, ns, ns, kmesh_coarse->nk);

    for (i = 0; i < kmesh_coarse->nk; ++i) {

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
            for (i = 0; i < kmesh_coarse->nk; ++i) {
                dymat_inout[is][js][i] = dymat_all[is][js][i];
            }
        }
    }
    deallocate(dymat_all);
}

void Scph::setup_eigvecs()
{
    const auto ns = dynamical->neval;

    if (mympi->my_rank == 0) {
        std::cout << std::endl
                  << " Diagonalizing dynamical matrices for all k points ... ";
    }

    allocate(evec_harmonic, kmesh_dense->nk, ns, ns);
    allocate(omega2_harmonic, kmesh_dense->nk, ns);

    // Calculate phonon eigenvalues and eigenvectors for all k-points for scph

#pragma omp parallel for
    for (int ik = 0; ik < kmesh_dense->nk; ++ik) {

        dynamical->eval_k(kmesh_dense->xk[ik],
                          kmesh_dense->kvec_na[ik],
                          fcs_phonon->fc2_ext,
                          omega2_harmonic[ik],
                          evec_harmonic[ik], true);
    }

    if (mympi->my_rank == 0) {
        std::cout << "done !" << std::endl;
    }
}

void Scph::setup_pp_interaction()
{
    // Prepare information for calculating ph-ph interaction coefficients.

    if (mympi->my_rank == 0) {
        if (relax_coordinate || bubble > 0) {
            std::cout << " Preparing for calculating V3 & V4  ...";
        } else {
            std::cout << " Preparing for calculating V4  ...";
        }
    }

    if (anharmonic_core->quartic_mode != 1) {
        exit("setup_pp_interaction",
             "quartic_mode should be 1 for SCPH");
    }

    // Setup for V3 if relax_coordinate = True.
    if (relax_coordinate || bubble > 0) {
        allocate(phi3_reciprocal, anharmonic_core->get_ngroup_fcs(3));
    }
    allocate(phi4_reciprocal, anharmonic_core->get_ngroup_fcs(4));

    phase_factor_scph = new PhaseFactorStorage(kmesh_dense->nk_i);
    phase_factor_scph->create(true);

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
    const auto nk = kmesh_coarse->nk;
    const auto nat = system->natmin;
    unsigned int iat;

    int **shift_cell, **shift_cell_super;
    double **xf_p;
    double ****x_all;

    const int nkx = static_cast<int>(kmesh_coarse->nk_i[0]); // This should be int (must not be unsigned int).
    const int nky = static_cast<int>(kmesh_coarse->nk_i[1]); // same as above
    const int nkz = static_cast<int>(kmesh_coarse->nk_i[2]); // same as above

    const auto ncell = nk;
    const auto ncell_s = 27;

    allocate(shift_cell, ncell, 3);
    allocate(shift_cell_super, ncell_s, 3);
    allocate(xf_p, nat, 3);
    allocate(x_all, ncell_s, ncell, nat, 3);

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

    allocate(mindist_list_scph, nat, nat, ncell);

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

    deallocate(shift_cell);
    deallocate(shift_cell_super);
    deallocate(xf_p);
    deallocate(x_all);
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

    allocate(mat_tmp, ns, ns);
    allocate(eval_real, ns);
    allocate(mat_harmonic, ns, ns);

    if (dynamical->nonanalytic) {
        allocate(mat_harmonic_na, ns, ns);
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

    deallocate(eval_real);
    deallocate(mat_tmp);
    deallocate(mat_harmonic);

    if (dynamical->nonanalytic) {
        deallocate(mat_harmonic_na);
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
                for (const auto &it: mindist_list_scph[iat][jat][icell].shift) {

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
    allocate(RWORK, 3 * ns - 2);
    allocate(WORK, LWORK);

    if (require_evec) {
        JOBZ = 'V';
    } else {
        JOBZ = 'N';
    }

    char UPLO = 'U';

    allocate(amat, ns * ns);

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

    deallocate(amat);
    deallocate(WORK);
    deallocate(RWORK);
}

void Scph::find_degeneracy(std::vector<int> *degeneracy_out,
                           const unsigned int nk_in,
                           double **eval_in) const
{
    // eval is omega^2 in atomic unit

    const auto ns = dynamical->neval;
    const auto tol_omega = 1.0e-7;

    for (unsigned int ik = 0; ik < nk_in; ++ik) {

        degeneracy_out[ik].clear();

        auto omega_prev = eval_in[ik][0];
        auto ideg = 1;

        for (unsigned int is = 1; is < ns; ++is) {
            const auto omega_now = eval_in[ik][is];

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

    allocate(polarization_matrix, ns2);
    allocate(mat_tmp, ns2);
    allocate(eigval_matrix, ns2);
    allocate(beta, ns);
    allocate(dmat, ns2);
    allocate(dymat_q, ns, ns, kmesh_coarse->nk);
    allocate(dymat_harmonic, ns, ns);

    for (is = 0; is < ns; ++is) beta[is] = std::complex<double>(0.0, 0.0);

    for (ik = 0; ik < kmesh_coarse->nk; ++ik) {

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
        dynamical->calc_analytic_k(kmesh_coarse->xk[ik],
                                   fcs_phonon->fc2_ext,
                                   dymat_harmonic);

        for (is = 0; is < ns; ++is) {
            for (js = 0; js < ns; ++js) {
                dymat_q[is][js][ik] -= dymat_harmonic[is][js];
            }
        }
    }

    deallocate(polarization_matrix);
    deallocate(mat_tmp);
    deallocate(eigval_matrix);
    deallocate(beta);
    deallocate(dmat);
    deallocate(dymat_harmonic);

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

                for (ik = 0; ik < kmesh_coarse->nk; ++ik) {

                    duplicate_xk_boundary(kmesh_coarse->xk[ik], xk_dup);

                    auto cexp_phase = std::complex<double>(0.0, 0.0);

                    for (const auto &i: xk_dup) {

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
                        dymat_out[is][js][icell] /= static_cast<double>(kmesh_coarse->nk);
                    }
                }

                ++icell;
            }
        }
    }

    deallocate(dymat_q);
}

void Scph::compute_anharmonic_frequency(std::complex<double> ***v4_array_all,
                                        double **omega2_out,
                                        std::complex<double> ***evec_anharm_scph,
                                        const double temp,
                                        bool &flag_converged,
                                        std::complex<double> ***cmat_convert,
                                        const bool offdiag,
                                        const unsigned int verbosity)
{
    // This is the main function of the SCPH equation.
    // The detailed algorithm can be found in PRB 92, 054301 (2015).
    // Eigen3 library is used for the compact notation of matrix-matrix products.

    using namespace Eigen;

    int ik, jk;
    unsigned int i;
    unsigned int is, js, ks;
    unsigned int kk;
    const auto nk = kmesh_dense->nk;
    const auto ns = dynamical->neval;
    unsigned int knum, knum_interpolate;
    const auto nk_interpolate = kmesh_coarse->nk;
    const auto nk_irred_interpolate = kmesh_coarse->nk_irred;
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

    allocate(mat_omega2_harmonic, nk_interpolate, ns, ns);
    allocate(eval_interpolate, nk, ns);
    allocate(evec_initial, nk, ns, ns);
    allocate(evec_new, nk, ns, ns);
    allocate(dmat_convert, nk, ns, ns);
    allocate(dmat_convert_old, nk, ns, ns);
    allocate(dymat_new, ns, ns, nk_interpolate);
    allocate(dymat_q, ns, ns, nk_interpolate);
    allocate(dymat_harmonic, nk_interpolate, ns, ns);
    allocate(Fmat0, nk_irred_interpolate, ns, ns);

    const auto T_in = temp;

    std::cout << " Temperature = " << T_in << " K" << std::endl;

    // Set initial values

    for (ik = 0; ik < nk; ++ik) {
        for (is = 0; is < ns; ++is) {

            if (flag_converged) {
                if (omega2_out[ik][is] < 0.0 && std::abs(omega2_out[ik][is]) > 1.0e-16 && verbosity > 0) {
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

        dynamical->calc_analytic_k(kmesh_coarse->xk[ik],
                                   fcs_phonon->fc2_ext,
                                   dymat_harmonic[ik]);
    }

    for (ik = 0; ik < nk_irred_interpolate; ++ik) {

        knum_interpolate = kmesh_coarse->kpoint_irred_all[ik][0].knum;

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

            knum_interpolate = kmesh_coarse->kpoint_irred_all[ik][0].knum;
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

#pragma omp parallel for private(jk, kk, ks, ctmp), reduction(+:re_tmp, im_tmp)
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

#pragma omp parallel for private(jk, kk, ks, ctmp), reduction(+:re_tmp, im_tmp)
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

                    if (verbosity > 1) {
                        std::cout << " Detect imaginary : ";
                        std::cout << "  knum = " << knum + 1 << " is = " << is + 1 << '\n';
                        for (int j = 0; j < 3; ++j) {
                            std::cout << "  xk = " << std::setw(15) << kmesh_dense->xk[knum][j];
                        }
                        std::cout << '\n';
                    }

                    if (v4_array_all[nk * ik + knum][(ns + 1) * is][(ns + 1) * is].real() > 0.0) {
                        if (verbosity > 1) {
                            std::cout << "  onsite V4 is positive\n\n";
                        }

                        if (flag_converged) {
                            ++icount;
                            eval_tmp(is) = omega2_out[knum][is] * std::pow(0.99, icount);
                        } else {
                            ++icount;
                            eval_tmp(is) = -eval_tmp(is) * std::pow(0.99, icount);
                        }
                    } else {
                        if (verbosity > 1) {
                            std::cout << "  onsite V4 is negative\n\n";
                        }
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

#ifdef _DEBUG2
            Dymat_sym = Dymat;
            symmetrize_dynamical_matrix(ik, Dymat_sym);
            std::complex<double> **dymat_exact;
            allocate(dymat_exact, ns, ns);
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
            deallocate(dymat_exact);

#endif
            symmetrize_dynamical_matrix(ik, Dymat);

            for (is = 0; is < ns; ++is) {
                for (js = 0; js < ns; ++js) {
                    dymat_q[is][js][knum_interpolate] = Dymat(is, js);
                }
            }
        } // close loop ik

        replicate_dymat_for_all_kpoints(dymat_q);

#ifdef _DEBUG2
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
                                                  reinterpret_cast<fftw_complex *>(dymat_q[is][js]),
                                                  reinterpret_cast<fftw_complex *>(dymat_new[is][js]),
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
                           kmesh_dense->xk,
                           kmesh_dense->kvec_na,
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
            if (verbosity > 0) {
                std::cout << "  SCPH ITER " << std::setw(5) << iloop + 1 << " :  DIFF = N/A" << std::endl;
            }

            for (ik = 0; ik < nk; ++ik) {
                for (is = 0; is < ns; ++is) {
                    omega_old(ik, is) = omega_now(ik, is);
                }
            }

        } else {

            diff = 0.0;

            for (ik = 0; ik < nk_interpolate; ++ik) {
                knum = kmap_interpolate_to_scph[ik];
                for (is = 0; is < ns; ++is) {
                    diff += std::pow(omega_now(knum, is) - omega_old(knum, is), 2.0);
                }
            }
            diff /= static_cast<double>(nk_interpolate * ns);
            if (verbosity > 0) {
                std::cout << "  SCPH ITER " << std::setw(5) << iloop + 1 << " : ";
                std::cout << " DIFF = " << std::setw(15) << std::sqrt(diff) << std::endl;
            }
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
                    if (verbosity > 0) std::cout << "  DIFF < SCPH_TOL : break SCPH loop\n";
                    break;
                }
                if (verbosity > 0) std::cout << "  DIFF < SCPH_TOL but a negative frequency is detected.\n";
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
        if (verbosity > 0) {
            std::cout << " Temp = " << T_in;
            std::cout << " : convergence achieved in " << std::setw(5)
                      << iloop + 1 << " iterations." << std::endl;
        }
        flag_converged = true;
    } else {
        if (verbosity > 0) {
            std::cout << "Temp = " << T_in;
            std::cout << " : not converged." << std::endl;
        }
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

    if (verbosity > 1) {
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
    }

    deallocate(mat_omega2_harmonic);
    deallocate(eval_interpolate);
    deallocate(evec_initial);
    deallocate(dmat_convert);
    deallocate(dmat_convert_old);
    deallocate(evec_new);
    deallocate(dymat_new);
    deallocate(dymat_q);
    deallocate(dymat_harmonic);
    deallocate(Fmat0);
}

void Scph::compute_free_energy_bubble_SCPH(const unsigned int kmesh[3],
                                           std::complex<double> ****delta_dymat_scph)
{
    const auto NT = static_cast<unsigned int>((system->Tmax - system->Tmin) / system->dT) + 1;
    const auto nk_ref = dos->kmesh_dos->nk;
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
                     "  the number of temperature grids.\n\n";
    }

    allocate(thermodynamics->FE_bubble, NT);
    allocate(eval, NT, nk_ref, ns);
    allocate(evec, NT, nk_ref, ns, ns); // This requires lots of RAM

    for (auto iT = 0; iT < NT; ++iT) {
        exec_interpolation(kmesh,
                           delta_dymat_scph[iT],
                           nk_ref,
                           dos->kmesh_dos->xk,
                           dos->kmesh_dos->kvec_na,
                           eval[iT],
                           evec[iT]);
    }

    thermodynamics->compute_FE_bubble_SCPH(eval, evec, thermodynamics->FE_bubble);

    deallocate(eval);
    deallocate(evec);

    if (mympi->my_rank == 0) {
        std::cout << " done!" << std::endl << std::endl;
    }
}

void Scph::bubble_correction(std::complex<double> ****delta_dymat_scph,
                             std::complex<double> ****delta_dymat_scph_plus_bubble)
{
    const auto NT = static_cast<unsigned int>((system->Tmax - system->Tmin) / system->dT) + 1;
    const auto ns = dynamical->neval;

    auto epsilon = integration->epsilon;
    const auto nk_irred_interpolate = kmesh_coarse->nk_irred;
    const auto nk_scph = kmesh_dense->nk;

    double **eval = nullptr;
    double ***eval_bubble = nullptr;
    std::complex<double> ***evec;
    double *real_self = nullptr;
    std::vector<std::complex<double>> omegalist;

    if (mympi->my_rank == 0) {
        std::cout << std::endl;
        std::cout << " -----------------------------------------------------------------"
                  << std::endl;
        std::cout << " Calculating the bubble self-energy " << std::endl;
        std::cout << " on top of the SCPH calculation." << std::endl;
        std::cout << '\n';
    }

    allocate(eval, nk_scph, ns);
    allocate(evec, nk_scph, ns, ns);

    if (mympi->my_rank == 0) {
        allocate(eval_bubble, NT, nk_scph, ns);
        for (auto iT = 0; iT < NT; ++iT) {
            for (auto ik = 0; ik < nk_scph; ++ik) {
                for (auto is = 0; is < ns; ++is) {
                    eval_bubble[iT][ik][is] = 0.0;
                }
            }
        }
        allocate(real_self, ns);
    }

    std::vector<int> *degeneracy_at_k;
    allocate(degeneracy_at_k, nk_scph);

    for (auto iT = 0; iT < NT; ++iT) {
        const auto temp = system->Tmin + system->dT * float(iT);

        exec_interpolation(kmesh_interpolate,
                           delta_dymat_scph[iT],
                           nk_scph,
                           kmesh_dense->xk,
                           kmesh_dense->kvec_na,
                           eval,
                           evec);

        find_degeneracy(degeneracy_at_k,
                        nk_scph,
                        eval);

        if (mympi->my_rank == 0) std::cout << " Temperature (K) : " << std::setw(6) << temp << '\n';

        for (auto ik = 0; ik < nk_irred_interpolate; ++ik) {

            auto knum_interpolate = kmesh_coarse->kpoint_irred_all[ik][0].knum;
            auto knum = kmap_interpolate_to_scph[knum_interpolate];

            if (mympi->my_rank == 0) {
                std::cout << "  Irred. k: " << std::setw(5) << ik + 1 << " (";
                for (auto m = 0; m < 3; ++m) std::cout << std::setw(15) << kmesh_dense->xk[knum][m];
                std::cout << ")\n";
            }

            for (unsigned int snum = 0; snum < ns; ++snum) {

                if (eval[knum][snum] < eps8) {
                    if (mympi->my_rank == 0) real_self[snum] = 0.0;
                } else {
                    omegalist.clear();

                    if (bubble == 1) {

                        omegalist.push_back(im * epsilon);

                        auto se_bubble = get_bubble_selfenergy(kmesh_dense,
                                                               ns,
                                                               eval,
                                                               evec,
                                                               knum,
                                                               snum,
                                                               temp,
                                                               omegalist);

                        if (mympi->my_rank == 0) real_self[snum] = se_bubble[0].real();

                    } else if (bubble == 2) {

                        omegalist.push_back(eval[knum][snum] + im * epsilon);

                        auto se_bubble = get_bubble_selfenergy(kmesh_dense,
                                                               ns,
                                                               eval,
                                                               evec,
                                                               knum,
                                                               snum,
                                                               temp,
                                                               omegalist);

                        if (mympi->my_rank == 0) real_self[snum] = se_bubble[0].real();

                    } else if (bubble == 3) {

                        auto maxfreq = eval[knum][snum] + 50.0 * time_ry / Hz_to_kayser;
                        auto minfreq = eval[knum][snum] - 50.0 * time_ry / Hz_to_kayser;

                        if (minfreq < 0.0) minfreq = 0.0;

                        const auto domega = 0.1 * time_ry / Hz_to_kayser;
                        auto nomega = static_cast<unsigned int>((maxfreq - minfreq) / domega) + 1;

                        for (auto iomega = 0; iomega < nomega; ++iomega) {
                            omegalist.push_back(minfreq + static_cast<double>(iomega) * domega + im * epsilon);
                        }

                        auto se_bubble = get_bubble_selfenergy(kmesh_dense,
                                                               ns,
                                                               eval,
                                                               evec,
                                                               knum,
                                                               snum,
                                                               temp,
                                                               omegalist);

                        if (mympi->my_rank == 0) {

                            std::vector<double> nonlinear_func(nomega);
                            for (auto iomega = 0; iomega < nomega; ++iomega) {
                                nonlinear_func[iomega] = omegalist[iomega].real() * omegalist[iomega].real()
                                                         - eval[knum][snum] * eval[knum][snum]
                                                         + 2.0 * eval[knum][snum] * se_bubble[iomega].real();
                            }

                            // find a root of nonlinear_func = 0 from the sign change.
                            int count_root = 0;
                            std::vector<unsigned int> root_index;

                            for (auto iomega = 0; iomega < nomega - 1; ++iomega) {
                                if (nonlinear_func[iomega] * nonlinear_func[iomega + 1] < 0.0) {
                                    ++count_root;
                                    root_index.push_back(iomega);
                                }
                            }

                            if (count_root == 0) {
                                warn("bubble_correction",
                                     "Could not find a root in the nonlinear equation at this temperature. "
                                     "Use the w=0 component.");

                                real_self[snum] = se_bubble[0].real();

                            } else {
                                if (count_root > 1) {
                                    warn("bubble_correction",
                                         "Multiple roots were found in the nonlinear equation at this temperature. "
                                         "Use the lowest-frequency solution");
                                    std::cout << "   solution found at the following frequencies:\n";
                                    for (auto iroot = 0; iroot < count_root; ++iroot) {
                                        std::cout << std::setw(15)
                                                  << writes->in_kayser(omegalist[root_index[iroot]].real());
                                    }
                                    std::cout << '\n';
                                }

                                // Instead of performing a linear interpolation (secant method) of nonlinear_func,
                                // we interpolate the bubble self-energy. Since the frequency grid is dense (0.1 cm^-1 step),
                                // this approximation should not make any problems (hopefully).

                                double omega_solution = omegalist[root_index[0] + 1].real()
                                                        - nonlinear_func[root_index[0] + 1]
                                                          * domega / (nonlinear_func[root_index[0] + 1] -
                                                                      nonlinear_func[root_index[0]]);

                                real_self[snum] = (se_bubble[root_index[0] + 1].real()
                                                   - se_bubble[root_index[0]].real())
                                                  * (omega_solution - omegalist[root_index[0] + 1].real()) / domega
                                                  + se_bubble[root_index[0] + 1].real();
                            }
                        }
                    }
                }
                if (mympi->my_rank == 0) {
                    std::cout << "   branch : " << std::setw(5) << snum + 1;
                    std::cout << " omega (SC1) = " << std::setw(15) << writes->in_kayser(eval[knum][snum])
                              << " (cm^-1); ";
                    std::cout << " Re[Self] = " << std::setw(15) << writes->in_kayser(real_self[snum]) << " (cm^-1)\n";
                }
            }

            if (mympi->my_rank == 0) {
                // average self energy of degenerate modes
                int ishift = 0;
                double real_self_avg = 0.0;

                for (const auto &it: degeneracy_at_k[knum]) {
                    for (auto m = 0; m < it; ++m) {
                        real_self_avg += real_self[m + ishift];
                    }
                    real_self_avg /= static_cast<double>(it);

                    for (auto m = 0; m < it; ++m) {
                        real_self[m + ishift] = real_self_avg;
                    }
                    real_self_avg = 0.0;
                    ishift += it;
                }

                for (unsigned int snum = 0; snum < ns; ++snum) {
                    eval_bubble[iT][knum][snum] = eval[knum][snum] * eval[knum][snum]
                                                  - 2.0 * eval[knum][snum] * real_self[snum];
                    for (auto jk = 1; jk < kmesh_coarse->kpoint_irred_all[ik].size(); ++jk) {
                        auto knum2 = kmap_interpolate_to_scph[kmesh_coarse->kpoint_irred_all[ik][jk].knum];
                        eval_bubble[iT][knum2][snum] = eval_bubble[iT][knum][snum];
                    }
                }

                std::cout << '\n';
            }
        }

        if (mympi->my_rank == 0) {
            calc_new_dymat_with_evec(delta_dymat_scph_plus_bubble[iT],
                                     eval_bubble[iT],
                                     evec);
        }
    }

    deallocate(eval);
    deallocate(evec);
    deallocate(degeneracy_at_k);

    if (eval_bubble) deallocate(eval_bubble);

    if (mympi->my_rank == 0) {
        std::cout << " done!" << std::endl << std::endl;
    }
}

std::vector<std::complex<double>> Scph::get_bubble_selfenergy(const KpointMeshUniform *kmesh_in,
                                                              const unsigned int ns_in,
                                                              const double *const *eval_in,
                                                              const std::complex<double> *const *const *evec_in,
                                                              const unsigned int knum,
                                                              const unsigned int snum,
                                                              const double temp_in,
                                                              const std::vector<std::complex<double>> &omegalist)
{
    unsigned int arr_cubic[3];
    double xk_tmp[3];
    std::complex<double> omega_sum[2];

    double factor = 1.0 / (static_cast<double>(kmesh_in->nk) * std::pow(2.0, 4));
    const auto ns2 = ns_in * ns_in;
    const auto nks = kmesh_in->nk * ns2;

    double n1, n2;
    double f1, f2;

    auto knum_minus = kmesh_in->kindex_minus_xk[knum];
    arr_cubic[0] = ns_in * knum_minus + snum;

    std::vector<std::complex<double>> se_bubble(omegalist.size());

    const auto nomega = omegalist.size();

    std::complex<double> *ret_sum, *ret_mpi;
    allocate(ret_sum, nomega);
    allocate(ret_mpi, nomega);

    for (auto iomega = 0; iomega < nomega; ++iomega) {
        ret_sum[iomega] = std::complex<double>(0.0, 0.0);
        ret_mpi[iomega] = std::complex<double>(0.0, 0.0);
    }

    for (auto iks = mympi->my_rank; iks < nks; iks += mympi->nprocs) {

        auto ik1 = iks / ns2;
        auto is1 = (iks % ns2) / ns_in;
        auto is2 = iks % ns_in;

        for (auto m = 0; m < 3; ++m) xk_tmp[m] = kmesh_in->xk[knum][m] - kmesh_in->xk[ik1][m];
        auto ik2 = kmesh_in->get_knum(xk_tmp);

        double omega1 = eval_in[ik1][is1];
        double omega2 = eval_in[ik2][is2];

        arr_cubic[1] = ns_in * ik1 + is1;
        arr_cubic[2] = ns_in * ik2 + is2;

        double v3_tmp = std::norm(anharmonic_core->V3(arr_cubic,
                                                      kmesh_in->xk,
                                                      eval_in,
                                                      evec_in,
                                                      phase_factor_scph));

        if (thermodynamics->classical) {
            n1 = thermodynamics->fC(omega1, temp_in);
            n2 = thermodynamics->fC(omega2, temp_in);
            f1 = n1 + n2;
            f2 = n2 - n1;
        } else {
            n1 = thermodynamics->fB(omega1, temp_in);
            n2 = thermodynamics->fB(omega2, temp_in);
            f1 = n1 + n2 + 1.0;
            f2 = n2 - n1;
        }
        for (auto iomega = 0; iomega < nomega; ++iomega) {
            omega_sum[0] = 1.0 / (omegalist[iomega] + omega1 + omega2) - 1.0 / (omegalist[iomega] - omega1 - omega2);
            omega_sum[1] = 1.0 / (omegalist[iomega] + omega1 - omega2) - 1.0 / (omegalist[iomega] - omega1 + omega2);
            ret_mpi[iomega] += v3_tmp * (f1 * omega_sum[0] + f2 * omega_sum[1]);
        }
    }
    for (auto iomega = 0; iomega < nomega; ++iomega) {
        ret_mpi[iomega] *= factor;
    }

    MPI_Reduce(&ret_mpi[0], &ret_sum[0], nomega, MPI_COMPLEX16, MPI_SUM, 0, MPI_COMM_WORLD);

    for (auto iomega = 0; iomega < nomega; ++iomega) {
        se_bubble[iomega] = ret_sum[iomega];
    }

    deallocate(ret_mpi);
    deallocate(ret_sum);

    return se_bubble;
}

double Scph::distance(double *x1,
                      double *x2)
{
    auto dist = std::pow(x1[0] - x2[0], 2)
                + std::pow(x1[1] - x2[1], 2)
                + std::pow(x1[2] - x2[2], 2);
    dist = std::sqrt(dist);

    return dist;
}

void Scph::duplicate_xk_boundary(double *xk_in,
                                 std::vector<std::vector<double>> &vec_xk)
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
                                           const unsigned int NT,
                                           const int type)
{
    unsigned int i, j;
    const auto Tmin = system->Tmin;
    const auto dT = system->dT;
    double ***delta_fc2;
    double **xtmp;
    const auto ns = dynamical->neval;
    unsigned int is, js, icell;
    unsigned int iat, jat;

    std::string file_fc2;
    std::ofstream ofs_fc2;

    if (type == 0) {
        file_fc2 = input->job_title + ".scph_dfc2";
    } else if (type == 1) {
        file_fc2 = input->job_title + ".scph+bubble(0)_dfc2";
    } else if (type == 2) {
        file_fc2 = input->job_title + ".scph+bubble(w)_dfc2";
    } else if (type == 3) {
        file_fc2 = input->job_title + ".scph+bubble(wQP)_dfc2";
    }

    ofs_fc2.open(file_fc2.c_str(), std::ios::out);
    if (!ofs_fc2)
        exit("write_anharmonic_correction_fc2",
             "Cannot open file_fc2");

    const auto ncell = kmesh_interpolate[0] * kmesh_interpolate[1] * kmesh_interpolate[2];

    allocate(delta_fc2, ns, ns, ncell);

    allocate(xtmp, system->natmin, 3);

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
            ofs_fc2 << std::setw(20) << xtmp[i][j];
        }
        ofs_fc2 << std::setw(5) << system->kd[system->map_p2s[i][0]] + 1 << std::endl;
    }

    deallocate(xtmp);

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

    deallocate(delta_fc2);

    ofs_fc2.close();
    std::cout << "  " << std::setw(input->job_title.length() + 12) << std::left << file_fc2;

    if (type == 0) {
        std::cout << " : Anharmonic corrections to the second-order IFCs (SCPH)" << std::endl;
    } else if (type == 1) {
        std::cout << " : Anharmonic corrections to the second-order IFCs (SCPH+Bubble(0))" << std::endl;
    } else if (type == 2) {
        std::cout << " : Anharmonic corrections to the second-order IFCs (SCPH+Bubble(w))" << std::endl;
    } else if (type == 3) {
        std::cout << " : Anharmonic corrections to the second-order IFCs (SCPH+Bubble(wQP))" << std::endl;
    }
}

void Scph::mpi_bcast_complex(std::complex<double> ****data,
                             const unsigned int NT,
                             const unsigned int nk,
                             const unsigned int ns)
{
    const int _NT = static_cast<int>(NT);
    const int _nk = static_cast<int>(nk);
    const int _ns = static_cast<int>(ns);

#ifdef MPI_CXX_DOUBLE_COMPLEX
    MPI_Bcast(&data[0][0][0][0], _NT * _nk * _ns * _ns,
              MPI_CXX_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
#elif defined MPI_DOUBLE_COMPLEX
    MPI_Bcast(&data[0][0][0][0], _NT * _nk * _ns * _ns, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);
#else
    MPI_Bcast(&data[0][0][0][0], _NT * _nk * _ns * _ns, MPI_COMPLEX16, 0, MPI_COMM_WORLD);
#endif
}

void Scph::get_derivative_central_diff(const double delta_t,
                                       const unsigned int nk,
                                       double **omega0,
                                       double **omega2,
                                       double **domega_dt)
{
    const auto ns = dynamical->neval;
    const auto inv_dt = 1.0 / (2.0 * delta_t);
    for (auto ik = 0; ik < nk; ++ik) {
        for (auto is = 0; is < ns; ++is) {
            domega_dt[ik][is] = (omega2[ik][is] - omega0[ik][is]) * inv_dt;
            //    std::cout << "domega_dt = " << domega_dt[ik][is] << '\n';
        }
    }
}
