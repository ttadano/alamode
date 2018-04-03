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
};


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
    eval_harmonic = nullptr;
    evec_harmonic = nullptr;
    fcs_group = nullptr;
    fcs_group2 = nullptr;
    knum_minus_scph = nullptr;
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
    if (knum_minus_scph) {
        memory->deallocate(knum_minus_scph);
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
    if (eval_harmonic) {
        memory->deallocate(eval_harmonic);
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
    unsigned int nk_ref = kpoint->nk;
    unsigned int ns = dynamical->neval;
    double Tmin = system->Tmin;
    double Tmax = system->Tmax;
    double dT = system->dT;

    double ***eval_anharm = nullptr;
    std::complex<double> ****evec_anharm = nullptr;
    std::complex<double> ****delta_dymat_scph = nullptr;

    auto NT = static_cast<unsigned int>((Tmax - Tmin) / dT) + 1;

    MPI_Bcast(&restart_scph, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD);
    MPI_Bcast(&selfenergy_offdiagonal, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD);
    MPI_Bcast(&ialgo, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

    //   if (mympi->my_rank == 0) {
    memory->allocate(delta_dymat_scph, NT, ns, ns, nk_interpolate);
    //   }

    if (restart_scph) {

        // Read anharmonic correction to the dynamical matrix from the existing file
        load_scph_dymat_from_file(delta_dymat_scph);

    } else {

        // Solve the SCPH equation and obtain the correction to the dynamical matrix
        exec_scph_main(delta_dymat_scph);

        if (mympi->my_rank == 0) {
            store_scph_dymat_to_file(delta_dymat_scph);
            write_anharmonic_correction_fc2(delta_dymat_scph, NT);
        }
    }

    //   if (mympi->my_rank == 0) {
    memory->allocate(eval_anharm, NT, nk_ref, ns);
    memory->allocate(evec_anharm, NT, nk_ref, ns, ns); // This requires lots of RAM

    for (auto iT = 0; iT < NT; ++iT) {
        exec_interpolation(delta_dymat_scph[iT],
                           eval_anharm[iT],
                           evec_anharm[iT]);
    }
    //    }
    //   if (delta_dymat_scph) {
    memory->deallocate(delta_dymat_scph);
    //   }

    if (kpoint->kpoint_mode == 2) {
        if (thermodynamics->calc_FE_bubble) {

            /*
            if (mympi->my_rank > 0) {
                memory->allocate(eval_anharm, NT, nk_ref, ns);
                memory->allocate(evec_anharm, NT, nk_ref, ns, ns); // Memory intensive
            }

            MPI_Bcast(&eval_anharm[0][0][0], NT * nk_ref * ns, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            mpi_bcast_complex(evec_anharm, NT, nk_ref, ns);
            */
            thermodynamics->compute_free_energy_bubble_SCPH(eval_anharm,
                                                            evec_anharm);
        }
    }

    if (mympi->my_rank == 0) {
        if (kpoint->kpoint_mode == 0) {
            write_scph_energy(eval_anharm);
        } else if (kpoint->kpoint_mode == 1) {
            write_scph_bands(eval_anharm);
        } else if (kpoint->kpoint_mode == 2) {
            if (dos->compute_dos) {
                write_scph_dos(eval_anharm);
            }
            write_scph_thermodynamics(eval_anharm, evec_anharm);
            if (writes->print_msd) {
                write_scph_msd(eval_anharm, evec_anharm);
            }
        }
    }

    if (eval_anharm) {
        memory->deallocate(eval_anharm);
    }
    if (evec_anharm) {
        memory->deallocate(evec_anharm);
    }
}


void Scph::load_scph_dymat_from_file(std::complex<double> ****dymat_out)
{
    unsigned int ns = dynamical->neval;
    unsigned int NT;
    double Tmin = system->Tmin;
    double Tmax = system->Tmax;
    double dT = system->dT;
    NT = static_cast<unsigned int>((Tmax - Tmin) / dT) + 1;
    std::vector<double> Temp_array(NT);

    for (int i = 0; i < NT; ++i) {
        Temp_array[i] = Tmin + dT * static_cast<double>(i);
    }

    if (mympi->my_rank == 0) {

        int iT, ik, is, js;
        unsigned int nk = nk_scph;
        unsigned int nk_ref = kpoint->nk;
        bool consider_offdiagonal = selfenergy_offdiagonal;
        double temp;
        std::ifstream ifs_dymat;
        std::string file_dymat = input->job_title + ".scph_dymat";
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
            error->exit("load_scph_dymat_from_file",
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
        for (iT = 0; iT < NT_ref; ++iT) {
            ifs_dymat >> str_dummy >> temp;
            for (is = 0; is < ns; ++is) {
                for (js = 0; js < ns; ++js) {
                    for (ik = 0; ik < nk_interpolate; ++ik) {
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
    int iT, ik, is, js;
    unsigned int ns = dynamical->neval;
    double Tmin = system->Tmin;
    double Tmax = system->Tmax;
    double dT = system->dT;
    double temp;
    unsigned int NT;
    std::ofstream ofs_dymat;
    std::string file_dymat = input->job_title + ".scph_dymat";

    NT = static_cast<unsigned int>((Tmax - Tmin) / dT) + 1;

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

    for (iT = 0; iT < NT; ++iT) {
        temp = Tmin + static_cast<double>(iT) * dT;
        ofs_dymat << "# " << temp << std::endl;
        for (is = 0; is < ns; ++is) {
            for (js = 0; js < ns; ++js) {
                for (ik = 0; ik < nk_interpolate; ++ik) {
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
}

void Scph::exec_scph_main(std::complex<double> ****dymat_anharm)
{
    int i, iT, ik, is, js;
    unsigned int nk = nk_scph;
    unsigned int ns = dynamical->neval;
    unsigned int NT;
    unsigned int nk_reduced_scph = kp_irred_scph.size();
    unsigned int nk_irred_interpolate = kp_irred_interpolate.size();
    double Tmin = system->Tmin;
    double Tmax = system->Tmax;
    double dT = system->dT;
    double temp;
    double ***omega2_anharm;
    std::complex<double> ***evec_anharm_tmp;
    std::complex<double> ***v3_array_all;
    std::complex<double> ***v4_array_all;
    bool converged_prev;

    std::vector<double> vec_temp;

    NT = static_cast<unsigned int>((Tmax - Tmin) / dT) + 1;

    // Find the degeneracy at each irreducible k points.

    std::vector<int> *degeneracy_at_k;
    memory->allocate(degeneracy_at_k, nk_reduced_scph);
    find_degeneracy(degeneracy_at_k, nk_reduced_scph,
                    kp_irred_scph, eval_harmonic);

    // Compute matrix element of 4-phonon interaction

    memory->allocate(omega2_anharm, NT, nk, ns);
    memory->allocate(evec_anharm_tmp, nk, ns, ns);
    memory->allocate(v4_array_all, nk_irred_interpolate * nk_scph,
                     ns * ns, ns * ns);

    // Calculate v4 array. 
    // This operation is the most expensive part of the calculation.
    if (selfenergy_offdiagonal & (ialgo == 1)) {
        compute_V4_array_all2(v4_array_all, evec_harmonic, selfenergy_offdiagonal);
    } else {
        compute_V4_array_all(v4_array_all, evec_harmonic,
                             selfenergy_offdiagonal, relax_coordinate);
    }

    if (relax_coordinate) {
        memory->allocate(v3_array_all, nk, ns, ns * ns);
        compute_V3_array_all(v3_array_all, evec_harmonic, selfenergy_offdiagonal);
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

        converged_prev = false;

        for (i = 0; i < vec_temp.size(); ++i) {
            temp = vec_temp[i];

            iT = static_cast<unsigned int>((temp - Tmin) / dT);

            // Initialize phonon eigenvectors with harmonic values

            for (ik = 0; ik < nk; ++ik) {
                for (is = 0; is < ns; ++is) {
                    for (js = 0; js < ns; ++js) {
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

            compute_anharmonic_frequency(v4_array_all, omega2_anharm[iT],
                                         evec_anharm_tmp, temp, degeneracy_at_k,
                                         converged_prev, cmat_convert,
                                         selfenergy_offdiagonal);
            calc_new_dymat_with_evec(dymat_anharm[iT], omega2_anharm[iT], evec_anharm_tmp);

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


void Scph::compute_V3_array_all(std::complex<double> ***v3_out,
                                std::complex<double> ***evec_in,
                                const bool self_offdiag)
{
    // Calculate the matrix elements of quartic terms in reciprocal space.
    // This is the most expensive part of the SCPH calculation.

    int ik_mpi;
    unsigned int nk_reduced_scph = kp_irred_scph.size();
    unsigned int nk_reduced_interpolate = kp_irred_interpolate.size();
    unsigned int ns = dynamical->neval;
    unsigned int ns2 = ns * ns;
    unsigned int ns3 = ns * ns * ns;
    unsigned int ns4 = ns * ns * ns * ns;
    unsigned int ik, jk, is, js, ks, ls;
    unsigned int knum, knum_minus, jk_minus;
    unsigned int nsize = fcs_phonon->force_constant_with_cell[1].size();
    double phase, phase3[3];
    int loc, loc3[3];
    unsigned int **ind;
    unsigned int i, j, k, ielem;
    std::complex<double> sum_tmp;
    std::complex<double> ret;
    long int ii;

    double inv2pi = 1.0 / (2.0 * pi);
    double dnk_represent = static_cast<double>(nk_represent);

    double factor = std::pow(0.5, 2) / static_cast<double>(nk_scph);
    static std::complex<double> complex_zero = std::complex<double>(0.0, 0.0);
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

    for (ik = mympi->my_rank; ik < nk_scph; ik += mympi->nprocs) {

        for (is = 0; is < ngroup; ++is) v3_array_at_kpair[is] = complex_zero;

        ielem = 0;

        for (i = 0; i < ngroup; ++i) {

            sum_tmp = std::complex<double>(0.0, 0.0);
            for (j = 0; j < 3; ++j) ind[i][j] = evec_index3[ielem][j];

            if (tune_type == 0) {
                for (j = 0; j < fcs_group[i].size(); ++j) {
                    phase
                        = xk_scph[ik][0] * (vec_for_v3[ielem][0][0] - vec_for_v3[ielem][1][0])
                        + xk_scph[ik][1] * (vec_for_v3[ielem][0][1] - vec_for_v3[ielem][1][1])
                        + xk_scph[ik][2] * (vec_for_v3[ielem][0][2] - vec_for_v3[ielem][1][2]);

                    loc = nint(phase * dnk_represent * inv2pi) % nk_represent + nk_represent - 1;

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
                        * evec_in[ik][js][ind[i][1]] * std::conj(evec_in[ik][ks][ind[i][2]]);
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
                            * evec_in[ik][js][ind[i][1]] * std::conj(evec_in[ik][ks][ind[i][2]]);
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
                            * evec_in[ik][js][ind[i][1]] * std::conj(evec_in[ik][js][ind[i][2]]);
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

    bool *is_acoustic;
    memory->allocate(is_acoustic, ns);
    int nacoustic;
    auto threshould = eps12;

    do {
        nacoustic = 0;
        for (is = 0; is < ns; ++is) {
            if (std::abs(eval_harmonic[0][is]) < threshould) {
                is_acoustic[is] = true;
                ++nacoustic;
            } else {
                is_acoustic[is] = false;
            }
        }
        if (nacoustic > 3) {
            error->exit("compute_V4_array_all",
                        "Could not assign acoustic modes at Gamma.");
        }
        threshould *= 2.0;
    } while (nacoustic < 3);

    // Set V3 to zeros so as to avoid mixing with gamma acoustic modes
    // jk = 0;
    for (is = 0; is < ns; ++is) {
        for (ks = 0; ks < ns; ++ks) {
            for (ls = 0; ls < ns; ++ls) {
                if (is_acoustic[ks] || is_acoustic[ls]) {
                    v3_out[0][is][ns * ks + ls] = complex_zero;
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
                        v3_out[jk][is][ns * ks + ls] = complex_zero;
                    }
                }
            }
        }
    }

    memory->deallocate(is_acoustic);

    if (mympi->my_rank == 0) {
        std::cout << " done !" << std::endl;
        timer->print_elapsed();
    }
}

void Scph::compute_V4_array_all(std::complex<double> ***v4_out,
                                std::complex<double> ***evec_in,
                                const bool self_offdiag,
                                const bool relax)
{
    // Calculate the matrix elements of quartic terms in reciprocal space.
    // This is the most expensive part of the SCPH calculation.

    int nk2_prod;
    unsigned int nk_reduced_scph = kp_irred_scph.size();
    unsigned int nk_reduced_interpolate = kp_irred_interpolate.size();
    unsigned int ns = dynamical->neval;
    unsigned int ns2 = ns * ns;
    unsigned int ns3 = ns * ns * ns;
    unsigned int ns4 = ns * ns * ns * ns;
    unsigned int ik, jk, is, js, ks, ls;
    unsigned int knum, knum_minus, jk_minus;
    unsigned int nsize = fcs_phonon->force_constant_with_cell[2].size();
    double phase, phase3[3];
    int loc, loc3[3];
    unsigned int **ind;
    unsigned int i, j, k, ielem;
    std::complex<double> sum_tmp;
    std::complex<double> ret;
    long int ii;

    double inv2pi = 1.0 / (2.0 * pi);
    double dnk_represent = static_cast<double>(nk_represent);

    double factor = std::pow(0.5, 2) / static_cast<double>(nk_scph);
    static std::complex<double> complex_zero = std::complex<double>(0.0, 0.0);
    std::complex<double> *v4_array_at_kpair;
    std::complex<double> ***v4_mpi;

    nk2_prod = nk_reduced_interpolate * nk_scph;

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
        ik = ik_prod / nk_scph;
        jk = ik_prod % nk_scph;

        knum = kmap_interpolate_to_scph[kp_irred_interpolate[ik][0].knum];

        knum_minus = knum_minus_scph[knum];
        jk_minus = knum_minus_scph[jk];

        for (is = 0; is < ngroup2; ++is) v4_array_at_kpair[is] = complex_zero;

        ielem = 0;

        for (i = 0; i < ngroup2; ++i) {

            sum_tmp = std::complex<double>(0.0, 0.0);
            for (j = 0; j < 4; ++j) ind[i][j] = evec_index4[ielem][j];

            if (tune_type == 0) {
                for (j = 0; j < fcs_group2[i].size(); ++j) {
                    phase
                        = - xk_scph[knum][0] * vec_for_v4[ielem][0][0]
                        - xk_scph[knum][1] * vec_for_v4[ielem][0][1]
                        - xk_scph[knum][2] * vec_for_v4[ielem][0][2]
                        + xk_scph[jk][0] * (vec_for_v4[ielem][1][0] - vec_for_v4[ielem][2][0])
                        + xk_scph[jk][1] * (vec_for_v4[ielem][1][1] - vec_for_v4[ielem][2][1])
                        + xk_scph[jk][2] * (vec_for_v4[ielem][1][2] - vec_for_v4[ielem][2][2]);

                    loc = nint(phase * dnk_represent * inv2pi) % nk_represent + nk_represent - 1;

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
                        * evec_in[knum][is][ind[i][0]] * std::conj(evec_in[knum][js][ind[i][1]])
                        * evec_in[jk][ks][ind[i][2]] * std::conj(evec_in[jk][ls][ind[i][3]]);
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
                            * evec_in[knum][is][ind[i][0]] * std::conj(evec_in[knum][js][ind[i][1]])
                            * evec_in[jk][ks][ind[i][2]] * std::conj(evec_in[jk][ls][ind[i][3]]);
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
                            * evec_in[knum][is][ind[i][0]] * std::conj(evec_in[knum][is][ind[i][1]])
                            * evec_in[jk][js][ind[i][2]] * std::conj(evec_in[jk][js][ind[i][3]]);
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

    bool *is_acoustic;
    memory->allocate(is_acoustic, ns);
    int nacoustic;
    auto threshould = eps12;

    do {
        nacoustic = 0;
        for (is = 0; is < ns; ++is) {
            if (std::abs(eval_harmonic[0][is]) < threshould) {
                is_acoustic[is] = true;
                ++nacoustic;
            } else {
                is_acoustic[is] = false;
            }
        }
        if (nacoustic > 3) {
            error->exit("compute_V4_array_all",
                        "Could not assign acoustic modes at Gamma.");
        }
        threshould *= 2.0;
    } while (nacoustic < 3);

    // Set V4 to zeros so as to avoid mixing with gamma acoustic modes
    // jk = 0;
    for (ik = 0; ik < nk_reduced_interpolate; ++ik) {
        for (is = 0; is < ns; ++is) {
            for (js = 0; js < ns; ++js) {
                for (ks = 0; ks < ns; ++ks) {
                    for (ls = 0; ls < ns; ++ls) {
                        if (is_acoustic[ks] || is_acoustic[ls]) {
                            v4_out[nk_scph * ik][ns * is + js][ns * ks + ls] = complex_zero;
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
                            v4_out[jk][ns * is + js][ns * ks + ls] = complex_zero;
                        }
                    }
                }
            }
        }
    }

    memory->deallocate(is_acoustic);

    if (mympi->my_rank == 0) {
        std::cout << " done !" << std::endl;
        timer->print_elapsed();
    }
}

void Scph::compute_V4_array_all2(std::complex<double> ***v4_out,
                                 std::complex<double> ***evec_in,
                                 const bool self_offdiag)
{
    // Calculate the matrix elements of quartic terms in reciprocal space.
    // This is the most expensive part of the SCPH calculation.

    int ik_prod;
    unsigned int nk_reduced_interpolate = kp_irred_interpolate.size();
    unsigned int ns = dynamical->neval;
    unsigned int ns2 = ns * ns;
    unsigned int ns4 = ns * ns * ns * ns;
    int is, js;
    unsigned int ks, ls;
    unsigned int knum;
    double phase, phase3[3];
    int loc, loc3[3];
    unsigned int **ind;
    unsigned int i, j, k, ielem;
    std::complex<double> sum_tmp;
    std::complex<double> ret;
    long int *nset_mpi;

    int ik_now, jk_now;
    int is_now, js_now, is_prod;

    double inv2pi = 1.0 / (2.0 * pi);
    double dnk_represent = static_cast<double>(nk_represent);

    double factor = std::pow(0.5, 2) / static_cast<double>(nk_scph);
    static std::complex<double> complex_zero = std::complex<double>(0.0, 0.0);
    std::complex<double> *v4_array_at_kpair;
    std::complex<double> ***v4_mpi;

    std::vector<int> ik_vec, jk_vec, is_vec, js_vec;

    int nk2_prod = nk_reduced_interpolate * nk_scph;

    if (mympi->my_rank == 0) {
        if (self_offdiag) {
            std::cout << " IALGO = 1 : Use different algorithm efficient when ns >> nk" << std::endl;
            std::cout << " SELF_OFFDIAG = 1: Calculating all components of v4_array ... ";
        } else {
            error->exit("compute_V4_array_all2",
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

        ik_now = ik_vec[ii];
        jk_now = jk_vec[ii];
        is_now = is_vec[ii];
        js_now = js_vec[ii];

        if (!((ik_now == ik_old) && (jk_now == jk_old))) {

            // Update v4_array_at_kpair and ind

            knum = kmap_interpolate_to_scph[kp_irred_interpolate[ik_now][0].knum];

            for (is = 0; is < ngroup2; ++is) v4_array_at_kpair[is] = complex_zero;

            ielem = 0;

            for (i = 0; i < ngroup2; ++i) {

                sum_tmp = std::complex<double>(0.0, 0.0);
                for (j = 0; j < 4; ++j) ind[i][j] = evec_index4[ielem][j];

                if (tune_type == 0) {
                    for (j = 0; j < fcs_group2[i].size(); ++j) {
                        phase
                            = - xk_scph[knum][0] * vec_for_v4[ielem][0][0]
                            - xk_scph[knum][1] * vec_for_v4[ielem][0][1]
                            - xk_scph[knum][2] * vec_for_v4[ielem][0][2]
                            + xk_scph[jk_now][0] * (vec_for_v4[ielem][1][0] - vec_for_v4[ielem][2][0])
                            + xk_scph[jk_now][1] * (vec_for_v4[ielem][1][1] - vec_for_v4[ielem][2][1])
                            + xk_scph[jk_now][2] * (vec_for_v4[ielem][1][2] - vec_for_v4[ielem][2][2]);

                        loc = nint(phase * dnk_represent * inv2pi) % nk_represent + nk_represent - 1;

                        sum_tmp += fcs_group2[i][j] * invmass_for_v4[ielem] * exp_phase[loc];

                        ++ielem;
                    }
                } else if (tune_type == 1) {
                    for (j = 0; j < fcs_group2[i].size(); ++j) {

                        for (k = 0; k < 3; ++k) {
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
        is_prod = ns * is_now + js_now;

#pragma omp parallel for private (ks, ls, ret, i)
        for (js = 0; js < ns2; ++js) {

            ks = js / ns;
            ls = js % ns;

            ret = std::complex<double>(0.0, 0.0);

            for (i = 0; i < ngroup2; ++i) {

                ret += v4_array_at_kpair[i]
                    * evec_in[knum][is_now][ind[i][0]] * std::conj(evec_in[knum][js_now][ind[i][1]])
                    * evec_in[jk_now][ks][ind[i][2]] * std::conj(evec_in[jk_now][ls][ind[i][3]]);
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

    bool *is_acoustic;
    memory->allocate(is_acoustic, ns);
    int nacoustic;
    auto threshould = eps12;

    do {
        nacoustic = 0;
        for (is = 0; is < ns; ++is) {
            if (std::abs(eval_harmonic[0][is]) < threshould) {
                is_acoustic[is] = true;
                ++nacoustic;
            } else {
                is_acoustic[is] = false;
            }
        }
        if (nacoustic > 3) {
            error->exit("compute_V4_array_all2",
                        "Could not assign acoustic modes at Gamma.");
        }
        threshould *= 2.0;
    } while (nacoustic < 3);

    // Set V4 to zeros so as to avoid mixing with gamma acoustic modes
    // jk = 0;
    for (unsigned int ik = 0; ik < nk_reduced_interpolate; ++ik) {
        for (is = 0; is < ns; ++is) {
            for (js = 0; js < ns; ++js) {
                for (ks = 0; ks < ns; ++ks) {
                    for (ls = 0; ls < ns; ++ls) {
                        if (is_acoustic[ks] || is_acoustic[ls]) {
                            v4_out[nk_scph * ik][ns * is + js][ns * ks + ls] = complex_zero;
                        }
                    }
                }
            }
        }
    }
    // ik = 0;
    for (unsigned int jk = 0; jk < nk_scph; ++jk) {
        for (is = 0; is < ns; ++is) {
            for (js = 0; js < ns; ++js) {
                if (is_acoustic[is] || is_acoustic[js]) {
                    for (ks = 0; ks < ns; ++ks) {
                        for (ls = 0; ls < ns; ++ls) {
                            v4_out[jk][ns * is + js][ns * ks + ls] = complex_zero;
                        }
                    }
                }
            }
        }
    }

    memory->deallocate(is_acoustic);

    if (mympi->my_rank == 0) {
        std::cout << " done !" << std::endl;
        timer->print_elapsed();
    }
}

void Scph::setup_kmesh()
{
    unsigned int ik;
    unsigned int i;
    unsigned int ns = dynamical->neval;
    int ik_minus, loc;
    double xtmp[3];
    double norm;

    // Set up k points for SCPH equation
    MPI_Bcast(&kmesh_scph[0], 3, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    MPI_Bcast(&kmesh_interpolate[0], 3, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

    nk_scph = kmesh_scph[0] * kmesh_scph[1] * kmesh_scph[2];
    memory->allocate(xk_scph, nk_scph, 3);
    memory->allocate(kvec_na_scph, nk_scph, 3);

    kpoint->gen_kmesh(true, kmesh_scph, xk_scph, kp_irred_scph);

    for (ik = 0; ik < nk_scph; ++ik) {
        for (i = 0; i < 3; ++i) {
            kvec_na_scph[ik][i] = 0.0;
        }
    }

    for (ik = 0; ik < nk_scph; ++ik) {

        for (i = 0; i < 3; ++i) xtmp[i] = xk_scph[ik][i];
        rotvec(xtmp, xtmp, system->rlavec_p, 'T');
        norm = xtmp[0] * xtmp[0] + xtmp[1] * xtmp[1] + xtmp[2] * xtmp[2];

        if (norm > eps) {
            for (i = 0; i < 3; ++i) kvec_na_scph[ik][i] = xk_scph[ik][i] / std::sqrt(norm);
        }
    }

    memory->allocate(knum_minus_scph, nk_scph);

    for (ik = 0; ik < nk_scph; ++ik) {
        for (i = 0; i < 3; ++i) xtmp[i] = -xk_scph[ik][i];

        ik_minus = kpoint->get_knum(xtmp, kmesh_scph);

        if (ik_minus == -1) error->exit("setup_kmesh", "Could not find the k point.");
        if (ik_minus < ik) continue;

        knum_minus_scph[ik] = ik_minus;
        knum_minus_scph[ik_minus] = ik;
    }

    // Set up k points for Fourier interpolation

    nk_interpolate = kmesh_interpolate[0] * kmesh_interpolate[1] * kmesh_interpolate[2];
    memory->allocate(xk_interpolate, nk_interpolate, 3);
    kpoint->gen_kmesh(true, kmesh_interpolate,
                      xk_interpolate, kp_irred_interpolate);

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

        loc = kpoint->get_knum(xtmp, kmesh_scph);

        if (loc == -1)
            error->exit("setup_kmesh",
                        "KMESH_INTERPOLATE should be a integral multiple of KMESH_SCPH");
        kmap_interpolate_to_scph[ik] = loc;
    }
}


void Scph::setup_transform_symmetry()
{
    unsigned int i;
    unsigned int ik;
    unsigned int is, js;
    unsigned int iat, jat, icrd, jcrd;
    unsigned int knum;
    unsigned int isym;
    unsigned int natmin = system->natmin;
    unsigned int ns = dynamical->neval;
    double phase;
    double x1[3], x2[3], k[3], k_minus[3], Sk[3], xtmp[3];
    double S_cart[3][3], S_frac[3][3], S_frac_inv[3][3];
    double S_recip[3][3];
    std::complex<double> im(0.0, 1.0);
    int knum_sym, knum_minus;

    std::complex<double> **gamma_tmp;
    std::complex<double> **dymat, **dymat_sym;
    Eigen::MatrixXcd Dymat(ns, ns), Gamma(ns, ns), Dymat2(ns, ns);
    Eigen::MatrixXcd Mat(ns, ns);
    bool *flag;

    memory->allocate(gamma_tmp, ns, ns);
    memory->allocate(dymat, ns, ns);
    memory->allocate(dymat_sym, ns, ns);
    memory->allocate(mat_transform_sym, kp_irred_interpolate.size(),
                     symmetry->nsym, ns, ns);
    memory->allocate(small_group_at_k, kp_irred_interpolate.size());
    memory->allocate(symop_minus_at_k, kp_irred_interpolate.size());
    memory->allocate(kpoint_map_symmetry, nk_interpolate);
    memory->allocate(flag, nk_interpolate);

    for (ik = 0; ik < nk_interpolate; ++ik) {
        flag[ik] = false;
    }

    for (ik = 0; ik < kp_irred_interpolate.size(); ++ik) {
        knum = kp_irred_interpolate[ik][0].knum;

        for (icrd = 0; icrd < 3; ++icrd) {
            k[icrd] = xk_interpolate[knum][icrd];
            k_minus[icrd] = -k[icrd];
        }
        knum_minus = kpoint->get_knum(k_minus, kmesh_interpolate);

        dynamical->calc_analytic_k(k, fcs_phonon->fc2_ext, dymat);

        isym = 0;

        small_group_at_k[ik].clear();
        symop_minus_at_k[ik].clear();

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
            for (i = 0; i < 3; ++i) Sk[i] = Sk[i] - nint(Sk[i]);

            knum_sym = kpoint->get_knum(Sk, kmesh_interpolate);
            if (knum_sym == -1)
                error->exit("modify_eigenvectors_sym",
                            "kpoint not found");
            if (knum_sym == knum) {
                small_group_at_k[ik].push_back(isym);
            }
            if (knum_sym == knum_minus) {
                symop_minus_at_k[ik].push_back(isym);
            }

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

            for (jat = 0; jat < natmin; ++jat) {
                iat = it.mapping[jat];

                // Fractional coordinates of x1 and x2
                for (icrd = 0; icrd < 3; ++icrd) {
                    x1[icrd] = system->xr_p[system->map_p2s[iat][0]][icrd];
                    x2[icrd] = system->xr_p[system->map_p2s[jat][0]][icrd];
                }

                rotvec(xtmp, x1, S_frac_inv);
                for (icrd = 0; icrd < 3; ++icrd) {
                    xtmp[icrd] = xtmp[icrd] - x2[icrd];
                }

                phase = 2.0 * pi * (k[0] * xtmp[0] + k[1] * xtmp[1] + k[2] * xtmp[2]);

                for (icrd = 0; icrd < 3; ++icrd) {
                    for (jcrd = 0; jcrd < 3; ++jcrd) {
                        gamma_tmp[3 * iat + icrd][3 * jat + jcrd]
                            = S_cart[icrd][jcrd] * std::exp(im * phase);
                    }
                }
            }

#ifdef _DEBUG
            Eigen::MatrixXcd Dymat3(ns,ns);
            dynamical->calc_analytic_k(Sk, fcs_phonon->fc2_ext, dymat_sym);

            for (is = 0; is < ns; ++is) {
                for (js = 0; js < ns; ++js) {
                    Dymat(is,js) = dymat[is][js];
                    Gamma(is,js) = gamma_tmp[is][js];
                }
            }
            Dymat2 = Gamma * Dymat * (Gamma.transpose()).conjugate();

            double diff = 0.0;
            for (is = 0; is < ns; ++is) {
                for (js = 0; js < ns; ++js) {
                    diff += std::norm(Dymat2(is,js) - dymat_sym[is][js]);
                }
            }
            std::cout << "DIFF = " << std::sqrt(diff) << std::endl;
            if (std::sqrt(diff) > eps12) {
                for (is = 0; is < ns; ++is) {
                    for (js = 0; js < ns; ++js) {
                        Dymat3(is,js) = dymat_sym[is][js];
                    }
                }
                std::cout << "isym = " << std::setw(4) << isym + 1 << std::endl;
                std::cout << " k = ";
                for (icrd = 0; icrd < 3; ++icrd) {
                    std::cout << std::setw(15) << k[icrd];
                }
                std::cout << std::endl;
                std::cout << "Sk = ";
                for (icrd = 0; icrd < 3; ++icrd) {
                    std::cout << std::setw(15) << Sk[icrd];
                }
                std::cout << std::endl;

                std::cout << " D(k):" << std::endl;
                std::cout << Dymat << std::endl;
                std::cout << " D(Sk) transform : " << std::endl;
                std::cout << Dymat2 << std::endl;
                std::cout << " D(Sk) exact : " << std::endl;
                std::cout << Dymat3 << std::endl;
                std::cout << std::endl;
            }
#endif

            for (is = 0; is < ns; ++is) {
                for (js = 0; js < ns; ++js) {
                    mat_transform_sym[ik][isym][is][js] = gamma_tmp[is][js];
                }
            }

            ++isym;
        }
    }

    memory->deallocate(dymat);
    memory->deallocate(dymat_sym);
    memory->deallocate(gamma_tmp);
}

void Scph::symmetrize_dynamical_matrix(const unsigned int ik,
                                       Eigen::MatrixXcd &dymat)
{
    using namespace Eigen;
    unsigned int i, isym;
    unsigned int is, js;
    unsigned int ns = dynamical->neval;
    MatrixXcd dymat_sym = MatrixXcd::Zero(ns, ns);
    MatrixXcd dymat_tmp(ns, ns), gamma(ns, ns);

    unsigned int nsym_small = small_group_at_k[ik].size();
    unsigned int nsym_minus = symop_minus_at_k[ik].size();

    for (i = 0; i < nsym_minus; ++i) {
        isym = symop_minus_at_k[ik][i];

        for (is = 0; is < ns; ++is) {
            for (js = 0; js < ns; ++js) {
                gamma(is, js) = mat_transform_sym[ik][isym][is][js];
            }
        }

        dymat_tmp = gamma * dymat * (gamma.transpose()).conjugate();
        dymat_sym += dymat_tmp.conjugate();
    }

    for (i = 0; i < nsym_small; ++i) {
        isym = small_group_at_k[ik][i];

        for (is = 0; is < ns; ++is) {
            for (js = 0; js < ns; ++js) {
                gamma(is, js) = mat_transform_sym[ik][isym][is][js];
            }
        }

        dymat_tmp = gamma * dymat * (gamma.transpose()).conjugate();
        dymat_sym += dymat_tmp;
    }


    dymat = dymat_sym / static_cast<double>(nsym_small + nsym_minus);
}

void Scph::replicate_dymat_for_all_kpoints(std::complex<double> ***dymat_inout)
{
    using namespace Eigen;
    unsigned int i, isym;
    unsigned int is, js;
    unsigned int ik_irred, ik_orig;
    unsigned int ns = dynamical->neval;
    MatrixXcd dymat_tmp(ns, ns), gamma(ns, ns), dymat(ns, ns);

    std::complex<double> ***dymat_all;

    memory->allocate(dymat_all, ns, ns, nk_interpolate);

    for (i = 0; i < nk_interpolate; ++i) {

        ik_irred = kpoint_map_symmetry[i].knum_irred_orig;
        ik_orig = kpoint_map_symmetry[i].knum_orig;
        isym = kpoint_map_symmetry[i].symmetry_op;

        for (is = 0; is < ns; ++is) {
            for (js = 0; js < ns; ++js) {
                gamma(is, js) = mat_transform_sym[ik_irred][isym][is][js];
                dymat(is, js) = dymat_inout[is][js][ik_orig];
            }
        }
        dymat_tmp = gamma * dymat * (gamma.transpose()).conjugate();

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
    int ik;
    unsigned int is;
    unsigned int nk = kpoint->nk;
    unsigned int ns = dynamical->neval;

    if (mympi->my_rank == 0) {
        std::cout << std::endl
            << " Diagonalizing dynamical matrices for all k points ... ";
    }

    memory->allocate(eval_harmonic, nk_scph, ns);
    memory->allocate(evec_harmonic, nk_scph, ns, ns);
    memory->allocate(omega2_harmonic, nk_scph, ns);

    // Calculate phonon eigenvalues and eigenvectors for all k-points for scph

#pragma omp parallel for private (is)
    for (ik = 0; ik < nk_scph; ++ik) {

        dynamical->eval_k(xk_scph[ik], kvec_na_scph[ik],
                          fcs_phonon->fc2_ext, eval_harmonic[ik],
                          evec_harmonic[ik], true);

        // Phonon energy is the square-root of the eigenvalue 
        for (is = 0; is < ns; ++is) {
            omega2_harmonic[ik][is] = eval_harmonic[ik][is];
            eval_harmonic[ik][is] = dynamical->freq(eval_harmonic[ik][is]);
        }
    }

    if (mympi->my_rank == 0) {
        std::cout << "done !" << std::endl;
    }
}

void Scph::setup_pp_interaction()
{
    unsigned int nk = nk_scph;
    unsigned int ns = dynamical->neval;

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

    int ii, jj, kk;

    if (tune_type == 0) {

        double phase;

        memory->allocate(exp_phase, 2 * nk_represent - 1);
        for (ii = 0; ii < 2 * nk_represent - 1; ++ii) {
            phase = 2.0 * pi * static_cast<double>(ii - nk_represent + 1)
                / static_cast<double>(nk_represent);
            exp_phase[ii] = std::exp(im * phase);
        }

    } else if (tune_type == 1) {
        double phase[3];

        tune_type = 1;
        memory->allocate(exp_phase3, 2 * nk_grid[0] - 1, 2 * nk_grid[1] - 1, 2 * nk_grid[2] - 1);

        for (ii = 0; ii < 2 * nk_grid[0] - 1; ++ii) {
            phase[0] = 2.0 * pi * static_cast<double>(ii - nk_grid[0] + 1) / dnk[0];
            for (jj = 0; jj < 2 * nk_grid[1] - 1; ++jj) {
                phase[1] = 2.0 * pi * static_cast<double>(jj - nk_grid[1] + 1) / dnk[1];
                for (kk = 0; kk < 2 * nk_grid[2] - 1; ++kk) {
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
    int i, j;
    int ix, iy, iz;
    unsigned int nk = nk_interpolate;
    unsigned int nat = system->natmin;
    unsigned int iat, jat;

    int **shift_cell, **shift_cell_super;
    double **xf_p;
    double ****x_all;

    int nkx = kmesh_interpolate[0];
    int nky = kmesh_interpolate[1];
    int nkz = kmesh_interpolate[2];

    unsigned int ncell = nk;
    unsigned int ncell_s = 27;

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

    double dist, dist_min;
    std::vector<DistList> dist_tmp;
    ShiftCell shift_tmp;
    std::vector<int> vec_tmp;

    memory->allocate(mindist_list_scph, nat, nat, ncell);

    for (iat = 0; iat < nat; ++iat) {
        for (jat = 0; jat < nat; ++jat) {
            for (icell = 0; icell < ncell; ++icell) {

                dist_tmp.clear();
                for (i = 0; i < ncell_s; ++i) {
                    dist = distance(x_all[0][0][iat], x_all[i][icell][jat]);
                    dist_tmp.emplace_back(i, dist);
                }
                std::sort(dist_tmp.begin(), dist_tmp.end());

                dist_min = dist_tmp[0].dist;
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


void Scph::exec_interpolation(std::complex<double> ***dymat_r,
                              double **eval_out,
                              std::complex<double> ***evec_out)
{
    unsigned int i, j, is;
    unsigned int ns = dynamical->neval;
    unsigned int nk = nk_interpolate;
    unsigned int nk1 = kmesh_interpolate[0];
    unsigned int nk2 = kmesh_interpolate[1];
    unsigned int nk3 = kmesh_interpolate[2];
    unsigned int nk_ref = kpoint->nk;

    double eval_tmp;
    double *eval_real;

    std::complex<double> **mat_tmp;
    std::complex<double> **mat_harmonic, **mat_harmonic_na;
    std::vector<double> eval_vec;

    memory->allocate(mat_tmp, ns, ns);
    memory->allocate(eval_real, ns);
    memory->allocate(mat_harmonic, ns, ns);

    if (dynamical->nonanalytic) {
        memory->allocate(mat_harmonic_na, ns, ns);
    }

    for (unsigned int ik = 0; ik < nk_ref; ++ik) {
        dynamical->calc_analytic_k(kpoint->xk[ik], fcs_phonon->fc2_ext, mat_harmonic);
        r2q(kpoint->xk[ik], nk1, nk2, nk3, ns, dymat_r, mat_tmp);

        for (i = 0; i < ns; ++i) {
            for (j = 0; j < ns; ++j) {
                mat_tmp[i][j] += mat_harmonic[i][j];
            }
        }

        if (dynamical->nonanalytic) {

            if (dynamical->nonanalytic == 1) {
                dynamical->calc_nonanalytic_k(kpoint->xk[ik],
                                              kpoint->kvec_na[ik],
                                              mat_harmonic_na);
            } else if (dynamical->nonanalytic == 2) {
                dynamical->calc_nonanalytic_k2(kpoint->xk[ik],
                                               kpoint->kvec_na[ik],
                                               mat_harmonic_na);
            }

            for (i = 0; i < ns; ++i) {
                for (j = 0; j < ns; ++j) {
                    mat_tmp[i][j] += mat_harmonic_na[i][j];
                }
            }
        }

        diagonalize_interpolated_matrix(mat_tmp, eval_real, evec_out[ik], true);

        eval_vec.clear();

        for (is = 0; is < ns; ++is) {
            eval_tmp = eval_real[is];

            if (eval_tmp < 0.0) {
                eval_vec.push_back(-std::sqrt(-eval_tmp));
            } else {
                eval_vec.push_back(std::sqrt(eval_tmp));
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

void Scph::exec_interpolation2(std::complex<double> ***dymat_r,
                               double **eval_out,
                               std::complex<double> ***evec_out)
{
    unsigned int i, j, is;
    unsigned int ns = dynamical->neval;
    unsigned int nk = nk_interpolate;
    unsigned int nk1 = kmesh_interpolate[0];
    unsigned int nk2 = kmesh_interpolate[1];
    unsigned int nk3 = kmesh_interpolate[2];
    unsigned int nk_ref = scph->nk_scph;

    double eval_tmp;
    double *eval_real;

    std::complex<double> **mat_tmp;
    std::complex<double> **mat_harmonic, **mat_harmonic_na;
    std::vector<double> eval_vec;

    memory->allocate(mat_tmp, ns, ns);
    memory->allocate(eval_real, ns);
    memory->allocate(mat_harmonic, ns, ns);
    if (dynamical->nonanalytic) {
        memory->allocate(mat_harmonic_na, ns, ns);
    }

    for (unsigned int ik = 0; ik < nk_ref; ++ik) {
        dynamical->calc_analytic_k(xk_scph[ik], fcs_phonon->fc2_ext, mat_harmonic);
        r2q(xk_scph[ik], nk1, nk2, nk3, ns, dymat_r, mat_tmp);

        for (i = 0; i < ns; ++i) {
            for (j = 0; j < ns; ++j) {
                mat_tmp[i][j] += mat_harmonic[i][j];
            }
        }

        if (dynamical->nonanalytic) {

            if (dynamical->nonanalytic == 1) {
                dynamical->calc_nonanalytic_k(scph->xk_scph[ik],
                                              scph->kvec_na_scph[ik],
                                              mat_harmonic_na);
            } else if (dynamical->nonanalytic == 2) {
                dynamical->calc_nonanalytic_k2(scph->xk_scph[ik],
                                               scph->kvec_na_scph[ik],
                                               mat_harmonic_na);
            }

            for (i = 0; i < ns; ++i) {
                for (j = 0; j < ns; ++j) {
                    mat_tmp[i][j] += mat_harmonic_na[i][j];
                }
            }
        }
        diagonalize_interpolated_matrix(mat_tmp, eval_real, evec_out[ik], true);
        eval_vec.clear();

        for (is = 0; is < ns; ++is) {
            eval_tmp = eval_real[is];

            if (eval_tmp < 0.0) {
                eval_vec.push_back(-std::sqrt(-eval_tmp));
            } else {
                eval_vec.push_back(std::sqrt(eval_tmp));
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
               std::complex<double> **dymat_k_out)
{
    unsigned int icell;
    unsigned int i, j;
    double phase;
    unsigned int iat, jat;
    std::complex<double> im(0.0, 1.0);
    std::complex<double> exp_phase;

    unsigned int ncell = nx * ny * nz;

    for (i = 0; i < ns; ++i) {

        iat = i / 3;

        for (j = 0; j < ns; ++j) {

            jat = j / 3;

            dymat_k_out[i][j] = std::complex<double>(0.0, 0.0);

            for (icell = 0; icell < ncell; ++icell) {

                exp_phase = std::complex<double>(0.0, 0.0);

                // This operation is necessary for the Hermiticity of the dynamical matrix.
                for (auto it = mindist_list_scph[iat][jat][icell].shift.cbegin();
                     it != mindist_list_scph[iat][jat][icell].shift.cend(); ++it) {

                    phase = 2.0 * pi
                        * (static_cast<double>((*it).sx) * xk_in[0]
                            + static_cast<double>((*it).sy) * xk_in[1]
                            + static_cast<double>((*it).sz) * xk_in[2]);

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
                                           const bool require_evec)
{
    unsigned int i, j, k;
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

    k = 0;
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
                           double **eval)
{
    unsigned int knum;
    unsigned int ns = dynamical->neval;
    double omega_prev, omega_now;
    int ideg;
    double tol_omega = 1.0e-5;

    for (unsigned int ik = 0; ik < nk_irred; ++ik) {
        knum = kp_info[ik][0].knum;

        degeneracy_out[ik].clear();

        omega_prev = eval[knum][0];
        ideg = 1;

        for (unsigned int is = 1; is < ns; ++is) {
            omega_now = eval[knum][is];

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
    int knum;
    int ns = dynamical->neval;

    unsigned int ns2 = ns * ns;
    unsigned int m;

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

        knum = kmap_interpolate_to_scph[ik];

        // create eigval matrix

        for (is = 0; is < ns2; ++is) eigval_matrix[is] = std::complex<double>(0.0, 0.0);

        m = 0;
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
        dynamical->calc_analytic_k(xk_interpolate[ik], fcs_phonon->fc2_ext,
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

    unsigned int nk1 = kmesh_interpolate[0];
    unsigned int nk2 = kmesh_interpolate[1];
    unsigned int nk3 = kmesh_interpolate[2];
    unsigned int i, j;


    /*
        fftw_plan plan;
    
        for (i = 0; i < ns; ++i) {
            for (j = 0; j < ns; ++j) {
                plan = fftw_plan_dft_3d(nk1, nk2, nk3, reinterpret_cast<fftw_complex*>(dymat_q[i][j]), 
                    reinterpret_cast<fftw_complex*>(dymat_out[i][j]), FFTW_FORWARD, FFTW_ESTIMATE);
                fftw_execute(plan);
                fftw_destroy_plan(plan);
    
                for (ik = 0; ik < nk_interpolate; ++ik) dymat_out[i][j][ik] /= static_cast<double>(nk_interpolate);
            }
        }
    */


    double phase;
    std::complex<double> cexp_phase;
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

                    cexp_phase = std::complex<double>(0.0, 0.0);

                    for (i = 0; i < xk_dup.size(); ++i) {

                        phase = 2.0 * pi * (xk_dup[i][0] * static_cast<double>(ix)
                            + xk_dup[i][1] * static_cast<double>(iy)
                            + xk_dup[i][2] * static_cast<double>(iz));
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
    unsigned int is, js, ks, ls;
    unsigned int kk;
    unsigned int nk = nk_scph;
    unsigned int ns = dynamical->neval;
    unsigned int knum, knum_interpolate;
    unsigned int nk_irred_interpolate = kp_irred_interpolate.size();
    unsigned int nk1 = kmesh_interpolate[0];
    unsigned int nk2 = kmesh_interpolate[1];
    unsigned int nk3 = kmesh_interpolate[2];
    int icount;
    int iloop;

    MatrixXd mat_evec(ns, ns);
    MatrixXd omega_now(nk, ns), omega_old(nk, ns);
    MatrixXd omega2_harmonic(nk, ns);
    MatrixXcd mat_tmp(ns, ns), evec_tmp(ns, ns);

    VectorXd eval_tmp(ns), eval_orig(ns);
    MatrixXcd Dymat(ns, ns), Dymat_sym(ns, ns);
    MatrixXcd Fmat(ns, ns);
    MatrixXcd Qmat = MatrixXcd::Zero(ns, ns);
    MatrixXcd Cmat(ns, ns), Dmat(ns, ns);

    double omega_tmp, diff;
    double omega1, n1;
    double conv_tol = tolerance_scph;
    double alpha = mixalpha;

    double **eval_interpolate;
    double omega2_tmp;
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

    static std::complex<double> complex_one = std::complex<double>(1.0, 0.0);
    static std::complex<double> complex_zero = std::complex<double>(0.0, 0.0);

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

    double T_in = temp;

    std::cout << " Temperature = " << T_in << " K" << std::endl;

    // Set initial values

    for (ik = 0; ik < nk; ++ik) {
        for (is = 0; is < ns; ++is) {
            omega_tmp = eval_harmonic[ik][is]; // This is harmonic frequency

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
                omega_now(ik, is) = std::abs(omega_tmp);
            }

            if (omega_tmp < 0.0) {
                omega2_harmonic(ik, is) = -std::pow(omega_tmp, 2.0);
            } else {
                omega2_harmonic(ik, is) = std::pow(omega_tmp, 2.0);
            }

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
            mat_omega2_harmonic[ik][is][is] = std::complex<double>(omega2_harmonic(knum, is), 0.0);
            eval_interpolate[knum][is] = omega2_harmonic(knum, is);
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

    icount = 0;

    // Main loop
    for (iloop = 0; iloop < maxiter; ++iloop) {

        for (ik = 0; ik < nk; ++ik) {
            for (is = 0; is < ns; ++is) {
                omega1 = omega_now(ik, is);
                if (std::abs(omega1) < eps8) {
                    Qmat(is, is) = complex_zero;
                } else {
                    // Note that the missing factor 2 in the denominator of Qmat is 
                    // already considered in the v4_array_all.
                    if (thermodynamics->classical) {
                        Qmat(is, is) = std::complex<double>(2.0 * T_in * thermodynamics->T_to_Ryd / (omega1 * omega1),
                                                            0.0);
                    } else {
                        n1 = thermodynamics->fB(omega1, T_in);
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

#pragma omp parallel for private(jk,kk,ks,ls), reduction(+:re_tmp, im_tmp)
                    for (jk = 0; jk < nk; ++jk) {

                        kk = nk * ik + jk;

                        for (ks = 0; ks < ns; ++ks) {
                            ctmp = v4_array_all[kk][i][(ns + 1) * ks] * dmat_convert[jk][ks][ks];
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

#pragma omp parallel for private(jk,kk,ks,ls), reduction(+:re_tmp, im_tmp)
                        for (jk = 0; jk < nk; ++jk) {

                            kk = nk * ik + jk;

                            for (ks = 0; ks < ns; ++ks) {
                                for (ls = 0; ls < ns; ++ls) {
                                    ctmp = v4_array_all[kk][i][ns * ks + ls] * dmat_convert[jk][ks][ls];
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

                omega2_tmp = eval_tmp(is);

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

        // Inverse Fourier transform of delta Dymat.
        fftw_plan plan;

        for (is = 0; is < ns; ++is) {
            for (js = 0; js < ns; ++js) {
                plan = fftw_plan_dft_3d(nk1, nk2, nk3,
                                        reinterpret_cast<fftw_complex*>(dymat_q[is][js]),
                                        reinterpret_cast<fftw_complex*>(dymat_new[is][js]),
                                        FFTW_FORWARD, FFTW_ESTIMATE);
                fftw_execute(plan);
                fftw_destroy_plan(plan);

                for (ik = 0; ik < nk_interpolate; ++ik)
                    dymat_new[is][js][ik] /= static_cast<double>(nk_interpolate);
            }
        }

        exec_interpolation2(dymat_new, eval_interpolate, evec_new);

#ifdef _DEBUG
        for (ik = 0; ik < nk; ++ik) {
            std::cout << " ik = " << std::setw(4) << ik + 1 << " : ";
            for (i = 0; i < 3; ++i)  std::cout << std::setw(15) << xk_scph[ik][i];
            std::cout << std::endl;
            for (is = 0; is < ns; ++is) {
                std::cout << std::setw(15) << eval_harmonic[ik][is];
                std::cout << std::setw(15) << eval_interpolate[ik][is];
                std::cout << std::setw(15) << eval_harmonic[ik][is] - eval_interpolate[ik][is] << std::endl;
            }
        }
#endif

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


void Scph::write_scph_energy(double ***eval)
{
    unsigned int nk = kpoint->nk;
    unsigned int ns = dynamical->neval;
    double Tmin = system->Tmin;
    double Tmax = system->Tmax;
    double dT = system->dT;
    double temp;
    unsigned int NT = static_cast<unsigned int>((Tmax - Tmin) / dT) + 1;

    std::ofstream ofs_energy;
    std::string file_energy = input->job_title + ".scph_eval";

    ofs_energy.open(file_energy.c_str(), std::ios::out);
    if (!ofs_energy) error->exit("write_scph_energy", "cannot open file_energy");

    ofs_energy << "# K point, mode, Temperature [K], Eigenvalues [cm^-1]" << std::endl;


    for (unsigned int ik = 0; ik < nk; ++ik) {
        for (unsigned int is = 0; is < ns; ++is) {
            for (unsigned int iT = 0; iT < NT; ++iT) {
                temp = Tmin + static_cast<double>(iT) * dT;

                ofs_energy << std::setw(5) << ik + 1;
                ofs_energy << std::setw(5) << is + 1;
                ofs_energy << std::setw(8) << temp;
                ofs_energy << std::setw(15) << writes->in_kayser(eval[iT][ik][is]);
                ofs_energy << std::endl;
            }
            ofs_energy << std::endl;
        }
        ofs_energy << std::endl;
    }

    ofs_energy.close();
}

void Scph::write_scph_bands(double ***eval)
{
    std::ofstream ofs_bands;
    std::string file_bands = input->job_title + ".scph_bands";

    ofs_bands.open(file_bands.c_str(), std::ios::out);
    if (!ofs_bands) error->exit("write_scph_bands", "cannot open file_bands");

    unsigned int i, j;
    unsigned int nk = kpoint->nk;

    double *kaxis = kpoint->kaxis;
    double temp;
    double Tmin = system->Tmin;
    double Tmax = system->Tmax;
    double dT = system->dT;
    unsigned int NT = static_cast<unsigned int>((Tmax - Tmin) / dT) + 1;
    unsigned int ns = dynamical->neval;
    int kcount = 0;

    std::string str_tmp = "NONE";
    std::string str_kpath = "";
    std::string str_kval = "";

    for (i = 0; i < kpoint->kpInp.size(); ++i) {
        if (str_tmp != kpoint->kpInp[i].kpelem[0]) {
            str_tmp = kpoint->kpInp[i].kpelem[0];
            str_kpath += " " + str_tmp;

            std::ostringstream ss;
            ss << std::fixed << std::setprecision(6) << kpoint->kaxis[kcount];
            str_kval += " " + ss.str();
        }
        kcount += std::atoi(kpoint->kpInp[i].kpelem[8].c_str());

        if (str_tmp != kpoint->kpInp[i].kpelem[4]) {
            str_tmp = kpoint->kpInp[i].kpelem[4];
            str_kpath += " " + str_tmp;

            std::ostringstream ss;
            ss << std::fixed << std::setprecision(6) << kpoint->kaxis[kcount - 1];
            str_kval += " " + ss.str();
        }
    }

    ofs_bands << "# " << str_kpath << std::endl;
    ofs_bands << "#" << str_kval << std::endl;
    ofs_bands << "# Temperature [K], k-axis, Eigenvalues [cm^-1]" << std::endl;

    for (unsigned int iT = 0; iT < NT; ++iT) {
        temp = Tmin + static_cast<double>(iT) * dT;

        for (i = 0; i < nk; ++i) {
            ofs_bands << std::setw(15) << std::fixed << temp;
            ofs_bands << std::setw(15) << std::fixed << kaxis[i];
            for (j = 0; j < ns; ++j) {
                ofs_bands << std::setw(15) << std::scientific << writes->in_kayser(eval[iT][i][j]);
            }
            ofs_bands << std::endl;
        }
        ofs_bands << std::endl;
    }

    ofs_bands.close();

    std::cout << "  " << std::setw(input->job_title.length() + 12) << std::left << file_bands;
    std::cout << " : SCPH band structure" << std::endl;
}

void Scph::write_scph_dos(double ***eval)
{
    unsigned int iT;
    double Tmin = system->Tmin;
    double Tmax = system->Tmax;
    double dT = system->dT;
    unsigned int NT = static_cast<unsigned int>((Tmax - Tmin) / dT) + 1;

    double **dos_scph;

    memory->allocate(dos_scph, NT, dos->n_energy);
    dos->calc_dos_scph(eval, dos_scph);

    std::ofstream ofs_dos;
    std::string file_dos = input->job_title + ".scph_dos";

    ofs_dos.open(file_dos.c_str(), std::ios::out);
    if (!ofs_dos) error->exit("write_scph_dos", "cannot open file_dos");

    ofs_dos << "# ";

    for (iT = 0; iT < NT; ++iT) {
        ofs_dos << std::setw(15) << Tmin + static_cast<double>(iT) * dT;
    }
    ofs_dos << std::endl;

    for (unsigned int j = 0; j < dos->n_energy; ++j) {
        ofs_dos << std::setw(15) << dos->energy_dos[j];

        for (iT = 0; iT < NT; ++iT) {
            ofs_dos << std::setw(15) << dos_scph[iT][j];
        }
        ofs_dos << std::endl;
    }

    ofs_dos << std::endl;

    ofs_dos.close();

    memory->deallocate(dos_scph);
}

void Scph::write_scph_thermodynamics(double ***eval,
                                     std::complex<double> ****evec)
{
    int i;
    unsigned int ik, is;
    unsigned int nk = kpoint->nk;
    unsigned int ns = dynamical->neval;
    double Tmin = system->Tmin;
    double Tmax = system->Tmax;
    double dT = system->dT;
    unsigned int NT = static_cast<unsigned int>((Tmax - Tmin) / dT) + 1;
    double temp;
    double tmp1, tmp2, tmp3;
    double omega, x;
    double T_to_Ryd = thermodynamics->T_to_Ryd;

    int N = nk * ns;

    std::ofstream ofs_thermo;
    std::string file_thermo = input->job_title + ".scph_thermo";
    ofs_thermo.open(file_thermo.c_str(), std::ios::out);
    if (!ofs_thermo)
        error->exit("write_scph_thermodynamics",
                    "cannot open file_thermo");

    if (thermodynamics->calc_FE_bubble) {
        ofs_thermo << "# The bubble free-energy calculated on top of the SCPH wavefunction is also shown." << std::endl;
        ofs_thermo <<
            "# Temperature [K], Cv [in kB unit], F_{vib} (QHA term) [Ry], F_{vib} (SCPH) [Ry], F_{vib} (Bubble) [Ry]"
            << std::endl;
    } else {
        ofs_thermo << "# Temperature [K], Cv [in kB unit], F_{vib} (QHA term) [Ry], F_{vib} (SCPH) [Ry]"
            << std::endl;
    }

    for (unsigned int iT = 0; iT < NT; ++iT) {

        temp = Tmin + static_cast<double>(iT) * dT;

        tmp1 = 0.0;
        tmp2 = 0.0;

#pragma omp parallel for private(ik, is, omega, x), reduction(+:tmp1,tmp2)
        for (i = 0; i < N; ++i) {
            ik = i / ns;
            is = i % ns;
            omega = eval[iT][ik][is];

            if (omega <= eps8) continue;

            tmp1 += thermodynamics->Cv(omega, temp);

            if (std::abs(temp) < eps) {
                tmp2 += 0.5 * omega;
            } else {
                x = omega / (temp * T_to_Ryd);
                tmp2 += 0.5 * x + std::log(1.0 - std::exp(-x));
            }
        }

        tmp1 /= static_cast<double>(nk);
        if (std::abs(temp) < eps) {
            tmp2 /= static_cast<double>(nk);
        } else {
            tmp2 *= temp * T_to_Ryd / static_cast<double>(nk);
        }
        tmp3 = tmp2 + FE_scph(iT, eval[iT], evec[iT]);

        ofs_thermo << std::setw(16) << std::fixed << temp;
        ofs_thermo << std::setw(18) << std::scientific << tmp1 / k_Boltzmann;
        ofs_thermo << std::setw(18) << tmp2;
        ofs_thermo << std::setw(18) << tmp3;

        if (thermodynamics->calc_FE_bubble) {
            ofs_thermo << std::setw(18) << thermodynamics->FE_bubble[iT];
        }
        ofs_thermo << std::endl;
    }

    ofs_thermo.close();
}

double Scph::FE_scph(unsigned int iT,
                     double **eval,
                     std::complex<double> ***evec)
{
    // This function calculates the additional free energy by the SCPH approximation
    int ik, is, js, ks, ls;
    int i;
    unsigned int nk = kpoint->nk;
    unsigned int ns = dynamical->neval;
    double ret;
    double omega, omega2_harm;
    std::complex<double> tmp_c;
    double temp = system->Tmin + static_cast<double>(iT) * system->dT;
    int N = nk * ns;

    ret = 0.0;

#pragma omp parallel for private(ik, is, js, ks, ls, omega, tmp_c), reduction(+ : ret)
    for (i = 0; i < N; ++i) {
        ik = i / ns;
        is = i % ns;
        omega = eval[ik][is];
        if (std::abs(omega) < eps6) continue;

        tmp_c = std::complex<double>(0.0, 0.0);

        for (js = 0; js < ns; ++js) {
            omega2_harm = dynamical->eval_phonon[ik][js];
            if (omega2_harm >= 0.0) {
                omega2_harm = std::pow(omega2_harm, 2);
            } else {
                omega2_harm = -std::pow(omega2_harm, 2);
            }

            for (ks = 0; ks < ns; ++ks) {
                for (ls = 0; ls < ns; ++ls) {
                    tmp_c += omega2_harm
                        * dynamical->evec_phonon[ik][js][ks] * std::conj(dynamical->evec_phonon[ik][js][ls])
                        * std::conj(evec[ik][is][ks]) * evec[ik][is][ls];
                }
            }
        }

        ret += (tmp_c.real() - omega * omega) * (1.0 + 2.0 * thermodynamics->fB(omega, temp)) / (8.0 * omega);
    }

    return ret / static_cast<double>(nk);
}


void Scph::write_scph_msd(double ***eval,
                          std::complex<double> ****evec)
{
    unsigned int i, iT;
    unsigned int ik, is;
    unsigned int nk = kpoint->nk;
    unsigned int ns = dynamical->neval;
    double Tmin = system->Tmin;
    double Tmax = system->Tmax;
    double dT = system->dT;
    unsigned int NT = static_cast<unsigned int>((Tmax - Tmin) / dT) + 1;
    double *Cv;
    double temp;
    double tmp;
    double omega;
    double **msd;

    memory->allocate(Cv, NT);
    memory->allocate(msd, NT, ns);

    for (iT = 0; iT < NT; ++iT) {

        temp = Tmin + static_cast<double>(iT) * dT;

        for (i = 0; i < ns; ++i) {
            tmp = 0.0;

            for (ik = 0; ik < nk; ++ik) {
                for (is = 0; is < ns; ++is) {
                    omega = eval[iT][ik][is];

                    if (omega < eps8) continue;
                    if (thermodynamics->classical) {
                        tmp += std::norm(evec[iT][ik][is][i])
                            * thermodynamics->fC(omega, temp) / omega;
                    } else {
                        tmp += std::norm(evec[iT][ik][is][i])
                            * (thermodynamics->fB(omega, temp) + 0.5) / omega;
                    }
                }
            }
            msd[iT][i] = tmp * std::pow(Bohr_in_Angstrom, 2.0)
                / (static_cast<double>(nk) * system->mass[system->map_p2s[i / 3][0]]);
        }
    }

    std::ofstream ofs_msd;
    std::string file_msd = input->job_title + ".scph_msd";
    ofs_msd.open(file_msd.c_str(), std::ios::out);
    if (!ofs_msd) error->exit("write_scph_msd", "cannot open file_thermo");
    ofs_msd << "# Mean Square Displacements at a function of temperature." << std::endl;
    ofs_msd << "# Temperature [K], <(u_{1}^{x})^{2}>, <(u_{1}^{y})^{2}>, <(u_{1}^{z})^{2}>, .... [Angstrom^2]" << std::
        endl;

    for (iT = 0; iT < NT; ++iT) {
        temp = Tmin + static_cast<double>(iT) * dT;

        ofs_msd << std::setw(15) << temp;
        for (i = 0; i < ns; ++i) {
            ofs_msd << std::setw(15) << msd[iT][i];
        }
        ofs_msd << std::endl;
    }

    ofs_msd.close();

    memory->deallocate(msd);
}

double Scph::distance(double *x1,
                      double *x2)
{
    double dist = std::pow(x1[0] - x2[0], 2) + std::pow(x1[1] - x2[1], 2) + std::pow(x1[2] - x2[2], 2);
    dist = std::sqrt(dist);

    return dist;
}

void Scph::duplicate_xk_boundary(double *xk_in,
                                 std::vector<std::vector<double>> &vec_xk)
{
    int i, j, k, l;
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
        for (j = 0; j < n[1]; ++j) {
            sign[1] = 1.0 - 2.0 * static_cast<double>(j);
            for (k = 0; k < n[2]; ++k) {
                sign[2] = 1.0 - 2.0 * static_cast<double>(k);

                vec_tmp.clear();
                for (l = 0; l < 3; ++l) {
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
    double Tmin = system->Tmin;
    double Tmax = system->Tmax;
    double dT = system->dT;
    double temp;
    double ***delta_fc2;
    double **xtmp;
    unsigned int ns = dynamical->neval;
    unsigned int is, js, icell;
    unsigned int iat, jat, icrd, jcrd;
    unsigned int nmulti;

    std::ofstream ofs_fc2;
    std::string file_fc2 = input->job_title + ".scph_fc2_correction";

    ofs_fc2.open(file_fc2.c_str(), std::ios::out);
    if (!ofs_fc2)
        error->exit("write_anharmonic_correction_fc2",
                    "Cannot open file_fc2");

    unsigned int ncell = kmesh_interpolate[0] * kmesh_interpolate[1] * kmesh_interpolate[2];

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
            ofs_fc2 << std::setw(15) << system->lavec_p[j][i];
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
        temp = Tmin + dT * static_cast<double>(iT);

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
                icrd = is % 3;

                for (js = 0; js < ns; ++js) {
                    jat = js / 3;
                    jcrd = js % 3;

                    nmulti = mindist_list_scph[iat][jat][icell].shift.size();

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
}

void Scph::mpi_bcast_complex(std::complex<double> ****data,
                             const int NT,
                             const int nk,
                             const int ns)
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
