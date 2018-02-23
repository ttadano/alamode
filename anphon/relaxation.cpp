/*
relaxation.cpp

Copyright (c) 2014, 2015, 2016 Terumasa Tadano

This file is distributed under the terms of the MIT license.
Please see the file 'LICENCE.txt' in the root directory
or http://opensource.org/licenses/mit-license.php for information.
*/

#include "mpi_common.h"
#include "relaxation.h"
#include "constants.h"
#include "dynamical.h"
#include "error.h"
#include "fcs_phonon.h"
#include "integration.h"
#include "kpoint.h"
#include "mathfunctions.h"
#include "memory.h"
#include "parsephon.h"
#include "phonon_dos.h"
#include "selfenergy.h"
#include "symmetry_core.h"
#include "system.h"
#include "thermodynamics.h"
#include "write_phonons.h"
#include <boost/lexical_cast.hpp>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace PHON_NS;

Relaxation::Relaxation(PHON *phon) : Pointers(phon)
{
    set_default_variables();
}

Relaxation::~Relaxation()
{
    deallocate_variables();
};

void Relaxation::set_default_variables()
{
    im = std::complex<double>(0.0, 1.0);
    quartic_mode = 0;
    ks_analyze_mode = false;
    atom_project_mode = false;
    calc_realpart = false;
    calc_fstate_omega = false;
    calc_fstate_k = false;
    print_V3 = false;
    use_tuned_ver = true;
    spectral_func = false;
    use_triplet_symmetry = true;
    is_imaginary = nullptr;
    vec_for_v3 = nullptr;
    vec_for_v4 = nullptr;
    invmass_for_v3 = nullptr;
    invmass_for_v4 = nullptr;
    evec_index = nullptr;
    evec_index4 = nullptr;
    fcs_group = nullptr;
    fcs_group2 = nullptr;
    exp_phase = nullptr;
    exp_phase3 = nullptr;
}

void Relaxation::deallocate_variables()
{
    if (is_imaginary) {
        memory->deallocate(is_imaginary);
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
    if (evec_index) {
        memory->deallocate(evec_index);
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
}


void Relaxation::setup_relaxation()
{
    nk = kpoint->nk;
    ns = dynamical->neval;
    nks = ns * nk;
    int nk_tmp[3];

    if (mympi->my_rank == 0) {
        std::cout << std::endl;
        std::cout << " -----------------------------------------------------------------" << std::endl << std::endl;
        std::cout << " Now, move on to phonon lifetime calculations." << std::endl;
    }

    setup_mode_analysis();
    setup_cubic();
    sym_permutation = true;

    if (ks_analyze_mode) {

        if (kpoint->kpoint_mode == 2 && use_triplet_symmetry) {
            use_triplet_symmetry = false;
            if (mympi->my_rank == 0) {
                std::cout << std::endl;
                std::cout << " TRISYM was automatically set to 0." << std::endl;
                std::cout << std::endl;
            }
        }

        if (calc_fstate_omega) sym_permutation = false;

        if (quartic_mode > 0) {

            // This is for quartic vertexes.

            if (mympi->my_rank == 0) {
                std::cout << " QUARTIC = 1 : Frequency shift due to the loop diagram associated with" << std::endl;
                std::cout << "               quartic anharmonicity will be calculated." << std::endl;
                std::cout << "               Please check the accuracy of the quartic IFCs " << std::endl;
                std::cout << "               before doing serious calculations." << std::endl;
                std::cout << std::endl;
            }

            setup_quartic();
        }

        if (calc_realpart && integration->ismear != 0) {
            error->exit("setup_relaxation",
                        "Sorry. REALPART = 1 can be used only with ISMEAR = 0");
        }

        if (spectral_func && integration->ismear != -1) {
            error->exit("setup_relaxation",
                        "Sorry. SELF_W = 1 can be used only with the tetrahedron method (ISMEAR = -1).");
        }

        dynamical->modify_eigenvectors();
    }

    if (calc_fstate_k) {
        use_tuned_ver = false;
    } else {
        use_tuned_ver = true;
        nk_tmp[0] = kpoint->nkx;
        nk_tmp[1] = kpoint->nky;
        nk_tmp[2] = kpoint->nkz;
        store_exponential_for_acceleration(nk_tmp, nk_represent,
                                           exp_phase, exp_phase3);
    }

    if (phon->mode == "RTA") {
        detect_imaginary_branches(dynamical->eval_phonon);
    }
}


void Relaxation::detect_imaginary_branches(double **eval)
{
    int ik, is;
    auto nk = kpoint->nk;
    auto ns = dynamical->neval;
    auto nks = ns * nk;
    int knum;
    double omega;
    int ndup;

    memory->allocate(is_imaginary, kpoint->nk_irred, ns);

    bool is_anyof_imaginary = false;

    for (ik = 0; ik < kpoint->nk_irred; ++ik) {
        for (is = 0; is < ns; ++is) {
            knum = kpoint->kpoint_irred_all[ik][0].knum;
            omega = eval[knum][is];

            if (omega < 0.0) {
                is_imaginary[ik][is] = true;
                is_anyof_imaginary = true;
            } else {
                is_imaginary[ik][is] = false;
            }
        }
    }
    if (mympi->my_rank == 0) {

        if (is_anyof_imaginary) {
            int count = 0;
            std::cout << std::endl;
            std::cout << " WARNING: Imaginary frequency detected at the following branches:" << std::endl;
            for (ik = 0; ik < kpoint->nk_irred; ++ik) {
                for (is = 0; is < ns; ++is) {
                    if (is_imaginary[ik][is]) {
                        ndup = kpoint->kpoint_irred_all[ik].size();
                        count += ndup;
                        for (int i = 0; i < ndup; ++i) {
                            knum = kpoint->kpoint_irred_all[ik][i].knum;
                            omega = eval[knum][is];
                            for (int j = 0; j < 3; ++j) {
                                std::cout << std::setw(15) << kpoint->xk[knum][j];
                            }
                            std::cout << std::setw(4) << is + 1 << " :"
                                << std::setw(10) << std::fixed
                                << writes->in_kayser(omega) << " (cm^-1)" << std::endl;
                            std::cout << std::scientific;
                        }
                    }
                }
            }
            std::cout << std::setw(5) << count << " imaginary branches out of "
                << std::setw(5) << nks << " total branches." << std::endl;
            std::cout << std::endl;
            std::cout << " Phonon-phonon scattering rate and thermal conductivity involving these" << std::endl;
            std::cout << " imaginary branches will be treated as zero in the following calculations." << std::endl;
            std::cout << " If imaginary branches are acoustic phonons at Gamma point (0, 0, 0), " << std::endl;
            std::cout << " you can safely ignore this message." << std::endl << std::endl << std::flush;
        }
    }
}

void Relaxation::prepare_relative_vector(const std::vector<FcsArrayWithCell> &fcs_in,
                                         const unsigned int N,
                                         double ***vec_out)
{
    int i, j, k;
    int ix, iy, iz;

    double vec[3];
    double **xshift_s;

    std::vector<unsigned int> atm_super, atm_prim;
    std::vector<unsigned int> xyz;
    std::vector<unsigned int> cells;

    double mat_convert[3][3];

    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            mat_convert[i][j] = 0.0;
            for (k = 0; k < 3; ++k) {
                mat_convert[i][j] += system->rlavec_p[i][k] * system->lavec_s_anharm[k][j];
            }
        }
    }

    memory->allocate(xshift_s, 27, 3);

    for (i = 0; i < 3; ++i) xshift_s[0][i] = 0.0;

    unsigned int icell = 0;

    for (ix = -1; ix <= 1; ++ix) {
        for (iy = -1; iy <= 1; ++iy) {
            for (iz = -1; iz <= 1; ++iz) {
                if (ix == 0 && iy == 0 && iz == 0) continue;

                ++icell;

                xshift_s[icell][0] = static_cast<double>(ix);
                xshift_s[icell][1] = static_cast<double>(iy);
                xshift_s[icell][2] = static_cast<double>(iz);
            }
        }
    }

    unsigned int atm_p, atm_s;
    unsigned int tran_tmp;
    unsigned int icount = 0;

    for (const auto &it : fcs_in) {

        atm_super.clear();
        atm_prim.clear();
        xyz.clear();
        cells.clear();

        for (i = 0; i < it.pairs.size(); ++i) {
            atm_p = it.pairs[i].index / 3;
            tran_tmp = it.pairs[i].tran;
            atm_s = system->map_p2s_anharm[atm_p][tran_tmp];

            atm_prim.push_back(atm_p);
            atm_super.push_back(atm_s);
            cells.push_back(it.pairs[i].cell_s);
        }


        for (i = 0; i < N - 1; ++i) {

            for (j = 0; j < 3; ++j) {
                vec[j] = system->xr_s_anharm[atm_super[i + 1]][j] + xshift_s[cells[i + 1]][j]
                    - system->xr_s_anharm[system->map_p2s_anharm[atm_prim[i + 1]][0]][j];
            }

            rotvec(vec, vec, mat_convert);

            for (j = 0; j < 3; ++j) {
                //   vec_out[j][i][icount] = vec[j];
                vec_out[icount][i][j] = vec[j];
            }
        }
        ++icount;
    }
    memory->deallocate(xshift_s);
}

void Relaxation::prepare_group_of_force_constants(const std::vector<FcsArrayWithCell> &fcs_in,
                                                  const unsigned int N,
                                                  int &number_of_groups,
                                                  std::vector<double> *&fcs_group_out)
{
    // Find the number of groups which has different evecs.

    unsigned int i;
    std::vector<int> arr_old, arr_tmp;

    number_of_groups = 0;

    arr_old.clear();
    for (i = 0; i < N; ++i) {
        arr_old.push_back(-1);
    }

    for (const auto &it : fcs_in) {

        arr_tmp.clear();

        for (i = 0; i < it.pairs.size(); ++i) {
            arr_tmp.push_back(it.pairs[i].index);
        }

        if (arr_tmp != arr_old) {
            ++number_of_groups;
            arr_old.clear();
            arr_old.reserve(arr_tmp.size());
            std::copy(arr_tmp.begin(), arr_tmp.end(), std::back_inserter(arr_old));
        }
    }

    memory->allocate(fcs_group_out, number_of_groups);

    int igroup = -1;

    arr_old.clear();
    for (i = 0; i < N; ++i) {
        arr_old.push_back(-1);
    }

    for (const auto &it : fcs_in) {

        arr_tmp.clear();

        for (i = 0; i < it.pairs.size(); ++i) {
            arr_tmp.push_back(it.pairs[i].index);
        }

        if (arr_tmp != arr_old) {
            ++igroup;
            arr_old.clear();
            arr_old.reserve(arr_tmp.size());
            std::copy(arr_tmp.begin(), arr_tmp.end(), std::back_inserter(arr_old));
        }

        fcs_group_out[igroup].push_back(it.fcs_val);
    }
}

void Relaxation::setup_mode_analysis()
{
    // Judge if ks_analyze_mode should be turned on or not.

    unsigned int i;

    if (mympi->my_rank == 0) {
        if (!ks_input.empty()) {
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
            int knum_tmp;

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
                    knum_tmp = kpoint->get_knum(ktmp[0], ktmp[1], ktmp[2]);

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

            ks_analyze_mode = true;
            ifs_ks.close();

        } else {

            ks_analyze_mode = false;

        }
    }
    MPI_Bcast(&ks_analyze_mode, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD);
    MPI_Bcast(&calc_realpart, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD);
    MPI_Bcast(&atom_project_mode, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD);
    MPI_Bcast(&calc_fstate_k, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD);
    MPI_Bcast(&calc_fstate_omega, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD);

    unsigned int *kslist_arr;
    unsigned int nlist;
    double **vec_tmp;
    unsigned int *mode_tmp;

    if (kpoint->kpoint_mode == 3) {
        int j;

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
}

std::complex<double> Relaxation::V3(const unsigned int ks[3])
{
    unsigned int i, j;
    unsigned int kn[3], sn[3];
    unsigned int nsize_group;

    double phase, omega[3];
    double phase3[3];

    std::complex<double> ret = std::complex<double>(0.0, 0.0);
    std::complex<double> ret_in, vec_tmp;

    int iloc, loc[3];
    int ii;
    const auto inv2pi = 1.0 / (2.0 * pi);
    double dnk_represent = static_cast<double>(nk_represent);

    dnk_represent *= inv2pi;

    for (i = 0; i < 3; ++i) {
        kn[i] = ks[i] / ns;
        sn[i] = ks[i] % ns;
        omega[i] = dynamical->eval_phonon[kn[i]][sn[i]];
    }
    // Return zero if any of the involving phonon has imaginary frequency
    if (omega[0] < 0.0 || omega[1] < 0.0 || omega[2] < 0.0) return 0.0;

    unsigned int ielem = 0;

    if (use_tuned_ver) {
        if (tune_type == 0) {

            // Tuned version used when nk1=nk2=nk3.
            for (i = 0; i < ngroup; ++i) {

                vec_tmp
                    = dynamical->evec_phonon[kn[0]][sn[0]][evec_index[i][0]]
                    * dynamical->evec_phonon[kn[1]][sn[1]][evec_index[i][1]]
                    * dynamical->evec_phonon[kn[2]][sn[2]][evec_index[i][2]]
                    * invmass_for_v3[i];

                ret_in = std::complex<double>(0.0, 0.0);

                nsize_group = fcs_group[i].size();

                for (j = 0; j < nsize_group; ++j) {

                    phase
                        = vec_for_v3[ielem][0][0] * kpoint->xk[kn[1]][0]
                        + vec_for_v3[ielem][0][1] * kpoint->xk[kn[1]][1]
                        + vec_for_v3[ielem][0][2] * kpoint->xk[kn[1]][2]
                        + vec_for_v3[ielem][1][0] * kpoint->xk[kn[2]][0]
                        + vec_for_v3[ielem][1][1] * kpoint->xk[kn[2]][1]
                        + vec_for_v3[ielem][1][2] * kpoint->xk[kn[2]][2];

                    iloc = nint(phase * dnk_represent) % nk_represent + nk_represent - 1;
                    ret_in += fcs_group[i][j] * exp_phase[iloc];

                    ++ielem;
                }
                ret += ret_in * vec_tmp;
            }
        } else if (tune_type == 1) {

            // Tuned version used when nk1=nk2=nk3 is not met.
            for (i = 0; i < ngroup; ++i) {

                vec_tmp
                    = dynamical->evec_phonon[kn[0]][sn[0]][evec_index[i][0]]
                    * dynamical->evec_phonon[kn[1]][sn[1]][evec_index[i][1]]
                    * dynamical->evec_phonon[kn[2]][sn[2]][evec_index[i][2]]
                    * invmass_for_v3[i];

                ret_in = std::complex<double>(0.0, 0.0);

                nsize_group = fcs_group[i].size();

                for (j = 0; j < nsize_group; ++j) {

                    for (ii = 0; ii < 3; ++ii) {
                        phase3[ii]
                            = vec_for_v3[ielem][0][ii] * kpoint->xk[kn[1]][ii]
                            + vec_for_v3[ielem][1][ii] * kpoint->xk[kn[2]][ii];

                        loc[ii] = nint(phase3[ii] * dnk[ii] * inv2pi) % nk_grid[ii] + nk_grid[ii] - 1;
                    }

                    ret_in += fcs_group[i][j] * exp_phase3[loc[0]][loc[1]][loc[2]];

                    ++ielem;
                }
                ret += ret_in * vec_tmp;
            }
        }

    } else {
        // Original version
        for (i = 0; i < ngroup; ++i) {

            vec_tmp
                = dynamical->evec_phonon[kn[0]][sn[0]][evec_index[i][0]]
                * dynamical->evec_phonon[kn[1]][sn[1]][evec_index[i][1]]
                * dynamical->evec_phonon[kn[2]][sn[2]][evec_index[i][2]]
                * invmass_for_v3[i];

            ret_in = std::complex<double>(0.0, 0.0);

            nsize_group = fcs_group[i].size();

            for (j = 0; j < nsize_group; ++j) {

                phase
                    = vec_for_v3[ielem][0][0] * kpoint->xk[kn[1]][0]
                    + vec_for_v3[ielem][0][1] * kpoint->xk[kn[1]][1]
                    + vec_for_v3[ielem][0][2] * kpoint->xk[kn[1]][2]
                    + vec_for_v3[ielem][1][0] * kpoint->xk[kn[2]][0]
                    + vec_for_v3[ielem][1][1] * kpoint->xk[kn[2]][1]
                    + vec_for_v3[ielem][1][2] * kpoint->xk[kn[2]][2];

                ret_in += fcs_group[i][j] * std::exp(im * phase);

                ++ielem;
            }
            ret += ret_in * vec_tmp;
        }
    }

    return ret / std::sqrt(omega[0] * omega[1] * omega[2]);
}


std::complex<double> Relaxation::V4(const unsigned int ks[4])
{
    int ii;
    unsigned int i, j;
    unsigned int kn[4], sn[4];

    double phase, phase3[3];
    double omega[4];

    std::complex<double> ctmp, ret_in, vec_tmp;
    std::complex<double> ret = std::complex<double>(0.0, 0.0);

    int iloc, loc[3];
    const auto inv2pi = 1.0 / (2.0 * pi);
    double dnk_represent = static_cast<double>(nk_represent);

    dnk_represent *= inv2pi;

    for (i = 0; i < 4; ++i) {
        kn[i] = ks[i] / ns;
        sn[i] = ks[i] % ns;
        omega[i] = dynamical->eval_phonon[kn[i]][sn[i]];
    }

    unsigned int ielem = 0;

    if (use_tuned_ver) {
        if (tune_type == 0) {
            for (i = 0; i < ngroup2; ++i) {

                vec_tmp
                    = dynamical->evec_phonon[kn[0]][sn[0]][evec_index4[i][0]]
                    * dynamical->evec_phonon[kn[1]][sn[1]][evec_index4[i][1]]
                    * dynamical->evec_phonon[kn[2]][sn[2]][evec_index4[i][2]]
                    * dynamical->evec_phonon[kn[3]][sn[3]][evec_index4[i][3]]
                    * invmass_for_v4[i];

                ret_in = std::complex<double>(0.0, 0.0);

                for (j = 0; j < fcs_group2[i].size(); ++j) {
                    phase
                        = vec_for_v4[ielem][0][0] * kpoint->xk[kn[1]][0]
                        + vec_for_v4[ielem][0][1] * kpoint->xk[kn[1]][1]
                        + vec_for_v4[ielem][0][2] * kpoint->xk[kn[1]][2]
                        + vec_for_v4[ielem][1][0] * kpoint->xk[kn[2]][0]
                        + vec_for_v4[ielem][1][1] * kpoint->xk[kn[2]][1]
                        + vec_for_v4[ielem][1][2] * kpoint->xk[kn[2]][2]
                        + vec_for_v4[ielem][2][0] * kpoint->xk[kn[3]][0]
                        + vec_for_v4[ielem][2][1] * kpoint->xk[kn[3]][1]
                        + vec_for_v4[ielem][2][2] * kpoint->xk[kn[3]][2];

                    iloc = nint(phase * dnk_represent) % nk_represent + nk_represent - 1;

                    ctmp = fcs_group2[i][j] * exp_phase[iloc];
                    ret_in += ctmp;

                    ++ielem;
                }
                ret += ret_in * vec_tmp;
            }
        } else if (tune_type == 1) {
            for (i = 0; i < ngroup2; ++i) {

                vec_tmp
                    = dynamical->evec_phonon[kn[0]][sn[0]][evec_index4[i][0]]
                    * dynamical->evec_phonon[kn[1]][sn[1]][evec_index4[i][1]]
                    * dynamical->evec_phonon[kn[2]][sn[2]][evec_index4[i][2]]
                    * dynamical->evec_phonon[kn[3]][sn[3]][evec_index4[i][3]]
                    * invmass_for_v4[i];

                ret_in = std::complex<double>(0.0, 0.0);

                for (j = 0; j < fcs_group2[i].size(); ++j) {

                    for (ii = 0; ii < 3; ++ii) {
                        phase3[ii]
                            = vec_for_v4[ielem][0][ii] * kpoint->xk[kn[1]][ii]
                            + vec_for_v4[ielem][1][ii] * kpoint->xk[kn[2]][ii]
                            + vec_for_v4[ielem][2][ii] * kpoint->xk[kn[3]][ii];

                        loc[ii] = nint(phase3[ii] * dnk[ii] * inv2pi) % nk_grid[ii] + nk_grid[ii] - 1;
                    }

                    ctmp = fcs_group2[i][j] * exp_phase3[loc[0]][loc[1]][loc[2]];
                    ret_in += ctmp;

                    ++ielem;
                }
                ret += ret_in * vec_tmp;
            }
        }

    } else {
        for (i = 0; i < ngroup2; ++i) {

            vec_tmp
                = dynamical->evec_phonon[kn[0]][sn[0]][evec_index4[i][0]]
                * dynamical->evec_phonon[kn[1]][sn[1]][evec_index4[i][1]]
                * dynamical->evec_phonon[kn[2]][sn[2]][evec_index4[i][2]]
                * dynamical->evec_phonon[kn[3]][sn[3]][evec_index4[i][3]]
                * invmass_for_v4[i];

            ret_in = std::complex<double>(0.0, 0.0);

            for (j = 0; j < fcs_group2[i].size(); ++j) {
                phase
                    = vec_for_v4[ielem][0][0] * kpoint->xk[kn[1]][0]
                    + vec_for_v4[ielem][0][1] * kpoint->xk[kn[1]][1]
                    + vec_for_v4[ielem][0][2] * kpoint->xk[kn[1]][2]
                    + vec_for_v4[ielem][1][0] * kpoint->xk[kn[2]][0]
                    + vec_for_v4[ielem][1][1] * kpoint->xk[kn[2]][1]
                    + vec_for_v4[ielem][1][2] * kpoint->xk[kn[2]][2]
                    + vec_for_v4[ielem][2][0] * kpoint->xk[kn[3]][0]
                    + vec_for_v4[ielem][2][1] * kpoint->xk[kn[3]][1]
                    + vec_for_v4[ielem][2][2] * kpoint->xk[kn[3]][2];

                ctmp = fcs_group2[i][j] * std::exp(im * phase);
                ret_in += ctmp;

                ++ielem;
            }
            ret += ret_in * vec_tmp;
        }
    }


    return ret / std::sqrt(omega[0] * omega[1] * omega[2] * omega[3]);
}


std::complex<double> Relaxation::V3_mode(int mode,
                                         double *xk2,
                                         double *xk3,
                                         int is,
                                         int js,
                                         double **eval,
                                         std::complex<double> ***evec)
{
    int i, j;
    int nsize_group;

    double phase;
    std::complex<double> ctmp = std::complex<double>(0.0, 0.0);
    std::complex<double> vec_tmp, ret_in;

    // Return zero if any of the involving phonon has imaginary frequency
    if (eval[0][mode] < 0.0 || eval[1][is] < 0.0 || eval[2][js] < 0.0) return 0.0;

    unsigned int ielem = 0;

    for (i = 0; i < ngroup; ++i) {

        vec_tmp
            = evec[0][mode][evec_index[i][0]]
            * evec[1][is][evec_index[i][1]]
            * evec[2][js][evec_index[i][2]]
            * invmass_for_v3[i];

        ret_in = std::complex<double>(0.0, 0.0);

        nsize_group = fcs_group[i].size();

        for (j = 0; j < nsize_group; ++j) {

            phase
                = vec_for_v3[ielem][0][0] * xk2[0]
                + vec_for_v3[ielem][0][1] * xk2[1]
                + vec_for_v3[ielem][0][2] * xk2[2]
                + vec_for_v3[ielem][1][0] * xk3[0]
                + vec_for_v3[ielem][1][1] * xk3[1]
                + vec_for_v3[ielem][1][2] * xk3[2];

            ret_in += fcs_group[i][j] * std::exp(im * phase);

            ++ielem;
        }
        ctmp += ret_in * vec_tmp;
    }

    return ctmp / std::sqrt(eval[0][mode] * eval[1][is] * eval[2][js]);
}


void Relaxation::calc_damping_smearing(const unsigned int N,
                                       double *T,
                                       const double omega,
                                       const unsigned int ik_in,
                                       const unsigned int snum,
                                       double *ret)
{
    // This function returns the imaginary part of phonon self-energy 
    // for the given frequency omega.
    // Lorentzian or Gaussian smearing will be used.
    // This version employs the crystal symmetry to reduce the computational cost

    unsigned int i;
    int ik;
    unsigned int is, js;
    unsigned int arr[3];

    int k1, k2;

    double T_tmp;
    double n1, n2;
    double omega_inner[2];

    int knum, knum_minus;
    double multi;

    for (i = 0; i < N; ++i) ret[i] = 0.0;

    double **v3_arr;
    double ***delta_arr;
    double ret_tmp;

    double f1, f2;

    double epsilon = integration->epsilon;

    std::vector<KsListGroup> triplet;

    get_unique_triplet_k(ik_in,
                         use_triplet_symmetry,
                         sym_permutation,
                         triplet);

    int npair_uniq = triplet.size();

    memory->allocate(v3_arr, npair_uniq, ns * ns);
    memory->allocate(delta_arr, npair_uniq, ns * ns, 2);

    knum = kpoint->kpoint_irred_all[ik_in][0].knum;
    knum_minus = kpoint->knum_minus[knum];
#ifdef _OPENMP
#pragma omp parallel for private(multi, arr, k1, k2, is, js, omega_inner)
#endif
    for (ik = 0; ik < npair_uniq; ++ik) {
        multi = static_cast<double>(triplet[ik].group.size());

        arr[0] = ns * knum_minus + snum;

        k1 = triplet[ik].group[0].ks[0];
        k2 = triplet[ik].group[0].ks[1];

        for (is = 0; is < ns; ++is) {
            arr[1] = ns * k1 + is;
            omega_inner[0] = dynamical->eval_phonon[k1][is];

            for (js = 0; js < ns; ++js) {
                arr[2] = ns * k2 + js;
                omega_inner[1] = dynamical->eval_phonon[k2][js];

                v3_arr[ik][ns * is + js] = std::norm(V3(arr)) * multi;

                if (integration->ismear == 0) {
                    delta_arr[ik][ns * is + js][0]
                        = delta_lorentz(omega - omega_inner[0] - omega_inner[1], epsilon)
                        - delta_lorentz(omega + omega_inner[0] + omega_inner[1], epsilon);
                    delta_arr[ik][ns * is + js][1]
                        = delta_lorentz(omega - omega_inner[0] + omega_inner[1], epsilon)
                        - delta_lorentz(omega + omega_inner[0] - omega_inner[1], epsilon);
                } else if (integration->ismear == 1) {
                    delta_arr[ik][ns * is + js][0]
                        = delta_gauss(omega - omega_inner[0] - omega_inner[1], epsilon)
                        - delta_gauss(omega + omega_inner[0] + omega_inner[1], epsilon);
                    delta_arr[ik][ns * is + js][1]
                        = delta_gauss(omega - omega_inner[0] + omega_inner[1], epsilon)
                        - delta_gauss(omega + omega_inner[0] - omega_inner[1], epsilon);
                }
            }
        }
    }

    for (i = 0; i < N; ++i) {
        T_tmp = T[i];
        ret_tmp = 0.0;
#ifdef _OPENMP
#pragma omp parallel for private(k1, k2, is, js, omega_inner, n1, n2, f1, f2), reduction(+:ret_tmp)
#endif
        for (ik = 0; ik < npair_uniq; ++ik) {

            k1 = triplet[ik].group[0].ks[0];
            k2 = triplet[ik].group[0].ks[1];

            for (is = 0; is < ns; ++is) {

                omega_inner[0] = dynamical->eval_phonon[k1][is];

                for (js = 0; js < ns; ++js) {

                    omega_inner[1] = dynamical->eval_phonon[k2][js];

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

                    ret_tmp += v3_arr[ik][ns * is + js]
                        * (n1 * delta_arr[ik][ns * is + js][0]
                            - n2 * delta_arr[ik][ns * is + js][1]);
                }
            }
        }
        ret[i] = ret_tmp;
    }

    memory->deallocate(v3_arr);
    memory->deallocate(delta_arr);
    triplet.clear();

    for (i = 0; i < N; ++i) ret[i] *= pi * std::pow(0.5, 4) / static_cast<double>(nk);
}

void Relaxation::calc_damping_tetrahedron(const unsigned int N,
                                          double *T,
                                          const double omega,
                                          const unsigned int ik_in,
                                          const unsigned int snum,
                                          double *ret)
{
    // This function returns the imaginary part of phonon self-energy 
    // for the given frequency omega.
    // Tetrahedron method will be used.
    // This version employs the crystal symmetry to reduce the computational cost

    int ik, ib;
    int ns2 = ns * ns;

    unsigned int i;
    unsigned int jk;
    unsigned int is, js;
    unsigned int k1, k2;
    unsigned int arr[3];

    double T_tmp;
    double n1, n2;
    double f1, f2;

    double xk_tmp[3];
    double omega_inner[2];

    double ret_tmp;

    int *kmap_identity;
    double **energy_tmp;
    double **weight_tetra;
    double **v3_arr;
    double ***delta_arr;

    std::vector<KsListGroup> triplet;

    for (i = 0; i < N; ++i) ret[i] = 0.0;

    get_unique_triplet_k(ik_in,
                         use_triplet_symmetry,
                         sym_permutation,
                         triplet);

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
                integration->calc_weight_tetrahedron(nk, kmap_identity,
                                                     weight_tetra[i], energy_tmp[i], omega);
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

                    v3_arr[ik][ib] = std::norm(V3(arr));
                } else {
                    v3_arr[ik][ib] = 0.0;
                }
            }
        }

        memory->deallocate(energy_tmp);
        memory->deallocate(weight_tetra);
    }

    for (i = 0; i < N; ++i) {
        T_tmp = T[i];
        ret_tmp = 0.0;
#ifdef _OPENMP
#pragma omp parallel for private(k1, k2, is, js, omega_inner, n1, n2, f1, f2), reduction(+:ret_tmp)
#endif
        for (ik = 0; ik < npair_uniq; ++ik) {

            k1 = triplet[ik].group[0].ks[0];
            k2 = triplet[ik].group[0].ks[1];

            for (is = 0; is < ns; ++is) {

                omega_inner[0] = dynamical->eval_phonon[k1][is];

                for (js = 0; js < ns; ++js) {

                    omega_inner[1] = dynamical->eval_phonon[k2][js];

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

                    ret_tmp += v3_arr[ik][ns * is + js]
                        * (n1 * delta_arr[ik][ns * is + js][0]
                            - n2 * delta_arr[ik][ns * is + js][1]);
                }
            }
        }
        ret[i] = ret_tmp;
    }

    memory->deallocate(v3_arr);
    memory->deallocate(delta_arr);
    memory->deallocate(kmap_identity);

    for (i = 0; i < N; ++i) ret[i] *= pi * std::pow(0.5, 4);
}


void Relaxation::calc_frequency_resolved_final_state(const unsigned int N,
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

    double epsilon = integration->epsilon;

    std::vector<KsListGroup> triplet;

    get_unique_triplet_k(ik_in,
                         use_triplet_symmetry,
                         sym_permutation,
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

                v3_tmp = std::norm(V3(arr));

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


void Relaxation::calc_frequency_resolved_final_state_tetrahedron(const unsigned int N,
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

    int ns2 = ns * ns;

    int *kmap_identity;
    double **energy_tmp;
    double **weight_tetra;
    double **v3_arr;
    double ***delta_arr;
    double prod_tmp[2];

    double epsilon = integration->epsilon;

    std::vector<KsListGroup> triplet;

    get_unique_triplet_k(ik_in,
                         use_triplet_symmetry,
                         sym_permutation,
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

                    v3_arr[ik][ib] = std::norm(V3(arr));
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

void Relaxation::perform_mode_analysis()
{
    unsigned int i, j;
    unsigned int NT;
    unsigned int knum, snum;

    double Tmax = system->Tmax;
    double Tmin = system->Tmin;
    double dT = system->dT;
    double omega;
    double *T_arr;

    std::ofstream ofs_linewidth, ofs_shift, ofs_fstate_w;
    std::string file_linewidth, file_shift, file_fstate_w;


    NT = static_cast<unsigned int>((Tmax - Tmin) / dT) + 1;
    memory->allocate(T_arr, NT);
    for (i = 0; i < NT; ++i) T_arr[i] = Tmin + static_cast<double>(i) * dT;

    double epsilon = integration->epsilon;

    if (print_V3) {

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

            get_unique_triplet_k(ik_irred,
                                 use_triplet_symmetry,
                                 sym_permutation,
                                 triplet);
            nk_size = triplet.size();

            memory->allocate(v3norm, nk_size, ns * ns);

            calc_V3norm2(ik_irred, snum, v3norm);

            if (mympi->my_rank == 0) {

                file_V3 = input->job_title + ".V3." + std::to_string(i + 1);
                ofs_V3.open(file_V3.c_str(), std::ios::out);
                if (!ofs_V3)
                    error->exit("perform_mode_analysis",
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

    } else if (calc_fstate_k) {

        // Momentum-resolved final state amplitude
        print_momentum_resolved_final_state(NT, T_arr, epsilon);

    } else if (calc_fstate_omega) {

        print_frequency_resolved_final_state(NT, T_arr);

    } else if (spectral_func) {

        int ik_irred, iomega, iT;
        double **self3_imag, **self3_real;
        std::string file_self;
        std::ofstream ofs_self;
        double *omega_array;
        double Omega_min = dos->emin;
        double Omega_max = dos->emax;
        double delta_omega = dos->delta_e;
        double T_now, omega2[2];

        int nomega = static_cast<unsigned int>((Omega_max - Omega_min) / delta_omega) + 1;

        memory->allocate(omega_array, nomega);
        memory->allocate(self3_imag, NT, nomega);
        memory->allocate(self3_real, NT, nomega);

        for (i = 0; i < nomega; ++i) {
            omega_array[i] = Omega_min + delta_omega * static_cast<double>(i);
            omega_array[i] *= time_ry / Hz_to_kayser;
        }

        for (i = 0; i < relaxation->kslist.size(); ++i) {
            knum = relaxation->kslist[i] / ns;
            snum = relaxation->kslist[i] % ns;
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
                if (!ofs_self) error->exit("perform_mode_analysis", "Cannot open file file_shift");

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

                calc_self3omega_tetrahedron(T_now,
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

    } else {

        double *damping_a;
        double omega_shift;
        std::complex<double> *self_tadpole;
        std::complex<double> *self_a, *self_b, *self_c, *self_d, *self_e;
        std::complex<double> *self_f, *self_g, *self_h, *self_i, *self_j;

        if (mympi->my_rank == 0) {
            std::cout << std::endl;
            std::cout << " Calculate the line width (FWHM) of phonons" << std::endl;
            std::cout << " due to 3-phonon interactions for given "
                << kslist.size() << " modes." << std::endl;

            if (calc_realpart) {
                if (quartic_mode == 1) {
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

            if (quartic_mode == 2) {
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

        if (quartic_mode == 2) {
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
                calc_damping_tetrahedron(NT, T_arr, omega, ik_irred, snum, damping_a);
            } else {
                selfenergy->selfenergy_a(NT, T_arr, omega, knum, snum, self_a);
                for (j = 0; j < NT; ++j) damping_a[j] = self_a[j].imag();
            }
            if (quartic_mode == 2) {
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
                    error->exit("perform_mode_analysis",
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
                if (quartic_mode == 2) ofs_linewidth << ", 2*Gamma4(cm^-1) <-- specific diagram only";
                ofs_linewidth << std::endl;

                for (j = 0; j < NT; ++j) {
                    ofs_linewidth << std::setw(10) << T_arr[j]
                        << std::setw(15) << writes->in_kayser(2.0 * damping_a[j]);

                    if (quartic_mode == 2) {
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

                if (quartic_mode == 1) {
                    selfenergy->selfenergy_b(NT, T_arr, omega, knum, snum, self_b);
                }

                if (mympi->my_rank == 0) {

                    file_shift = input->job_title + ".Shift." + std::to_string(i + 1);
                    ofs_shift.open(file_shift.c_str(), std::ios::out);
                    if (!ofs_shift)
                        error->exit("perform_mode_analysis",
                                    "Cannot open file file_shift");

                    ofs_shift << "# xk = ";

                    for (j = 0; j < 3; ++j) {
                        ofs_shift << std::setw(15) << kpoint->xk[knum][j];
                    }
                    ofs_shift << std::endl;
                    ofs_shift << "# mode = " << snum + 1 << std::endl;
                    ofs_shift << "# Frequency = " << writes->in_kayser(omega) << std::endl;
                    ofs_shift << "## T[K], Shift3 (cm^-1) (tadpole), Shift3 (cm^-1) (bubble)";
                    if (quartic_mode == 1) ofs_shift << ", Shift4 (cm^-1) (loop)";
                    ofs_shift << ", Shifted frequency (cm^-1)";
                    ofs_shift << std::endl;


                    for (j = 0; j < NT; ++j) {
                        ofs_shift << std::setw(10) << T_arr[j];
                        ofs_shift << std::setw(15) << writes->in_kayser(-self_tadpole[j].real());
                        ofs_shift << std::setw(15) << writes->in_kayser(-self_a[j].real());

                        omega_shift = omega - self_tadpole[j].real() - self_a[j].real();

                        if (quartic_mode == 1) {
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

        memory->deallocate(damping_a);
        memory->deallocate(self_a);
        memory->deallocate(self_b);
        memory->deallocate(self_tadpole);

        if (quartic_mode == 2) {
            memory->deallocate(self_c);
            memory->deallocate(self_d);
            memory->deallocate(self_e);
            memory->deallocate(self_f);
            memory->deallocate(self_g);
            memory->deallocate(self_h);
            memory->deallocate(self_i);
            memory->deallocate(self_j);
        }

    }
    memory->deallocate(T_arr);
}

void Relaxation::print_frequency_resolved_final_state(const unsigned int NT,
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

void Relaxation::print_momentum_resolved_final_state(const unsigned int NT,
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
            error->exit("setup_relaxation", "Invalid ksum_mode");
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
            dynamical->eval_k(xk_plane[j], kvec_plane[j],
                              fcs_phonon->fc2_ext, eval[j], evec[0], false);
            for (k = 0; k < ns; ++k) {
                eval[j][k] = dynamical->freq(eval[j][k]);
            }
        }

        // Loop over k points to analyze the final state amplitude

        for (j = 0; j < kslist_fstate_k.size(); ++j) {

            for (k = 0; k < 3; ++k) xk1[k] = -kslist_fstate_k[j].xk[k];
            mode = kslist_fstate_k[j].nmode;

            for (k = 0; k < 3; ++k) kvec[k] = dynamical->fold(xk1[k]);
            rotvec(kvec, kvec, system->rlavec_p, 'T');
            norm = std::sqrt(kvec[0] * kvec[0]
                + kvec[1] * kvec[1]
                + kvec[2] * kvec[2]);

            if (norm > eps) {
                for (k = 0; k < 3; ++k) kvec[k] /= norm;
            }

            for (k = 0; k < 3; ++k) xk1[k] = dynamical->fold(xk1[k]);

            dynamical->eval_k(xk1, kvec, fcs_phonon->fc2_ext, eval_tmp, evec[0], false);
            for (k = 0; k < ns; ++k) eval_tmp[k] = dynamical->freq(eval_tmp[k]);

            // Calculate xk's for the third index that satisfy the momentum conservation

            for (k = 0; k < nk_plane; ++k) {
                for (l = 0; l < 3; ++l) {
                    xk_plane2[k][l] = dynamical->fold(-xk1[l] - xk_plane[k][l]);
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
                dynamical->eval_k(xk_plane2[k], kvec_plane[k],
                                  fcs_phonon->fc2_ext, eval2[k], evec[0], false);
                for (l = 0; l < ns; ++l) {
                    eval2[k][l] = dynamical->freq(eval2[k][l]);
                }
            }

            // Find a list of k points which satisfy the energy conservation

            for (const auto &it : kpoint->kp_planes_tri[i]) {

                // K point indexes for each triangle
                for (k = 0; k < 3; ++k) knum_triangle[k] = it.knum[k];

                for (is = 0; is < ns; ++is) {
                    for (js = 0; js < ns; ++js) {

                        // The case of delta(w1 - w2 - w3) 

                        for (k = 0; k < 3; ++k) {
                            omega_sum[k] = eval_tmp[mode]
                                - eval[knum_triangle[k]][is]
                                - eval2[knum_triangle[k]][js];
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

            for (j = 0; j < 3; ++j) xk_sym[j] = xk_sym[j] - nint(xk_sym[j]);

            diff = 0.0;
            for (j = 0; j < 3; ++j) diff += std::pow(xk_sym[j] - xk1[j], 2);

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

                        V3norm = std::norm(V3_mode(mode, xk_sym, xk3, is, js, eval, evec));

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


bool Relaxation::is_proper(const int isym)
{
    double S[3][3];
    bool ret;

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            S[i][j] = static_cast<double>(symmetry->SymmList[isym].rot[i][j]);
        }
    }

    double det = S[0][0] * (S[1][1] * S[2][2] - S[2][1] * S[1][2])
        - S[1][0] * (S[0][1] * S[2][2] - S[2][1] * S[0][2])
        + S[2][0] * (S[0][1] * S[1][2] - S[1][1] * S[0][2]);

    if (std::abs(det - 1.0) < eps12) {
        ret = true;
    } else if (std::abs(det + 1.0) < eps12) {
        ret = false;
    } else {
        error->exit("is_proper", "This cannot happen.");
    }

    return ret;
}

bool Relaxation::is_symmorphic(const int isym)
{
    double tran[3];

    for (int i = 0; i < 3; ++i) tran[i] = symmetry->SymmList[isym].tran[i];

    bool ret = (std::abs(tran[0]) < eps) && (std::abs(tran[1]) < eps) && (std::abs(tran[2]) < eps);
    return ret;
}

void Relaxation::get_unique_triplet_k(const int ik,
                                      const bool use_triplet_symmetry,
                                      const bool use_permutation_symmetry,
                                      std::vector<KsListGroup> &triplet)
{
    int i, ik1, ik2, isym;
    int num_group_k, tmp;
    int ks_in[2];
    int knum = kpoint->kpoint_irred_all[ik][0].knum;
    bool *flag_found;
    std::vector<KsList> kslist;
    double xk[3], xk1[3], xk2[3];

    memory->allocate(flag_found, kpoint->nk);

    if (use_triplet_symmetry) {
        num_group_k = kpoint->small_group_of_k[ik].size();
    } else {
        num_group_k = 1;
    }

    for (i = 0; i < 3; ++i) xk[i] = kpoint->xk[knum][i];
    for (i = 0; i < kpoint->nk; ++i) flag_found[i] = false;

    triplet.clear();

    for (ik1 = 0; ik1 < nk; ++ik1) {

        for (i = 0; i < 3; ++i) xk1[i] = kpoint->xk[ik1][i];
        for (i = 0; i < 3; ++i) xk2[i] = xk[i] - xk1[i];

        ik2 = kpoint->get_knum(xk2[0], xk2[1], xk2[2]);

        kslist.clear();

        if (ik1 > ik2 && use_permutation_symmetry) continue;

        // Add symmety-connected triplets to kslist
        for (isym = 0; isym < num_group_k; ++isym) {

            ks_in[0] = kpoint->knum_sym(ik1, kpoint->small_group_of_k[ik][isym]);
            ks_in[1] = kpoint->knum_sym(ik2, kpoint->small_group_of_k[ik][isym]);

            if (!flag_found[ks_in[0]]) {
                kslist.emplace_back(2, ks_in, kpoint->small_group_of_k[ik][isym]);
                flag_found[ks_in[0]] = true;
            }

            if (ks_in[0] != ks_in[1] && use_permutation_symmetry && (!flag_found[ks_in[1]])) {
                tmp = ks_in[0];
                ks_in[0] = ks_in[1];
                ks_in[1] = tmp;

                kslist.emplace_back(2, ks_in, kpoint->small_group_of_k[ik][isym]);
                flag_found[ks_in[0]] = true;
            }
        }
        if (!kslist.empty()) {
            triplet.emplace_back(kslist);
        }
    }

    memory->deallocate(flag_found);
}


void Relaxation::calc_V3norm2(const unsigned int ik_in,
                              const unsigned int snum,
                              double **ret)
{
    int ib;
    unsigned int ik;
    unsigned int is, js;
    unsigned int k1, k2;
    unsigned int arr[3];
    unsigned int knum, knum_minus;

    int ns2 = ns * ns;

    double factor = std::pow(0.5, 3) * std::pow(Hz_to_kayser / time_ry, 2);
    std::vector<KsListGroup> triplet;

    knum = kpoint->kpoint_irred_all[ik_in][0].knum;
    knum_minus = kpoint->knum_minus[knum];

    get_unique_triplet_k(ik_in,
                         use_triplet_symmetry,
                         sym_permutation,
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

            ret[ik][ib] = std::norm(V3(arr)) * factor;
        }
    }
}

void Relaxation::setup_cubic()
{
    int i, j, k;
    double *invsqrt_mass_p;

    // Sort force_constant[1] using the operator defined in fcs_phonons.h
    // This sorting is necessary.
    std::sort(fcs_phonon->force_constant_with_cell[1].begin(),
              fcs_phonon->force_constant_with_cell[1].end());
    prepare_group_of_force_constants(fcs_phonon->force_constant_with_cell[1],
                                     3, ngroup, fcs_group);

    memory->allocate(vec_for_v3, fcs_phonon->force_constant_with_cell[1].size(), 2, 3);

    memory->allocate(invmass_for_v3, ngroup);
    memory->allocate(evec_index, ngroup, 3);

    prepare_relative_vector(fcs_phonon->force_constant_with_cell[1], 3, vec_for_v3);

    memory->allocate(invsqrt_mass_p, system->natmin);

    for (i = 0; i < system->natmin; ++i) {
        invsqrt_mass_p[i] = std::sqrt(1.0 / system->mass[system->map_p2s[i][0]]);
    }

    k = 0;
    for (i = 0; i < ngroup; ++i) {
        for (j = 0; j < 3; ++j) {
            evec_index[i][j] = fcs_phonon->force_constant_with_cell[1][k].pairs[j].index;
        }
        invmass_for_v3[i]
            = invsqrt_mass_p[evec_index[i][0] / 3]
            * invsqrt_mass_p[evec_index[i][1] / 3]
            * invsqrt_mass_p[evec_index[i][2] / 3];
        k += fcs_group[i].size();
    }

    memory->deallocate(invsqrt_mass_p);
}

void Relaxation::setup_quartic()
{
    int i, j, k;
    double *invsqrt_mass_p;
    std::sort(fcs_phonon->force_constant_with_cell[2].begin(),
              fcs_phonon->force_constant_with_cell[2].end());
    prepare_group_of_force_constants(fcs_phonon->force_constant_with_cell[2],
                                     4, ngroup2, fcs_group2);

    memory->allocate(vec_for_v4, fcs_phonon->force_constant_with_cell[2].size(), 3, 3);

    memory->allocate(invmass_for_v4, ngroup2);
    memory->allocate(evec_index4, ngroup2, 4);

    prepare_relative_vector(fcs_phonon->force_constant_with_cell[2], 4, vec_for_v4);

    memory->allocate(invsqrt_mass_p, system->natmin);

    for (i = 0; i < system->natmin; ++i) {
        invsqrt_mass_p[i] = std::sqrt(1.0 / system->mass[system->map_p2s[i][0]]);
    }

    k = 0;
    for (i = 0; i < ngroup2; ++i) {
        for (j = 0; j < 4; ++j) {
            evec_index4[i][j] = fcs_phonon->force_constant_with_cell[2][k].pairs[j].index;
        }
        invmass_for_v4[i]
            = invsqrt_mass_p[evec_index4[i][0] / 3]
            * invsqrt_mass_p[evec_index4[i][1] / 3]
            * invsqrt_mass_p[evec_index4[i][2] / 3]
            * invsqrt_mass_p[evec_index4[i][3] / 3];
        k += fcs_group2[i].size();
    }

    memory->deallocate(invsqrt_mass_p);
}

void Relaxation::store_exponential_for_acceleration(const int nk_in[3],
                                                    int &nkrep_out,
                                                    std::complex<double> *exp_out,
                                                    std::complex<double> ***exp3_out)
{
    // For accelerating function V3 and V4 by avoiding continual call of std::exp.

    int i;

    MPI_Bcast(&use_tuned_ver, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD);

    if (use_tuned_ver) {

        nk_grid[0] = nk_in[0];
        nk_grid[1] = nk_in[1];
        nk_grid[2] = nk_in[2];

        for (i = 0; i < 3; ++i) dnk[i] = static_cast<double>(nk_grid[i]);

        if (nk_grid[0] == nk_grid[1] && nk_grid[1] == nk_grid[2]) {
            nkrep_out = nk_grid[0];
            tune_type = 0;

        } else if (nk_grid[0] == nk_grid[1] && nk_grid[2] == 1) {
            nkrep_out = nk_grid[0];
            tune_type = 0;

        } else if (nk_grid[1] == nk_grid[2] && nk_grid[0] == 1) {
            nkrep_out = nk_grid[1];
            tune_type = 0;

        } else if (nk_grid[2] == nk_grid[0] && nk_grid[1] == 1) {
            nkrep_out = nk_grid[2];
            tune_type = 0;

        } else if (nk_grid[0] == 1 && nk_grid[1] == 1) {
            nkrep_out = nk_grid[2];
            tune_type = 0;

        } else if (nk_grid[1] == 1 && nk_grid[2] == 1) {
            nkrep_out = nk_grid[0];
            tune_type = 0;

        } else if (nk_grid[2] == 1 && nk_grid[0] == 1) {
            nkrep_out = nk_grid[1];
            tune_type = 0;

        } else {
            tune_type = 1;
        }

        int ii, jj, kk;

        if (tune_type == 0) {

            double phase;

            memory->allocate(exp_phase, 2 * nkrep_out - 1);
#ifdef _OPENMP
#pragma omp parallel for private(phase)
#endif
            for (ii = 0; ii < 2 * nkrep_out - 1; ++ii) {
                phase = 2.0 * pi * static_cast<double>(ii - nkrep_out + 1)
                    / static_cast<double>(nkrep_out);
                exp_phase[ii] = std::exp(im * phase);
            }

        } else if (tune_type == 1) {

            double phase[3];

            memory->allocate(exp_phase3,
                             2 * nk_grid[0] - 1,
                             2 * nk_grid[1] - 1,
                             2 * nk_grid[2] - 1);
#ifdef _OPENMP
#pragma omp parallel for private(phase, jj, kk)
#endif
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
    }
}


void Relaxation::calc_self3omega_tetrahedron(const double Temp,
                                             double **eval,
                                             std::complex<double> ***evec,
                                             const unsigned int ik_in,
                                             const unsigned int snum,
                                             const unsigned int nomega,
                                             double *omega,
                                             double *ret)
{
    // This function returns the imaginary part of phonon self-energy 
    // for the given frequency range of omega, phonon frequency (eval) and phonon eigenvectors (evec).
    // The tetrahedron method will be used.
    // This version employs the crystal symmetry to reduce the computational cost
    // In addition, both MPI and OpenMP parallizations are used in a hybrid way inside this function.

    int ik, ib, iomega;
    int ns2 = ns * ns;

    unsigned int i;
    unsigned int is, js;
    unsigned int k1, k2;
    unsigned int arr[3];
    unsigned int npair_uniq;
    unsigned int nk_tmp;

    int ik_now;

    double n1, n2;
    double f1, f2;
    double omega_inner[2];

    int *kmap_identity, **kpairs;
    double **energy_tmp;
    double **weight_tetra;
    double **v3_arr, *v3_arr_loc;
    double *ret_private;

    std::vector<KsListGroup> triplet;
    std::vector<int> vk_l;

    int knum = kpoint->kpoint_irred_all[ik_in][0].knum;
    int knum_minus = kpoint->knum_minus[knum];
    double omega0 = eval[knum_minus][snum];

    get_unique_triplet_k(ik_in,
                         false,
                         false,
                         triplet);

    npair_uniq = triplet.size();

    if (npair_uniq != nk) {
        error->exit("hoge", "Something is wrong.");
    }

    memory->allocate(kpairs, nk, 2);
    memory->allocate(kmap_identity, nk);

    for (i = 0; i < nk; ++i) kmap_identity[i] = i;

    for (iomega = 0; iomega < nomega; ++iomega) ret[iomega] = 0.0;

    for (ik = 0; ik < npair_uniq; ++ik) {
        kpairs[ik][0] = triplet[ik].group[0].ks[0];
        kpairs[ik][1] = triplet[ik].group[0].ks[1];
    }

    if (nk % mympi->nprocs != 0) {
        nk_tmp = nk / mympi->nprocs + 1;
    } else {
        nk_tmp = nk / mympi->nprocs;
    }

    vk_l.clear();

    for (ik = 0; ik < nk; ++ik) {
        if (ik % mympi->nprocs == mympi->my_rank) {
            vk_l.push_back(ik);
        }
    }

    if (vk_l.size() < nk_tmp) {
        vk_l.push_back(-1);
    }

    memory->allocate(v3_arr_loc, ns2);
    memory->allocate(v3_arr, nk_tmp * mympi->nprocs, ns2);

    for (ik = 0; ik < nk_tmp; ++ik) {

        ik_now = vk_l[ik];

        if (ik_now == -1) {

            for (ib = 0; ib < ns2; ++ib) v3_arr_loc[ib] = 0.0; // do nothing

        } else {
#ifdef _OPENMP
#pragma omp parallel for private(is, js, arr)
#endif
            for (ib = 0; ib < ns2; ++ib) {

                is = ib / ns;
                js = ib % ns;

                arr[0] = ns * knum_minus + snum;
                arr[1] = ns * kpairs[ik_now][0] + is;
                arr[2] = ns * kpairs[ik_now][1] + js;

                v3_arr_loc[ib] = std::norm(V3(arr));
            }
        }
        MPI_Gather(&v3_arr_loc[0], ns2, MPI_DOUBLE,
                   v3_arr[ik * mympi->nprocs], ns2,
                   MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
    memory->deallocate(v3_arr_loc);

    if (mympi->my_rank == 0) {

#ifdef _OPENMP
#pragma omp parallel private(is, js, k1, k2, energy_tmp, i, \
                             iomega, weight_tetra, ik, \
                             omega_inner, f1, f2, n1, n2) 
#endif
        {
            memory->allocate(energy_tmp, 2, nk);
            memory->allocate(weight_tetra, 2, nk);
#ifdef _OPENMP
        const int nthreads = omp_get_num_threads();
        const int ithread = omp_get_thread_num();
#else
            const int nthreads = 1;
            const int ithread = 0;
#endif

#ifdef _OPENMP
#pragma omp single
#endif
            {
                memory->allocate(ret_private, nthreads * nomega);
                for (i = 0; i < nthreads * nomega; ++i) ret_private[i] = 0.0;
            }
#ifdef _OPENMP
#pragma omp for
#endif
            for (ib = 0; ib < ns2; ++ib) {

                is = ib / ns;
                js = ib % ns;

                for (ik = 0; ik < nk; ++ik) {
                    k1 = kpairs[ik][0];
                    k2 = kpairs[ik][1];

                    energy_tmp[0][ik] = eval[k1][is] + eval[k2][js];
                    energy_tmp[1][ik] = eval[k1][is] - eval[k2][js];
                }
                for (iomega = 0; iomega < nomega; ++iomega) {
                    for (i = 0; i < 2; ++i) {
                        integration->calc_weight_tetrahedron(nk,
                                                             kmap_identity,
                                                             weight_tetra[i],
                                                             energy_tmp[i],
                                                             omega[iomega]);
                    }

                    for (ik = 0; ik < nk; ++ik) {
                        k1 = kpairs[ik][0];
                        k2 = kpairs[ik][1];

                        omega_inner[0] = eval[k1][is];
                        omega_inner[1] = eval[k2][js];
                        if (thermodynamics->classical) {
                            f1 = thermodynamics->fC(omega_inner[0], Temp);
                            f2 = thermodynamics->fC(omega_inner[1], Temp);
                            n1 = f1 + f2;
                            n2 = f1 - f2;
                        } else {
                            f1 = thermodynamics->fB(omega_inner[0], Temp);
                            f2 = thermodynamics->fB(omega_inner[1], Temp);
                            n1 = f1 + f2 + 1.0;
                            n2 = f1 - f2;
                        }

                        //#pragma omp critical
                        ret_private[nomega * ithread + iomega]
                            += v3_arr[ik][ib] * (n1 * weight_tetra[0][ik] - 2.0 * n2 * weight_tetra[1][ik]);
                    }
                }
            }
#ifdef _OPENMP
#pragma omp for
#endif
            for (iomega = 0; iomega < nomega; ++iomega) {
                for (int t = 0; t < nthreads; t++) {
                    ret[iomega] += ret_private[nomega * t + iomega];
                }
            }
            memory->deallocate(energy_tmp);
            memory->deallocate(weight_tetra);
        }
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (iomega = 0; iomega < nomega; ++iomega) {
            ret[iomega] *= pi * std::pow(0.5, 4);
        }
        memory->deallocate(ret_private);
    }

    memory->deallocate(v3_arr);
    memory->deallocate(kmap_identity);
    memory->deallocate(kpairs);
}
