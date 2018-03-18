/*
anharmonic_core.cpp

Copyright (c) 2014, 2015, 2016 Terumasa Tadano

This file is distributed under the terms of the MIT license.
Please see the file 'LICENCE.txt' in the root directory
or http://opensource.org/licenses/mit-license.php for information.
*/

#include "mpi_common.h"
#include "anharmonic_core.h"
#include "constants.h"
#include "dynamical.h"
#include "error.h"
#include "fcs_phonon.h"
#include "integration.h"
#include "kpoint.h"
#include "mathfunctions.h"
#include "memory.h"
#include "mode_analysis.h"
#include "system.h"
#include "thermodynamics.h"
#include "write_phonons.h"
#include <boost/lexical_cast.hpp>
#include <iomanip>
#include <algorithm>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace PHON_NS;

AnharmonicCore::AnharmonicCore(PHON *phon) : Pointers(phon)
{
    set_default_variables();
}

AnharmonicCore::~AnharmonicCore()
{
    deallocate_variables();
};

void AnharmonicCore::set_default_variables()
{
    im = std::complex<double>(0.0, 1.0);
    quartic_mode = 0;
    use_tuned_ver = true;
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

void AnharmonicCore::deallocate_variables()
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


void AnharmonicCore::setup()
{
    if (mympi->my_rank == 0) {
        std::cout << std::endl;
        std::cout << " -----------------------------------------------------------------" << std::endl << std::endl;
        std::cout << " Now, move on to phonon lifetime calculations." << std::endl;
    }

    setup_cubic();
    if (quartic_mode > 0) setup_quartic();

    sym_permutation = true;
    use_tuned_ver = true;

    if (!mode_analysis->calc_fstate_k) {
        int nk_tmp[3];
        nk_tmp[0] = kpoint->nkx;
        nk_tmp[1] = kpoint->nky;
        nk_tmp[2] = kpoint->nkz;
        store_exponential_for_acceleration(nk_tmp,
                                           nk_represent,
                                           exp_phase,
                                           exp_phase3);
    }

    if (phon->mode == "RTA" && !mode_analysis->calc_fstate_k) {
        detect_imaginary_branches(dynamical->eval_phonon);
    }
}


void AnharmonicCore::detect_imaginary_branches(double **eval)
{
    int ik, is;
    auto nk = kpoint->nk;
    auto ns = dynamical->neval;
    auto nks = ns * nk;
    int knum;
    double omega;
    int ndup;

    bool is_anyof_imaginary = false;
    if (mympi->my_rank == 0) {

        memory->allocate(is_imaginary, kpoint->nk_irred, ns);

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

void AnharmonicCore::prepare_relative_vector(const std::vector<FcsArrayWithCell> &fcs_in,
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

void AnharmonicCore::prepare_group_of_force_constants(const std::vector<FcsArrayWithCell> &fcs_in,
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


std::complex<double> AnharmonicCore::V3(const unsigned int ks[3])
{
    unsigned int i, j;
    unsigned int kn[3], sn[3];
    unsigned int nsize_group;
    int ns = dynamical->neval;

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
    if (omega[0] < eps8 || omega[1] < eps8 || omega[2] < eps8) return 0.0;

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


std::complex<double> AnharmonicCore::V4(const unsigned int ks[4])
{
    int ii;
    unsigned int i, j;
    unsigned int kn[4], sn[4];
    int ns = dynamical->neval;

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
    // Return zero if any of the involving phonon has imaginary frequency
    if (omega[0] < eps8 || omega[1] < eps8 || omega[2] < eps8 || omega[3] < eps8) return 0.0;

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


std::complex<double> AnharmonicCore::V3_mode(int mode,
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
    if (eval[0][mode] < eps8 || eval[1][is] < eps8 || eval[2][js] < eps8) return 0.0;

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


void AnharmonicCore::calc_damping_smearing(const unsigned int N,
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

    int nk = kpoint->nk;
    int ns = dynamical->neval;
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

void AnharmonicCore::calc_damping_tetrahedron(const unsigned int N,
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

    int nk = kpoint->nk;
    int ns = dynamical->neval;

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


void AnharmonicCore::get_unique_triplet_k(const int ik,
                                          const bool use_triplet_symmetry,
                                          const bool use_permutation_symmetry,
                                          std::vector<KsListGroup> &triplet)
{
    int nk = kpoint->nk;

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


void AnharmonicCore::setup_cubic()
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

void AnharmonicCore::setup_quartic()
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

void AnharmonicCore::store_exponential_for_acceleration(const int nk_in[3],
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


void AnharmonicCore::calc_self3omega_tetrahedron(const double Temp,
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

    int nk = kpoint->nk;
    int ns = dynamical->neval;

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
