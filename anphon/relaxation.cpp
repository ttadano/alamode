/*
relaxation.cpp

Copyright (c) 2014 Terumasa Tadano

This file is distributed under the terms of the MIT license.
Please see the file 'LICENCE.txt' in the root directory 
or http://opensource.org/licenses/mit-license.php for information.
*/

#include "mpi_common.h"
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif 
#include "conductivity.h"
#include "dynamical.h"
#include "error.h"
#include "fcs_phonon.h"
#include "integration.h"
#include "kpoint.h"
#include "memory.h"
#include "parsephon.h"
#include "phonon_dos.h"
#include "phonon_thermodynamics.h"
#include "phonon_velocity.h"
#include "relaxation.h"
#include "selfenergy.h"
#include "symmetry_core.h"
#include "system.h"
#include "write_phonons.h"
#include "timer.h"
#include "constants.h"
#include "mathfunctions.h"
#include <set>
#include "integration.h"
#include <boost/lexical_cast.hpp>

using namespace PHON_NS;

Relaxation::Relaxation(PHON *phon): Pointers(phon) {
    im = std::complex<double>(0.0, 1.0);
}

Relaxation::~Relaxation(){};

void Relaxation::setup_relaxation()
{

    nk = kpoint->nk;
    ns = dynamical->neval;
    nks = ns*nk;

    unsigned int i, j, k;
    double *invsqrt_mass_p;


    if (mympi->my_rank == 0) {
        std::cout << " Setting up the relaxation time calculation ...";

        if (calc_realpart && integration->ismear == -1) {
            error->exit("setup_relaxation", "Sorry. REALPART = 1 can be used only with ISMEAR = 0");
        }
    }

    // Sort force_constant[1] using the operator defined in fcs_phonons.h
    // This sorting is necessarily.
    std::sort(fcs_phonon->force_constant_with_cell[1].begin(), fcs_phonon->force_constant_with_cell[1].end());
    prepare_group_of_force_constants(fcs_phonon->force_constant_with_cell[1], 3, ngroup, fcs_group);

    memory->allocate(v3_arr, nk, ns*ns);
    memory->allocate(delta_arr, nk, ns*ns, 4);
    memory->allocate(vec_for_v3, 3, 2, fcs_phonon->force_constant_with_cell[1].size());
    memory->allocate(invmass_for_v3, fcs_phonon->force_constant_with_cell[1].size());
    memory->allocate(evec_index, fcs_phonon->force_constant_with_cell[1].size(), 3);

    memory->allocate(invsqrt_mass_p, system->natmin);

    for (i = 0; i < system->natmin; ++i){
        invsqrt_mass_p[i] = std::sqrt(1.0 / system->mass[system->map_p2s[i][0]]);
    }
    j = 0;
    for (std::vector<FcsArrayWithCell>::const_iterator it  = fcs_phonon->force_constant_with_cell[1].begin();
        it != fcs_phonon->force_constant_with_cell[1].end(); ++it) {
            invmass_for_v3[j] 
            = invsqrt_mass_p[(*it).pairs[0].index / 3]
            * invsqrt_mass_p[(*it).pairs[1].index / 3]
            * invsqrt_mass_p[(*it).pairs[2].index / 3];

            ++j;
    }

    prepare_relative_vector(fcs_phonon->force_constant_with_cell[1], 3, vec_for_v3);

    for (i = 0; i < fcs_phonon->force_constant_with_cell[1].size(); ++i) {
        for (j = 0; j < 3; ++j) {
            evec_index[i][j] = fcs_phonon->force_constant_with_cell[1][i].pairs[j].index;
        }
    }

    if (quartic_mode > 0) {

        // This is for quartic vertexes.

        if (mympi->my_rank == 0) {
            std::cout << std::endl << std::endl;
            std::cout << " **********************************************************" << std::endl;
            std::cout << "     QUARTIC = 1: quartic_mode is on !                     " << std::endl;
            std::cout << "     Be careful! This mode is still under test.            " << std::endl;
            std::cout << "     There can be bugs and the computation is very heavy   " << std::endl;
            std::cout << " **********************************************************" << std::endl;
            std::cout << std::endl;
        }

        std::sort(fcs_phonon->force_constant_with_cell[2].begin(), fcs_phonon->force_constant_with_cell[2].end());
        prepare_group_of_force_constants(fcs_phonon->force_constant_with_cell[2], 4, ngroup2, fcs_group2);

        memory->allocate(vec_for_v4, 3, 3, fcs_phonon->force_constant_with_cell[2].size());
        memory->allocate(invmass_for_v4, fcs_phonon->force_constant_with_cell[2].size());
        memory->allocate(evec_index4, fcs_phonon->force_constant_with_cell[2].size(), 4);

        j = 0;
        for (std::vector<FcsArrayWithCell>::const_iterator it  = fcs_phonon->force_constant_with_cell[2].begin(); 
            it != fcs_phonon->force_constant_with_cell[2].end(); ++it) {
                invmass_for_v4[j] 
                = invsqrt_mass_p[(*it).pairs[0].index / 3] 
                * invsqrt_mass_p[(*it).pairs[1].index / 3] 
                * invsqrt_mass_p[(*it).pairs[2].index / 3] 
                * invsqrt_mass_p[(*it).pairs[3].index / 3];

                ++j;     
        }
        prepare_relative_vector(fcs_phonon->force_constant_with_cell[2], 4, vec_for_v4);

        for (i = 0; i < fcs_phonon->force_constant_with_cell[2].size(); ++i) {
            for (j = 0; j < 4; ++j) {
                evec_index4[i][j] = fcs_phonon->force_constant_with_cell[2][i].pairs[j].index;
            }
        }

        dynamical->modify_eigenvectors();
    }
    memory->deallocate(invsqrt_mass_p);

    MPI_Bcast(&calc_realpart, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD);
    MPI_Bcast(&atom_project_mode, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD);
    MPI_Bcast(&calc_fstate_k, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD);
    MPI_Bcast(&calc_fstate_omega, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD);

    if (mympi->my_rank == 0) {
        //     if (kpoint->kpoint_mode == 2) print_minimum_energy_diff();
        if (calc_fstate_omega) sym_permutation = false;
    }
    MPI_Bcast(&sym_permutation, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD);

    if (kpoint->kpoint_mode == 2) {

        if (ks_analyze_mode && use_triplet_symmetry) {
            use_triplet_symmetry = false;
            if (mympi->my_rank == 0) {
                std::cout << std::endl;
                std::cout << " TRISYM was automatically set to 0." << std::endl;
                std::cout << std::endl;
            }
        }

        generate_triplet_k(use_triplet_symmetry, sym_permutation);
    }

    if (mympi->my_rank == 0) {
        std::cout << " done!" << std::endl;
    }
}

void Relaxation::finish_relaxation()
{
    memory->deallocate(vec_for_v3);
    memory->deallocate(invmass_for_v3);
    memory->deallocate(evec_index);
    memory->deallocate(fcs_group);
    memory->deallocate(v3_arr);
    memory->deallocate(delta_arr);

    if (quartic_mode > 0) {
        memory->deallocate(vec_for_v4);
        memory->deallocate(invmass_for_v4);
        memory->deallocate(evec_index4);
        memory->deallocate(fcs_group2);
    }
}

void Relaxation::print_minimum_energy_diff()
{
    int i, j;
    unsigned int nk_near = 0;
    double domega_min;
    double dist_k_min, dist_k;
    double xk_tmp[3], xk_tmp2[3];
    int ik;

    domega_min = 0.0;

    if (nk > 1) {

        for (i = 0; i < 3; ++i) {
            xk_tmp[i] = 0.5;
        }
        rotvec(xk_tmp2, xk_tmp, system->rlavec_p, 'T');
        dist_k_min = std::sqrt(xk_tmp2[0]*xk_tmp2[0] + xk_tmp2[1]*xk_tmp2[1] + xk_tmp2[2]*xk_tmp2[2]);

        for (ik = 1; ik < nk; ++ik) {
            for (j = 0; j < 3; ++j) {
                xk_tmp[j] = kpoint->xk[ik][j];
            }
            rotvec(xk_tmp2, xk_tmp, system->rlavec_p, 'T');

            dist_k = std::sqrt(xk_tmp2[0]*xk_tmp2[0] + xk_tmp2[1]*xk_tmp2[1] + xk_tmp2[2]*xk_tmp2[2]);

            if (dist_k <= dist_k_min) {
                dist_k_min = dist_k;
                nk_near = ik;
            }
        }
        domega_min =  writes->in_kayser(dynamical->eval_phonon[nk_near][0]);	
    } else {
        std::cout << "There is only 1 reciprocal point." << std::endl;
    }

    std::cout << std::endl;
    std::cout << " Estimated minimum energy difference (cm^-1) = " << domega_min << std::endl;
    std::cout << std::endl;
}

void Relaxation::prepare_relative_vector(std::vector<FcsArrayWithCell> fcs_in, const unsigned int N, double ***vec_out)
{

    int i, j, k;
    int ix, iy, iz;

    double vec[3];
    double **xshift_s;

    unsigned int atm1_s, atm2_s;
    unsigned int atm1_p, atm2_p;
    unsigned int xyz1, xyz2;
    unsigned int icell;

    std::vector<unsigned int> atm_super, atm_prim;
    std::vector<unsigned int> xyz;
    std::vector<unsigned int> cells;

    double mat_convert[3][3];

    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j){
            mat_convert[i][j] = 0.0;
            for (k = 0; k < 3; ++k){
                mat_convert[i][j] += system->rlavec_p[i][k] * system->lavec_s[k][j]; 
            }
        }
    }

    memory->allocate(xshift_s, 27, 3);

    for (i = 0; i < 3; ++i) xshift_s[0][i] = 0.0;

    icell =0;

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

    for (std::vector<FcsArrayWithCell>::const_iterator it = fcs_in.begin(); it != fcs_in.end(); ++it) {

        atm_super.clear();
        atm_prim.clear();
        xyz.clear();
        cells.clear();

        for (i = 0; i < (*it).pairs.size(); ++i) {
            atm_p = (*it).pairs[i].index / 3;
            tran_tmp = (*it).pairs[i].tran;
            atm_s = system->map_p2s[atm_p][tran_tmp];

            atm_prim.push_back(atm_p);
            atm_super.push_back(atm_s);
            cells.push_back((*it).pairs[i].cell_s);
        }


        for (i = 0; i < N - 1; ++i) {

            for (j = 0; j < 3; ++j) {
                vec[j] = system->xr_s[atm_super[i + 1]][j] + xshift_s[cells[i + 1]][j] 
                - system->xr_s[system->map_p2s[atm_prim[i + 1]][0]][j];
            }

            rotvec(vec, vec, mat_convert);

            for (j = 0; j < 3; ++j) {
                vec_out[j][i][icount] = -vec[j];
            }
        }
        ++icount;
    }
    memory->deallocate(xshift_s);
}

void Relaxation::prepare_group_of_force_constants(std::vector<FcsArrayWithCell> fcs_in, const unsigned int N, 
                                                  int &number_of_groups, std::vector<double> *&fcs_group_out) 
{
    // Find the number of groups which has different evecs.

    unsigned int i;
    number_of_groups = 0;


    std::vector<int> arr_old, arr_tmp;

    arr_old.clear();
    for (i = 0; i < N; ++i) {
        arr_old.push_back(-1);
    }

    for (std::vector<FcsArrayWithCell>::const_iterator it = fcs_in.begin(); it != fcs_in.end(); ++it) {

        arr_tmp.clear();

        for (i = 0; i < (*it).pairs.size(); ++i) {
            arr_tmp.push_back((*it).pairs[i].index);
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

    for (std::vector<FcsArrayWithCell>::const_iterator it = fcs_in.begin(); it != fcs_in.end(); ++it) {

        arr_tmp.clear();

        for (i = 0; i < (*it).pairs.size(); ++i) {
            arr_tmp.push_back((*it).pairs[i].index);
        }

        if (arr_tmp != arr_old) {
            ++igroup;
            arr_old.clear();
            arr_old.reserve(arr_tmp.size());
            std::copy(arr_tmp.begin(), arr_tmp.end(), std::back_inserter(arr_old));
        }

        fcs_group_out[igroup].push_back((*it).fcs_val);  
    }
}

void Relaxation::setup_mode_analysis()
{
    // Judge if ks_analyze_mode should be turned on or not.

    unsigned int i;

    if (mympi->my_rank == 0) {
        if (!ks_input.empty()) {
            std::cout << std::endl;
            std::cout << " KS_INPUT is given." << std::endl;
            std::cout << " Analysis on specific k points will be performed " << std::endl;
            std::cout << " instead of thermal conductivity calculations." << std::endl;
            std::cout << std::endl;


            std::ifstream ifs_ks;
            ifs_ks.open(ks_input.c_str(), std::ios::in);
            if (!ifs_ks) error->exit("setup_mode_analysis", "Cannot open file KS_INPUT");

            unsigned int nlist;
            double ktmp[3];
            unsigned int snum_tmp;
            int knum_tmp;

            ifs_ks >> nlist;

            if (nlist <= 0) error->exit("setup_mode_analysis", 
                "First line in KS_INPUT files should be a positive integer.");

            if (calc_fstate_k) {
                kslist_fstate_k.clear();

                for (i = 0; i < nlist; ++i) {

                    ifs_ks >> ktmp[0] >> ktmp[1] >> ktmp[2] >> snum_tmp;

                    if (snum_tmp <= 0 || snum_tmp > dynamical->neval) {
                        error->exit("setup_mode_analysis", "Mode index out of range.");
                    }

                    kslist_fstate_k.push_back(KsListMode(ktmp, snum_tmp - 1));
                }
                std::cout << " The number of entries = " << kslist_fstate_k.size() << std::endl;

            } else {
                kslist.clear();
                for (i = 0; i < nlist; ++i) {

                    ifs_ks >> ktmp[0] >> ktmp[1] >> ktmp[2] >> snum_tmp;
                    knum_tmp = kpoint->get_knum(ktmp[0], ktmp[1], ktmp[2]);

                    if (knum_tmp == -1) error->exit("setup_mode_analysis", 
                        "Given kpoint is not exist in given k-point grid.");
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
                kslist_fstate_k.push_back(KsListMode(vec_tmp[i], mode_tmp[i]));
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

    unsigned int i, ielem;
    unsigned int kn[3], sn[3];

    int pos = -1;

    double phase;
    double omega[3];

    std::complex<double> ctmp;
    std::complex<double> ret = std::complex<double>(0.0, 0.0);
    std::complex<double> ret_in, vec_tmp;


    for (i = 0; i < 3; ++i){
        kn[i] = ks[i] / ns;
        sn[i] = ks[i] % ns;
        omega[i] = dynamical->eval_phonon[kn[i]][sn[i]];
    }

    ielem = 0;

    for (i = 0; i < ngroup; ++i) {

        vec_tmp = dynamical->evec_phonon[kn[0]][sn[0]][evec_index[ielem][0]] 
        * dynamical->evec_phonon[kn[1]][sn[1]][evec_index[ielem][1]]
        * dynamical->evec_phonon[kn[2]][sn[2]][evec_index[ielem][2]];

        ret_in = std::complex<double>(0.0, 0.0);

        for (int j = 0; j < fcs_group[i].size(); ++j) {

            phase = vec_for_v3[0][0][ielem] * kpoint->xk[kn[1]][0] 
            + vec_for_v3[1][0][ielem] * kpoint->xk[kn[1]][1]
            + vec_for_v3[2][0][ielem] * kpoint->xk[kn[1]][2]
            + vec_for_v3[0][1][ielem] * kpoint->xk[kn[2]][0] 
            + vec_for_v3[1][1][ielem] * kpoint->xk[kn[2]][1] 
            + vec_for_v3[2][1][ielem] * kpoint->xk[kn[2]][2];

            ctmp = fcs_group[i][j] * invmass_for_v3[ielem] * std::exp(im*phase);
            ret_in += ctmp;

            ++ielem;
        }
        ret += ret_in * vec_tmp;
    }

    return ret / std::sqrt(omega[0] * omega[1] * omega[2]);
}

std::complex<double> Relaxation::V4(const unsigned int ks[4]) 
{

    unsigned int i, j, ielem;
    unsigned int kn[4], sn[4];

    double phase;
    double omega[4];

    std::complex<double> ctmp, ret_in, vec_tmp;
    std::complex<double> ret = std::complex<double>(0.0, 0.0);

    for (i = 0; i < 4; ++i){
        kn[i] = ks[i] / ns;
        sn[i] = ks[i] % ns;
        omega[i] = dynamical->eval_phonon[kn[i]][sn[i]];
    }

    ielem = 0;

    for (i = 0; i < ngroup2; ++i) {

        vec_tmp = dynamical->evec_phonon[kn[0]][sn[0]][evec_index4[ielem][0]] 
        * dynamical->evec_phonon[kn[1]][sn[1]][evec_index4[ielem][1]]
        * dynamical->evec_phonon[kn[2]][sn[2]][evec_index4[ielem][2]]
        * dynamical->evec_phonon[kn[3]][sn[3]][evec_index4[ielem][3]];

        ret_in = std::complex<double>(0.0, 0.0);

        for (j = 0; j < fcs_group2[i].size(); ++j) {
            phase =
              vec_for_v4[0][0][ielem] * kpoint->xk[kn[1]][0] 
            + vec_for_v4[1][0][ielem] * kpoint->xk[kn[1]][1] 
            + vec_for_v4[2][0][ielem] * kpoint->xk[kn[1]][2]
            + vec_for_v4[0][1][ielem] * kpoint->xk[kn[2]][0] 
            + vec_for_v4[1][1][ielem] * kpoint->xk[kn[2]][1] 
            + vec_for_v4[2][1][ielem] * kpoint->xk[kn[2]][2]
            + vec_for_v4[0][2][ielem] * kpoint->xk[kn[3]][0] 
            + vec_for_v4[1][2][ielem] * kpoint->xk[kn[3]][1] 
            + vec_for_v4[2][2][ielem] * kpoint->xk[kn[3]][2];

            ctmp = fcs_group2[i][j] * invmass_for_v4[ielem] * std::exp(im * phase);
            ret_in += ctmp;

            ++ielem;
        }
        ret += ret_in * vec_tmp;
    }

    return ret / std::sqrt(omega[0] * omega[1] * omega[2] * omega[3]);
}

std::complex<double> Relaxation::V3_mode(int mode, double *xk2, double *xk3, 
                                         int is, int js, double **eval, std::complex<double> ***evec)
{
    int ielem;

    double phase;
    std::complex<double> ctmp = std::complex<double>(0.0, 0.0);

    for (ielem = 0; ielem < fcs_phonon->force_constant_with_cell[1].size(); ++ielem) {

        phase = vec_for_v3[0][0][ielem] * xk2[0] 
        + vec_for_v3[1][0][ielem] * xk2[1]
        + vec_for_v3[2][0][ielem] * xk2[2] 
        + vec_for_v3[0][1][ielem] * xk3[0]
        + vec_for_v3[1][1][ielem] * xk3[1]
        + vec_for_v3[2][1][ielem] * xk3[2];


        ctmp += fcs_phonon->force_constant_with_cell[1][ielem].fcs_val * invmass_for_v3[ielem] * std::exp(im * phase)
            * evec[0][mode][evec_index[ielem][0]] * evec[1][is][evec_index[ielem][1]] * evec[2][js][evec_index[ielem][2]];
    }

    return ctmp / std::sqrt(eval[0][mode] * eval[1][is] * eval[2][js]);
}

// 
// void Relaxation::calc_realpart_V4(const unsigned int N, double *T, const double omega, 
//                                   const unsigned int knum, const unsigned int snum, double *ret)
// {
//     unsigned int i, ik, is;
//     unsigned int arr[4];
//     double n1, omega1;
//     double v4_tmp, T_tmp;
// 
//     for (i = 0; i < N; ++i) ret[i] = 0.0;
// 
//     arr[0] = ns * kpoint->knum_minus[knum] + snum;
//     arr[1] = ns * knum + snum;
// 
//     for (ik = 0; ik < nk; ++ik) {
//         for (is = 0; is < ns; ++is) {
// 
//             arr[2] = ns * ik + is;
//             arr[3] = ns * kpoint->knum_minus[ik] + is;
// 
//             v4_tmp = V4(arr).real();
// 
//             omega1 = dynamical->eval_phonon[ik][is];
// 
//             for (i = 0; i < N; ++i) {
//                 T_tmp = T[i];
//                 n1 = phonon_thermodynamics->fB(omega1, T_tmp);
// 
//                 ret[i] += v4_tmp * (2.0 * n1 + 1.0);
//             }
//         }
//     }
// 
//     for (i = 0; i < N; ++i) ret[i] *= - 1.0 / (8.0 * static_cast<double>(nk));
// }


void Relaxation::calc_damping_smearing(const unsigned int N, double *T, const double omega, 
                                       const unsigned int ik_in, const unsigned int snum, double *ret)
{
    // This function returns the imaginary part of phonon self-energy 
    // for the given frequency omega.
    // Lorentzian or Gaussian smearing will be used.
    // This version employs the crystal symmetry to reduce the computational cost

    unsigned int i;
    int ik;
    unsigned int jk;
    unsigned int is, js; 
    unsigned int arr[3];

    int k1, k2;

    double T_tmp;
    double n1, n2;
    double v3_tmp;
    double xk_tmp[3];
    double omega_inner[2];

    int knum, knum_minus;
    double multi;

    for (i = 0; i < N; ++i) ret[i] = 0.0;

    int iloc, jloc, kloc;

    double **v3_arr;
    double ***delta_arr;
    double ret_tmp;

    double f1, f2;

    double epsilon = integration->epsilon;

    memory->allocate(v3_arr, pair_uniq[ik_in].size(), ns * ns);
    memory->allocate(delta_arr, pair_uniq[ik_in].size(), ns * ns, 2);

    knum = kpoint->kpoint_irred_all[ik_in][0].knum;
    knum_minus = kpoint->knum_minus[knum];

#pragma omp parallel for private(multi, arr, k1, k2, is, js, omega_inner)
    for (ik = 0; ik < pair_uniq[ik_in].size(); ++ik) {
        multi = static_cast<double>(pair_uniq[ik_in][ik].group.size());

        arr[0] = ns * knum_minus + snum;

        k1 = pair_uniq[ik_in][ik].group[0].ks[0];
        k2 = pair_uniq[ik_in][ik].group[0].ks[1];

        for (is = 0; is < ns; ++is) {
            arr[1] = ns * k1 + is;
            omega_inner[0] = dynamical->eval_phonon[k1][is];

            for (js = 0; js < ns; ++js) {
                arr[2] = ns * k2 + js;
                omega_inner[1] = dynamical->eval_phonon[k2][js];

                v3_arr[ik][ns * is + js] = std::norm(V3(arr)) * multi;

                if (integration->ismear == 0) {
                    delta_arr[ik][ns * is + js][0] = delta_lorentz(omega - omega_inner[0] - omega_inner[1], epsilon)
                        - delta_lorentz(omega + omega_inner[0] + omega_inner[1], epsilon);
                    delta_arr[ik][ns * is + js][1] = delta_lorentz(omega - omega_inner[0] + omega_inner[1], epsilon)
                        - delta_lorentz(omega + omega_inner[0] - omega_inner[1], epsilon);
                } else if (integration->ismear == 1) {
                    delta_arr[ik][ns * is + js][0] = delta_gauss(omega - omega_inner[0] - omega_inner[1], epsilon)
                        - delta_gauss(omega + omega_inner[0] + omega_inner[1], epsilon);
                    delta_arr[ik][ns * is + js][1] = delta_gauss(omega - omega_inner[0] + omega_inner[1], epsilon)
                        - delta_gauss(omega + omega_inner[0] - omega_inner[1], epsilon);
                }
            }
        }   
    }

    for (i = 0; i < N; ++i) {
        T_tmp = T[i];
        ret_tmp = 0.0;

#pragma omp parallel for private(k1, k2, is, js, omega_inner, n1, n2, f1, f2), reduction(+:ret_tmp)
        for (ik = 0; ik < pair_uniq[ik_in].size(); ++ik) {

            k1 = pair_uniq[ik_in][ik].group[0].ks[0];
            k2 = pair_uniq[ik_in][ik].group[0].ks[1];

            for (is = 0; is < ns; ++is){
                for (js = 0; js < ns; ++js) {

                    omega_inner[0] = dynamical->eval_phonon[k1][is];
                    omega_inner[1] = dynamical->eval_phonon[k2][js];

                    f1 = phonon_thermodynamics->fB(omega_inner[0], T_tmp);
                    f2 = phonon_thermodynamics->fB(omega_inner[1], T_tmp);
                    n1 =  f1 + f2 + 1.0;
                    n2 =  f1 - f2;

                    ret_tmp += v3_arr[ik][ns * is + js] 
                    * (n1 * delta_arr[ik][ns * is + js][0] - n2 * delta_arr[ik][ns * is + js][1]);
                }
            }
        }
        ret[i] = ret_tmp;
    }

    memory->deallocate(v3_arr);
    memory->deallocate(delta_arr);

    for (i = 0; i < N; ++i) ret[i] *=  pi * std::pow(0.5, 4) / static_cast<double>(nk);
}

void Relaxation::calc_damping_tetrahedron(const unsigned int N, double *T, const double omega, 
                                          const unsigned int ik_in, const unsigned int snum, double *ret)
{
    // This function returns the imaginary part of phonon self-energy 
    // for the given frequency omega.
    // Tetrahedron method will be used.
    // This version employs the crystal symmetry to reduce the computational cost

    int ik, ib;

    unsigned int i, j;
    unsigned int jk;
    unsigned int is, js; 
    unsigned int k1, k2;
    unsigned int arr[3];

    int knum, knum_minus;

    double T_tmp;
    double n1, n2;
    double f1, f2;

    double xk_tmp[3];
    double omega_inner[2];
    double multi;

    double ret_tmp;
    double epsilon = integration->epsilon;

    int *kmap_identity;
    double **energy_tmp;
    double **weight_tetra;
    double **v3_arr;
    double ***delta_arr;


    for (i = 0; i < N; ++i) ret[i] = 0.0;

    memory->allocate(v3_arr, pair_uniq[ik_in].size(), ns * ns);
    memory->allocate(delta_arr, pair_uniq[ik_in].size(), ns * ns, 2);

    knum = kpoint->kpoint_irred_all[ik_in][0].knum;
    knum_minus = kpoint->knum_minus[knum];

    memory->allocate(kmap_identity, nk);

    for (i = 0; i < nk; ++i) kmap_identity[i] = i;


#pragma omp parallel private(is, js, k1, k2, xk_tmp, energy_tmp, i, weight_tetra, ik, jk, multi, arr)
    {
        memory->allocate(energy_tmp, 3, nk);
        memory->allocate(weight_tetra, 3, nk);

#pragma omp for
        for (ib = 0; ib < ns * ns; ++ib) {
            is = ib / ns;
            js = ib % ns;

            for (k1 = 0; k1 < nk; ++k1) {

                // Prepare two-phonon frequency for tetrahedron method

                for (i = 0; i < 3; ++i) xk_tmp[i] = kpoint->xk[knum][i] - kpoint->xk[k1][i];

                k2 = kpoint->get_knum(xk_tmp[0], xk_tmp[1], xk_tmp[2]);

                energy_tmp[0][k1] = dynamical->eval_phonon[k1][is] + dynamical->eval_phonon[k2][js];
                energy_tmp[1][k1] = dynamical->eval_phonon[k1][is] - dynamical->eval_phonon[k2][js];
                energy_tmp[2][k1] = -energy_tmp[1][k1];
            }

            for (i = 0; i < 3; ++i) {
                integration->calc_weight_tetrahedron(nk, kmap_identity, weight_tetra[i], energy_tmp[i], omega);
            }

            for (ik = 0; ik < pair_uniq[ik_in].size(); ++ik) {

                multi = static_cast<double>(pair_uniq[ik_in][ik].group.size());

                k1 = pair_uniq[ik_in][ik].group[0].ks[0];
                k2 = pair_uniq[ik_in][ik].group[0].ks[1];

                arr[0] = ns * knum_minus + snum;
                arr[1] = ns * k1 + is;
                arr[2] = ns * k2 + js;

                v3_arr[ik][ib] = std::norm(V3(arr));

                delta_arr[ik][ib][0] = 0.0;
                delta_arr[ik][ib][1] = 0.0;

                for (i = 0; i < pair_uniq[ik_in][ik].group.size(); ++i) {
                    jk = pair_uniq[ik_in][ik].group[i].ks[0];
                    delta_arr[ik][ib][0] += weight_tetra[0][jk];
                    delta_arr[ik][ib][1] += weight_tetra[1][jk] - weight_tetra[2][jk];
                }
            }
        }

        memory->deallocate(energy_tmp);
        memory->deallocate(weight_tetra);
    }



    for (i = 0; i < N; ++i) {
        T_tmp = T[i];
        ret_tmp = 0.0;

#pragma omp parallel for private(k1, k2, is, js, omega_inner, n1, n2, f1, f2), reduction(+:ret_tmp)
        for (ik = 0; ik < pair_uniq[ik_in].size(); ++ik) {

            k1 = pair_uniq[ik_in][ik].group[0].ks[0];
            k2 = pair_uniq[ik_in][ik].group[0].ks[1];

            for (is = 0; is < ns; ++is){
                for (js = 0; js < ns; ++js) {

                    omega_inner[0] = dynamical->eval_phonon[k1][is];
                    omega_inner[1] = dynamical->eval_phonon[k2][js];

                    f1 = phonon_thermodynamics->fB(omega_inner[0], T_tmp);
                    f2 = phonon_thermodynamics->fB(omega_inner[1], T_tmp);
                    n1 =  f1 + f2 + 1.0;
                    n2 =  f1 - f2;

                    ret_tmp += v3_arr[ik][ns * is + js] 
                    * (n1 * delta_arr[ik][ns * is + js][0] - n2 * delta_arr[ik][ns * is + js][1]);
                }
            }
        }
        ret[i] = ret_tmp;
    }

    memory->deallocate(v3_arr);
    memory->deallocate(delta_arr);
    memory->deallocate(kmap_identity);

    for (i = 0; i < N; ++i) ret[i] *=  pi * std::pow(0.5, 4);

}

void Relaxation::calc_frequency_resolved_final_state(const unsigned int N, double *T, const double omega0, 
                                                     const unsigned int M, const double *omega, const unsigned int ik_in, const unsigned int snum, double **ret)
{
    int i, j, ik;

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
    double **ret_mpi;

    double epsilon = integration->epsilon;

    memory->allocate(ret_mpi, N, M);

    for (i = 0; i < N; ++i) {
        for (j = 0; j < M; ++j) {
            ret_mpi[i][j] = 0.0;
        }
    }

    for (ik = mympi->my_rank; ik < pair_uniq[ik_in].size(); ik += mympi->nprocs) {

        multi = static_cast<double>(pair_uniq[ik_in][ik].group.size());
        knum = kpoint->kpoint_irred_all[ik_in][0].knum;
        knum_minus = kpoint->knum_minus[knum];

        arr[0] = ns * knum_minus + snum;

        k1 = pair_uniq[ik_in][ik].group[0].ks[0];
        k2 = pair_uniq[ik_in][ik].group[0].ks[1];

        for (is = 0; is < ns; ++is) {
            for (js = 0; js < ns; ++js) {
                arr[1] = ns * k1 + is;
                arr[2] = ns * k2 + js;

                omega_inner[0] = dynamical->eval_phonon[k1][is];
                omega_inner[1] = dynamical->eval_phonon[k2][js];

                v3_tmp = std::norm(V3(arr));

                for (i = 0; i < N; ++i) {
                    T_tmp = T[i];

                    f1 = phonon_thermodynamics->fB(omega_inner[0], T_tmp);
                    f2 = phonon_thermodynamics->fB(omega_inner[1], T_tmp);
                    n1 = f1 + f2 + 1.0;
                    n2 = f1 - f2;

                    if (integration->ismear == 0) {
                        prod_tmp[0] = n1 * (delta_lorentz(omega0 - omega_inner[0] - omega_inner[1], epsilon) 
                            - delta_lorentz(omega0 + omega_inner[0] + omega_inner[1], epsilon));
                        prod_tmp[1] = n2 * (delta_lorentz(omega0 + omega_inner[0] - omega_inner[1], epsilon)
                            - delta_lorentz(omega0 - omega_inner[0] + omega_inner[1], epsilon));

                        for (j = 0; j < M; ++j) {
                            ret_mpi[i][j] += v3_tmp * multi * delta_lorentz(omega[j] - omega_inner[0], epsilon)
                                * (prod_tmp[0] + prod_tmp[1]);
                        }
                    } else if (integration->ismear == 1) {
                        prod_tmp[0] = n1 * (delta_gauss(omega0 - omega_inner[0] - omega_inner[1], epsilon) 
                            - delta_gauss(omega0 + omega_inner[0] + omega_inner[1], epsilon));
                        prod_tmp[1] = n2 * (delta_gauss(omega0 + omega_inner[0] - omega_inner[1], epsilon)
                            - delta_gauss(omega0 - omega_inner[0] + omega_inner[1], epsilon));

                        for (j = 0; j < M; ++j) {
                            ret_mpi[i][j] += v3_tmp * multi * delta_gauss(omega[j] - omega_inner[0], epsilon)
                                * (prod_tmp[0] + prod_tmp[1]);
                        }
                    }
                    
                }
            }
        }
    }
    for (i = 0; i < N; ++i) {
        for (j = 0; j < M; ++j) {
            ret_mpi[i][j] *=  pi * std::pow(0.5, 4) / static_cast<double>(nk);
        }
    }

    MPI_Reduce(&ret_mpi[0][0], &ret[0][0], N*M, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    memory->deallocate(ret_mpi);
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


    NT = static_cast<unsigned int>((Tmax - Tmin) / dT);
    memory->allocate(T_arr, NT);
    for (i = 0; i < NT; ++i) T_arr[i] = Tmin + static_cast<double>(i)*dT;

    double epsilon = integration->epsilon;

    if (calc_fstate_k) {

        // Momentum-resolved final state amplitude
        print_momentum_resolved_final_state(NT, T_arr, epsilon);

    } else if (calc_fstate_omega) {

        print_frequency_resolved_final_state(NT, T_arr);

    } else {

        double *damping_a;
        double omega_shift;
        std::complex<double> *self_a, *self_b, *self_c, *self_d, *self_e;
        std::complex<double> *self_f, *self_g, *self_h, *self_i, *self_j;

        if (mympi->my_rank == 0) {
            std::cout << std::endl;
            std::cout << " Calculate the line width (FWHM) of phonons" << std::endl;
            std::cout << " due to 3-phonon interactions for given " << kslist.size() << " modes." << std::endl;

            if (calc_realpart) {
                if (quartic_mode == 1) {
                    std::cout << " REALPART = 1 and " << std::endl;
                    std::cout << " QUARTIC = 1      : Additionally, frequency shift of phonons due to 3-phonon" << std::endl;
                    std::cout << "                    and 4-phonon interactions will be calculated." << std::endl;
                } else {
                    std::cout << " REALPART = 1 : Additionally, frequency shift of phonons due to 3-phonon" << std::endl;
                    std::cout << "                interactions will be calculated." << std::endl;
                }
            }

            if (quartic_mode == 2) {
                std::cout << std::endl;
                std::cout << " QUARTIC = 2 : Additionally, phonon line width due to 4-phonon" << std::endl;
                std::cout << "               interactions will be calculated." << std::endl;
                std::cout << " WARNING: This is very very expensive." << std::endl;
            }
        }

        memory->allocate(damping_a, NT);
        memory->allocate(self_a, NT);
        memory->allocate(self_b, NT);
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
                std::cout << "  Frequency (cm^-1) : " << std::setw(15) << writes->in_kayser(omega) << std::endl;
            }

            int ik_irred = kpoint->kmap_to_irreducible[knum];

            if (integration->ismear == -1) {
                calc_damping_tetrahedron(NT, T_arr, omega, ik_irred, snum, damping_a);
            } else {
                calc_damping_smearing(NT, T_arr, omega, ik_irred, snum, damping_a);
            }
            if (quartic_mode == 2) {
                selfenergy->selfenergy_c(NT, T_arr, omega, knum, snum, self_c);
                selfenergy->selfenergy_d(NT, T_arr, omega, knum, snum, self_d);
                selfenergy->selfenergy_e(NT, T_arr, omega, knum, snum, self_e);
                selfenergy->selfenergy_f(NT, T_arr, omega, knum, snum, self_f);
                selfenergy->selfenergy_g(NT, T_arr, omega, knum, snum, self_g);
                selfenergy->selfenergy_h(NT, T_arr, omega, knum, snum, self_h);
                selfenergy->selfenergy_i(NT, T_arr, omega, knum, snum, self_i);
                selfenergy->selfenergy_j(NT, T_arr, omega, knum, snum, self_j);
            }

            if (mympi->my_rank == 0) {
                file_linewidth = input->job_title + ".Gamma." + boost::lexical_cast<std::string>(i + 1);
                ofs_linewidth.open(file_linewidth.c_str(), std::ios::out);
                if (!ofs_linewidth) error->exit("perform_mode_analysis", "Cannot open file file_linewidth");

                ofs_linewidth << "# xk = ";

                for (j = 0; j < 3; ++j) {
                    ofs_linewidth << std::setw(15) << kpoint->xk[knum][j];
                }
                ofs_linewidth << std::endl;
                ofs_linewidth << "# mode = " << snum + 1<< std::endl;
                ofs_linewidth << "# Frequency = " << writes->in_kayser(omega) << std::endl;
                ofs_linewidth << "## Temperature dependence of Gamma for given mode" << std::endl;
                ofs_linewidth << "## T[K], Gamma3 (cm^-1)";
                if (quartic_mode == 2) ofs_linewidth << ", Gamma4(cm^-1) <-- specific diagram only";
                ofs_linewidth << std::endl;

                for (j = 0; j < NT; ++j) {
                    ofs_linewidth << std::setw(10) << T_arr[j] << std::setw(15) << writes->in_kayser(2.0 * damping_a[j]);

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

                selfenergy->selfenergy_a(NT, T_arr, omega, knum, snum, self_a);
                
                if (quartic_mode == 1) {
                    selfenergy->selfenergy_b(NT, T_arr, omega, knum, snum, self_b);
                }

                if (mympi->my_rank == 0) {

                    file_shift = input->job_title + ".Shift." + boost::lexical_cast<std::string>(i + 1);
                    ofs_shift.open(file_shift.c_str(), std::ios::out);
                    if (!ofs_shift) error->exit("perform_mode_analysis", "Cannot open file file_shift");

                    ofs_shift << "# xk = ";

                    for (j = 0; j < 3; ++j) {
                        ofs_shift << std::setw(15) << kpoint->xk[knum][j];
                    }
                    ofs_shift << std::endl;
                    ofs_shift << "# mode = " << snum + 1<< std::endl;
                    ofs_shift << "# Frequency = " << writes->in_kayser(omega) << std::endl;
                    ofs_shift << "## T[K], Shift3 (cm^-1)";
                    if (quartic_mode == 1) ofs_shift << ", Shift4 (cm^-1) <-- linear term in lambda";
                    ofs_shift << ", Shifted frequency (cm^-1)";
                    ofs_shift << std::endl;

                    for (j = 0; j < NT; ++j) {
                        ofs_shift << std::setw(10) << T_arr[j] << std::setw(15) << writes->in_kayser(-self_a[j].real());

                        omega_shift = omega - self_a[j].real();

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
void Relaxation::print_frequency_resolved_final_state(const unsigned int NT, double *T_arr)
{
    int i, j;
    unsigned int knum, snum;
    double omega, omega0;
    double **gamma_final;
    double *freq_array;
    int ienergy;
    std::ofstream ofs_omega;
    std::string file_omega;

    memory->allocate(gamma_final, NT, dos->n_energy);
    memory->allocate(freq_array, dos->n_energy);

    for (i = 0; i < dos->n_energy; ++i) {
        freq_array[i] = dos->energy_dos[i] * time_ry / Hz_to_kayser;
    }

    if (mympi->my_rank == 0) {

        std::cout << std::endl;
        std::cout << " FSTATE_W = 1 : Calculate the frequency-resolved final state amplitude" << std::endl;
        std::cout << "                due to 3-phonon interactions." << std::endl;

        if (integration->ismear == -1) {
            error->exit("print_frequency_resolved_final_state",
                "Sorry, ISMEAR=-1 cannot be used with FSTATE_W = 1");
        }
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
            std::cout << "  Frequency (cm^-1) : " << std::setw(15) << writes->in_kayser(omega0) << std::endl;
        }

        calc_frequency_resolved_final_state(NT, T_arr, omega0, dos->n_energy, 
            freq_array, kpoint->kmap_to_irreducible[knum], snum, gamma_final);

        if (mympi->my_rank == 0) {

            file_omega = input->job_title + ".fw." + boost::lexical_cast<std::string>(i + 1);
            ofs_omega.open(file_omega.c_str(), std::ios::out);
            if (!ofs_omega) error->exit("print_frequency_resolved_final_state", "Cannot open file file_omega");

            ofs_omega << "# xk = ";

            for (j = 0; j < 3; ++j) {
                ofs_omega << std::setw(15) << kpoint->xk[knum][j];
            }
            ofs_omega << std::endl;
            ofs_omega << "# mode = " << snum << std::endl;
            ofs_omega << "# Frequency = " << writes->in_kayser(omega0) << std::endl;

            ofs_omega<< "## Frequency-resolved final state amplitude for given modes" << std::endl;
            ofs_omega << "## Gamma[omega][temperature] in cm^-1";
            ofs_omega << std::endl;

            ofs_omega << "## ";
            for (j = 0; j < NT; ++j) {
                ofs_omega << std::setw(10) << T_arr[j];
            }
            ofs_omega << std::endl;
            for (ienergy = 0; ienergy < dos->n_energy; ++ienergy) {
                omega = dos->energy_dos[ienergy];

                ofs_omega << std::setw(10) << omega;
                for (j = 0; j < NT; ++j) ofs_omega << std::setw(15) << writes->in_kayser(gamma_final[j][ienergy]);
                ofs_omega << std::endl;
            }
            ofs_omega.close();
            std::cout << "  Frequency-resolved final state amplitude is printed in " << file_omega << std::endl;
        }
    }

    memory->deallocate(freq_array);
    memory->deallocate(gamma_final);
}

void Relaxation::print_momentum_resolved_final_state(const unsigned int NT, double *T_arr, double epsilon)
{

    int i, j, k, l;
    int iT;
    int is, js;
    int nklist;
    int mode;
    double xk1[3], xk2[3], xk3[3];
    double kvec[3];
    double f1, f2, n1, n2;
    double norm, T_tmp;
    double V3norm;
    double **eval;
    double **gamma_k, **gamma_k_mpi;
    double delta_tmp[2];

    std::complex<double> ***evec;

    std::ofstream ofs_mode_tau;
    std::string file_mode_tau;

    if (mympi->my_rank == 0) {
        std::cout << std::endl;
        if (integration->ismear == -1) {
            std::cout << " ISMEAR = -1: Tetrahedron method will be used." << std::endl;
            std::cout << " Sorry. ISMEAR = -1 is not used with FSTATE_K = 1.";
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
        std::cout << "                due to 3-phonon interactions for given " << kslist_fstate_k.size() 
            << " entries." << std::endl;
        std::cout << std::endl;
    }

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

        if (mympi->my_rank == 0) {
            std::cout << " Number : " << std::setw(5) << i + 1 << std::endl;
            std::cout << "  Phonon at k = (";
            for (j = 0; j < 3; ++j) {
                std::cout << std::setw(10) << std::fixed << kslist_fstate_k[i].xk[j];
                if (j < 2) std::cout << ",";
            }
            std::cout << ")" << std::endl;
            std::cout << "  Mode index = " << std::setw(5) << mode + 1 << std::endl;
            std::cout << "  Frequency (cm^-1) : " << std::setw(15) << writes->in_kayser(eval[0][mode]) << std::endl;
        }

        for (j = 0; j < kpoint->nplanes; ++j) {

            nklist = kpoint->kp_planes[j].size();

            memory->allocate(gamma_k, nklist, NT);
            memory->allocate(gamma_k_mpi, nklist, NT);

            for (k = 0; k < nklist; ++k) {
                for (l = 0; l < NT; ++l) {
                    gamma_k[k][l] = 0.0;
                    gamma_k_mpi[k][l] = 0.0;
                }
            }

            for (k = mympi->my_rank; k < nklist; k += mympi->nprocs) {

                for (l = 0; l < 3; ++l) xk2[l] = dynamical->fold(kpoint->kp_planes[j][k].k[l]);
                for (l = 0; l < 3; ++l) xk3[l] = dynamical->fold(-xk1[l]-xk2[l]);

                for (l = 0; l < 3; ++l) kvec[l] = xk2[l];
                rotvec(kvec, kvec, system->rlavec_p, 'T');
                norm = std::sqrt(kvec[0] * kvec[0] + kvec[1] * kvec[1] + kvec[2] * kvec[2]);

                if (norm > eps) for (l = 0; l < 3; ++l) kvec[l] /= norm;

                dynamical->eval_k(xk2, kvec, fcs_phonon->fc2_ext, eval[1], evec[1], true);

                for (l = 0; l < 3; ++l) kvec[l] = xk3[l];
                rotvec(kvec, kvec, system->rlavec_p, 'T');
                norm = std::sqrt(kvec[0] * kvec[0] + kvec[1] * kvec[1] + kvec[2] * kvec[2]);

                if (norm > eps) for (l = 0; l < 3; ++l) kvec[l] /= norm;

                dynamical->eval_k(xk3, kvec, fcs_phonon->fc2_ext, eval[2], evec[2], true);

                for (l = 0; l < ns; ++l) {
                    eval[1][l] = dynamical->freq(eval[1][l]);
                    eval[2][l] = dynamical->freq(eval[2][l]);
                }

                for (is = 0; is < ns; ++is) {
                    for (js = 0; js < ns; ++js) {
                        V3norm = std::norm(V3_mode(mode, xk2, xk3, is, js, eval, evec));

                        if (integration->ismear == 0) {
                            delta_tmp[0] = delta_lorentz(eval[0][mode] - eval[1][is] - eval[2][js], epsilon)
                                - delta_lorentz(eval[0][mode] + eval[1][is] + eval[2][js], epsilon);
                            delta_tmp[1] = delta_lorentz(eval[0][mode] + eval[1][is] - eval[2][js], epsilon)
                                - delta_lorentz(eval[0][mode] - eval[1][is] + eval[2][js], epsilon);
                        } else {
                            delta_tmp[0] = delta_gauss(eval[0][mode] - eval[1][is] - eval[2][js], epsilon)
                                - delta_gauss(eval[0][mode] + eval[1][is] + eval[2][js], epsilon);
                            delta_tmp[1] = delta_gauss(eval[0][mode] + eval[1][is] - eval[2][js], epsilon)
                                - delta_gauss(eval[0][mode] - eval[1][is] + eval[2][js], epsilon);
                        }


                        for (iT = 0; iT < NT; ++iT) {
                            T_tmp = T_arr[iT];

                            f1 = phonon_thermodynamics->fB(eval[1][is], T_tmp);
                            f2 = phonon_thermodynamics->fB(eval[2][js], T_tmp);
                            n1 = f1 + f2 + 1.0;
                            n2 = f1 - f2;

                            gamma_k_mpi[k][iT] += V3norm * (n1 * delta_tmp[0] + n2 * delta_tmp[1]);
                        }
                    }
                }

                for (iT = 0; iT < NT; ++iT) gamma_k_mpi[k][iT] *= pi * std::pow(0.5, 4);
            }

            MPI_Reduce(&gamma_k_mpi[0][0], &gamma_k[0][0], NT*nklist, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

            if (mympi->my_rank == 0) {

                file_mode_tau = input->job_title + ".fk." + boost::lexical_cast<std::string>(i + 1);
                ofs_mode_tau.open(file_mode_tau.c_str(), std::ios::out);
                if (!ofs_mode_tau) error->exit("compute_mode_tau", "Cannot open file file_mode_tau");

                ofs_mode_tau << "## Momentum-resolved final state amplitude" << std::endl;

                ofs_mode_tau << "# " << "Gamma at ";
                for (l = 0; l < 3; ++l) ofs_mode_tau << std::setw(10) << -xk1[l];
                ofs_mode_tau << " , mode = " << mode << std::endl;

                for (iT = 0; iT < NT; ++iT) {
                    ofs_mode_tau << "# T = " << std::setw(10) << T_arr[iT] << std::endl;
                    ofs_mode_tau << "# Plane = " << std::setw(5) << j + 1 << std::endl;
                    for (k = 0; k < nklist; ++k) {
                        ofs_mode_tau << std::setw(5) << kpoint->kp_planes[j][k].n[0];
                        ofs_mode_tau << std::setw(5) << kpoint->kp_planes[j][k].n[1];
                        ofs_mode_tau << std::setw(15) << gamma_k[k][iT] << std::endl;
                    }
                }

                ofs_mode_tau.close();
                std::cout << "  The result is saved in " << file_mode_tau << std::endl;
                std::cout << std::endl;
            }

            memory->deallocate(gamma_k);
            memory->deallocate(gamma_k_mpi);
        }
    }
    memory->deallocate(eval);
    memory->deallocate(evec);
}

int Relaxation::knum_sym(const int nk_in, const int symop_num) {

    int i, j;
    double srot[3][3];
    double srot_inv[3][3], srot_inv_t[3][3];
    double xk_orig[3], xk_sym[3];

    if (symop_num < 0 || symop_num >= symmetry->nsym) {
        error->exit("knum_sym", "Invalid symop_num");
    }

    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            srot[i][j] = static_cast<double>(symmetry->SymmList[symop_num].symop[3 * i + j]);
        }
    }

    invmat3(srot_inv, srot);
    transpose3(srot_inv_t, srot_inv);

    for (i = 0; i < 3; ++i) xk_orig[i] = kpoint->xk[nk_in][i];

    rotvec(xk_sym, xk_orig, srot_inv_t);
    for (i = 0; i < 3; ++i){
        xk_sym[i] = xk_sym[i] - nint(xk_sym[i]);
    }

    int ret = kpoint->get_knum(xk_sym[0], xk_sym[1], xk_sym[2]);

    return ret;
}

bool Relaxation::is_proper(const int isym)
{
    int i, j;
    double det;
    double S[3][3];
    bool ret;

    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            S[i][j] = static_cast<double>(symmetry->SymmList[isym].symop[3 * i + j]);
        }
    }

    det = S[0][0] * (S[1][1] * S[2][2] - S[2][1] * S[1][2])
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
    int i;
    int tran[3];
    bool ret;

    for (i = 0; i < 3; ++i) tran[i] = symmetry->SymmList[isym].symop[9 + i];

    if (tran[0] == 0 && tran[1] == 0 && tran[2] == 0) {
        ret = true;
    } else {
        ret = false;
    }
    return ret;
}

void Relaxation::generate_triplet_k(const bool use_triplet_symmetry, const bool use_permutation_symmetry)
{
    int i, j;
    int *num_group_k;
    int **symmetry_group_k;

    int knum, isym, ksym;
    int ik1, ik2;
    int ks_in[2], tmp;
    double xk[3], xk1[3], xk2[3];

    bool *flag_found;

    std::vector<KsList> kslist;

    memory->allocate(num_group_k, kpoint->nk_reduced);
    memory->allocate(symmetry_group_k, kpoint->nk_reduced, symmetry->nsym);
    memory->allocate(pair_uniq, kpoint->nk_reduced);
    memory->allocate(flag_found, kpoint->nk);

    for (i = 0; i < kpoint->nk_reduced; ++i) {

        knum = kpoint->kpoint_irred_all[i][0].knum;

        if (use_triplet_symmetry) {

            num_group_k[i] = 0;
            j = 0;

            for (isym = 0; isym < symmetry->nsym; ++isym) {

                ksym = knum_sym(knum, isym);
                if (ksym == knum) {
                    num_group_k[i] += 1;
                    symmetry_group_k[i][j++] = isym;
                }
            }
        } else {
            num_group_k[i] = 1;
            symmetry_group_k[i][0] = 0; // Identity matrix
        }

        for (j = 0; j < 3; ++j) xk[j] = kpoint->xk[knum][j];

        for (j = 0; j < kpoint->nk; ++j) flag_found[j] = false;

        pair_uniq[i].clear();

        for (ik1 = 0; ik1 < nk; ++ik1) {

            for (j = 0; j < 3; ++j) xk1[j] = kpoint->xk[ik1][j];
            for (j = 0; j < 3; ++j) xk2[j] = xk[j] - xk1[j];

            ik2 = kpoint->get_knum(xk2[0], xk2[1], xk2[2]);

            kslist.clear();

            if (ik1 > ik2 && use_permutation_symmetry) continue;

            for (isym = 0; isym < num_group_k[i]; ++isym) {

                ks_in[0] = knum_sym(ik1, symmetry_group_k[i][isym]);
                ks_in[1] = knum_sym(ik2, symmetry_group_k[i][isym]);

                if (!flag_found[ks_in[0]]) {
                    kslist.push_back(KsList(2, ks_in, symmetry_group_k[i][isym]));
                    flag_found[ks_in[0]] = true;
                }

                if (ks_in[0] != ks_in[1] && use_permutation_symmetry) {
                    tmp = ks_in[0];
                    ks_in[0] = ks_in[1];
                    ks_in[1] = tmp;

                    if (!flag_found[ks_in[0]]) {
                        kslist.push_back(KsList(2, ks_in, symmetry_group_k[i][isym]));
                        flag_found[ks_in[0]] = true;
                    }
                }
            }

            if (kslist.size() > 0) {
                pair_uniq[i].push_back(kslist);
            }
        }
    }

    memory->deallocate(num_group_k);
    memory->deallocate(symmetry_group_k);
    memory->deallocate(flag_found);
}
