#include "mpi_common.h"
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <vector>
#include <omp.h>
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

using namespace PHON_NS;

Relaxation::Relaxation(PHON *phon): Pointers(phon) {
    im = std::complex<double>(0.0, 1.0);
}

Relaxation::~Relaxation(){};

void Relaxation::setup_relaxation()
{
    if (mympi->my_rank == 0) {
        std::cout << "Setting up the relaxation time calculation ...";

        if (calc_realpart && ksum_mode == -1) {
            error->exit("setup_relaxation", "Sorry. REALPART = 1 can be used only with ISMEAR = 0");
        }
    }

    nk = kpoint->nk;
    ns = dynamical->neval;
    nks = ns*nk;

    unsigned int i, j, k;
    double mat_convert[3][3];
    double ***relvec;
    double *invsqrt_mass_p;

    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j){
            mat_convert[i][j] = 0.0;
            for (k = 0; k < 3; ++k){
                mat_convert[i][j] += system->rlavec_p[i][k] * system->lavec_s[k][j]; 
            }
        }
    }

    memory->allocate(relvec, system->nat, system->nat, 3);
    memory->allocate(invsqrt_mass_p, system->natmin);

    if (mympi->my_rank == 0) {

        double vec[3];

        for (i = 0; i < system->nat; ++i){
            for (j = 0; j < system->nat; ++j){

                for (k = 0; k < 3; ++k){

                    vec[k] = system->xr_s[i][k] - system->xr_s[j][k];
                    vec[k] = dynamical->fold(vec[k]);

                    vec[k] += system->xr_s[system->map_p2s[system->map_s2p[j].atom_num][0]][k]
                    -system->xr_s[system->map_p2s[system->map_s2p[i].atom_num][0]][k];


                    relvec[i][j][k] = -vec[k];
                }
                rotvec(relvec[i][j], relvec[i][j], mat_convert);
            }
        }
    }
    MPI_Bcast(&relvec[0][0][0], 3*system->nat*system->nat, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    for (i = 0; i < system->natmin; ++i){
        invsqrt_mass_p[i] = std::sqrt(1.0 / system->mass[system->map_p2s[i][0]]);
    }

    // Sort force_constant[1] using the operator defined in fcs_phonons.h
    std::sort(fcs_phonon->force_constant[1].begin(), fcs_phonon->force_constant[1].end());


    // Find the number of groups which has different evecs.
    ngroup = 0;

    int arr_old[3] = {-1, -1, -1};
    int arr_tmp[3];
    for (std::vector<FcsClass>::const_iterator it = fcs_phonon->force_constant[1].begin(); it != fcs_phonon->force_constant[1].end(); ++it) {
        for (i = 0; i < 3; ++i) arr_tmp[i] = 3 * (*it).elems[i].atom + (*it).elems[i].xyz;

        if (arr_tmp[0] != arr_old[0] || arr_tmp[1] != arr_old[1] || arr_tmp[2] != arr_old[2]) {
            ++ngroup;
            for (i = 0; i < 3; ++i) arr_old[i] = arr_tmp[i];
        }
    }

    memory->allocate(fcs_group, ngroup);

    for (i = 0; i < 3; ++i) arr_old[i] = -1;
    int igroup = -1;
    for (std::vector<FcsClass>::const_iterator it = fcs_phonon->force_constant[1].begin(); it != fcs_phonon->force_constant[1].end(); ++it) {
        for (i = 0; i < 3; ++i) arr_tmp[i] = 3 * (*it).elems[i].atom + (*it).elems[i].xyz;

        if (arr_tmp[0] != arr_old[0] || arr_tmp[1] != arr_old[1] || arr_tmp[2] != arr_old[2]) {
            ++igroup;
            for (i = 0; i < 3; ++i) arr_old[i] = arr_tmp[i];
        }
        fcs_group[igroup].push_back(*it);  
    }

    memory->allocate(v3_arr, nk, ns*ns);
    memory->allocate(delta_arr, nk, ns*ns, 4);

    memory->allocate(vec_for_v3, 3, 2, fcs_phonon->force_constant[1].size());
    memory->allocate(invmass_for_v3, fcs_phonon->force_constant[1].size());

    j = 0;
    unsigned int atom_num[3];

    for (std::vector<FcsClass>::const_iterator it = fcs_phonon->force_constant[1].begin(); it != fcs_phonon->force_constant[1].end(); ++it) {

        for (i = 0; i < 3; ++i) atom_num[i] = system->map_p2s[(*it).elems[i].atom][(*it).elems[i].cell];

        for (i = 0; i < 3; ++i) {
            vec_for_v3[i][0][j] = relvec[atom_num[1]][atom_num[0]][i];
            vec_for_v3[i][1][j] = relvec[atom_num[2]][atom_num[0]][i];
        }

        invmass_for_v3[j] 
        = invsqrt_mass_p[(*it).elems[0].atom] 
        * invsqrt_mass_p[(*it).elems[1].atom] 
        * invsqrt_mass_p[(*it).elems[2].atom];

        ++j;     
    }

    memory->allocate(evec_index, fcs_phonon->force_constant[1].size(), 3);

    for (i = 0; i < fcs_phonon->force_constant[1].size(); ++i) {
        for (j = 0; j < 3; ++j) {
            evec_index[i][j] 
            = 3 * fcs_phonon->force_constant[1][i].elems[j].atom 
                + fcs_phonon->force_constant[1][i].elems[j].xyz;
        }
    }

    if (quartic_mode) {

        // This is for quartic vertexes.

        if (mympi->my_rank == 0) {
            std::cout << std::endl << std::endl;
            std::cout << "**********************************************************" << std::endl;
            std::cout << "    QUARTIC = 1: quartic_mode is on !                     " << std::endl;
            std::cout << "    Be careful! This mode is still under test.            " << std::endl;
            std::cout << "    There can be bugs and the computation is very heavy   " << std::endl;
            std::cout << "**********************************************************" << std::endl;
            std::cout << std::endl;
        }

        memory->allocate(vec_for_v4, fcs_phonon->force_constant[2].size(), 3, 3);
        memory->allocate(invmass_for_v4, fcs_phonon->force_constant[2].size());

        j = 0;
        unsigned int atom_num4[4];

        for (std::vector<FcsClass>::const_iterator it = fcs_phonon->force_constant[2].begin(); 
            it != fcs_phonon->force_constant[2].end(); ++it) {

                for (i = 0; i < 4; ++i) atom_num4[i] = system->map_p2s[(*it).elems[i].atom][(*it).elems[i].cell];

                for (i = 0; i < 3; ++i) {
                    vec_for_v4[j][i][0] = relvec[atom_num4[1]][atom_num4[0]][i];
                    vec_for_v4[j][i][1] = relvec[atom_num4[2]][atom_num4[0]][i];
                    vec_for_v4[j][i][2] = relvec[atom_num4[3]][atom_num4[0]][i];
                }

                invmass_for_v4[j] 
                = invsqrt_mass_p[(*it).elems[0].atom] 
                * invsqrt_mass_p[(*it).elems[1].atom] 
                * invsqrt_mass_p[(*it).elems[2].atom] 
                * invsqrt_mass_p[(*it).elems[3].atom];

                ++j;     
        }

        memory->allocate(evec_index4, fcs_phonon->force_constant[2].size(), 4);

        for (i = 0; i < fcs_phonon->force_constant[2].size(); ++i) {
            for (j = 0; j < 4; ++j) {
                evec_index4[i][j] 
                = 3 * fcs_phonon->force_constant[2][i].elems[j].atom 
                    + fcs_phonon->force_constant[2][i].elems[j].xyz;
            }
        }
    }

    MPI_Bcast(&ksum_mode, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&calc_realpart, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD);
    MPI_Bcast(&atom_project_mode, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD);
    MPI_Bcast(&calc_fstate_k, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD);
    MPI_Bcast(&calc_fstate_omega, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD);

    // For tetrahedron method
    if (ksum_mode == -1) {
        memory->allocate(e_tmp, 4, nk);
        memory->allocate(f_tmp, 4, nk);
    }

    if (mympi->my_rank == 0) {

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
                for (int j = 0; j < 3; ++j) {
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
        }

        std::cout << std::endl;
        std::cout << "Estimated minimum energy difference (cm^-1) = " << domega_min << std::endl;
        std::cout << "Given epsilon (cm^-1) = " << epsilon << std::endl << std::endl;

        if (ksum_mode == 0) {
            std::cout << "Lorentzian broadening will be used." << std::endl;
        } else if (ksum_mode == 1) {
            std::cout << "Gaussian broadening will be used." << std::endl;
        } else if (ksum_mode == -1) {
            std::cout << "Tetrahedron method will be used." << std::endl;
        } else {
            error->exit("setup_relaxation", "Invalid ksum_mode");
        }
        std::cout << std::endl;
    }

    memory->deallocate(relvec);
    memory->deallocate(invsqrt_mass_p);

    epsilon *= time_ry / Hz_to_kayser;
    MPI_Bcast(&epsilon, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (mympi->my_rank == 0) {
        if (calc_fstate_omega) sym_permutation = false;
    }
    MPI_Bcast(&sym_permutation, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD);

    generate_triplet_k(use_triplet_symmetry, sym_permutation);
    // gen_pair_uniq();
    if (mympi->my_rank == 1) gensym_kpairs();

    if (mympi->my_rank == 0) {
        std::cout << " done!" << std::endl;
    }
}

void Relaxation::setup_mode_analysis()
{
    // Judge if ks_analyze_mode should be turned on or not.

    unsigned int i;

    if (mympi->my_rank == 0) {
        if (!ks_input.empty()) {
            std::cout << std::endl;
            std::cout << "KS_INPUT is given." << std::endl;
            std::cout << "Analysis on specific k points will be performed instead of thermal conductivity calculations." << std::endl;
            std::cout << std::endl;


            std::ifstream ifs_ks;
            ifs_ks.open(ks_input.c_str(), std::ios::in);
            if (!ifs_ks) error->exit("setup_relaxation", "Cannot open file KS_INPUT");

            unsigned int nlist;
            double ktmp[3];
            unsigned int snum_tmp;
            int knum_tmp;

            ifs_ks >> nlist;

            if (nlist <= 0) error->exit("setup_relaxation", "First line in KS_INPUT files should be a positive integer.");

            if (calc_fstate_k) {
                kslist_fstate_k.clear();

                for (i = 0; i < nlist; ++i) {

                    ifs_ks >> ktmp[0] >> ktmp[1] >> ktmp[2] >> snum_tmp;

                    kslist_fstate_k.push_back(KsListMode(ktmp, snum_tmp));
                }
                std::cout << "The number of entries = " << kslist_fstate_k.size() << std::endl;

            } else {
                kslist.clear();
                for (i = 0; i < nlist; ++i) {

                    ifs_ks >> ktmp[0] >> ktmp[1] >> ktmp[2] >> snum_tmp;
                    knum_tmp = kpoint->get_knum(ktmp[0], ktmp[1], ktmp[2]);

                    if (knum_tmp == -1) error->exit("setup_relaxation", "Given kpoint is not exist in given k-point grid.");
                    kslist.push_back(knum_tmp * dynamical->neval + snum_tmp);
                }
                std::cout << "The number of entries = " << kslist.size() << std::endl;
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

void Relaxation::finish_relaxation()
{
    memory->deallocate(vec_for_v3);
    memory->deallocate(invmass_for_v3);
    memory->deallocate(evec_index);
    memory->deallocate(fcs_group);
    memory->deallocate(v3_arr);
    memory->deallocate(delta_arr);

    if (ksum_mode == -1) {
        memory->deallocate(e_tmp);
        memory->deallocate(f_tmp);
    }
}

std::complex<double> Relaxation::V3(const unsigned int ks[3])
{
    /* 
    This version requires massive RAM to store cexp_phase
    */

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

            ctmp = fcs_group[i][j].fcs_val * invmass_for_v3[ielem] * std::exp(im*phase);
            ret_in += ctmp;

            ++ielem;
        }
        ret += ret_in * vec_tmp;
    }

    return ret / std::sqrt(omega[0] * omega[1] * omega[2]);
}


// std::complex<double> Relaxation::V3(const unsigned int ks[3])
// {
// 	/* 
// 	This version requires massive RAM to store cexp_phase
// 	*/
// 
// 	unsigned int i;
// 	unsigned int kn[3];
// 	unsigned int sn[3];
// 
// 	int ielem;
// 
// 	double omega[3];
// 
// 	double ret_re = 0.0;
// 	double ret_im = 0.0;
// 
// 	for (i = 0; i < 3; ++i){
// 		kn[i] = ks[i] / ns;
// 		sn[i] = ks[i] % ns;
// 		omega[i] = dynamical->eval_phonon[kn[i]][sn[i]];
// 	}
// 
// #pragma omp parallel 
// 	{
// 		std::complex<double> ctmp;
// 		double phase;
// 		std::complex<double> evecs_tmp;
// 
// #pragma omp for reduction(+: ret_re, ret_im)
// 		for (ielem = 0; ielem < fcs_phonon->force_constant[1].size(); ++ielem) {
// 
// 			phase = vec_for_v3[ielem][0][0] * kpoint->xk[kn[1]][0] 
// 			+ vec_for_v3[ielem][1][0] * kpoint->xk[kn[1]][1]
// 			+ vec_for_v3[ielem][2][0] * kpoint->xk[kn[1]][2]
// 			+ vec_for_v3[ielem][0][1] * kpoint->xk[kn[2]][0] 
// 			+ vec_for_v3[ielem][1][1] * kpoint->xk[kn[2]][1] 
// 			+ vec_for_v3[ielem][2][1] * kpoint->xk[kn[2]][2];
// 
// 			ctmp = fcs_phonon->force_constant[1][ielem].fcs_val * invmass_for_v3[ielem] * std::exp(im*phase)
// 				* dynamical->evec_phonon[kn[0]][sn[0]][evec_index[ielem][0]]
// 			* dynamical->evec_phonon[kn[1]][sn[1]][evec_index[ielem][1]]
// 			* dynamical->evec_phonon[kn[2]][sn[2]][evec_index[ielem][2]];
// 
// 			ret_re += ctmp.real();
// 			ret_im += ctmp.imag();
// 		}
// 	}
// 
// 	return (ret_re + im * ret_im) / std::sqrt(omega[0] * omega[1] * omega[2]);
// }

std::complex<double> Relaxation::V4(const unsigned int ks[4]) 
{

    unsigned int i, ielem;
    unsigned int kn[4], sn[4];

    double phase;
    double omega[4];

    std::complex<double> ctmp;
    std::complex<double> ret = std::complex<double>(0.0, 0.0);

    for (i = 0; i < 4; ++i){
        kn[i] = ks[i] / ns;
        sn[i] = ks[i] % ns;
        omega[i] = dynamical->eval_phonon[kn[i]][sn[i]];
    }

    for (ielem = 0; ielem < fcs_phonon->force_constant[2].size(); ++ielem) {

        phase = 
            vec_for_v4[ielem][0][0] * kpoint->xk[kn[1]][0] 
        + vec_for_v4[ielem][1][0] * kpoint->xk[kn[1]][1] 
        + vec_for_v4[ielem][2][0] * kpoint->xk[kn[1]][2]
        + vec_for_v4[ielem][0][1] * kpoint->xk[kn[2]][0] 
        + vec_for_v4[ielem][1][1] * kpoint->xk[kn[2]][1] 
        + vec_for_v4[ielem][2][1] * kpoint->xk[kn[2]][2]
        + vec_for_v4[ielem][0][2] * kpoint->xk[kn[3]][0] 
        + vec_for_v4[ielem][1][2] * kpoint->xk[kn[3]][1] 
        + vec_for_v4[ielem][2][2] * kpoint->xk[kn[3]][2];

        ctmp = fcs_phonon->force_constant[2][ielem].fcs_val * invmass_for_v4[ielem] * std::exp(im*phase)
            * dynamical->evec_phonon[kn[0]][sn[0]][evec_index4[ielem][0]] 
        * dynamical->evec_phonon[kn[1]][sn[1]][evec_index4[ielem][1]] 
        * dynamical->evec_phonon[kn[2]][sn[2]][evec_index4[ielem][2]] 
        * dynamical->evec_phonon[kn[3]][sn[3]][evec_index4[ielem][3]];

        ret += ctmp;
    }

    return ret / std::sqrt(omega[0] * omega[1] * omega[2] * omega[3]);
}

std::complex<double> Relaxation::V3_mode(int mode, double *xk2, double *xk3, int is, int js, double **eval, std::complex<double> ***evec)
{
    /* 
    This version requires massive RAM to store cexp_phase
    */
    int ielem;

    double phase;
    std::complex<double> ctmp = std::complex<double>(0.0, 0.0);

    for (ielem = 0; ielem < fcs_phonon->force_constant[1].size(); ++ielem) {

        phase = vec_for_v3[ielem][0][0] * xk2[0] 
        + vec_for_v3[ielem][1][0] * xk2[1]
        + vec_for_v3[ielem][2][0] * xk2[2] 
        + vec_for_v3[ielem][0][1] * xk3[0]
        + vec_for_v3[ielem][1][1] * xk3[1]
        + vec_for_v3[ielem][2][1] * xk3[2];


        ctmp += fcs_phonon->force_constant[1][ielem].fcs_val * invmass_for_v3[ielem] * std::exp(im * phase)
            * evec[0][mode][evec_index[ielem][0]] * evec[1][is][evec_index[ielem][1]] * evec[2][js][evec_index[ielem][2]];
    }

    return ctmp / std::sqrt(eval[0][mode] * eval[1][is] * eval[2][js]);
}

// void Relaxation::v3_test() {
// 
// 	int i;
// 	unsigned int stmp[3], kstmp[3];
// 	int nkplus, nkminus;
// 	double k_tmp1[3], k_tmp2[3], k_tmp3[3];
// 
// 	nkplus = 1;
// 	nkminus= kpoint->knum_minus[nkplus];
// 
// 
// 	stmp[0] = 3;
// 	stmp[1] = 1;
// 	stmp[2] = 2;
// 
// 	k_tmp1[0] = 0.4; k_tmp1[1] = 0.0; k_tmp1[2] = 0.0;
// 	k_tmp2[0] = -0.4; k_tmp2[1] = 0.4; k_tmp2[2] = 0.0;
// 	k_tmp3[0] = 0.0; k_tmp3[1] = -0.4; k_tmp3[2] = 0.0;
// 
// 	kstmp[0] = kpoint->get_knum(k_tmp1[0], k_tmp1[1], k_tmp1[2]);
// 	kstmp[1] = kpoint->get_knum(k_tmp2[0], k_tmp2[1], k_tmp2[2]);
// 	kstmp[2] = kpoint->get_knum(k_tmp3[0], k_tmp3[1], k_tmp3[2]);
// 
// 	for (i = 0; i < 3; ++i) {
// 		std::cout << std::setw(15) << kpoint->xk[kstmp[0]][i];
// 	}
// 	std::cout << std::endl;
// 	for (i = 0; i < 3; ++i) {
// 		std::cout << std::setw(15) << kpoint->xk[kstmp[1]][i];
// 	}
// 	std::cout << std::endl;
// 	for (i = 0; i < 3; ++i) {
// 		std::cout << std::setw(15) << kpoint->xk[kstmp[2]][i];
// 	}
// 	std::cout << std::endl;
// 
// 	kstmp[0] = ns * kpoint->get_knum(k_tmp1[0], k_tmp1[1], k_tmp1[2]) + stmp[0];
// 	kstmp[1] = ns * kpoint->get_knum(k_tmp2[0], k_tmp2[1], k_tmp2[2]) + stmp[1];
// 	kstmp[2] = ns * kpoint->get_knum(k_tmp3[0], k_tmp3[1], k_tmp3[2]) + stmp[2];
// 
// 
// 	double t1, t2;
// 	t1 = timer->elapsed();
// 	std::cout << std::norm(V3(kstmp[0], kstmp[1], kstmp[2]));
// 	t2 = timer->elapsed();
// 	std::cout << std::setw(15) << t2 - t1 << std::endl;
// 
// 	t1 = timer->elapsed();
// 	std::cout << std::norm(V3(kstmp));
// 	t2 = timer->elapsed();
// 	std::cout << std::setw(15) << t2 - t1 << std::endl;
// 
// 	t1 = timer->elapsed();
// 	std::cout << std::norm(V32(kstmp));
// 	t2 = timer->elapsed();
// 	std::cout << std::setw(15) << t2 - t1 << std::endl;
// 
// // 	t1 = timer->elapsed();
// // 	std::cout << std::norm(V33(kstmp));
// // 	t2 = timer->elapsed();
// // 	std::cout << std::setw(15) << t2 - t1 << std::endl;
// 
// 
// 	error->exit("v3_test", "finished!");
// 
// 
// 	unsigned int ik0, ik1;
// 	int ik2;
// 	unsigned int is0, is1, is2;
// 	double omega1, omega2;
// 	double norm;
// 	double xk_tmp[3];
// 
// 	ik0 = 2;
// 	is0 = 1;
// 
// 	for (ik1 = 0; ik1 < kpoint->nk; ++ik1) {
// 
// 		xk_tmp[0] = - (kpoint->xk[ik0][0] + kpoint->xk[ik1][0]);
// 		xk_tmp[1] = - (kpoint->xk[ik0][1] + kpoint->xk[ik1][1]);
// 		xk_tmp[2] = - (kpoint->xk[ik0][2] + kpoint->xk[ik1][2]);
// 
// 		ik2 = kpoint->get_knum(xk_tmp[0], xk_tmp[1], xk_tmp[2]);
// 
// 		if (ik2 == -1) {
// 			error->exit("hoge","hoge");
// 		}
// 
// 		for (is1 = 0; is1 < ns; ++is1) {
// 			for (is2 = 0; is2 < ns; ++is2) {
// 				omega1 = dynamical->eval_phonon[ik1][is1];
// 				omega2 = dynamical->eval_phonon[ik2][is2];
// 
// 				std::cout << std::setw(5) << ik1 << std::setw(5) << is1;
// 				std::cout << std::setw(15) << writes->in_kayser(omega1);
// 				std::cout << std::setw(5) << ik2 << std::setw(5) << is2;
// 				std::cout << std::setw(15) << writes->in_kayser(omega2);
// 
// 
// 				norm =   std::pow(std::fmod(kpoint->xk[ik0][0] + kpoint->xk[ik1][0] + kpoint->xk[ik2][0], 1.0), 2)
// 					+ std::pow(std::fmod(kpoint->xk[ik0][1] + kpoint->xk[ik1][1] + kpoint->xk[ik2][1], 1.0), 2)
// 					+ std::pow(std::fmod(kpoint->xk[ik0][2] + kpoint->xk[ik1][2] + kpoint->xk[ik2][2], 1.0), 2);
// 				std::cout << std::setw(15) << norm;
// 
// 				kstmp[0] = ns * ik0 + is0;
// 				kstmp[1] = ns * ik1 + is1;
// 				kstmp[2] = ns * ik2 + is2;
// 				std::cout << std::setw(15) << std::norm(V3(kstmp[0],kstmp[1],kstmp[2])) << std::endl;
// 			}
// 		}
// 	}
// 
// 
// }
// 
// 
// void Relaxation::v4_test() {
// 
// 	int i;
// 	unsigned int stmp[4], kstmp[4];
// 	int nkplus, nkminus;
// 
// 	nkplus = 2;
// 	nkminus= kpoint->knum_minus[nkplus];
// 
// 	stmp[0] = 0;
// 	stmp[1] = 1;
// 	stmp[2] = 2;
// 	stmp[3] = 0;
// 
// 	for (i = 0; i < 4; ++i) {
// 		std::cout << std::setw(15) << kpoint->xk[nkplus][i];
// 	}
// 	std::cout << std::endl;
// 
// 	for (i = 0; i < 4; ++i) {
// 		kstmp[i] = dynamical->neval * nkplus + stmp[i];
// 	}
// 
// 
// 	std::cout << V4(kstmp) << std::endl;
// 
// 	for (i = 0; i < 4; ++i) {
// 		kstmp[i] = dynamical->neval * nkminus + stmp[i];
// 	}
// 	for (i = 0; i < 4; ++i) {
// 		std::cout << std::setw(15) << kpoint->xk[nkminus][i];
// 	}
// 	std::cout << std::endl;
// 
// 	std::cout << V4(kstmp) << std::endl;
// 
// 
// 	error->exit("v4_test", "finished!");
// }


void Relaxation::calc_realpart_V4(const unsigned int N, double *T, const double omega, const unsigned int knum, const unsigned int snum, double *ret)
{
    unsigned int i, ik, is;
    unsigned int arr[4];
    double n1, omega1;
    double v4_tmp, T_tmp;

    for (i = 0; i < N; ++i) ret[i] = 0.0;

    arr[0] = ns * kpoint->knum_minus[knum] + snum;
    arr[1] = ns * knum + snum;

    for (ik = 0; ik < nk; ++ik) {
        for (is = 0; is < ns; ++is) {

            arr[2] = ns * ik + is;
            arr[3] = ns * kpoint->knum_minus[ik] + is;

            v4_tmp = V4(arr).real();

            omega1 = dynamical->eval_phonon[ik][is];

            for (i = 0; i < N; ++i) {
                T_tmp = T[i];
                n1 = phonon_thermodynamics->fB(omega1, T_tmp);

                ret[i] += v4_tmp * (2.0 * n1 + 1.0);
            }
        }
    }

    for (i = 0; i < N; ++i) ret[i] *= - 1.0 / (8.0 * static_cast<double>(nk));
}

void Relaxation::calc_damping(const unsigned int N, double *T, const double omega, 
                              const unsigned int knum, const unsigned int snum, double *ret)
{
    unsigned int i;
    unsigned int ik, jk;
    unsigned int is, js; 
    unsigned int arr[3];

    double T_tmp;
    double n1, n2;
    double v3_tmp;
    double xk_tmp[3];
    double omega_inner[2];

    double multi;
    double f1, f2;

    for (i = 0; i < N; ++i) ret[i] = 0.0;

    arr[0] = ns * kpoint->knum_minus[knum] + snum;

    unsigned int nkx = kpoint->nkx;
    unsigned int nky = kpoint->nky;
    unsigned int nkz = kpoint->nkz;

    int iloc, jloc, kloc;

    for (ik = 0; ik < nk; ++ik) {

        xk_tmp[0] = kpoint->xk[knum][0] - kpoint->xk[ik][0];
        xk_tmp[1] = kpoint->xk[knum][1] - kpoint->xk[ik][1];
        xk_tmp[2] = kpoint->xk[knum][2] - kpoint->xk[ik][2];

        iloc = (nint(xk_tmp[0]*static_cast<double>(nkx) + static_cast<double>(2*nkx))) % nkx;
        jloc = (nint(xk_tmp[1]*static_cast<double>(nky) + static_cast<double>(2*nky))) % nky;
        kloc = (nint(xk_tmp[2]*static_cast<double>(nkz) + static_cast<double>(2*nkz))) % nkz;

        jk = kloc + nkz * jloc + nky * nkz * iloc;

        for (is = 0; is < ns; ++is){
            for (js = 0; js < ns; ++js){

                arr[1] = ns * ik + is;
                arr[2] = ns * jk + js;

                if (arr[1] > arr[2]) continue;

                if (arr[1] == arr[2]) {
                    multi = 1.0;
                } else {
                    multi = 2.0;
                }

                v3_tmp = std::norm(V3(arr));

                omega_inner[0] = dynamical->eval_phonon[ik][is];
                omega_inner[1] = dynamical->eval_phonon[jk][js];

                for (i = 0; i < N; ++i) {
                    T_tmp = T[i];

                    f1 = phonon_thermodynamics->fB(omega_inner[0], T_tmp);
                    f2 = phonon_thermodynamics->fB(omega_inner[1], T_tmp);
                    n1 =  f1 + f2 + 1.0;
                    n2 =  f1 - f2;

                    if (ksum_mode == 0) {
                        ret[i] += v3_tmp *multi
                            * ( - n1 * delta_lorentz(omega + omega_inner[0] + omega_inner[1])
                            + n1 * delta_lorentz(omega - omega_inner[0] - omega_inner[1])
                            - n2 * delta_lorentz(omega - omega_inner[0] + omega_inner[1])
                            + n2 * delta_lorentz(omega + omega_inner[0] - omega_inner[1]));
                    } else if (ksum_mode == 1) {
                        ret[i] += v3_tmp * multi
                            * ( - n1 * delta_gauss(omega + omega_inner[0] + omega_inner[1])
                            + n1 * delta_gauss(omega - omega_inner[0] - omega_inner[1])
                            - n2 * delta_gauss(omega - omega_inner[0] + omega_inner[1])
                            + n2 * delta_gauss(omega + omega_inner[0] - omega_inner[1]));
                    }
                }
            }
        }
    }

    for (i = 0; i < N; ++i) ret[i] *=  pi * std::pow(0.5, 4) / static_cast<double>(nk);
}

void Relaxation::calc_damping_tune(const unsigned int N, double *T, const double omega,
                                   const unsigned int knum, const unsigned int snum, double *ret)
{
    unsigned int i;
    int ik;
    unsigned int jk;
    unsigned int is, js;
    unsigned int arr[3];

    double T_tmp;
    double n1, n2;
    double v3_tmp;
    double xk_tmp[3];
    double omega_inner[2];

    double multi;
    double f1, f2;
    double ret_tmp;

    for (i = 0; i < N; ++i) ret[i] = 0.0;

    arr[0] = ns * kpoint->knum_minus[knum] + snum;

    unsigned int nkx = kpoint->nkx;
    unsigned int nky = kpoint->nky;
    unsigned int nkz = kpoint->nkz;

    int iloc, jloc, kloc;

    //	double **xks = kpoint->xk;
    //	double **eval= dynamical->eval_phonon;



    //#pragma offload target(mic) in(nkx, nky, nkz, knum, snum) inout(v3_arr:length(nk*ns*ns)) inout(delta_arr:length(4*nk*ns*ns)) inout(xks:length(3*nk)) inout(eval:length(nk*ns))
    {
#pragma omp parallel for private(xk_tmp, iloc, jloc, kloc, jk, is, js, arr, multi, omega_inner)
        for (ik = 0; ik < nk; ++ik) {

            xk_tmp[0] = kpoint->xk[knum][0] - kpoint->xk[ik][0];
            xk_tmp[1] = kpoint->xk[knum][1] - kpoint->xk[ik][1];
            xk_tmp[2] = kpoint->xk[knum][2] - kpoint->xk[ik][2];

            // 			xk_tmp[0] = xks[knum][0] - xks[ik][0];
            // 			xk_tmp[1] = xks[knum][1] - xks[ik][1];
            // 			xk_tmp[2] = xks[knum][2] - xks[ik][2];

            iloc = (nint(xk_tmp[0]*static_cast<double>(nkx) + static_cast<double>(2*nkx))) % nkx;
            jloc = (nint(xk_tmp[1]*static_cast<double>(nky) + static_cast<double>(2*nky))) % nky;
            kloc = (nint(xk_tmp[2]*static_cast<double>(nkz) + static_cast<double>(2*nkz))) % nkz;

            jk = kloc + nkz * jloc + nky * nkz * iloc;

            for (is = 0; is < ns; ++is){
                for (js = 0; js < ns; ++js){

                    arr[1] = ns * ik + is;
                    arr[2] = ns * jk + js;

                    if (arr[1] > arr[2]) continue;

                    if (arr[1] == arr[2]) {
                        multi = 1.0;
                    } else {
                        multi = 2.0;
                    }

                    arr[0] = ns * kpoint->knum_minus[knum] + snum;
                    v3_arr[ik][ns * is + js] = std::norm(V3(arr)) * multi;

                    omega_inner[0] = dynamical->eval_phonon[ik][is];
                    omega_inner[1] = dynamical->eval_phonon[jk][js];

                    delta_arr[ik][ns * is + js][0] = delta_lorentz(omega + omega_inner[0] + omega_inner[1]);
                    delta_arr[ik][ns * is + js][1] = delta_lorentz(omega - omega_inner[0] - omega_inner[1]);
                    delta_arr[ik][ns * is + js][2] = delta_lorentz(omega - omega_inner[0] + omega_inner[1]);
                    delta_arr[ik][ns * is + js][3] = delta_lorentz(omega + omega_inner[0] - omega_inner[1]);
                }
            }
        }
    }

    for (i = 0; i < N; ++i) {
        T_tmp = T[i];
        ret_tmp = 0.0;

#pragma omp parallel for private(xk_tmp, iloc, jloc, kloc, jk, is, js, arr, omega_inner, n1, n2, f1, f2), reduction(+:ret_tmp)
        for (ik = 0; ik < nk; ++ik) {
            xk_tmp[0] = kpoint->xk[knum][0] - kpoint->xk[ik][0];
            xk_tmp[1] = kpoint->xk[knum][1] - kpoint->xk[ik][1];
            xk_tmp[2] = kpoint->xk[knum][2] - kpoint->xk[ik][2];

            iloc = (nint(xk_tmp[0]*static_cast<double>(nkx) + static_cast<double>(2*nkx))) % nkx;
            jloc = (nint(xk_tmp[1]*static_cast<double>(nky) + static_cast<double>(2*nky))) % nky;
            kloc = (nint(xk_tmp[2]*static_cast<double>(nkz) + static_cast<double>(2*nkz))) % nkz;

            jk = kloc + nkz * jloc + nky * nkz * iloc;

            for (is = 0; is < ns; ++is){
                for (js = 0; js < ns; ++js) {

                    arr[1] = ns * ik + is;
                    arr[2] = ns * jk + js;

                    if (arr[1] > arr[2]) continue;

                    omega_inner[0] = dynamical->eval_phonon[ik][is];
                    omega_inner[1] = dynamical->eval_phonon[jk][js];

                    f1 = phonon_thermodynamics->fB(omega_inner[0], T_tmp);
                    f2 = phonon_thermodynamics->fB(omega_inner[1], T_tmp);
                    n1 =  f1 + f2 + 1.0;
                    n2 =  f1 - f2;

                    ret_tmp += v3_arr[ik][ns * is + js]
                    * ( -n1 * delta_arr[ik][ns * is + js][0] + n1 * delta_arr[ik][ns * is + js][1]
                    -n2 * delta_arr[ik][ns * is + js][2] + n2 * delta_arr[ik][ns * is + js][3]);

                }
            }
        }
        ret[i] = ret_tmp;
    }

    for (i = 0; i < N; ++i) ret[i] *=  pi * std::pow(0.5, 4) / static_cast<double>(nk);
}


void Relaxation::calc_damping2(const unsigned int N, double *T, const double omega, 
                               const unsigned int ik_in, const unsigned int snum, double *ret)
{
    // This is under test.

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

    memory->allocate(v3_arr, pair_uniq[ik_in].size(), ns * ns);
    memory->allocate(delta_arr, pair_uniq[ik_in].size(), ns * ns, 2);


#pragma omp parallel for private(multi, knum, knum_minus, arr, k1, k2, is, js, omega_inner)
    for (ik = 0; ik < pair_uniq[ik_in].size(); ++ik) {
        multi = static_cast<double>(pair_uniq[ik_in][ik].group.size());
        knum = kpoint->k_reduced[ik_in][0];
        knum_minus = kpoint->knum_minus[knum];

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

                delta_arr[ik][ns * is + js][0] = delta_lorentz(omega - omega_inner[0] - omega_inner[1])
                    - delta_lorentz(omega + omega_inner[0] + omega_inner[1]);
                delta_arr[ik][ns * is + js][1] = delta_lorentz(omega - omega_inner[0] + omega_inner[1])
                    - delta_lorentz(omega + omega_inner[0] - omega_inner[1]);
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

    //     for (ik = 0; ik < pair_uniq[ik_in].size(); ++ik) {
    // 
    //         multi = static_cast<double>(pair_uniq[ik_in][ik].group.size());
    //         knum = kpoint->k_reduced[ik_in][0];
    //         knum_minus = kpoint->knum_minus[knum];
    // 
    //         arr[0] = ns * knum_minus + snum;
    // 
    //         k1 = pair_uniq[ik_in][ik].group[0].ks[0];
    //         k2 = pair_uniq[ik_in][ik].group[0].ks[1];
    // 
    //         for (is = 0; is < ns; ++is) {
    //             for (js = 0; js < ns; ++js) {
    //                 arr[1] = ns * k1 + is;
    //                 arr[2] = ns * k2 + js;
    // 
    //                 omega_inner[0] = dynamical->eval_phonon[k1][is];
    //                 omega_inner[1] = dynamical->eval_phonon[k2][js];
    // 
    //                 v3_tmp = std::norm(V3(arr));
    // 
    //                 for (i = 0; i < N; ++i) {
    //                     T_tmp = T[i];
    // 
    //                     if (conductivity->use_classical_Cv == 0) {
    //                         n1 = phonon_thermodynamics->fB(omega_inner[0], T_tmp) + phonon_thermodynamics->fB(omega_inner[1], T_tmp) + 1.0;
    //                         n2 = phonon_thermodynamics->fB(omega_inner[0], T_tmp) - phonon_thermodynamics->fB(omega_inner[1], T_tmp);
    //                     } else if (conductivity->use_classical_Cv == 1) {
    //                         n1 = phonon_thermodynamics->fC(omega_inner[0], T_tmp) + phonon_thermodynamics->fC(omega_inner[1], T_tmp) + 1.0;
    //                         n2 = phonon_thermodynamics->fC(omega_inner[0], T_tmp) - phonon_thermodynamics->fC(omega_inner[1], T_tmp);
    //                     }
    // 
    //                     if (ksum_mode == 0) {
    //                         ret[i] += v3_tmp * multi
    //                             * ( - n1 * delta_lorentz(omega + omega_inner[0] + omega_inner[1])
    //                             + n1 * delta_lorentz(omega - omega_inner[0] - omega_inner[1])
    //                             - n2 * delta_lorentz(omega - omega_inner[0] + omega_inner[1])
    //                             + n2 * delta_lorentz(omega + omega_inner[0] - omega_inner[1]));
    //                     } else if (ksum_mode == 1) {
    //                         ret[i] += v3_tmp * multi
    //                             * ( - n1 * delta_gauss(omega + omega_inner[0] + omega_inner[1])
    //                             + n1 * delta_gauss(omega - omega_inner[0] - omega_inner[1])
    //                             - n2 * delta_gauss(omega - omega_inner[0] + omega_inner[1])
    //                             + n2 * delta_gauss(omega + omega_inner[0] - omega_inner[1]));
    //                     }
    //                 }
    //             }
    //         }
    //     }

    memory->deallocate(v3_arr);
    memory->deallocate(delta_arr);

    for (i = 0; i < N; ++i) ret[i] *=  pi * std::pow(0.5, 4) / static_cast<double>(nk);
}

void Relaxation::calc_damping_tetra(const unsigned int N, double *T, const double omega,
                                    const unsigned int knum, const unsigned int snum, double *ret)
{
    unsigned int i, j;
    unsigned int is, js, ik, jk;
    unsigned int ks_tmp[3];

    double xk_tmp[3];
    double n1, n2;
    double *v3_tmp;
    double **omega_inner;

    unsigned int nkx = kpoint->nkx;
    unsigned int nky = kpoint->nky;
    unsigned int nkz = kpoint->nkz;

    int iloc, jloc, kloc;

    for (i = 0; i < N; ++i) ret[i] = 0.0;

    ks_tmp[0] = ns * kpoint->knum_minus[knum] + snum;

    memory->allocate(v3_tmp, nk);
    memory->allocate(omega_inner, nk, 2);

    for (is = 0; is < ns; ++is){
        for (js = 0; js < ns; ++js){

            for (ik = 0; ik < nk; ++ik) {

                xk_tmp[0] = kpoint->xk[knum][0] - kpoint->xk[ik][0];
                xk_tmp[1] = kpoint->xk[knum][1] - kpoint->xk[ik][1];
                xk_tmp[2] = kpoint->xk[knum][2] - kpoint->xk[ik][2];

                iloc = (nint(xk_tmp[0]*static_cast<double>(nkx) + static_cast<double>(2*nkx))) % nkx;
                jloc = (nint(xk_tmp[1]*static_cast<double>(nky) + static_cast<double>(2*nky))) % nky;
                kloc = (nint(xk_tmp[2]*static_cast<double>(nkz) + static_cast<double>(2*nkz))) % nkz;

                jk = kloc + nkz * jloc + nky * nkz * iloc;

                ks_tmp[1] = ik * ns + is;
                ks_tmp[2] = jk * ns + js;

                omega_inner[ik][0] = dynamical->eval_phonon[ik][is];
                omega_inner[ik][1] = dynamical->eval_phonon[jk][js];

                v3_tmp[ik] = std::norm(V3(ks_tmp));

                // e_tmp[0][kcount] = -omega_inner[kcount][0] - omega_inner[kcount][1];
                e_tmp[1][ik] = omega_inner[ik][0] + omega_inner[ik][1];
                e_tmp[2][ik] = omega_inner[ik][0] - omega_inner[ik][1];
                e_tmp[3][ik] = -omega_inner[ik][0] + omega_inner[ik][1];
            }

            for (j = 0; j < N; ++j){
                for (i = 0; i < nk; ++i){

                    if (conductivity->use_classical_Cv == 0) {
                        n1 = phonon_thermodynamics->fB(omega_inner[i][0], T[j]) + phonon_thermodynamics->fB(omega_inner[i][1], T[j]) + 1.0;
                        n2 = phonon_thermodynamics->fB(omega_inner[i][0], T[j]) - phonon_thermodynamics->fB(omega_inner[i][1], T[j]);
                    } else if (conductivity->use_classical_Cv == 1) {
                        n1 = phonon_thermodynamics->fC(omega_inner[i][0], T[j]) + phonon_thermodynamics->fC(omega_inner[i][1], T[j]) + 1.0;
                        n2 = phonon_thermodynamics->fC(omega_inner[i][0], T[j]) - phonon_thermodynamics->fC(omega_inner[i][1], T[j]);
                    }

                    // f_tmp[0][i] = -v3_tmp[i] * n1;
                    f_tmp[1][i] = v3_tmp[i] * n1;
                    f_tmp[2][i] = -v3_tmp[i] * n2;
                    f_tmp[3][i] = v3_tmp[i] * n2;
                }

                for (i = 1; i < 4; ++i) {
                    ret[j] += integration->do_tetrahedron(e_tmp[i], f_tmp[i], omega);
                }

            }
        }
    }

    for (i = 0; i < N; ++i) {
        ret[i] *=  pi * std::pow(0.5, 4);
    }

    memory->deallocate(v3_tmp);
    memory->deallocate(omega_inner);
}

void Relaxation::calc_damping_atom(const unsigned  int N, double *T, const double omega, 
                                   const unsigned int knum, const unsigned int snum, double ***ret)
{
    unsigned int i, j, iks;
    unsigned int ik, jk;
    unsigned int is, js;
    unsigned int iat, jat;
    unsigned int arr[3];

    double T_tmp;
    double n1, n2;
    double v3_tmp, v3_tmp2;
    double xk_tmp[3];
    double omega_inner[2];

    double proj1, proj2;

    unsigned int natmin = system->natmin;

    for (i = 0; i < N; ++i) {
        for (iat = 0; iat < natmin; ++iat) {
            for (jat = 0; jat < natmin; ++jat) {
                ret[i][iat][jat] = 0.0;
            }
        }
    }

    arr[0] = ns * kpoint->knum_minus[knum] + snum;

    unsigned int nkx = kpoint->nkx;
    unsigned int nky = kpoint->nky;
    unsigned int nkz = kpoint->nkz;

    int iloc, jloc, kloc;
    unsigned int nks2 = nk * ns * ns;

    for (iks = mympi->my_rank; iks < nks2; iks += mympi->nprocs) {

        ik = iks / (ns * ns);
        is = (iks - ik * ns * ns) / ns;
        js = iks - ik * ns * ns - is * ns;

        xk_tmp[0] = kpoint->xk[knum][0] - kpoint->xk[ik][0];
        xk_tmp[1] = kpoint->xk[knum][1] - kpoint->xk[ik][1];
        xk_tmp[2] = kpoint->xk[knum][2] - kpoint->xk[ik][2];

        iloc = (nint(xk_tmp[0]*static_cast<double>(nkx) + static_cast<double>(2*nkx))) % nkx;
        jloc = (nint(xk_tmp[1]*static_cast<double>(nky) + static_cast<double>(2*nky))) % nky;
        kloc = (nint(xk_tmp[2]*static_cast<double>(nkz) + static_cast<double>(2*nkz))) % nkz;

        jk = kloc + nkz * jloc + nky * nkz * iloc;

        arr[1] = ns * ik + is;
        arr[2] = ns * jk + js;

        omega_inner[0] = dynamical->eval_phonon[ik][is];
        omega_inner[1] = dynamical->eval_phonon[jk][js];

        v3_tmp = std::norm(V3(arr));

        for (i = 0; i < N; ++i) {
            T_tmp = T[i];

            n1 = phonon_thermodynamics->fB(omega_inner[0], T_tmp) + phonon_thermodynamics->fB(omega_inner[1], T_tmp) + 1.0;
            n2 = phonon_thermodynamics->fB(omega_inner[0], T_tmp) - phonon_thermodynamics->fB(omega_inner[1], T_tmp);

            if (ksum_mode == 0) {
                v3_tmp2 = v3_tmp 
                    * (- n1 * delta_lorentz(omega + omega_inner[0] + omega_inner[1])
                    + n1 * delta_lorentz(omega - omega_inner[0] - omega_inner[1])
                    - n2 * delta_lorentz(omega - omega_inner[0] + omega_inner[1])
                    + n2 * delta_lorentz(omega + omega_inner[0] - omega_inner[1]));
            } else if (ksum_mode == 1) {
                v3_tmp2 = v3_tmp
                    * (- n1 * delta_gauss(omega + omega_inner[0] + omega_inner[1])
                    + n1 * delta_gauss(omega - omega_inner[0] - omega_inner[1])
                    - n2 * delta_gauss(omega - omega_inner[0] + omega_inner[1])
                    + n2 * delta_gauss(omega + omega_inner[0] - omega_inner[1]));
            }

            for (iat = 0; iat < natmin; ++iat) {
                proj1 = 0.0;
                for (j = 0; j < 3; ++j) proj1 += std::norm(dynamical->evec_phonon[ik][is][3*iat + j]);

                for (jat = 0; jat < natmin; ++jat) {
                    proj2 = 0.0;
                    for (j = 0; j < 3; ++j) proj2 += std::norm(dynamical->evec_phonon[jk][js][3*jat + j]);

                    ret[i][iat][jat] += v3_tmp2 * proj1 * proj2;

                }
            }
        }


    }

    for (i = 0; i < N; ++i) {
        for (iat = 0; iat < natmin; ++iat) {
            for (jat = 0; jat < natmin; ++jat) {
                ret[i][iat][jat] *=  pi * std::pow(0.5, 4) / static_cast<double>(nk);
            }
        }
    }

}

void Relaxation::calc_damping_tetra_atom(const unsigned int N, double *T, const double omega,
                                         const unsigned int knum, const unsigned int snum, double ***ret)
{
    unsigned int i, j, k;
    unsigned int is, js, ik, jk;
    unsigned int iat, jat;
    unsigned int ks_tmp[3];

    double xk_tmp[3];
    double n1, n2;
    double *v3_tmp;
    double **omega_inner;

    double ****f_tmp_atom;
    double ***v3_tmp_proj;
    double proj1, proj2;

    unsigned int nkx = kpoint->nkx;
    unsigned int nky = kpoint->nky;
    unsigned int nkz = kpoint->nkz;

    int iloc, jloc, kloc;

    memory->allocate(f_tmp_atom, system->natmin, system->natmin, 4, nk);

    unsigned int natmin = system->natmin;

    for (iat = 0; iat < natmin; ++iat) {
        for (jat = 0; jat < natmin; ++jat) {
            for (i = 0; i < N; ++i) {
                ret[iat][jat][i] = 0.0;
            }
        }
    }

    ks_tmp[0] = ns * kpoint->knum_minus[knum] + snum;

    memory->allocate(v3_tmp, nk);
    memory->allocate(v3_tmp_proj, natmin, natmin, nk);

    memory->allocate(omega_inner, nk, 2);

    for (is = 0; is < ns; ++is) {
        for (js = 0; js < ns; ++js)	{

            for (ik = 0; ik < nk; ++ik) {

                xk_tmp[0] = kpoint->xk[knum][0] - kpoint->xk[ik][0];
                xk_tmp[1] = kpoint->xk[knum][1] - kpoint->xk[ik][1];
                xk_tmp[2] = kpoint->xk[knum][2] - kpoint->xk[ik][2];

                iloc = (nint(xk_tmp[0]*static_cast<double>(nkx) + static_cast<double>(2*nkx))) % nkx;
                jloc = (nint(xk_tmp[1]*static_cast<double>(nky) + static_cast<double>(2*nky))) % nky;
                kloc = (nint(xk_tmp[2]*static_cast<double>(nkz) + static_cast<double>(2*nkz))) % nkz;

                jk = kloc + nkz * jloc + nky * nkz * iloc;

                ks_tmp[1] = ik * ns + is;
                ks_tmp[2] = jk * ns + js;

                omega_inner[ik][0] = dynamical->eval_phonon[ik][is];
                omega_inner[ik][1] = dynamical->eval_phonon[jk][js];

                v3_tmp[ik] = std::norm(V3(ks_tmp));

                for (iat = 0; iat < natmin; ++iat) {
                    proj1 = 0.0;
                    for (k = 0; k < 3; ++k) proj1 += std::norm(dynamical->evec_phonon[ik][is][3*iat + k]);

                    for (jat = 0; jat < natmin; ++jat) {
                        proj2 = 0.0;
                        for (k = 0; k < 3; ++k) proj2 += std::norm(dynamical->evec_phonon[jk][js][3*jat + k]);

                        v3_tmp_proj[iat][jat][ik] = v3_tmp[ik]*proj1*proj2;
                    }
                }

                // e_tmp[0][kcount] = -omega_inner[kcount][0] - omega_inner[kcount][1];
                e_tmp[1][ik] = omega_inner[ik][0] + omega_inner[ik][1];
                e_tmp[2][ik] = omega_inner[ik][0] - omega_inner[ik][1];
                e_tmp[3][ik] = -omega_inner[ik][0] + omega_inner[ik][1];
            }

            for (j = 0; j < N; ++j){
                for (iat = 0; iat < natmin; ++iat) {
                    for (jat = 0; jat < natmin; ++jat) {

                        for (i = 0; i < nk; ++i){

                            if (conductivity->use_classical_Cv == 0) {
                                n1 = phonon_thermodynamics->fB(omega_inner[i][0], T[j]) + phonon_thermodynamics->fB(omega_inner[i][1], T[j]) + 1.0;
                                n2 = phonon_thermodynamics->fB(omega_inner[i][0], T[j]) - phonon_thermodynamics->fB(omega_inner[i][1], T[j]);
                            } else if (conductivity->use_classical_Cv == 1) {
                                n1 = phonon_thermodynamics->fC(omega_inner[i][0], T[j]) + phonon_thermodynamics->fC(omega_inner[i][1], T[j]) + 1.0;
                                n2 = phonon_thermodynamics->fC(omega_inner[i][0], T[j]) - phonon_thermodynamics->fC(omega_inner[i][1], T[j]);
                            }

                            // f_tmp[0][i] = -v3_tmp[i] * n1;
                            f_tmp_atom[iat][jat][1][i]=  v3_tmp_proj[iat][jat][i] * n1;
                            f_tmp_atom[iat][jat][2][i]= -v3_tmp_proj[iat][jat][i] * n2;
                            f_tmp_atom[iat][jat][3][i] = v3_tmp_proj[iat][jat][i] * n2;
                        }

                        for (i = 1; i < 4; ++i) {
                            ret[iat][jat][j] += integration->do_tetrahedron(e_tmp[i], f_tmp_atom[iat][jat][i], omega);
                        }
                    }
                }
            }

        }
    }

    for (iat = 0; iat < natmin; ++iat){
        for (jat = 0; jat < natmin; ++jat) {
            for (i = 0; i < N; ++i) {
                ret[iat][jat][i] *=  pi * std::pow(0.5, 4);
            }
        }
    }

    memory->deallocate(v3_tmp);
    memory->deallocate(omega_inner);
    memory->deallocate(f_tmp_atom);
    memory->deallocate(v3_tmp_proj);
}

void Relaxation::calc_frequency_resolved_final_state(const unsigned int N, double *T, const double omega0, 
                                                     const double omega, const unsigned int knum, const unsigned int snum, double *ret)
{
    unsigned int i;
    unsigned int ik, jk;
    unsigned int is, js; 
    unsigned int arr[3];

    double T_tmp;
    double n1, n2;
    double v3_tmp;
    double xk_tmp[3];
    double omega_inner[2];

    double multi;
    double f1, f2;

    double *ret_mpi;

    memory->allocate(ret_mpi, N);

    for (i = 0; i < N; ++i) ret_mpi[i] = 0.0;

    arr[0] = ns * kpoint->knum_minus[knum] + snum;

    unsigned int nkx = kpoint->nkx;
    unsigned int nky = kpoint->nky;
    unsigned int nkz = kpoint->nkz;

    int iloc, jloc, kloc;

    for (ik = mympi->my_rank; ik < nk; ik += mympi->nprocs) {

        xk_tmp[0] = kpoint->xk[knum][0] - kpoint->xk[ik][0];
        xk_tmp[1] = kpoint->xk[knum][1] - kpoint->xk[ik][1];
        xk_tmp[2] = kpoint->xk[knum][2] - kpoint->xk[ik][2];

        iloc = (nint(xk_tmp[0]*static_cast<double>(nkx) + static_cast<double>(2*nkx))) % nkx;
        jloc = (nint(xk_tmp[1]*static_cast<double>(nky) + static_cast<double>(2*nky))) % nky;
        kloc = (nint(xk_tmp[2]*static_cast<double>(nkz) + static_cast<double>(2*nkz))) % nkz;

        jk = kloc + nkz * jloc + nky * nkz * iloc;

        for (is = 0; is < ns; ++is){

            arr[1] = ns * ik + is;
            omega_inner[0] = dynamical->eval_phonon[ik][is];

            for (js = 0; js < ns; ++js){

                arr[2] = ns * jk + js;
                omega_inner[1] = dynamical->eval_phonon[jk][js];

                v3_tmp = std::norm(V3(arr));

                for (i = 0; i < N; ++i) {
                    T_tmp = T[i];

                    f1 = phonon_thermodynamics->fB(omega_inner[0], T_tmp);
                    f2 = phonon_thermodynamics->fB(omega_inner[1], T_tmp);
                    n1 =  f1 + f2 + 1.0;
                    n2 =  f1 - f2;

                    if (ksum_mode == 0) {
                        ret_mpi[i] += v3_tmp * delta_lorentz(omega - omega_inner[0])
                            * ( - n1 * delta_lorentz(omega0 + omega_inner[0] + omega_inner[1])
                            + n1 * delta_lorentz(omega0 - omega_inner[0] - omega_inner[1])
                            - n2 * delta_lorentz(omega0 - omega_inner[0] + omega_inner[1])
                            + n2 * delta_lorentz(omega0 + omega_inner[0] - omega_inner[1]));
                    } else if (ksum_mode == 1) {
                        ret_mpi[i] += v3_tmp * delta_gauss(omega - omega_inner[0])
                            * ( - n1 * delta_gauss(omega0 + omega_inner[0] + omega_inner[1])
                            + n1 * delta_gauss(omega0 - omega_inner[0] - omega_inner[1])
                            - n2 * delta_gauss(omega0 - omega_inner[0] + omega_inner[1])
                            + n2 * delta_gauss(omega0 + omega_inner[0] - omega_inner[1]));
                    }
                }
            }
        }
    }

    for (i = 0; i < N; ++i) ret_mpi[i] *=  pi * std::pow(0.5, 4) / static_cast<double>(nk);

    MPI_Reduce(ret_mpi, ret, N, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    memory->deallocate(ret_mpi);
}


void Relaxation::calc_frequency_resolved_final_state2(const unsigned int N, double *T, const double omega0, 
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

    for (i = 0; i < N; ++i) {
        for (j = 0; j < M; ++j) {
            ret[i][j] = 0.0;
        }
    }

    for (ik = 0; ik < pair_uniq[ik_in].size(); ++ik) {

        multi = static_cast<double>(pair_uniq[ik_in][ik].group.size());
        knum = kpoint->k_reduced[ik_in][0];
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

                    prod_tmp[0] = n1 * (delta_lorentz(omega0 - omega_inner[0] - omega_inner[1]) 
                        - delta_lorentz(omega0 + omega_inner[0] + omega_inner[1]));
                    prod_tmp[1] = n2 * (delta_lorentz(omega0 + omega_inner[0] - omega_inner[1])
                        - delta_lorentz(omega0 - omega_inner[0] + omega_inner[1]));

                    for (j = 0; j < M; ++j) {
                        ret[i][j] += v3_tmp * multi * delta_lorentz(omega[j] - omega_inner[0])
                            * (prod_tmp[0] + prod_tmp[1]);
                    }
                }
            }
        }
    }
    for (i = 0; i < N; ++i) {
        for (j = 0; j < M; ++j) {
            ret[i][j] *=  pi * std::pow(0.5, 4) / static_cast<double>(nk);
        }
    }

}

inline double Relaxation::delta_lorentz(const double omega)
{
    return epsilon / (omega*omega + epsilon*epsilon) / pi;
}

double Relaxation::delta_gauss(const double omega)
{
    return std::exp(- omega * omega / (epsilon * epsilon)) / (epsilon * std::sqrt(pi));
}

void Relaxation::compute_mode_tau()
{
    unsigned int i, j;
    unsigned int NT;
    unsigned int knum, snum;

    double Tmax = system->Tmax;
    double Tmin = system->Tmin;
    double dT = system->dT;
    double omega;
    double *T_arr;

    std::ofstream ofs_mode_tau;
    std::string file_mode_tau;


    NT = static_cast<unsigned int>((Tmax - Tmin) / dT);
    memory->allocate(T_arr, NT);
    for (i = 0; i < NT; ++i) T_arr[i] = Tmin + static_cast<double>(i)*dT;


    if (calc_fstate_k) {

        // Momentum-resolved final state amplitude

        int k, l;
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


        if (mympi->my_rank == 0) {
            file_mode_tau = input->job_title + ".fstate_k";
            ofs_mode_tau.open(file_mode_tau.c_str(), std::ios::out);
            if (!ofs_mode_tau) error->exit("compute_mode_tau", "Cannot open file file_mode_tau");

            ofs_mode_tau << "## Momentum-resolved final state amplitude for given modes" << std::endl;
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

            dynamical->eval_k(xk1, kvec, fcs_phonon->fc2, eval[0], evec[0], true);
            for (j = 0; j < ns; ++j) eval[0][j] = dynamical->freq(eval[0][j]);


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

                    dynamical->eval_k(xk2, kvec, fcs_phonon->fc2, eval[1], evec[1], true);

                    for (l = 0; l < 3; ++l) kvec[l] = xk3[l];
                    rotvec(kvec, kvec, system->rlavec_p, 'T');
                    norm = std::sqrt(kvec[0] * kvec[0] + kvec[1] * kvec[1] + kvec[2] * kvec[2]);

                    if (norm > eps) for (l = 0; l < 3; ++l) kvec[l] /= norm;

                    dynamical->eval_k(xk3, kvec, fcs_phonon->fc2, eval[2], evec[2], true);

                    for (l = 0; l < ns; ++l) {
                        eval[1][l] = dynamical->freq(eval[1][l]);
                        eval[2][l] = dynamical->freq(eval[2][l]);
                    }

                    for (is = 0; is < ns; ++is) {
                        for (js = 0; js < ns; ++js) {
                            V3norm = std::norm(V3_mode(mode, xk2, xk3, is, js, eval, evec));

                            delta_tmp[0] = delta_lorentz(eval[0][mode] - eval[1][is] - eval[2][js])
                                - delta_lorentz(eval[0][mode] + eval[1][is] + eval[2][js]);
                            delta_tmp[1] = delta_lorentz(eval[0][mode] + eval[1][is] - eval[2][js])
                                - delta_lorentz(eval[0][mode] - eval[1][is] + eval[2][js]);

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

                }

                memory->deallocate(gamma_k);
                memory->deallocate(gamma_k_mpi);
            }
        }
        memory->deallocate(eval);
        memory->deallocate(evec);

    } else if (calc_fstate_omega) {

        double omega0;
        double **gamma_final;
        double *freq_array;
        int ienergy;

        memory->allocate(gamma_final, NT, dos->n_energy);
        memory->allocate(freq_array, dos->n_energy);

        for (i = 0; i < dos->n_energy; ++i) {
            freq_array[i] = dos->energy_dos[i] * time_ry / Hz_to_kayser;
        }

        if (mympi->my_rank == 0) {
            file_mode_tau = input->job_title + ".fstate_omega";

            ofs_mode_tau.open(file_mode_tau.c_str(), std::ios::out);
            if (!ofs_mode_tau) error->exit("compute_mode_tau", "Cannot open file file_mode_tau");

            ofs_mode_tau << "## Frequency-resolved final state amplitude for given modes" << std::endl;
            ofs_mode_tau << "## Gamma[omega][temperature] in cm^-1";
            ofs_mode_tau << std::endl;

            ofs_mode_tau << "## ";
            for (i = 0; i < NT; ++i) {
                ofs_mode_tau << std::setw(10) << T_arr[i];
            }
            ofs_mode_tau << std::endl;
        }

        for (i = 0; i < kslist.size(); ++i) {
            knum = kslist[i] / ns;
            snum = kslist[i] % ns;

            omega0 = dynamical->eval_phonon[knum][snum];

            if (mympi->my_rank == 0) {
                ofs_mode_tau << "# xk = ";

                for (j = 0; j < 3; ++j) {
                    ofs_mode_tau << std::setw(15) << kpoint->xk[knum][j];
                }
                ofs_mode_tau << std::endl;
                ofs_mode_tau << "# mode = " << snum << std::endl;
                ofs_mode_tau << "# Frequency = " << writes->in_kayser(omega0) << std::endl;
            }
            calc_frequency_resolved_final_state2(NT, T_arr, omega0, dos->n_energy, freq_array, kpoint->kmap_to_irreducible[knum], snum, gamma_final);
            if (mympi->my_rank == 0) {

                for (ienergy = 0; ienergy < dos->n_energy; ++ienergy) {
                    omega = dos->energy_dos[ienergy];
                    // calc_frequency_resolved_final_state(NT, T_arr, omega0, omega * time_ry / Hz_to_kayser, knum, snum, gamma_final);

                    ofs_mode_tau << std::setw(10) << omega;
                    for (j = 0; j < NT; ++j) ofs_mode_tau << std::setw(15) << writes->in_kayser(gamma_final[j][ienergy]);
                    ofs_mode_tau << std::endl;
                }
            }
        }

        memory->deallocate(freq_array);
        memory->deallocate(gamma_final);

    } else if (atom_project_mode) {

        /* Atom projection mode. Same as above except that the self-energy is projected on each atomic elements.
        calc_realpart is not used here.  */

        unsigned int natmin = system->natmin;
        int iat, jat;
        double ***damp3_atom, ***damp3_atom_g;
        double damp_sum;

        if (mympi->my_rank == 0) {
            file_mode_tau = input->job_title + ".mode_tau_atom";

            ofs_mode_tau.open(file_mode_tau.c_str(), std::ios::out);
            if (!ofs_mode_tau) error->exit("compute_mode_tau", "Cannot open file file_mode_tau");
            ofs_mode_tau << "## Temperature dependence of atom-projected Gamma for given mode" << std::endl;
            ofs_mode_tau << "## T[K], Gamma3 (cm^-1) (total, atomproj[i][j], i,j = 1, natmin)" << std::endl;
        }

        memory->allocate(damp3_atom, NT, natmin, natmin);
        memory->allocate(damp3_atom_g, NT, natmin, natmin);

        for (i = 0; i < kslist.size(); ++i) {

            knum = kslist[i] / ns;
            snum = kslist[i] % ns;

            omega = dynamical->eval_phonon[knum][snum];


            if (mympi->my_rank == 0) {
                ofs_mode_tau << "# xk = ";

                for (j = 0; j < 3; ++j) {
                    ofs_mode_tau << std::setw(15) << kpoint->xk[knum][j];
                }
                ofs_mode_tau << std::endl;
                ofs_mode_tau << "# mode = " << snum << std::endl;
                ofs_mode_tau << "# Frequency = " << writes->in_kayser(omega) << std::endl;
            }

            if (ksum_mode == -1) {

                std::cout << "myrank = " << mympi->my_rank << std::endl;

                memory->allocate(damp3_atom, natmin, natmin, NT);

                calc_damping_tetra_atom(NT, T_arr, omega, knum, snum, damp3_atom);

                if (mympi->my_rank == 0) {
                    for (j = 0; j < NT; ++j) {
                        ofs_mode_tau << std::setw(10) << T_arr[j];

                        damp_sum = 0.0;

                        for (iat = 0; iat < natmin; ++iat) {
                            for (jat = 0; jat < natmin; ++jat) {
                                damp_sum += damp3_atom[iat][jat][j];
                            }
                        }

                        ofs_mode_tau << std::setw(15) << writes->in_kayser(damp_sum);

                        for (iat = 0; iat < natmin; ++iat) {
                            for (jat = 0; jat < natmin; ++jat) {
                                ofs_mode_tau << std::setw(15) << writes->in_kayser(damp3_atom[iat][jat][j]);
                            }
                        }
                        ofs_mode_tau << std::endl; 
                    }
                }
                memory->deallocate(damp3_atom);

            } else {

                calc_damping_atom(NT, T_arr, omega, knum, snum, damp3_atom);
                MPI_Reduce(&damp3_atom[0][0][0], &damp3_atom_g[0][0][0], NT*natmin*natmin, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

                if (mympi->my_rank == 0) {
                    for (j = 0; j < NT; ++j) {
                        ofs_mode_tau << std::setw(10) << T_arr[j];

                        damp_sum = 0.0;

                        for (iat = 0; iat < natmin; ++iat) {
                            for (jat = 0; jat < natmin; ++jat) {
                                damp_sum += damp3_atom_g[j][iat][jat];
                            }
                        }

                        ofs_mode_tau << std::setw(15) << writes->in_kayser(damp_sum);

                        for (iat = 0; iat < natmin; ++iat) {
                            for (jat = 0; jat < natmin; ++jat) {
                                ofs_mode_tau << std::setw(15) << writes->in_kayser(damp3_atom_g[j][iat][jat]);
                            }
                        }
                        ofs_mode_tau << std::endl; 
                    }
                }

            }
        }
        if (mympi->my_rank == 0) ofs_mode_tau.close();
        memory->deallocate(damp3_atom);
        memory->deallocate(damp3_atom_g);

    } else {
        if (mympi->my_rank == 0) {
            file_mode_tau = input->job_title + ".mode_tau";

            ofs_mode_tau.open(file_mode_tau.c_str(), std::ios::out);
            if (!ofs_mode_tau) error->exit("compute_mode_tau", "Cannot open file file_mode_tau");
        }

        if (calc_realpart) {

            /* Calculate both real and imaginary part of self-energy.
            If quartic_mode == true, then the frequency shift from O(H_{4}) is also computed. */

            std::complex<double> *self_a, *self_b;
            double omega_shift;

            if (mympi->my_rank == 0) {
                ofs_mode_tau << "## Temperature dependence of self-energies of given mode" << std::endl;
                ofs_mode_tau << "## T[K], Gamma3 (cm^-1), Shift3 (cm^-1)";
                if (quartic_mode) ofs_mode_tau << ", Shift4 (cm^-1) <-- linear term in lambda";
                ofs_mode_tau << ", Shifted frequency (cm^-1)";
                ofs_mode_tau << std::endl;
            }

            memory->allocate(self_a, NT);

            if (quartic_mode) {
                memory->allocate(self_b, NT);
            }

            for (i = 0; i < kslist.size(); ++i) {
                knum = kslist[i] / ns;
                snum = kslist[i] % ns;

                omega = dynamical->eval_phonon[knum][snum];

                if (mympi->my_rank == 0) {
                    ofs_mode_tau << "# xk = ";

                    for (j = 0; j < 3; ++j) {
                        ofs_mode_tau << std::setw(15) << kpoint->xk[knum][j];
                    }
                    ofs_mode_tau << std::endl;
                    ofs_mode_tau << "# mode = " << snum << std::endl;
                    ofs_mode_tau << "# Frequency = " << writes->in_kayser(omega) << std::endl;
                }

                selfenergy->selfenergy_a(NT, T_arr, omega, knum, snum, self_a);

                if (quartic_mode) {
                    selfenergy->selfenergy_b(NT, T_arr, omega, knum, snum, self_b);
                }

                if (mympi->my_rank == 0) {
                    for (j = 0; j < NT; ++j) {
                        ofs_mode_tau << std::setw(10) << T_arr[j] << std::setw(15) << writes->in_kayser(self_a[j].imag());
                        ofs_mode_tau << std::setw(15) << writes->in_kayser(-self_a[j].real());

                        omega_shift = omega - self_a[j].real();

                        if (quartic_mode) { 
                            ofs_mode_tau << std::setw(15) << writes->in_kayser(-self_b[j].real());
                            omega_shift -= self_b[j].real();
                        }
                        ofs_mode_tau << std::setw(15) << writes->in_kayser(omega_shift);
                        ofs_mode_tau << std::endl; 
                    }
                }
            }

            memory->deallocate(self_a);
            if(quartic_mode) memory->deallocate(self_b);

        } else  {

            double *damp3, *damp4;
            std::complex<double> *self_a, *self_b, *self_c, *self_d, *self_e;
            std::complex<double> *self_f, *self_g, *self_h, *self_i, *self_j;

            /* Calculate the imaginary part of self-energy. 
            If quartic_mode == true, self-energy of O(H_{4}^{2}) is also calculated. */


            if (mympi->my_rank == 0) {
                ofs_mode_tau << "## Temperature dependence of Gamma for given mode" << std::endl;
                ofs_mode_tau << "## T[K], Gamma3 (cm^-1)";
                if(quartic_mode) ofs_mode_tau << ", Gamma4(cm^-1) <-- specific diagram only";
                ofs_mode_tau << std::endl;
            }

            memory->allocate(damp3, NT);
            memory->allocate(self_a, NT);
            if (quartic_mode) {
                memory->allocate(damp4, NT);
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
                    ofs_mode_tau << "# xk = ";

                    for (j = 0; j < 3; ++j) {
                        ofs_mode_tau << std::setw(15) << kpoint->xk[knum][j];
                    }
                    ofs_mode_tau << std::endl;
                    ofs_mode_tau << "# mode = " << snum << std::endl;
                    ofs_mode_tau << "# Frequency = " << writes->in_kayser(omega) << std::endl;
                }

                if (ksum_mode == -1) {
                    calc_damping_tetra(NT, T_arr, omega, knum, snum, damp3);
                } else {
                    selfenergy->selfenergy_a(NT, T_arr, omega, knum, snum, self_a);
                }

                if (quartic_mode) {
                    if (ksum_mode == -1) {
                        error->exit("compute_mode_tau", "ISMEAR = -1 is not supported for QUARTIC = 1");
                    } else {
                        //					calc_damping4(NT, T_arr, omega, knum, snum, damp4);
                        selfenergy->selfenergy_c(NT, T_arr, omega, knum, snum, self_c);
                        selfenergy->selfenergy_d(NT, T_arr, omega, knum, snum, self_d);
                        selfenergy->selfenergy_e(NT, T_arr, omega, knum, snum, self_e);
                        selfenergy->selfenergy_f(NT, T_arr, omega, knum, snum, self_f);
                        selfenergy->selfenergy_g(NT, T_arr, omega, knum, snum, self_g);
                        selfenergy->selfenergy_h(NT, T_arr, omega, knum, snum, self_h);
                        selfenergy->selfenergy_i(NT, T_arr, omega, knum, snum, self_i);
                        selfenergy->selfenergy_j(NT, T_arr, omega, knum, snum, self_j);
                    }
                }

                if (mympi->my_rank == 0) {
                    for (j = 0; j < NT; ++j) {
                        if (ksum_mode == -1) {
                            ofs_mode_tau << std::setw(10) << T_arr[j] << std::setw(15) << writes->in_kayser(damp3[j]);
                        } else {
                            ofs_mode_tau << std::setw(10) << T_arr[j] << std::setw(15) << writes->in_kayser(self_a[j].imag());
                        }

                        if (quartic_mode) {
                            //							ofs_mode_tau << std::setw(15) << writes->in_kayser(damp4[j]);
                            ofs_mode_tau << std::setw(15) << writes->in_kayser(self_c[j].imag());
                            ofs_mode_tau << std::setw(15) << writes->in_kayser(self_d[j].imag());
                            ofs_mode_tau << std::setw(15) << writes->in_kayser(self_e[j].imag());
                            ofs_mode_tau << std::setw(15) << writes->in_kayser(self_f[j].imag());
                            ofs_mode_tau << std::setw(15) << writes->in_kayser(self_g[j].imag());
                            ofs_mode_tau << std::setw(15) << writes->in_kayser(self_h[j].imag());
                            ofs_mode_tau << std::setw(15) << writes->in_kayser(self_i[j].imag());
                            ofs_mode_tau << std::setw(15) << writes->in_kayser(self_j[j].imag());
                        }

                        ofs_mode_tau << std::endl; 
                    }
                }
            }

            memory->deallocate(damp3);
            memory->deallocate(self_a);

            if (quartic_mode) {
                memory->deallocate(damp4);
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

        if (mympi->my_rank == 0) ofs_mode_tau.close();
    }

    memory->deallocate(T_arr);
}

void Relaxation::gensym_kpairs() {

    int i, j, k;
    int k1, k2, k3, k0;
    int l1, l2, l3;
    double xk_tmp[3];

    unsigned int arr[3];

    std::complex<double> V3tmp1, V3tmp2;

    for (i = 0; i < kpoint->nk_reduced; ++i) {
        k0 = kpoint->k_reduced[i][0];
        k1 = kpoint->knum_minus[k0];

        std::cout << "#Irred. K " << std::setw(5) << i + 1;
        std::cout << " #Number of Unique Triplets : " << std::setw(5) << pair_uniq[i].size() << std::endl;

        //         for (j = 0; j < pair_uniq[i].size(); ++j) {
        //             std::cout << "group " << std::setw(5) << j + 1 << std::endl;
        // 
        //             for (k = 0; k < pair_uniq[i][j].group.size(); ++k) {
        //                 k2 = pair_uniq[i][j].group[k].ks[0];
        //                 k3 = pair_uniq[i][j].group[k].ks[1];
        // 
        //                 std::cout << std::setw(5) << k1 + 1;
        //                 std::cout << std::setw(5) << k2 + 1;
        //                 std::cout << std::setw(5) << k3 + 1;
        // 
        // 
        //                 arr[0] = ns * k1 + 3;
        //                 arr[1] = ns * k2 + 3;
        //                 arr[2] = ns * k3 + 3;
        //                 V3tmp1 = V3(arr);
        // 
        //                 // 				arr[0] = ns * kpoint->knum_minus[k1] + 3;
        //                 // 				arr[1] = ns * kpoint->knum_minus[k2] + 3;
        //                 // 				arr[2] = ns * kpoint->knum_minus[k3] + 3;
        //                 // 				V3tmp2 = V3(arr) * V3tmp1;
        // 
        //                 std::cout << std::setw(15) << std::norm(V3tmp1);
        //                 std::cout << std::setw(5) << pair_uniq[i][j].group[k].symnum;
        //                 std::cout << std::setw(3) << is_proper(pair_uniq[i][j].group[k].symnum);
        //                 std::cout << std::setw(3) << is_symmorphic(pair_uniq[i][j].group[k].symnum) << std::endl;
        //                 // 								std::cout << std::setw(15) << V3tmp2.real();
        //                 // 								std::cout << std::setw(15) << V3tmp2.imag() << std::endl;
        //             }
        //         }
        //         std::cout << std::endl;
    }
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

void Relaxation::gen_pair_uniq()
{
    int i, j;
    int ik1, ik2;
    int isym;
    int nks = kpoint->nk_reduced * ns;
    int knum, knum_minus;
    int ksym;

    double xk[3], xk1[3], xk2[3];
    int ks_in[2];

    std::vector<int> *pointgroup;
    std::vector<int> ks;

    std::set<KsList> visited;
    std::vector<KsList> kslist;

    memory->allocate(pointgroup, kpoint->nk_reduced);

    for (i = 0; i < kpoint->nk_reduced; ++i) {

        pointgroup[i].clear();

        knum = kpoint->k_reduced[i][0];
        knum_minus = kpoint->knum_minus[knum];

        for (isym = 0; isym < symmetry->nsym; ++isym) {

            ksym = knum_sym(knum_minus, isym);

            if (ksym == knum_minus && is_proper(isym)) {
                pointgroup[i].push_back(isym);
            }
        }
    }

    memory->allocate(pair_uniq, kpoint->nk_reduced);

    for (i = 0; i < kpoint->nk_reduced; ++i) {

        knum = kpoint->k_reduced[i][0];
        knum_minus = kpoint->knum_minus[knum];

        for (j = 0; j < 3; ++j) xk[j] = kpoint->xk[knum_minus][j];

        pair_uniq[i].clear();

        for (ik1 = 0; ik1 < nk; ++ik1) {
            for (j = 0; j < 3; ++j) xk1[j] = kpoint->xk[ik1][j];
            for (j = 0; j < 3; ++j) xk2[j] = -xk[j] - xk1[j];
            ik2 = kpoint->get_knum(xk2[0], xk2[1], xk2[2]);

            kslist.clear();

            for (isym = 0; isym < pointgroup[i].size(); ++isym) {
                ks_in[0] = knum_sym(ik1, pointgroup[i][isym]);
                ks_in[1] = knum_sym(ik2, pointgroup[i][isym]);

                if (visited.find(KsList(2, ks_in, isym)) == visited.end()) {
                    visited.insert(KsList(2, ks_in, isym));
                    kslist.push_back(KsList(2, ks_in, pointgroup[i][isym]));
                }
            }

            if (kslist.size() > 0) {
                pair_uniq[i].push_back(kslist);
            }
        }
    }
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

        knum = kpoint->k_reduced[i][0];

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
