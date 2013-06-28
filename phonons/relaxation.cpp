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
#include "symmetry_core.h"
#include "system.h"
#include "write_phonons.h"
#include "../alm_c++/constants.h"

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

	memory->allocate(V, 1);

	unsigned int i, j, k;

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
		memory->allocate(vec_s, system->ntran, 3);

		for (i = 0; i < system->ntran; ++i){
			for (j = 0; j < 3; ++j){
				vec_s[i][j] = system->xr_s[system->map_p2s[0][i]][j];
			}
		}

		for (i = 0; i < system->nat; ++i){
			for (j = 0; j < system->nat; ++j){

				for (k = 0; k < 3; ++k){

					if (system->cell_dimension[k] == 1) {

						vec[k] = system->xr_s[i][k] - system->xr_s[j][k];

						if (std::abs(vec[k]) < 0.5) {
							vec[k] = 0.0;
						} else {
							if (system->xr_s[i][k] < 0.5) {
								vec[k] = 1.0;
							} else {
								vec[k] = -1.0;
							}
						}

					} else if (system->cell_dimension[k] == 2) {

						vec[k] = system->xr_s[system->map_p2s[0][system->map_s2p[i].tran_num]][k] 
						- system->xr_s[system->map_p2s[0][system->map_s2p[j].tran_num]][k];  
						vec[k] = dynamical->fold(vec[k]);
						if (std::abs(system->xr_s[i][k] - system->xr_s[j][k]) > 0.5) vec[k] *= -1.0;

					} else {

						//                          vec[k] = system->xr_s[system->map_p2s[0][system->map_s2p[i].tran_num]][k] 
						//                          - system->xr_s[system->map_p2s[0][system->map_s2p[j].tran_num]][k];  
						// 
						vec[k] = system->xr_s[i][k] - system->xr_s[j][k];
						vec[k] = dynamical->fold(vec[k]);

						vec[k] += system->xr_s[system->map_p2s[system->map_s2p[j].atom_num][0]][k]
						-system->xr_s[system->map_p2s[system->map_s2p[i].atom_num][0]][k];
					}

					relvec[i][j][k] = vec[k];
				}
				system->rotvec(relvec[i][j], relvec[i][j], mat_convert);
			}
		}
		memory->deallocate(vec_s);
	}
	MPI_Bcast(&relvec[0][0][0], 3*system->nat*system->nat, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	for (i = 0; i < system->natmin; ++i){
		invsqrt_mass_p[i] = std::sqrt(1.0 / system->mass[system->map_p2s[i][0]]);
	}

	memory->allocate(vec_for_v3, fcs_phonon->force_constant[1].size(), 3, 2);
	memory->allocate(invmass_for_v3, fcs_phonon->force_constant[1].size());

	j = 0;
	unsigned int atom_num[3];

	for (std::vector<FcsClass>::const_iterator it = fcs_phonon->force_constant[1].begin(); it != fcs_phonon->force_constant[1].end(); ++it) {

		for (i = 0; i < 3; ++i) atom_num[i] = system->map_p2s[(*it).elems[i].atom][(*it).elems[i].cell];

		for (i = 0; i < 3; ++i) {
			vec_for_v3[j][i][0] = relvec[atom_num[1]][atom_num[0]][i];
			vec_for_v3[j][i][1] = relvec[atom_num[2]][atom_num[0]][i];
		}

		invmass_for_v3[j] =  invsqrt_mass_p[(*it).elems[0].atom] * invsqrt_mass_p[(*it).elems[1].atom] * invsqrt_mass_p[(*it).elems[2].atom];
		++j;     
	}

	memory->allocate(evec_index, fcs_phonon->force_constant[1].size(), 3);

	for (i = 0; i < fcs_phonon->force_constant[1].size(); ++i) {
		evec_index[i][0] = 3 * fcs_phonon->force_constant[1][i].elems[0].atom + fcs_phonon->force_constant[1][i].elems[0].xyz;
		evec_index[i][1] = 3 * fcs_phonon->force_constant[1][i].elems[1].atom + fcs_phonon->force_constant[1][i].elems[1].xyz;
		evec_index[i][2] = 3 * fcs_phonon->force_constant[1][i].elems[2].atom + fcs_phonon->force_constant[1][i].elems[2].xyz;
	}


	// Tentative modification for tuning
	// This requires large amount of RAM.
	// 	memory->allocate(cexp_phase, nk, fcs_phonon->force_constant[1].size(), 2);
	// 
	// 	for (i = 0; i < nk; ++i){
	// 		for (j = 0; j < fcs_phonon->force_constant[1].size(); ++j){
	// 			cexp_phase[i][j][0] =  std::exp(im * (vec_for_v3[j][0][0] * kpoint->xk[i][0] + vec_for_v3[j][1][0] * kpoint->xk[i][1] + vec_for_v3[j][2][0] * kpoint->xk[i][2]));
	// 			cexp_phase[i][j][1] =  std::exp(im * (vec_for_v3[j][0][1] * kpoint->xk[i][0] + vec_for_v3[j][1][1] * kpoint->xk[i][1] + vec_for_v3[j][2][1] * kpoint->xk[i][2]));
	// 		}
	// 	}

	// 	memory->allocate(cexp_phase2, nk, system->natmin);
	// 
	// 	for (i = 0; i < nk; ++i){
	// 		for (j = 0; j < system->natmin; ++j){
	// 			double xtmp[3];
	// 			for (unsigned int m = 0; m < 3; ++m){
	// 				xtmp[m] = system->xr_s[system->map_p2s[j][0]][m];
	// 			}
	// 			system->rotvec(xtmp, xtmp, mat_convert);
	// 
	// 			cexp_phase2[i][j] = std::exp(im * (xtmp[0] * kpoint->xk[i][0] + xtmp[1] * kpoint->xk[i][1] + xtmp[2] * kpoint->xk[i][2]));
	// 		}
	// 	}


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

		for (std::vector<FcsClass>::const_iterator it = fcs_phonon->force_constant[2].begin(); it != fcs_phonon->force_constant[2].end(); ++it) {

			for (i = 0; i < 4; ++i) atom_num4[i] = system->map_p2s[(*it).elems[i].atom][(*it).elems[i].cell];

			for (i = 0; i < 3; ++i) {
				vec_for_v4[j][i][0] = relvec[atom_num4[1]][atom_num4[0]][i];
				vec_for_v4[j][i][1] = relvec[atom_num4[2]][atom_num4[0]][i];
				vec_for_v4[j][i][2] = relvec[atom_num4[3]][atom_num4[0]][i];
			}

			invmass_for_v4[j] =  invsqrt_mass_p[(*it).elems[0].atom] * invsqrt_mass_p[(*it).elems[1].atom] * invsqrt_mass_p[(*it).elems[2].atom] * invsqrt_mass_p[(*it).elems[3].atom];
			++j;     
		}

		memory->allocate(evec_index4, fcs_phonon->force_constant[2].size(), 4);

		for (i = 0; i < fcs_phonon->force_constant[2].size(); ++i) {
			evec_index4[i][0] = 3 * fcs_phonon->force_constant[2][i].elems[0].atom + fcs_phonon->force_constant[2][i].elems[0].xyz;
			evec_index4[i][1] = 3 * fcs_phonon->force_constant[2][i].elems[1].atom + fcs_phonon->force_constant[2][i].elems[1].xyz;
			evec_index4[i][2] = 3 * fcs_phonon->force_constant[2][i].elems[2].atom + fcs_phonon->force_constant[2][i].elems[2].xyz;
			evec_index4[i][3] = 3 * fcs_phonon->force_constant[2][i].elems[3].atom + fcs_phonon->force_constant[2][i].elems[3].xyz;
		}
	}

	MPI_Bcast(&ksum_mode, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&calc_realpart, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD);
	MPI_Bcast(&atom_project_mode, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD);

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
			system->rotvec(xk_tmp2, xk_tmp, system->rlavec_p, 'T');
			dist_k_min = std::sqrt(xk_tmp2[0]*xk_tmp2[0] + xk_tmp2[1]*xk_tmp2[1] + xk_tmp2[2]*xk_tmp2[2]);

			for (ik = 1; ik < nk; ++ik) {
				for (int j = 0; j < 3; ++j) {
					xk_tmp[j] = kpoint->xk[ik][j];
				}
				system->rotvec(xk_tmp2, xk_tmp, system->rlavec_p, 'T');

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

#ifdef _DEBUG

	// Check if e_{ks}^{*} = e_{-ks} is satisfied.
	// Note that the violation of the above equation does not necessary 
	// mean bugs.

	int ik, is, js;
	int nk_inv;
	double res;

	for (ik = 0; ik < nk; ++ik) {

		nk_inv = kpoint->knum_minus[ik];

		res = 0.0;

		for (is = 0; is < ns; ++is) {
			for (js = 0; js < ns; ++js) {
				res += std::norm(dynamical->evec_phonon[ik][is][js] - std::conj(dynamical->evec_phonon[nk_inv][is][js]));
			}
		}

		if (std::sqrt(res)/static_cast<double>(std::pow(ns, 2)) > eps12) {
			std::cout << "ik = " << ik << std::endl;
			error->warn("setup_relaxation", "e_{ks}^{*} = e_{-ks} is not satisfied.");
			std::cout << "The above message doesn't necessary mean a problem." << std::endl;
		}
	}

#endif

	modify_eigenvectors();

	epsilon *= time_ry / Hz_to_kayser;

	MPI_Bcast(&epsilon, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

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

			kslist.clear();
			for (i = 0; i < nlist; ++i) {

				ifs_ks >> ktmp[0] >> ktmp[1] >> ktmp[2] >> snum_tmp;
				knum_tmp = kpoint->get_knum(ktmp[0], ktmp[1], ktmp[2]);

				if (knum_tmp == -1) error->exit("setup_relaxation", "Given kpoint is not exist in given k-point grid.");
				kslist.push_back(knum_tmp * dynamical->neval + snum_tmp);
			}
			std::cout << "The number of entries = " << kslist.size() << std::endl;

			ks_analyze_mode = true;

		} else {

			ks_analyze_mode = false;

		}
	}
	MPI_Bcast(&ks_analyze_mode, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD);

	unsigned int *kslist_arr;
	unsigned int nlist;

	nlist = kslist.size();

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

void Relaxation::finish_relaxation()
{
	V[0].clear();
	memory->deallocate(V);

	memory->deallocate(relvec);
	memory->deallocate(invsqrt_mass_p);
	memory->deallocate(vec_for_v3);
	memory->deallocate(invmass_for_v3);
	memory->deallocate(evec_index);

	if (ksum_mode == -1) {
		memory->deallocate(e_tmp);
		memory->deallocate(f_tmp);
	}
}

void Relaxation::calc_ReciprocalV()
{
	// Calculate V_{ks, k's', k''s''}^(3) for Self-energy calculation

	int ik;
	int k1, k2, k3;
	unsigned int i;
	unsigned int b1, b2, b3;
	unsigned int nkp;

	double xk_tmp[3], xk_norm;

	bool is_necessary;

	StructKS ks_tmp;
	std::vector<StructKS> kslist;

	std::cout << std::endl;
	std::cout << "Calculating force constants in reciprocal space .." << std::endl;

	nkp = 0;

	kslist.clear();

	for (k1 = 0; k1 < nk; ++k1){
		for (k2 = k1; k2 < nk; ++k2){
			for (k3 = k2; k3 < nk; ++k3){

				is_necessary 
					= (kpoint->kpset_uniq.find(kpoint->knum_minus[k1]) != kpoint->kpset_uniq.end()) 
					|| (kpoint->kpset_uniq.find(kpoint->knum_minus[k2]) != kpoint->kpset_uniq.end())
					|| (kpoint->kpset_uniq.find(kpoint->knum_minus[k3]) != kpoint->kpset_uniq.end());

				if(!is_necessary) continue;

				xk_tmp[0] = kpoint->xk[k1][0] + kpoint->xk[k2][0] + kpoint->xk[k3][0];
				xk_tmp[1] = kpoint->xk[k1][1] + kpoint->xk[k2][1] + kpoint->xk[k3][1];
				xk_tmp[2] = kpoint->xk[k1][2] + kpoint->xk[k2][2] + kpoint->xk[k3][2];

				for (i = 0; i < 3; ++i){
					xk_tmp[i] = std::fmod(xk_tmp[i], 1.0);
				}
				xk_norm = std::pow(xk_tmp[0], 2) + std::pow(xk_tmp[1], 2) + std::pow(xk_tmp[2], 2);

				// If the momentum-conservation is not satisfied, we skip the loop.

				if (std::sqrt(xk_norm) > eps15) continue;

				for (b1 = 0; b1 < ns; ++b1){
					for (b2 = 0; b2 < ns; ++b2){
						for (b3 = 0; b3 < ns; ++b3){

							ks_tmp.ks1 = ns * k1 + b1;
							ks_tmp.ks2 = ns * k2 + b2;
							ks_tmp.ks3 = ns * k3 + b3;

							if (ks_tmp.ks1 > ks_tmp.ks2 || ks_tmp.ks2 > ks_tmp.ks3) continue;

							kslist.push_back(ks_tmp);

						}
					}
				}
			}
		}
	}

	nkp = kslist.size();

#pragma omp parallel
	{
		unsigned int ks_arr[3];
		std::complex<double> prod;

#pragma omp for schedule(static)
		for (ik = 0; ik < nkp; ++ik){

			ks_arr[0] = kslist[ik].ks1;
			ks_arr[1] = kslist[ik].ks2;
			ks_arr[2] = kslist[ik].ks3;
			prod = V3(ks_arr[0], ks_arr[1], ks_arr[2]);

			// Add to list
			if (std::abs(prod) > eps) {
#pragma omp critical
				V[0].push_back(ReciprocalVs(prod, ks_arr, 3));
			}
		}
	}

	kslist.clear();
	std::cout << "Done !" << std::endl;
	std::cout << "Number of nonzero V's: " << V[0].size() << std::endl;
}

std::complex<double> Relaxation::V3(const unsigned int ks1, const unsigned int ks2, const unsigned int ks3)
{
	unsigned int i;
	unsigned int k1, k2, k3;
	unsigned int b1, b2, b3;
	int ielem;

	double omega_prod;
	double omega[3];

	//  std::vector<FcsClass>::iterator it;

	double ret_re = 0.0;
	double ret_im = 0.0;

	k1 = ks1 / ns;
	k2 = ks2 / ns;
	k3 = ks3 / ns;

	b1 = ks1 % ns;
	b2 = ks2 % ns;
	b3 = ks3 % ns;

	omega[0] = dynamical->eval_phonon[k1][b1];
	omega[1] = dynamical->eval_phonon[k2][b2];
	omega[2] = dynamical->eval_phonon[k3][b3];
	omega_prod = omega[0] * omega[1] * omega[2];

#pragma omp parallel
	{
		unsigned int atom_num[3];
		double vec1[3], vec2[3];
		double invsqrt_mass_prod, phase; 
		std::complex<double> tmp;
		tmp = std::complex<double>(0.0, 0.0);

#pragma omp parallel for reduction(+: ret_re, ret_im)
		for (ielem = 0; ielem < fcs_phonon->force_constant[1].size(); ++ielem){


			for (i = 0; i < 3; ++i) atom_num[i] = system->map_p2s[fcs_phonon->force_constant[1][ielem].elems[i].atom][fcs_phonon->force_constant[1][ielem].elems[i].cell];

			for (i = 0; i < 3; ++i){
				vec1[i] = relvec[atom_num[1]][atom_num[0]][i];
				vec2[i] = relvec[atom_num[2]][atom_num[0]][i];
			}

			phase = vec1[0] * kpoint->xk[k2][0] + vec1[1] * kpoint->xk[k2][1] + vec1[2] * kpoint->xk[k2][2]
			+ vec2[0] * kpoint->xk[k3][0] + vec2[1] * kpoint->xk[k3][1] + vec2[2] * kpoint->xk[k3][2];

			invsqrt_mass_prod = invsqrt_mass_p[fcs_phonon->force_constant[1][ielem].elems[0].atom]
			* invsqrt_mass_p[fcs_phonon->force_constant[1][ielem].elems[1].atom] 
			* invsqrt_mass_p[fcs_phonon->force_constant[1][ielem].elems[2].atom];

			tmp = fcs_phonon->force_constant[1][ielem].fcs_val * std::exp(im * phase) * invsqrt_mass_prod
				* dynamical->evec_phonon[k1][b1][3 * fcs_phonon->force_constant[1][ielem].elems[0].atom + fcs_phonon->force_constant[1][ielem].elems[0].xyz]
			* dynamical->evec_phonon[k2][b2][3 * fcs_phonon->force_constant[1][ielem].elems[1].atom + fcs_phonon->force_constant[1][ielem].elems[1].xyz]
			* dynamical->evec_phonon[k3][b3][3 * fcs_phonon->force_constant[1][ielem].elems[2].atom + fcs_phonon->force_constant[1][ielem].elems[2].xyz];
			ret_re += tmp.real();
			ret_im += tmp.imag();
		}
	}

	return (ret_re + im * ret_im) / std::sqrt(omega_prod);
}

std::complex<double> Relaxation::V3new(const unsigned int ks[3])
{
	/* 
	This version requires massive RAM to store cexp_phase
	*/

	unsigned int i;
	unsigned int kn[3];
	unsigned int sn[3];

	int ielem;

	double omega[3];

	double ret_re = 0.0;
	double ret_im = 0.0;

	for (i = 0; i < 3; ++i){
		kn[i] = ks[i] / ns;
		sn[i] = ks[i] % ns;
		omega[i] = dynamical->eval_phonon[kn[i]][sn[i]];
	}

#pragma omp parallel 
	{
		std::complex<double> ctmp;
		double phase;

#pragma omp for reduction(+: ret_re, ret_im)
		for (ielem = 0; ielem < fcs_phonon->force_constant[1].size(); ++ielem) {

			phase = (vec_for_v3[ielem][0][0]*kpoint->xk[kn[1]][0] + vec_for_v3[ielem][1][0]*kpoint->xk[kn[1]][1] + vec_for_v3[ielem][2][0]*kpoint->xk[kn[1]][2]
			+vec_for_v3[ielem][0][1]*kpoint->xk[kn[2]][0] + vec_for_v3[ielem][1][1]*kpoint->xk[kn[2]][1] + vec_for_v3[ielem][2][1]*kpoint->xk[kn[2]][2]);

			ctmp = fcs_phonon->force_constant[1][ielem].fcs_val * invmass_for_v3[ielem] * std::exp(im*phase)
				* dynamical->evec_phonon[kn[0]][sn[0]][evec_index[ielem][0]] * dynamical->evec_phonon[kn[1]][sn[1]][evec_index[ielem][1]] * dynamical->evec_phonon[kn[2]][sn[2]][evec_index[ielem][2]];

			ret_re += ctmp.real();
			ret_im += ctmp.imag();
		}
	}

	return (ret_re + im * ret_im) / std::sqrt(omega[0] * omega[1] * omega[2]);
}

std::complex<double> Relaxation::V3new2(const unsigned int ks[3])
{
	unsigned int i;
	unsigned int kn[3];
	unsigned int sn[3];

	int ielem;

	double omega[3];

	double ret_re = 0.0;
	double ret_im = 0.0;

	for (i = 0; i < 3; ++i){
		kn[i] = ks[i] / ns;
		sn[i] = ks[i] % ns;
		omega[i] = dynamical->eval_phonon[kn[i]][sn[i]];
	}

#pragma omp parallel 
	{
		double invsqrt_mass_prod;
		double phase;
		std::complex<double> ctmp;
		double vec1[3], vec2[3];
		unsigned int atom_num[3];

#pragma omp for reduction(+: ret_re, ret_im)
		for (ielem = 0; ielem < fcs_phonon->force_constant[1].size(); ++ielem) {

			for (i = 0; i < 3; ++i) atom_num[i] = system->map_p2s[fcs_phonon->force_constant[1][ielem].elems[i].atom][fcs_phonon->force_constant[1][ielem].elems[i].cell];
			for (i = 0; i < 3; ++i) {
				vec1[i] = relvec[atom_num[1]][atom_num[0]][i];
				vec2[i] = relvec[atom_num[2]][atom_num[0]][i];
			}

			phase = 0.0;
			ctmp = std::complex<double>(1.0, 0.0);
			invsqrt_mass_prod = 1.0;

			for (i = 0; i < 3; ++i){
				phase += vec1[i] * kpoint->xk[kn[1]][i] + vec2[i] * kpoint->xk[kn[2]][i];
				invsqrt_mass_prod *= invsqrt_mass_p[fcs_phonon->force_constant[1][ielem].elems[i].atom];
				ctmp *= dynamical->evec_phonon[kn[i]][sn[i]][3 * fcs_phonon->force_constant[1][ielem].elems[i].atom + fcs_phonon->force_constant[1][ielem].elems[i].xyz];
			}

			ctmp *=  fcs_phonon->force_constant[1][ielem].fcs_val * std::exp(im * phase) * invsqrt_mass_prod;

			ret_re += ctmp.real();
			ret_im += ctmp.imag();
		}
	}

	return (ret_re + im * ret_im) / std::sqrt(omega[0] * omega[1] * omega[2]);
}


std::complex<double> Relaxation::V4(const unsigned int ks[4]) 
{

	unsigned int i;
	unsigned int kn[4];
	unsigned int sn[4];

	int ielem;

	double omega[4];

	double ret_re = 0.0;
	double ret_im = 0.0;

	for (i = 0; i < 4; ++i){
		kn[i] = ks[i] / ns;
		sn[i] = ks[i] % ns;
		omega[i] = dynamical->eval_phonon[kn[i]][sn[i]];
	}

#pragma omp parallel 
	{
		std::complex<double> ctmp;
		double phase;

#pragma omp for reduction(+: ret_re, ret_im)
		for (ielem = 0; ielem < fcs_phonon->force_constant[2].size(); ++ielem) {

			phase = (vec_for_v4[ielem][0][0]*kpoint->xk[kn[1]][0] + vec_for_v4[ielem][1][0]*kpoint->xk[kn[1]][1] + vec_for_v4[ielem][2][0]*kpoint->xk[kn[1]][2]
			+vec_for_v4[ielem][0][1]*kpoint->xk[kn[2]][0] + vec_for_v4[ielem][1][1]*kpoint->xk[kn[2]][1] + vec_for_v4[ielem][2][1]*kpoint->xk[kn[2]][2]
			+vec_for_v4[ielem][0][2]*kpoint->xk[kn[3]][0] + vec_for_v4[ielem][1][2]*kpoint->xk[kn[3]][1] + vec_for_v4[ielem][2][2]*kpoint->xk[kn[3]][2]);

			ctmp = fcs_phonon->force_constant[2][ielem].fcs_val * invmass_for_v4[ielem] * std::exp(im*phase)
				* dynamical->evec_phonon[kn[0]][sn[0]][evec_index4[ielem][0]] * dynamical->evec_phonon[kn[1]][sn[1]][evec_index4[ielem][1]] 
			* dynamical->evec_phonon[kn[2]][sn[2]][evec_index4[ielem][2]] * dynamical->evec_phonon[kn[3]][sn[3]][evec_index4[ielem][3]];


			ret_re += ctmp.real();
			ret_im += ctmp.imag();
		}
	}

	return (ret_re + im * ret_im) / std::sqrt(omega[0] * omega[1] * omega[2] * omega[3]);
}

void Relaxation::calc_selfenergy_V3(const unsigned int N, double *T, const double omega, const unsigned int knum, const unsigned int snum, std::complex<double> *ret)
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

	for (i = 0; i < N; ++i) ret[i] = std::complex<double>(0.0, 0.0);

	arr[0] = ns * kpoint->knum_minus[knum] + snum;

	unsigned int nkx = kpoint->nkx;
	unsigned int nky = kpoint->nky;
	unsigned int nkz = kpoint->nkz;

	int iloc, jloc, kloc;

	for (ik = 0; ik < nk; ++ik) {

		xk_tmp[0] = kpoint->xk[knum][0] - kpoint->xk[ik][0];
		xk_tmp[1] = kpoint->xk[knum][1] - kpoint->xk[ik][1];
		xk_tmp[2] = kpoint->xk[knum][2] - kpoint->xk[ik][2];

		iloc = (kpoint->nint(xk_tmp[0]*static_cast<double>(nkx) + static_cast<double>(2*nkx))) % nkx;
		jloc = (kpoint->nint(xk_tmp[1]*static_cast<double>(nky) + static_cast<double>(2*nky))) % nky;
		kloc = (kpoint->nint(xk_tmp[2]*static_cast<double>(nkz) + static_cast<double>(2*nkz))) % nkz;

		jk = kloc + nkz * jloc + nky * nkz * iloc;

		for (is = 0; is < ns; ++is){
			for (js = 0; js < ns; ++js){

				arr[1] = ns * ik + is;
				arr[2] = ns * jk + js;

				omega_inner[0] = dynamical->eval_phonon[ik][is];
				omega_inner[1] = dynamical->eval_phonon[jk][js];

				v3_tmp = std::norm(V3new(arr));

				for (i = 0; i < N; ++i) {
					T_tmp = T[i];

					if (conductivity->use_classical_Cv == 0) {
						n1 = phonon_thermodynamics->fB(omega_inner[0], T_tmp) + phonon_thermodynamics->fB(omega_inner[1], T_tmp) + 1.0;
						n2 = phonon_thermodynamics->fB(omega_inner[0], T_tmp) - phonon_thermodynamics->fB(omega_inner[1], T_tmp);
					} else if (conductivity->use_classical_Cv == 1) {
						n1 = phonon_thermodynamics->fC(omega_inner[0], T_tmp) + phonon_thermodynamics->fC(omega_inner[1], T_tmp) + 1.0;
						n2 = phonon_thermodynamics->fC(omega_inner[0], T_tmp) - phonon_thermodynamics->fC(omega_inner[1], T_tmp);
					}

					ret[i] += v3_tmp 
						* ( n1 / (omega + omega_inner[0] + omega_inner[1] + im * epsilon)
						- n1 / (omega - omega_inner[0] - omega_inner[1] + im * epsilon) 
						+ n2 / (omega - omega_inner[0] + omega_inner[1] + im * epsilon)
						- n2 / (omega + omega_inner[0] - omega_inner[1] + im * epsilon));

				}
			}
		}
	}

	for (i = 0; i < N; ++i) ret[i] *=  std::pow(0.5, 4) / static_cast<double>(nk);
}

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

		iloc = (kpoint->nint(xk_tmp[0]*static_cast<double>(nkx) + static_cast<double>(2*nkx))) % nkx;
		jloc = (kpoint->nint(xk_tmp[1]*static_cast<double>(nky) + static_cast<double>(2*nky))) % nky;
		kloc = (kpoint->nint(xk_tmp[2]*static_cast<double>(nkz) + static_cast<double>(2*nkz))) % nkz;

		jk = kloc + nkz * jloc + nky * nkz * iloc;

		for (is = 0; is < ns; ++is){
			for (js = 0; js < ns; ++js){

				arr[1] = ns * ik + is;
				arr[2] = ns * jk + js;

				omega_inner[0] = dynamical->eval_phonon[ik][is];
				omega_inner[1] = dynamical->eval_phonon[jk][js];

				v3_tmp = std::norm(V3new(arr));

				for (i = 0; i < N; ++i) {
					T_tmp = T[i];

					if (conductivity->use_classical_Cv == 0) {
						n1 = phonon_thermodynamics->fB(omega_inner[0], T_tmp) + phonon_thermodynamics->fB(omega_inner[1], T_tmp) + 1.0;
						n2 = phonon_thermodynamics->fB(omega_inner[0], T_tmp) - phonon_thermodynamics->fB(omega_inner[1], T_tmp);
					} else if (conductivity->use_classical_Cv == 1) {
						n1 = phonon_thermodynamics->fC(omega_inner[0], T_tmp) + phonon_thermodynamics->fC(omega_inner[1], T_tmp) + 1.0;
						n2 = phonon_thermodynamics->fC(omega_inner[0], T_tmp) - phonon_thermodynamics->fC(omega_inner[1], T_tmp);
					}

					if (ksum_mode == 0) {
						ret[i] += v3_tmp 
							* (- n1 * delta_lorentz(omega + omega_inner[0] + omega_inner[1])
							+ n1 * delta_lorentz(omega - omega_inner[0] - omega_inner[1])
							- n2 * delta_lorentz(omega - omega_inner[0] + omega_inner[1])
							+ n2 * delta_lorentz(omega + omega_inner[0] - omega_inner[1]));
					} else if (ksum_mode == 1) {
						ret[i] += v3_tmp 
							* (- n1 * delta_gauss(omega + omega_inner[0] + omega_inner[1])
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

				iloc = (kpoint->nint(xk_tmp[0]*static_cast<double>(nkx) + static_cast<double>(2*nkx))) % nkx;
				jloc = (kpoint->nint(xk_tmp[1]*static_cast<double>(nky) + static_cast<double>(2*nky))) % nky;
				kloc = (kpoint->nint(xk_tmp[2]*static_cast<double>(nkz) + static_cast<double>(2*nkz))) % nkz;

				jk = kloc + nkz * jloc + nky * nkz * iloc;

				ks_tmp[1] = ik * ns + is;
				ks_tmp[2] = jk * ns + js;

				//    if (kcount == nk) error->exitall("calc_damping_tetra", "The number of element is larger than nk");

				omega_inner[ik][0] = dynamical->eval_phonon[ik][is];
				omega_inner[ik][1] = dynamical->eval_phonon[jk][js];

				v3_tmp[ik] = std::norm(V3new(ks_tmp));

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

void Relaxation::calc_damping4(const unsigned int N, double *T, const double omega,
							   const unsigned int knum, const unsigned int snum, double *ret)
{
	unsigned int i;
	unsigned int ik1, ik2, ik3;
	unsigned int is1, is2, is3;
	unsigned int arr[4];
	double n1, n2, n3;
	double n12, n23, n31;
	double xk_tmp[3];
	double omega_inner[3];
	double v4_tmp, T_tmp;

	int iloc, jloc, kloc;

	for (i = 0; i < N; ++i) ret[i] = 0.0;

	arr[0] = ns * kpoint->knum_minus[knum] + snum;

	unsigned int nkx = kpoint->nkx;
	unsigned int nky = kpoint->nky;
	unsigned int nkz = kpoint->nkz;

	for (ik1 = 0; ik1 < nk; ++ik1) {
		for (ik2 = 0; ik2 < nk; ++ik2) {

			xk_tmp[0] = kpoint->xk[knum][0] - kpoint->xk[ik1][0] - kpoint->xk[ik2][0];
			xk_tmp[1] = kpoint->xk[knum][1] - kpoint->xk[ik1][1] - kpoint->xk[ik2][1];
			xk_tmp[2] = kpoint->xk[knum][2] - kpoint->xk[ik1][2] - kpoint->xk[ik2][2];

			iloc = (kpoint->nint(xk_tmp[0]*static_cast<double>(nkx) + static_cast<double>(2*nkx))) % nkx;
			jloc = (kpoint->nint(xk_tmp[1]*static_cast<double>(nky) + static_cast<double>(2*nky))) % nky;
			kloc = (kpoint->nint(xk_tmp[2]*static_cast<double>(nkz) + static_cast<double>(2*nkz))) % nkz;

			ik3 = kloc + nkz * jloc + nky * nkz * iloc;

			for (is1 = 0; is1 < ns; ++is1){
				for (is2 = 0; is2 < ns; ++is2){
					for (is3 = 0; is3 < ns; ++is3) {


						arr[1] = ns * ik1 + is1;
						arr[2] = ns * ik2 + is2;
						arr[3] = ns * ik3 + is3;

						omega_inner[0] = dynamical->eval_phonon[ik1][is1];
						omega_inner[1] = dynamical->eval_phonon[ik2][is2];
						omega_inner[2] = dynamical->eval_phonon[ik3][is3];

						v4_tmp = std::norm(V4(arr));

						for (i = 0; i < N; ++i) {
							T_tmp = T[i];

							n1 = phonon_thermodynamics->fB(omega_inner[0], T_tmp);
							n2 = phonon_thermodynamics->fB(omega_inner[1], T_tmp);
							n3 = phonon_thermodynamics->fB(omega_inner[2], T_tmp);

							n12 = n1 * n2;
							n23 = n2 * n3;
							n31 = n3 * n1;

							ret[i] += v4_tmp 
								* ((n12 + n23 + n31 + n1 + n2 + n3 + 1.0) * delta_lorentz(omega - omega_inner[0] - omega_inner[1] - omega_inner[2])
								+ (n12 - n23 - n31 - n3) * (delta_lorentz(omega + omega_inner[0] + omega_inner[1] - omega_inner[2]) - delta_lorentz(omega - omega_inner[0] - omega_inner[1] + omega_inner[2]))
								+ (n23 - n12 - n31 - n1) * (delta_lorentz(omega - omega_inner[0] + omega_inner[1] + omega_inner[2]) - delta_lorentz(omega + omega_inner[0] - omega_inner[1] - omega_inner[2]))
								+ (n31 - n12 - n23 - n2) * (delta_lorentz(omega + omega_inner[0] - omega_inner[1] + omega_inner[2]) - delta_lorentz(omega - omega_inner[0] + omega_inner[1] - omega_inner[2])));

						}
					}
				}
			}
		}
	}

	for (i = 0; i < N; ++i) ret[i] *=  pi / (std::pow(static_cast<double>(nk), 2) * 3.0 * std::pow(2.0, 5));
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

		iloc = (kpoint->nint(xk_tmp[0]*static_cast<double>(nkx) + static_cast<double>(2*nkx))) % nkx;
		jloc = (kpoint->nint(xk_tmp[1]*static_cast<double>(nky) + static_cast<double>(2*nky))) % nky;
		kloc = (kpoint->nint(xk_tmp[2]*static_cast<double>(nkz) + static_cast<double>(2*nkz))) % nkz;

		jk = kloc + nkz * jloc + nky * nkz * iloc;

		arr[1] = ns * ik + is;
		arr[2] = ns * jk + js;

		omega_inner[0] = dynamical->eval_phonon[ik][is];
		omega_inner[1] = dynamical->eval_phonon[jk][js];

		v3_tmp = std::norm(V3new(arr));

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

				iloc = (kpoint->nint(xk_tmp[0]*static_cast<double>(nkx) + static_cast<double>(2*nkx))) % nkx;
				jloc = (kpoint->nint(xk_tmp[1]*static_cast<double>(nky) + static_cast<double>(2*nky))) % nky;
				kloc = (kpoint->nint(xk_tmp[2]*static_cast<double>(nkz) + static_cast<double>(2*nkz))) % nkz;

				jk = kloc + nkz * jloc + nky * nkz * iloc;

				ks_tmp[1] = ik * ns + is;
				ks_tmp[2] = jk * ns + js;

				omega_inner[ik][0] = dynamical->eval_phonon[ik][is];
				omega_inner[ik][1] = dynamical->eval_phonon[jk][js];

				v3_tmp[ik] = std::norm(V3new(ks_tmp));

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
	memory->deallocate(v3_tmp_proj);
}

void Relaxation::modify_eigenvectors()
{
	bool *flag_done;
	unsigned int ik;
	unsigned int is, js;
	unsigned int nk_inv;
	std::complex<double> *evec_tmp;

	if (mympi->my_rank == 0) {
		std::cout << "**********      NOTICE      **********" << std::endl;
		std::cout << "For the brevity of the calculation, " << std::endl;
		std::cout << "phonon eigenvectors will be modified" << std::endl;
		std::cout << "so that e_{-ks}^{mu} = (e_{ks}^{mu})^{*}. " << std::endl;
	}

	memory->allocate(flag_done, nk);
	memory->allocate(evec_tmp, ns);



	for (ik = 0; ik < nk; ++ik) flag_done[ik] = false;

	for (ik = 0; ik < nk; ++ik){

		if (!flag_done[ik]) {

			nk_inv = kpoint->knum_minus[ik];   

			for (is = 0; is < ns; ++is){
				for (js = 0; js < ns; ++js){
					evec_tmp[js] = dynamical->evec_phonon[ik][is][js];
				}

				for (js = 0; js < ns; ++js){
					dynamical->evec_phonon[nk_inv][is][js] = std::conj(evec_tmp[js]);
				}
			}

			flag_done[ik] = true;
			flag_done[nk_inv] = true;
		}
	}

	memory->deallocate(flag_done);
	memory->deallocate(evec_tmp);

	MPI_Barrier(MPI_COMM_WORLD);
	if (mympi->my_rank == 0) {
		std::cout << "done !" << std::endl;
		std::cout << "**************************************" << std::endl;
	}
}

double Relaxation::delta_lorentz(const double omega)
{
	return epsilon / (omega*omega + epsilon*epsilon) / pi;
}

double Relaxation::delta_gauss(const double omega)
{
	return std::exp(- omega * omega / (epsilon * epsilon)) / (epsilon * std::sqrt(pi));
}

void Relaxation::calc_selfenergy()
{
	double Tmin = system->Tmin;
	double Tmax = system->Tmax;
	double dT = system->dT;

	unsigned int NT = static_cast<unsigned int>((Tmax - Tmin) / dT);
	unsigned int i;

	double *T_arr, T;

	std::string file_selfenergy;
	std::ofstream ofs_selfenergy;

	std::string file_test;
	std::ofstream ofs_test;

	int knum, snum;
	double omega;
	double *damping;


	double k_tmp[3];

	memory->allocate(T_arr, NT);
	memory->allocate(damping, NT);

	for (i = 0; i < NT; ++i) T_arr[i] = Tmin + dT * static_cast<double>(i);

	if (mympi->my_rank == 0) {
		file_test = input->job_title + ".damp_T";
		ofs_test.open(file_test.c_str(), std::ios::out);
		if(!ofs_test) error->exit("write_selfenergy", "cannot open file_test");

		std::string file_KS;

		file_KS = "KS_INPUT";
		std::ifstream ifs_KS;
		ifs_KS.open(file_KS.c_str(), std::ios::in);
		if (!ifs_KS) error->exit("write_selfenergy", "cannot open KS_INPUT");
		ifs_KS >> k_tmp[0] >> k_tmp[1] >> k_tmp[2];
		ifs_KS >> snum;
		ifs_KS.close();

		std::cout << "Given kpoints: ";
		for (i = 0; i < 3; ++i) {
			std::cout << std::setw(15) << k_tmp[i];
		}
		std::cout << std::endl;
		std::cout << "Given branch: " << snum + 1 << std::endl;

		knum = kpoint->get_knum(k_tmp[0], k_tmp[1], k_tmp[2]);
		if (knum == -1) error->exit("calc_selfenergy", "Corresponding k-point does not exist");

		omega = dynamical->eval_phonon[knum][snum];

		ofs_test << "# Damping function [cm] of a phonon at xk = ";
		for (i = 0; i < 3; ++i) {
			ofs_test << std::setw(15) << kpoint->xk[knum][i];
		}
		ofs_test << std::endl;
		ofs_test << "# Branch = " << snum << std::endl;

		if (relaxation->ksum_mode == 0 || relaxation->ksum_mode == 1) {
			relaxation->calc_damping(NT, T_arr, omega, knum, snum, damping);
		} else if (relaxation->ksum_mode == -1) {
			relaxation->calc_damping_tetra(NT, T_arr, omega, knum, snum, damping);
		}

		for (i = 0; i < NT; ++i) {
			T = Tmin + dT * static_cast<double>(i);
			ofs_test  << std::setw(5) << T;
			ofs_test << std::setw(15) << 2.0 * damping[i]/time_ry*Hz_to_kayser;
			ofs_test << std::setw(15) << time_ry / (2.0 * damping[i]) * 1.0e+12;
			ofs_test << std::endl;
		}
		ofs_test.close();
	}
	error->exitall("hoge", "tomare!");
}


void Relaxation::v3_test() {

	int i;
	unsigned int stmp[3], kstmp[3];
	int nkplus, nkminus;

	nkplus = 1;
	nkminus= kpoint->knum_minus[nkplus];


	stmp[0] = 0;
	stmp[1] = 1;
	stmp[2] = 2;

	for (i = 0; i < 3; ++i) {
		std::cout << std::setw(15) << kpoint->xk[nkplus][i];
	}
	std::cout << std::endl;

	for (i = 0; i < 3; ++i) {
		kstmp[i] = dynamical->neval * nkplus + stmp[i];
	}
	std::cout << V3(kstmp[0], kstmp[1], kstmp[2]) << std::endl;
	std::cout << V3new(kstmp) << std::endl;
	std::cout << V3new2(kstmp) << std::endl;

	for (i = 0; i < 3; ++i) {
		kstmp[i] = dynamical->neval * nkminus + stmp[i];
	}
	for (i = 0; i < 3; ++i) {
		std::cout << std::setw(15) << kpoint->xk[nkminus][i];
	}
	std::cout << std::endl;

	std::cout << V3(kstmp[0], kstmp[1], kstmp[2]) << std::endl;
	std::cout << V3new(kstmp) << std::endl;
	std::cout << V3new2(kstmp) << std::endl;


	// error->exit("v3_test", "finished!");
}


void Relaxation::v4_test() {

	int i;
	unsigned int stmp[4], kstmp[4];
	int nkplus, nkminus;

	nkplus = 2;
	nkminus= kpoint->knum_minus[nkplus];


	stmp[0] = 0;
	stmp[1] = 1;
	stmp[2] = 2;
	stmp[3] = 0;

	for (i = 0; i < 4; ++i) {
		std::cout << std::setw(15) << kpoint->xk[nkplus][i];
	}
	std::cout << std::endl;

	for (i = 0; i < 4; ++i) {
		kstmp[i] = dynamical->neval * nkplus + stmp[i];
	}


	std::cout << V4(kstmp) << std::endl;

	for (i = 0; i < 4; ++i) {
		kstmp[i] = dynamical->neval * nkminus + stmp[i];
	}
	for (i = 0; i < 4; ++i) {
		std::cout << std::setw(15) << kpoint->xk[nkminus][i];
	}
	std::cout << std::endl;

	std::cout << V4(kstmp) << std::endl;


	error->exit("v4_test", "finished!");
}

void Relaxation::compute_mode_tau()
{
	unsigned int i, j;
	double *damp3, *damp4;

	double *T_arr;

	unsigned int NT;

	double Tmax = system->Tmax;
	double Tmin = system->Tmin;
	double dT = system->dT;
	double omega;

	unsigned int knum, snum;

	NT = static_cast<unsigned int>((Tmax - Tmin) / dT);
	memory->allocate(T_arr, NT);

	for (i = 0; i < NT; ++i) T_arr[i] = Tmin + static_cast<double>(i)*dT;

	std::ofstream ofs_mode_tau;
	std::string file_mode_tau;

	if (!atom_project_mode) {

		if (mympi->my_rank == 0) {
			file_mode_tau = input->job_title + ".mode_tau";

			ofs_mode_tau.open(file_mode_tau.c_str(), std::ios::out);
			if (!ofs_mode_tau) error->exit("compute_mode_tau", "Cannot open file file_mode_tau");
		}

		if (calc_realpart) {

			std::complex<double> *self3;
			double *shift4;
			double omega_shift;

			if (mympi->my_rank == 0) {
				ofs_mode_tau << "## Temperature dependence of self-energies of given mode" << std::endl;
				ofs_mode_tau << "## T[K], Gamma3 (cm^-1), Shift3 (cm^-1)";
				if (quartic_mode) ofs_mode_tau << ", Shift4 (cm^-1) <-- linear term in lambda";
				ofs_mode_tau << ", Shifted frequency (cm^-1)";

				ofs_mode_tau << std::endl;
			}

			memory->allocate(self3, NT);

			if (quartic_mode) {
				memory->allocate(shift4, NT);
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

				calc_selfenergy_V3(NT, T_arr, omega, knum, snum, self3);

				if (quartic_mode) {
					calc_realpart_V4(NT, T_arr, omega, knum, snum, shift4);
				}

				if (mympi->my_rank == 0) {
					for (j = 0; j < NT; ++j) {
						ofs_mode_tau << std::setw(10) << T_arr[j] << std::setw(15) << writes->in_kayser(self3[j].imag());
						ofs_mode_tau << std::setw(15) << writes->in_kayser(-self3[j].real());

						omega_shift = omega - self3[j].real();

						if (quartic_mode) { 
							ofs_mode_tau << std::setw(15) << writes->in_kayser(-shift4[j]);
							omega_shift -= shift4[j];
						}
						ofs_mode_tau << std::setw(15) << writes->in_kayser(omega_shift);
						ofs_mode_tau << std::endl; 
					}
				}
			}

			memory->deallocate(self3);
			if(quartic_mode) memory->deallocate(shift4);

		} else {


			if (mympi->my_rank == 0) {
				ofs_mode_tau << "## Temperature dependence of Gamma for given mode" << std::endl;
				ofs_mode_tau << "## T[K], Gamma3 (cm^-1)";
				if(quartic_mode) ofs_mode_tau << ", Gamma4(cm^-1) <-- specific diagram only";
				ofs_mode_tau << std::endl;
			}

			memory->allocate(damp3, NT);
			if (quartic_mode) {
				memory->allocate(damp4, NT);
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
					calc_damping(NT, T_arr, omega, knum, snum, damp3);
				}

				if (quartic_mode) {
					if (ksum_mode == -1) {
						error->exit("compute_mode_tau", "ISMEAR = -1 is not supported for QUARTIC = 1");
					} else {
						calc_damping4(NT, T_arr, omega, knum, snum, damp4);
					}
				}

				if (mympi->my_rank == 0) {
					for (j = 0; j < NT; ++j) {
						ofs_mode_tau << std::setw(10) << T_arr[j] << std::setw(15) << writes->in_kayser(damp3[j]);

						if (quartic_mode) {
							ofs_mode_tau << std::setw(15) << writes->in_kayser(damp4[j]);
						}

						ofs_mode_tau << std::endl; 
					}
				}
			}

			memory->deallocate(damp3);
			if (quartic_mode) memory->deallocate(damp4);
		}

	} else {

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

				memory->allocate(damp3_atom, NT, natmin, natmin);
				memory->allocate(damp3_atom_g, NT, natmin, natmin);

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
				memory->deallocate(damp3_atom);
				memory->deallocate(damp3_atom_g);
			}
		}
	}
	memory->deallocate(T_arr);
}

