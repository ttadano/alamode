#include "mpi_common.h" 
#include "kpoint.h"
#include "memory.h"
#include "error.h"
#include "system.h"
#include "phonon_dos.h"
#include "symmetry_core.h"
#include "dynamical.h"
#include "timer.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <set>
#include <numeric>
#include "parsephon.h"

#ifdef _USE_BOOST
#include <boost/lexical_cast.hpp>
#endif

#ifdef _USE_EIGEN
#include <Eigen/Core>
#include <Eigen/LU>
#endif

using namespace PHON_NS;

Kpoint::Kpoint(PHON *phon): Pointers(phon) {
	npath = 0;
}

Kpoint::~Kpoint() {
	memory->deallocate(xk);

// 	if (kpoint_mode == 1) {
// 		memory->deallocate(kpoint_direction);
// 		memory->deallocate(kaxis);
// 		memory->deallocate(kp_symbol);
// 		memory->deallocate(kp_bound);
// 	}
}

void Kpoint::kpoint_setups()
{
	symmetry->symmetry_flag = true;

	unsigned int i, j;
	std::string str_tmp;

	MPI_Bcast(&kpoint_mode, 1, MPI_INT, 0, MPI_COMM_WORLD);

	if (mympi->my_rank == 0) {

		std::cout << "**k-points**" << std::endl;

		switch (kpoint_mode){

		case 0:
			std::cout << " KPMODE = 0 : calculation on given k-points" << std::endl;

			nk = kpInp.size();
			memory->allocate(xk, nk, 3);
			j = 0;

			for (std::vector<KpointInp>::const_iterator it = kpInp.begin(); it != kpInp.end(); ++it) {
				for (i = 0; i < 3; ++i) {
#ifdef _USE_BOOST
					xk[j][i] = boost::lexical_cast<double>((*it).kpelem[i]);
#else
					xk[j][i] = std::atof(((*it).kpelem[i]).c_str());
#endif 
				}	
				++j;
			}

			std::cout << " Number of k-points: " << nk << std::endl << std::endl;

			break;

		case 1:
			std::cout << " KPMODE = 1: band structure calculation" << std::endl;

			npath = kpInp.size();

			std::cout << " Number of paths: " << npath << std::endl;
			memory->allocate(kp_symbol, npath, 2);
			memory->allocate(kp_bound, npath, 2, 3);
			memory->allocate(nkp, npath);

			nk = 0;

			j = 0;
			for (std::vector<KpointInp>::const_iterator it = kpInp.begin(); it != kpInp.end(); ++it) {
				kp_symbol[j][0] = (*it).kpelem[0];
				kp_symbol[j][1] = (*it).kpelem[4];
#ifdef _USE_BOOST
				for (i = 0; i < 3; ++i) {
					kp_bound[j][0][i] = boost::lexical_cast<double>((*it).kpelem[i + 1]);
					kp_bound[j][1][i] = boost::lexical_cast<double>((*it).kpelem[i + 5]);
				}
				nkp[j] = boost::lexical_cast<int>((*it).kpelem[8]);
#else
				for (i = 0; i < 3; ++i) {
					kp_bound[j][0][i] = std::atof(((*it).kpelem[i + 1]).c_str());
					kp_bound[j][1][i] = std::atof(((*it).kpelem[i + 5]).c_str());

				}
				nkp[j] = std::atoi(((*it).kpelem[8]).c_str());
#endif
				nk += nkp[j];
				++j;
			}

			memory->allocate(xk, nk, 3);
			memory->allocate(kpoint_direction, nk, 3);

			gen_kpoints_band();
			memory->deallocate(kp_bound);

			double **xk_c, tmp[3];
			memory->allocate(kaxis, nk);
			memory->allocate(xk_c, nk, 3);

			for (i = 0; i < nk; ++i){
				system->rotvec(xk_c[i], xk[i], system->rlavec_p, 'T');
			}

			unsigned int j, k;
			unsigned int ik;

			ik = 0;
			std::cout << std::endl << " ---------- kpval at the edges ----------" << std::endl;
			std::cout << " kval";
			std::cout.unsetf(std::ios::scientific);

			for (i = 0; i < npath; ++i){
				for (j = 0; j < nkp[i]; ++j){
					if (j == 0){
						if(ik == 0) {
							kaxis[ik] = 0.0;
						} else {
							kaxis[ik] = kaxis[ik - 1];
						}
						std::cout << std::setw(10) << kaxis[ik];
					} else {
						for (k = 0; k < 3; ++k){
							tmp[k] = xk_c[ik][k] - xk_c[ik - 1][k];
						}
						kaxis[ik] = kaxis[ik - 1] + std::sqrt(tmp[0] * tmp[0] + tmp[1] * tmp[1] + tmp[2] * tmp[2]);
					}
					++ik;
				}
			}

			std::cout << std::setw(10) << kaxis[ik - 1];
			std::cout << std::endl;
			std::cout << " ----------------------------------------" << std::endl;

			memory->deallocate(xk_c);
			memory->deallocate(nkp);

			std::cout.setf(std::ios::scientific);

			break;
		case 2:
			std::cout << " KPMODE = 2: Uniform grid" << std::endl;

#ifdef _USE_BOOST
			nkx = boost::lexical_cast<int>(kpInp[0].kpelem[0]);
			nky = boost::lexical_cast<int>(kpInp[0].kpelem[1]);
			nkz = boost::lexical_cast<int>(kpInp[0].kpelem[2]);
#else 
			nkx = std::atoi((kpInp[0].kpelem[0]).c_str());
			nky = std::atoi((kpInp[0].kpelem[1]).c_str());
			nkz = std::atoi((kpInp[0].kpelem[2]).c_str());	
#endif
			nk = nkx * nky * nkz;
			memory->allocate(xk, nk, 3);

			gen_kmesh(symmetry->symmetry_flag);
			std::cout << " nkx: " << std::setw(6) << nkx;
			std::cout << " nky: " << std::setw(6) << nky;
			std::cout << " nkz: " << std::setw(6) << nkz;
			std::cout << std::endl;
			std::cout << " Number of k-points: " << nk << std::endl << std::endl;

			memory->allocate(knum_minus, nk);
			gen_nkminus();

			std::cout << " Number of irreducible k-points: " << nk_equiv.size() << std::endl << std::endl;
			std::cout << "#k-points (in unit of reciprocal vector), weights" << std::endl;

			ik = 0;
			for (i = 0; i < nk_equiv.size(); ++i){
				std::cout << std::setw(8) << i + 1 << ":";
				for (j = 0; j < 3; ++j){
					std::cout << std::setw(15) << std::scientific << kpIBZ[ik].kval[j];
				}
				std::cout << std::setw(15) << std::fixed << weight_k[i] << std::endl;
				ik += nk_equiv[i];
			}
			std::cout << std::endl;

			std::cout.unsetf(std::ios::fixed);

			break;
		default:
			error->exit("read_kpoints", "This cannot happen.");
		}
	}

	MPI_Bcast(&nkx, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
	MPI_Bcast(&nky, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
	MPI_Bcast(&nkz, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
	MPI_Bcast(&nk , 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

	if (mympi->my_rank > 0) {
		memory->allocate(xk, nk, 3);
	}
	MPI_Bcast(&xk[0][0], 3*nk, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	if (kpoint_mode == 2) {
		if (mympi->my_rank > 0) {
			memory->allocate(knum_minus, nk);
		}
		MPI_Bcast(&knum_minus[0], nk, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
		MPI_Bcast(&nk_reduced, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
		MPI_Bcast(&nequiv_max, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
		if (mympi->my_rank > 0) {
			memory->allocate(nk_equiv_arr, nk_reduced);
			memory->allocate(k_reduced, nk_reduced, nequiv_max);
		}
		MPI_Bcast(&nk_equiv_arr[0], nk_reduced, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
		MPI_Bcast(&k_reduced[0][0], nk_reduced*nequiv_max, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
	}
}

void Kpoint::gen_kpoints_band()
{
	unsigned int i, j, k;
	double xk_s[3], xk_e[3];
	double xk_direction[3], norm;
	unsigned int ik = 0;

	for (i = 0; i < npath; ++i){
		for (j = 0; j < 3; ++j){
			xk_s[j] = kp_bound[i][0][j];
			xk_e[j] = kp_bound[i][1][j];
			xk_direction[j] = xk_e[j] - xk_s[j];
		}

		system->rotvec(xk_direction, xk_direction, system->rlavec_p, 'T');
		norm = std::pow(xk_direction[0], 2) + std::pow(xk_direction[1], 2) + std::pow(xk_direction[2], 2);
		norm = std::sqrt(norm);

		if (norm > eps) {
			for (j = 0; j < 3; ++j) xk_direction[j] /= norm;
		}

		for(j = 0; j < nkp[i]; ++j){
			for(k = 0; k < 3; ++k){
				xk[ik][k] = xk_s[k] + (xk_e[k] - xk_s[k]) * static_cast<double>(j) / static_cast<double>(nkp[i] - 1);
				kpoint_direction[ik][k] = xk_direction[k];
			}
			++ik;
		}
	}
}

void Kpoint::gen_kmesh(bool usesym)
{
	unsigned int ix, iy, iz;
	unsigned int i, ik;
	double **xkr;

	std::cout << "Generating uniform k-point grid ..." << std::endl;

	memory->allocate(xkr, nk, 3);

	for (ix = 0; ix < nkx; ++ix){
		for (iy = 0; iy < nky; ++iy){
			for (iz = 0; iz < nkz; ++iz){
				ik = iz + iy * nkz + ix * nkz * nky;
				xkr[ik][0] = static_cast<double>(ix) / static_cast<double>(nkx);
				xkr[ik][1] = static_cast<double>(iy) / static_cast<double>(nky);
				xkr[ik][2] = static_cast<double>(iz) / static_cast<double>(nkz);
			}
		}
	}

	if (usesym) {
		reduce_kpoints(xkr);
	}

	for (ik = 0; ik < nk; ++ik){
		for (i = 0; i < 3; ++i){
			xk[ik][i] = xkr[ik][i] - static_cast<double>(nint(xkr[ik][i]));
		}
	}
	memory->deallocate(xkr);

	timer->print_elapsed();
}

void Kpoint::gen_nkminus()
{
	unsigned int ik;
	double diff[3];
	std::vector<KpointList> ksets;
	std::vector<KpointList>::iterator it;
	std::vector<double> ktmp;
	bool *found_minus;

	ksets.clear();

	for (ik = 0; ik < nk; ++ik){

		ktmp.clear();

		ktmp.push_back(xk[ik][0]);
		ktmp.push_back(xk[ik][1]);
		ktmp.push_back(xk[ik][2]);

		ksets.push_back(KpointList(ik, ktmp));
	}

	memory->allocate(found_minus, nk);

	for (ik = 0; ik < nk; ++ik) found_minus[ik] = false;

	for (ik = 0; ik < nk; ++ik){

		if (found_minus[ik]) continue;
		ktmp.clear();
		ktmp.push_back(-xk[ik][0]);
		ktmp.push_back(-xk[ik][1]);
		ktmp.push_back(-xk[ik][2]);

		for (it = ksets.begin(); it != ksets.end(); ++it){
			diff[0] = std::fmod(ktmp[0]-(*it).kval[0], 1.0);
			diff[1] = std::fmod(ktmp[1]-(*it).kval[1], 1.0);
			diff[2] = std::fmod(ktmp[2]-(*it).kval[2], 1.0);


			if (std::sqrt(diff[0]*diff[0] + diff[1]*diff[1] + diff[2]*diff[2]) < eps12) {
				knum_minus[ik] = (*it).knum;
				knum_minus[(*it).knum] = ik;
				found_minus[ik] = true;
				found_minus[(*it).knum] = true;
				break;
			} 
		}

		if (!found_minus[ik]) {
			error->exit("gen_nkminus", "k-point of the inverse point is not found");
		}
	}

	memory->deallocate(found_minus);
	ksets.clear();
}

void Kpoint::reduce_kpoints(double **xkr)
{
	unsigned int ik;
	unsigned int i, j;
	unsigned int nsame;
	int iloc, jloc, kloc;
	int nloc;

	unsigned int *kequiv;

	double diff[3];

	std::set<KpointList> ksets;
	std::vector<double> ktmp;

#ifdef _USE_EIGEN
	Eigen::Matrix3d srot;
	Eigen::Vector3d xk_sym, xk_orig;
#else
	double srot[3][3];
	double xk_sym[3], xk_orig[3];
#endif

	std::cout << "Reducing the k-points by using the crystal symmetry ... " << std::endl;

	kpIBZ.clear();
	nk_equiv.clear();
	weight_k.clear();

	memory->allocate(kequiv, nk);

	ksets.clear();

	for (ik = 0; ik < nk; ++ik) kequiv[ik] = ik;

	for (ik = 0; ik < nk; ++ik){

		if(kequiv[ik] != ik) continue;

		nsame = 1;

		ktmp.clear();
		ktmp.push_back(xkr[ik][0] - static_cast<double>(nint(xkr[ik][0])));
		ktmp.push_back(xkr[ik][1] - static_cast<double>(nint(xkr[ik][1])));
		ktmp.push_back(xkr[ik][2] - static_cast<double>(nint(xkr[ik][2])));

		kpIBZ.push_back(KpointList(ik, ktmp));

#ifdef _USE_EIGEN
		for (i = 0; i < 3; ++i) xk_orig(i) = xkr[ik][i];
#else
		for (i = 0; i < 3; ++i) xk_orig[i] = xkr[ik][i];
#endif

		for (std::vector<SymmetryOperation>::iterator isym = symmetry->SymmList.begin(); isym != symmetry->SymmList.end(); ++isym){

			for (i = 0; i < 3; ++i){
				for (j = 0; j < 3; ++j){
#ifdef _USE_EIGEN
					srot(i,j) = static_cast<double>((*isym).symop[3 * i + j]);
#else
					srot[i][j] = static_cast<double>((*isym).symop[3 * i + j]);
#endif
				}
			}
#ifdef _USE_EIGEN
			Eigen::FullPivLU< Eigen::Matrix3d > lu(srot);
			Eigen::Matrix3d S = lu.inverse();

			xk_sym = S.transpose() * xk_orig;

			for (i = 0; i < 3; ++i){
				xk_sym(i) = xk_sym(i) - nint(xk_sym(i));
			}    

			diff[0] = static_cast<double>(nint(xk_sym(0)*nkx)) - xk_sym(0)*nkx;
			diff[1] = static_cast<double>(nint(xk_sym(1)*nky)) - xk_sym(1)*nky;
			diff[2] = static_cast<double>(nint(xk_sym(2)*nkz)) - xk_sym(2)*nkz;

			if(std::abs(std::pow(diff[0], 2) + std::pow(diff[1], 2) + std::pow(diff[2], 2)) < eps12) {

				iloc = (nint(xk_sym(0)*nkx + 2 * nkx)) % nkx;
				jloc = (nint(xk_sym(1)*nky + 2 * nky)) % nky;
				kloc = (nint(xk_sym(2)*nkz + 2 * nkz)) % nkz;

				nloc = kloc + nkz * jloc + nky * nkz * iloc;

				if(nloc > ik && kequiv[nloc] == nloc) {
					kequiv[nloc] = ik;
					ktmp.clear();
					ktmp.push_back(xk_sym(0));
					ktmp.push_back(xk_sym(1));
					ktmp.push_back(xk_sym(2));
					kpIBZ.push_back(KpointList(nloc, ktmp));
					++nsame;
				}
			}

			// Time-reversal symmetry

			if (!symmetry->time_reversal_sym) continue;

			for (i = 0; i < 3; ++i) xk_sym(i) *= -1.0; 

			diff[0] = static_cast<double>(nint(xk_sym(0)*nkx)) - xk_sym(0)*nkx;
			diff[1] = static_cast<double>(nint(xk_sym(1)*nky)) - xk_sym(1)*nky;
			diff[2] = static_cast<double>(nint(xk_sym(2)*nkz)) - xk_sym(2)*nkz;

			if(std::abs(std::pow(diff[0], 2) + std::pow(diff[1], 2) + std::pow(diff[2], 2)) < eps12) {

				iloc = (nint(xk_sym(0)*nkx + 2 * nkx)) % nkx;
				jloc = (nint(xk_sym(1)*nky + 2 * nky)) % nky;
				kloc = (nint(xk_sym(2)*nkz + 2 * nkz)) % nkz;

				nloc = kloc + nkz * jloc + nky * nkz * iloc;

				if(nloc > ik && kequiv[nloc] == nloc) {
					kequiv[nloc] = ik;
					ktmp.clear();
					ktmp.push_back(xk_sym(0));
					ktmp.push_back(xk_sym(1));
					ktmp.push_back(xk_sym(2));
					kpIBZ.push_back(KpointList(nloc, ktmp));
					++nsame;
				}
			}

#else
			double srot_inv[3][3], srot_inv_t[3][3];

			system->invmat3(srot_inv, srot);
			system->transpose3(srot_inv_t, srot_inv);
			system->rotvec(xk_sym, xk_orig, srot_inv_t);


			for (i = 0; i < 3; ++i){
				xk_sym[i] = xk_sym[i] - nint(xk_sym[i]);
			}    

			diff[0] = static_cast<double>(nint(xk_sym[0]*nkx)) - xk_sym[0]*nkx;
			diff[1] = static_cast<double>(nint(xk_sym[1]*nky)) - xk_sym[1]*nky;
			diff[2] = static_cast<double>(nint(xk_sym[2]*nkz)) - xk_sym[2]*nkz;

			if(std::sqrt(diff[0]*diff[0] + diff[1]*diff[1] + diff[2]*diff[2]) < eps12) {

				iloc = (nint(xk_sym[0]*nkx + 2 * nkx)) % nkx;
				jloc = (nint(xk_sym[1]*nky + 2 * nky)) % nky;
				kloc = (nint(xk_sym[2]*nkz + 2 * nkz)) % nkz;

				nloc = kloc + nkz * jloc + nky * nkz * iloc;

				if(nloc > ik && kequiv[nloc] == nloc) {
					kequiv[nloc] = ik;
					ktmp.clear();
					ktmp.push_back(xk_sym[0]);
					ktmp.push_back(xk_sym[1]);
					ktmp.push_back(xk_sym[2]);
					kpIBZ.push_back(KpointList(nloc, ktmp));
					++nsame;
				}
			}

			// Time-reversal symmetry

			if (!symmetry->time_reversal_sym) continue;

			for (i = 0; i < 3; ++i) xk_sym[i] *= -1.0; 

			diff[0] = static_cast<double>(nint(xk_sym[0]*nkx)) - xk_sym[0]*nkx;
			diff[1] = static_cast<double>(nint(xk_sym[1]*nky)) - xk_sym[1]*nky;
			diff[2] = static_cast<double>(nint(xk_sym[2]*nkz)) - xk_sym[2]*nkz;

			if(std::sqrt(diff[0]*diff[0] + diff[1]*diff[1] + diff[2]*diff[2]) < eps12) {

				iloc = (nint(xk_sym[0]*nkx + 2 * nkx)) % nkx;
				jloc = (nint(xk_sym[1]*nky + 2 * nky)) % nky;
				kloc = (nint(xk_sym[2]*nkz + 2 * nkz)) % nkz;

				nloc = kloc + nkz * jloc + nky * nkz * iloc;

				if(nloc > ik && kequiv[nloc] == nloc) {
					kequiv[nloc] = ik;
					ktmp.clear();
					ktmp.push_back(xk_sym[0]);
					ktmp.push_back(xk_sym[1]);
					ktmp.push_back(xk_sym[2]);
					kpIBZ.push_back(KpointList(nloc, ktmp));
					++nsame;
				}
			}

#endif
		}

		if (nsame > 0) {
			nk_equiv.push_back(nsame);
		}
	}

	for (std::vector<unsigned int>::iterator p = nk_equiv.begin(); p != nk_equiv.end(); ++p){
		weight_k.push_back(static_cast<double>(*p)/static_cast<double>(nk));
	}

	// Construct the independent k-point list

	memory->allocate(nk_equiv_arr, nk_equiv.size());
	kpset_uniq.clear();
	unsigned int jk = 0;
	for (ik = 0; ik < nk_equiv.size(); ++ik){
		kpset_uniq.insert(kpIBZ[jk].knum);
		nk_equiv_arr[ik] = nk_equiv[ik];
		jk += nk_equiv[ik];
	}
	memory->deallocate(kequiv);

	nk_reduced = nk_equiv.size();
	nequiv_max = 0;

	for (ik = 0; ik < nk_reduced; ++ik){
		nequiv_max = std::max<unsigned int>(nequiv_max, nk_equiv[ik]);
	}

	memory->allocate(k_reduced, nk_reduced, nequiv_max);

	j = 0;
	for (ik = 0; ik < nk_reduced; ++ik) {
		for (i = 0; i < nequiv_max; ++i){
			k_reduced[ik][i] = 0;
		}
		for (i = 0; i < nk_equiv[ik]; ++i){
			k_reduced[ik][i] = kpIBZ[j++].knum;
		}
	}
}


int Kpoint::get_knum(const double kx, const double ky, const double kz)
{

	double diff[3];
	double dkx = static_cast<double>(nkx);
	double dky = static_cast<double>(nky);
	double dkz = static_cast<double>(nkz);

	diff[0] = static_cast<double>(nint(kx*dkx)) - kx*dkx;
	diff[1] = static_cast<double>(nint(ky*dky)) - ky*dky;
	diff[2] = static_cast<double>(nint(kz*dkz)) - kz*dkz;

	double norm = std::sqrt(diff[0]*diff[0] + diff[1]*diff[1] + diff[2]*diff[2]);

	if (norm >= eps12) {

		return -1;

	} else {

		int iloc, jloc, kloc;

		iloc = (nint(kx*dkx + 2.0 * dkx)) % nkx;
		jloc = (nint(ky*dky + 2.0 * dky)) % nky;
		kloc = (nint(kz*dkz + 2.0 * dkz)) % nkz;

		return kloc + nkz * jloc + nky * nkz * iloc;
	}
}


int Kpoint::nint(const double x)
{
	return static_cast<int>(x + 0.5 - (x < 0.0));
}
