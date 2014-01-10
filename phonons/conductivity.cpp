#include "mpi_common.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include "conductivity.h"
#include "dynamical.h"
#include "error.h"
#include "integration.h"
#include "kpoint.h"
#include "memory.h"
#include "parsephon.h"
#include "phonon_thermodynamics.h"
#include "phonon_velocity.h"
#include "relaxation.h"
#include "system.h"
#include "write_phonons.h"
#include "../alm_c++/constants.h"
#include <set>
#include <vector>
#include "../alm_c++/mathfunctions.h"
#include "isotope.h"

using namespace PHON_NS;

Conductivity::Conductivity(PHON *phon): Pointers(phon) {}
Conductivity::~Conductivity(){};

void Conductivity::setup_kappa()
{
	unsigned int i, j, k;
	unsigned int nks_total, nks_each_thread, nrem;

	nk = kpoint->nk;
	ns = dynamical->neval;

	ntemp = static_cast<unsigned int>((system->Tmax - system->Tmin) / system->dT);
	memory->allocate(Temperature, ntemp);
	memory->allocate(tau_l, ntemp);

	for (i = 0; i < ntemp; ++i) Temperature[i] = system->Tmin + static_cast<double>(i)*system->dT;

	nks_total = kpoint->nk_reduced * ns;
	nks_each_thread = nks_total / mympi->nprocs;
	nrem = nks_total - nks_each_thread * mympi->nprocs;

	if (nrem > 0) {
		memory->allocate(tau, (nks_each_thread + 1)*mympi->nprocs, ntemp);
	} else {
		memory->allocate(tau, nks_total, ntemp);
	}

	if (mympi->my_rank == 0) {
		memory->allocate(vel, nk, ns, 3);

		for (i = 0; i < nk; ++i){
			phonon_velocity->phonon_vel_k(kpoint->xk[i], vel[i]);

			// Generate phonon velocity in Cartesian coordinate
			for (j = 0; j < ns; ++j){
				rotvec(vel[i][j], vel[i][j], system->lavec_p);
				for (k = 0; k < 3; ++k) vel[i][j][k] /= 2.0 * pi;
				for (k = 0; k < 3; ++k) vel[i][j][k] *= Bohr_in_Angstrom*1.0e-10/time_ry;
			}
		}

		std::cout.setf(std::ios::fixed);
		std::cout << std::endl;
		std::cout << " Tmin = " << std::setw(10) << system->Tmin; 
		std::cout << " Tmax = " << std::setw(10) << system->Tmax; 
		std::cout << " dT   = " << std::setw(10) << system->dT; 
		std::cout << std::endl;

		std::cout.unsetf(std::ios::fixed);

		if (use_classical_Cv == 1) {
			std::cout << "Heat capacity will be replaced by kB (classical limit)" << std::endl;
		}
	}

	MPI_Bcast(&use_classical_Cv, 1, MPI_INT, 0, MPI_COMM_WORLD);

	vks_job.clear();

	for (i = 0; i < kpoint->nk_reduced; ++i) {
		for (j = 0; j < ns; ++j) {
			vks_job.insert(i*ns + j);
		}
	}
}

void Conductivity::prepare_restart()
{
	// Write phonon frequency to result file

	int ik, is;
	int i;
	std::set<int>::iterator it_set;
	std::string line_tmp;
	unsigned int nk_tmp, ns_tmp, nks_tmp;
	unsigned int multiplicity;
	int nks_done, *arr_done;

	double vel_dummy[3];

	nshift_restart = 0;

	vks_done.clear();

	if (mympi->my_rank == 0) {
		if (!phon->restart_flag) {

			writes->fs_result << "##Phonon Frequency" << std::endl;
			writes->fs_result << "#K-point (Symmetrically reduced), Branch, Omega (cm^-1)" << std::endl;

			for (i = 0; i < kpoint->nk_reduced; ++i) {
				ik = kpoint->k_reduced[i][0];
				for (is = 0; is < dynamical->neval; ++is) {
					writes->fs_result << std::setw(6) << i + 1 << std::setw(6) << is + 1;
					writes->fs_result << std::setw(15) << writes->in_kayser(dynamical->eval_phonon[ik][is]) << std::endl;
				}
			}

			writes->fs_result << "##END Phonon Frequency" << std::endl << std::endl;
			writes->fs_result << "##Phonon Relaxation Time" << std::endl;

		} else {

			while (writes->fs_result >> line_tmp) {

				if (line_tmp == "#TAU_EACH") {

					writes->fs_result >> nk_tmp >> ns_tmp;
					writes->fs_result >> multiplicity;

					nks_tmp = (nk_tmp - 1) * ns + ns_tmp -1;

					for (i = 0; i < multiplicity; ++i) {
						writes->fs_result >> vel_dummy[0] >> vel_dummy[1] >> vel_dummy[2];
					}

					for (i = 0; i < ntemp; ++i) {
						writes->fs_result >> tau[nks_tmp][i];
					}
					vks_done.push_back(nks_tmp);
				}
			}
		}

		writes->fs_result.close();
		writes->fs_result.open(writes->file_result.c_str(), std::ios::app | std::ios::out);
	}


	// Add vks_done list here

	if (mympi->my_rank == 0) {
		nks_done = vks_done.size();
	}
	MPI_Bcast(&nks_done, 1, MPI_INT, 0, MPI_COMM_WORLD);
	memory->allocate(arr_done, nks_done);

	if (mympi->my_rank == 0) {
		for (i = 0; i < nks_done; ++i) {
			arr_done[i] = vks_done[i];
		}
	}
	MPI_Bcast(arr_done, nks_done, MPI_INT, 0, MPI_COMM_WORLD);

	nshift_restart = nks_done;

	// Remove vks_done elements from vks_job

	for (i = 0; i < nks_done; ++i) {

		it_set = vks_job.find(arr_done[i]);

		if (it_set == vks_job.end()) {
			error->exit("prepare_restart", "This cannot happen");
		} else {
			vks_job.erase(it_set);
		}
	}

	memory->deallocate(arr_done);
	vks_done.clear();

}

void Conductivity::finish_kappa()
{
	if (mympi->my_rank == 0) {
		memory->deallocate(vel);
		memory->deallocate(kappa);
	}
	memory->deallocate(tau);
	memory->deallocate(Temperature);
}

void Conductivity::calc_anharmonic_tau()
{
	unsigned int nks_g;
	unsigned int i, j, k;
	unsigned int knum, snum;
	unsigned int *nks_thread;
	int iks;
	unsigned int iks_g;
	unsigned int nk_equiv;
	unsigned int ktmp;

	double omega, tau_tmp;

	// Distribute (k,s) to individual MPI threads

	nks_g = vks_job.size();
	vks_l.clear();

	unsigned int icount = 0;

	for (std::set<int>::iterator it = vks_job.begin(); it != vks_job.end(); ++it) {
		if (icount % mympi->nprocs == mympi->my_rank) {
			vks_l.push_back(*it);
		}
		++icount;
	}

	if (mympi->my_rank == 0) {
		memory->allocate(nks_thread, mympi->nprocs);
	}

	unsigned int nks_tmp = vks_l.size();
	MPI_Gather(&nks_tmp, 1, MPI_UNSIGNED, &nks_thread[mympi->my_rank], 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

	if (mympi->my_rank == 0) {
		std::cout << std::endl;
		std::cout << "Total Number of (k, s) pairs to be calculated : " << nks_g << std::endl;
		std::cout << "Assigned number of (k, s) pairs for each MPI threads below" << std::endl;
		for (i = 0; i < mympi->nprocs; ++i) {
			std::cout << " RANK: " << std::setw(5) << i + 1;
			std::cout << std::setw(8) << "NKS: " << std::setw(5) << nks_thread[i] << std::endl;
		}
		std::cout << std::endl;

		memory->deallocate(nks_thread);
	}

	unsigned int nk_tmp = nks_g / mympi->nprocs + 1;

	if (nks_g % mympi->nprocs !=0) {
		nk_tmp = nks_g / mympi->nprocs + 1;
	} else {
		nk_tmp = nks_g/ mympi->nprocs;
	}

	if (vks_l.size() < nk_tmp) {
		vks_l.push_back(-1);
	}

	for (i = 0; i < nk_tmp; ++i) {

		iks = vks_l[i];

		if (iks == -1) {

			for (j = 0; j < ntemp; ++j) tau_l[j] = 0.0; // do nothing

		} else {

			knum = kpoint->k_reduced[iks / ns][0];
			snum = iks % ns;

			omega = dynamical->eval_phonon[knum][snum];

			if (relaxation->ksum_mode == 0 || relaxation->ksum_mode == 1) {
				relaxation->calc_damping(ntemp, Temperature, omega, knum, snum, tau_l);
		//		relaxation->calc_damping2(ntemp, Temperature, omega, iks/ns, snum, tau_l);
			} else if (relaxation->ksum_mode == -1) {
				relaxation->calc_damping_tetra(ntemp, Temperature, omega, knum, snum, tau_l);
			}
		}

		MPI_Gather(&tau_l[0], ntemp, MPI_DOUBLE, tau[nshift_restart + i*mympi->nprocs], ntemp, MPI_DOUBLE, 0, MPI_COMM_WORLD);	

		if (mympi->my_rank == 0) {
			for (j = 0; j < mympi->nprocs; ++j) {

				iks_g = i*mympi->nprocs + j + nshift_restart;

				if (iks_g >= kpoint->nk_reduced*ns) break;

				writes->fs_result << "#TAU_EACH" << std::endl;
				writes->fs_result << iks_g / ns + 1 << " " << iks_g % ns + 1 << std::endl;

				nk_equiv = kpoint->nk_equiv_arr[iks_g / ns];

				writes->fs_result << nk_equiv << std::endl;
				for (k = 0; k < nk_equiv; ++k) {
					ktmp = kpoint->k_reduced[iks_g/ ns][k];
					writes->fs_result << std::setw(15) << vel[ktmp][iks_g % ns][0];
					writes->fs_result << std::setw(15) << vel[ktmp][iks_g % ns][1];
					writes->fs_result << std::setw(15) << vel[ktmp][iks_g % ns][2] << std::endl;
				}

				for (k = 0; k < ntemp; ++k) {
					tau_tmp = tau[iks_g][k];
					if (std::abs(tau_tmp) < eps) {
                                            //	tau[iks_g][k] = 0.0; 
					} else {
						tau[iks_g][k] = time_ry * 1.0e+12 / (2.0 * tau_tmp);
					}
					writes->fs_result << std::setw(15) << tau[iks_g][k] << std::endl;
				}
				writes->fs_result << "#END TAU_EACH" << std::endl;
			}
			std::cout <<  "ELEMENT " << std::setw(5) << i + 1 << " done." << std::endl;
		}
	}
}

void Conductivity::compute_kappa()
{
	unsigned int i, j, k;
	unsigned int iks;
	unsigned int knum, snum;

	double omega;
	unsigned int nk_equiv;
	unsigned int ktmp;

	if (mympi->my_rank == 0) {

		std::string file_kl;
		std::ofstream ofs_kl;
		double vv_tmp;
		int ieq;

		if (isotope->include_isotope) {
			for (iks = 0; iks < kpoint->nk_reduced*ns; ++iks) {
				knum = kpoint->k_reduced[iks / ns][0];
				snum = iks % ns;

				for (i = 0; i < ntemp; ++i) {
					tau[iks][i] = 1.0 / (1.0 / tau[iks][i] + 2.0 * isotope->gamma_isotope[iks/ns][snum] * 1.0e-12 / time_ry);
				}
			}
		}

		memory->allocate(kappa, ntemp, 3, 3);

		for (i = 0; i < ntemp; ++i) {
			for (j = 0; j < 3; ++j) {
				for (k = 0; k < 3; ++k) {

					kappa[i][j][k] = 0.0;

					for (iks = 0; iks < kpoint->nk_reduced*ns; ++iks) {

						knum = kpoint->k_reduced[iks / ns][0];
						snum = iks % ns;


						omega = dynamical->eval_phonon[knum][snum];

						vv_tmp = 0.0;
						nk_equiv = kpoint->nk_equiv_arr[iks / ns];

						for (ieq = 0; ieq < nk_equiv; ++ieq) {
							ktmp = kpoint->k_reduced[iks / ns][ieq];
							vv_tmp += vel[ktmp][snum][j] * vel[ktmp][snum][k];
						}

						if (use_classical_Cv == 1) {
							kappa[i][j][k] +=  phonon_thermodynamics->Cv_classical(omega, Temperature[i]) * vv_tmp * tau[iks][i];
						} else {
							kappa[i][j][k] +=  phonon_thermodynamics->Cv(omega, Temperature[i]) * vv_tmp * tau[iks][i];
						}
					}
					// Convert to SI unit
					kappa[i][j][k] *=  1.0e+18 / (std::pow(Bohr_in_Angstrom, 3) * system->volume_p * static_cast<double>(nk));
				}
			}
		}
	}
}