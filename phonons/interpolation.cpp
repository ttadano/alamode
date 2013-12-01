#include "mpi_common.h"
#include "interpolation.h"
#include "dynamical.h"
#include "kpoint.h"
#include "memory.h"
#include "fcs_phonon.h"
#include "system.h"
#include "conductivity.h"
#include "dynamical.h"
#include "../alm_c++/constants.h"
#include "write_phonons.h"
#include "error.h"
#include "relaxation.h"
#include "symmetry_core.h"
#include <fftw3.h>
#include <complex>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include "../alm_c++/mathfunctions.h"

#if defined(WIN32) || defined(_WIN32)
#pragma comment(lib, "libfftw3-3.lib")
#pragma comment(lib, "libfftw3f-3.lib")
#pragma comment(lib, "libfftw3l-3.lib")
#endif

using namespace PHON_NS;

Interpolation::Interpolation(PHON *phon) : Pointers(phon) {}
Interpolation::~Interpolation(){};

void Interpolation::prepare_interpolation()
{
	unsigned int ik, i, j;
	unsigned int nk_tmp[3];
	std::vector<KpointList> kplist_ref;

	parse_self_energy();

	nk_ref = nk1 * nk2 * nk3;
	nk_tmp[0] = nk1;
	nk_tmp[1] = nk2;
	nk_tmp[2] = nk3;

	// Generate symmetrically reduced k points

	memory->allocate(xk_interpolate, nk_ref, 3);	

	kpoint->gen_kmesh(symmetry->symmetry_flag, nk_tmp, xk_interpolate, nk_equiv_ref, kplist_ref);

	nk_reduced_ref = nk_equiv_ref.size();
	if (nk_reduced_ref != nksym_ref) {
		error->exit("prepare_interpolation", "The number of symmetry-reduced k points are not the same.");
	}

	nequiv_max_ref = 0;
	for (ik = 0; ik < nk_reduced_ref; ++ik){
		nequiv_max_ref = std::max<unsigned int>(nequiv_max_ref, nk_equiv_ref[ik]);
	}

	memory->allocate(k_reduced_ref, nk_reduced_ref, nequiv_max_ref);

	j = 0;
	for (ik = 0; ik < nk_reduced_ref; ++ik) {
		for (i = 0; i < nequiv_max_ref; ++i){
			k_reduced_ref[ik][i] = 0;
		}
		for (i = 0; i < nk_equiv_ref[ik]; ++i){
			k_reduced_ref[ik][i] = kplist_ref[j++].knum;
		}
	}

	// Generate k vectors for the non-analytic correction.

	double norm;

	memory->allocate(kvec_na_interpolate, nk_ref, 3);
	for (i = 0; i < nk_ref; ++i) {
		for (j = 0; j < 3; ++j) {
			kvec_na_interpolate[i][j] = dynamical->fold(xk_interpolate[i][j]);
		}
		rotvec(kvec_na_interpolate[i], kvec_na_interpolate[i], system->rlavec_p, 'T');
		norm = std::sqrt(kvec_na_interpolate[i][0] * kvec_na_interpolate[i][0] + kvec_na_interpolate[i][1] * kvec_na_interpolate[i][1] + kvec_na_interpolate[i][2] * kvec_na_interpolate[i][2]);

		if (norm > eps) {
			for (j = 0; j < 3; ++j) kvec_na_interpolate[i][j] /= norm;
		}
	}

	// Generate eigval and eigvec for interpolation
	prepare_dymat_for_interpolation();
	prepare_self_energy_extend();

// 	double xk_minus[3], diff[3];
// 	int iloc, jloc, kloc;
// 	int nk_inv;
// 
// 	for (ik = 0; ik < nk_ref; ++ik) {
// 		for (i = 0; i < 3; ++i) xk_minus[i] = -xk_interpolate[ik][i];
// 
// 		diff[0] = static_cast<double>(kpoint->nint(xk_minus[0]*static_cast<double>(nk1))) - xk_minus[0]*static_cast<double>(nk1);
// 		diff[1] = static_cast<double>(kpoint->nint(xk_minus[1]*static_cast<double>(nk2))) - xk_minus[1]*static_cast<double>(nk2);
// 		diff[2] = static_cast<double>(kpoint->nint(xk_minus[2]*static_cast<double>(nk3))) - xk_minus[2]*static_cast<double>(nk3);
// 
// 		norm = std::sqrt(diff[0]*diff[0] + diff[1]*diff[1] + diff[2]*diff[2]);
// 		if (norm > eps12) error->exit("prepare_dymat_for_interpolation", "Cannot find a k point.");
// 
// 		iloc = (kpoint->nint(xk_minus[0]*static_cast<double>(nk1) + 2.0 * static_cast<double>(nk1))) % nk1;
// 		jloc = (kpoint->nint(xk_minus[1]*static_cast<double>(nk2) + 2.0 * static_cast<double>(nk2))) % nk2;
// 		kloc = (kpoint->nint(xk_minus[2]*static_cast<double>(nk3) + 2.0 * static_cast<double>(nk3))) % nk3;
// 
// 		nk_inv = kloc + nk3 * jloc + nk2 * nk3 * iloc;  
// 
// 		std::cout << "ik = " << ik;
// 		for (i = 0; i < 3; ++i) std::cout << std::setw(15) << xk_interpolate[ik][i];
// 		std::cout << std::endl;
// 		std::cout << "minus ik = " << nk_inv;
// 		for (i = 0; i < 3; ++i) std::cout<< std::setw(15) << xk_interpolate[nk_inv][i];
// 		std::cout << std::endl;
// 
// 		for (i = 0; i < dynamical->neval; ++i) {
// 			for (j = 0; j < ntemp_ref; ++j) {
// 				if (std::abs(self_energy_extend[j][ik][i].imag() - self_energy_extend[j][nk_inv][i].imag())> eps) {
// 					std::cout << std::setw(15) << self_energy_extend[j][ik][i];
// 					std::cout << std::setw(15) << self_energy_extend[j][nk_inv][i];
// 					std::cout << std::endl;
// 				}
// 				
// 			}
// 			std::cout << std::endl;
// 		}
// 		std::cout << std::endl;
// 	}

}


void Interpolation::exec_interpolation()
{
	unsigned int i, j, k;
	unsigned int ik;
	unsigned int ns = dynamical->neval;
	double xk_tmp[3];
	std::complex<double> **mat_tmp;
	std::complex<double> *eval_complex;
	std::complex<double> ***matrix_k;
	std::vector<std::complex<double> > eval_vector;

	memory->allocate(mat_r, ns, ns, nk_ref);
	memory->allocate(mat_tmp, ns, ns);
	memory->allocate(eval_complex, ns);
	memory->allocate(matrix_k, ns, ns, nk_ref);

	create_matrix_k(1, matrix_k);

	fftw_plan plan;

	for (i = 0; i < ns; ++i) {
		for (j = 0; j < ns; ++j) {
			plan = fftw_plan_dft_3d(nk1, nk2, nk3, reinterpret_cast<fftw_complex*>(matrix_k[i][j]), 
				reinterpret_cast<fftw_complex*>(mat_r[i][j]), FFTW_FORWARD, FFTW_ESTIMATE);
			fftw_execute(plan);
			fftw_destroy_plan(plan);

			for (ik = 0; ik < nk_ref; ++ik) mat_r[i][j][ik] /= static_cast<double>(nk_ref);
		}
	}


// 	for (i = 0; i < ns; ++i) {
// 		for (j = 0; j < ns; ++j) {
// 			for (ik = 0; ik < nk_ref; ++ik) mat_r[i][j][ik] *= std::sqrt(system->mass[system->map_p2s[i / 3][0]] * system->mass[system->map_p2s[j / 3][0]]);
// 		}
// 	}
// 
// 	double x_tmp[3];
// 	double dist;
// 	int nktmp[3];
// 
// 	for (i = 0; i < system->natmin; ++i) {
// 		for (int mu = 0; mu < 3; ++mu) {
// 			for (int nu = 0; nu < 3; ++nu) {
// 				
// 				for (j = 0; j < system->natmin; ++j) {
// 					for (ik = 0; ik < nk_ref; ++ik) {
// 
// 						for (k = 0; k < 3; ++k) x_tmp[k] = system->xr_p[i][k] - system->xr_p[j][k];
// 						nktmp[0] = ik / (nk2 * nk3);
// 						nktmp[1] = (ik % (nk2 * nk3)) / nk3;
// 						nktmp[2] = ik % nk3;
// 
// 						if (nktmp[0] > nk1/2) nktmp[0] -= nk1;
// 						if (nktmp[1] > nk2/2) nktmp[1] -= nk2;
// 						if (nktmp[2] > nk3/2) nktmp[2] -= nk3;
// 
// 						x_tmp[0] -= static_cast<double>(nktmp[0]);
// 						x_tmp[1] -= static_cast<double>(nktmp[1]);
// 						x_tmp[2] -= static_cast<double>(nktmp[2]);
// 
// 						rotvec(x_tmp, x_tmp, system->lavec_p);
// 						dist = std::sqrt(x_tmp[0] * x_tmp[0] + x_tmp[1] * x_tmp[1] + x_tmp[2] * x_tmp[2]);
// 
// 						std::cout << "i = " << std::setw(4) << i + 1;
// 						std::cout << " mu = " << std::setw(3) << mu + 1 << " nu = " << std::setw(3) << nu + 1 ;
// 						std::cout << " j = " <<  std::setw(4) << j + 1;
// 						std::cout << " cell = " << std::setw(10) << ik + 1;
// 						std::cout << std::setw(15) << dist;
// 						std::cout  << std::setw(15) << mat_r[3 * i + mu][3 * j + nu][ik].real();
// 						std::cout << std:: setw(15) << mat_r[3 * i + mu][3 * j + nu][ik].imag() << std::endl;
// 					}
// 				}
// 			}
// 		}
// 	}
// 
// 	error->exit("hoge","hoge");

	for (ik = 0; ik < kpoint->nk; ++ik) {

		for (k = 0; k < 3; ++k) {
			xk_tmp[k] = kpoint->xk[ik][k];
		}

		// 		dynamical->calc_analytic_k(xk_tmp, fcs_phonon->fc2, mat_tmp);
		// 
		// 		std::cout << "ik = " << std::setw(6) << ik << std::endl;
		// 		std::cout << "ORIGINAL" << std::endl;
		// 		for (i = 0; i < ns; ++i) {
		// 			for (j = 0; j < ns; ++j) {
		// 				std::cout << std::setw(15) << mat_tmp[i][j].real() << std::setw(15) << mat_tmp[i][j].imag();
		// 			}
		// 			std::cout << std::endl;
		// 		}
		// 		std::cout << std::endl;

		for (k = 0; k < 3; ++k) {
			std::cout << std::setw(15) << kpoint->xk[ik][k];
		}
		std::cout << std::endl;

		r2q(xk_tmp, nk_ref, nk1, nk2, nk3, ns, mat_r, mat_tmp);

		std::cout << "INTERPOLATE" << std::endl;
		for (i = 0; i < ns; ++i) {
			for (j = 0; j < ns; ++j) {
				std::cout << std::setw(15) << mat_tmp[i][j].real() << std::setw(15) << mat_tmp[i][j].imag();
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;

		diagonalize_interpolated_matrix(mat_tmp, eval_complex);

		eval_vector.clear();

		for (i = 0; i < ns; ++i) {
			eval_vector.push_back(eval_complex[i]);
		}
		std::sort(eval_vector.begin(), eval_vector.end(), &compare_real);


		for (i = 0; i < ns; ++i) {
			std::cout << std::setw(15) << eval_vector[i].real() << std::setw(15) << -eval_vector[i].imag();
			std::cout << std::setw(15) << dynamical->eval_phonon[ik][i] << std::endl;
			//	std::cout << std::setw(15) << self_energy_extend[1][ik][i].imag() << std::endl;
		}
		std::cout << std::endl;
	}

	eval_vector.clear();
	memory->deallocate(mat_r);
	memory->deallocate(mat_tmp);
	memory->deallocate(eval_complex);
	memory->deallocate(matrix_k);
}


void Interpolation::finish_interpolation()
{
	memory->deallocate(xk_interpolate);	
	memory->deallocate(kvec_na_interpolate);
	memory->deallocate(k_reduced_ref);
	memory->deallocate(self_energy_extend);
	memory->deallocate(self_energy);
	memory->deallocate(eigval);
	memory->deallocate(eigvec);

}
void Interpolation::diagonalize_interpolated_matrix(std::complex<double> **mat_in, std::complex<double> *eval_out)
{
	unsigned int i, j, k;
	std::complex<double> *amat;

	char JOBZ, UPLO;
	char JOBVL, JOBVR;
	int LDVL, LDVR;
	int INFO, LWORK;
	double *RWORK;
	std::complex<double> *WORK;
	std::complex<double> *evec_left, *evec_right;

	int ns = dynamical->neval;

	LWORK = (2 * ns - 1) * 10;
	// memory->allocate(RWORK, 3* ns - 2);
	memory->allocate(RWORK, 2 * ns);
	memory->allocate(WORK, LWORK);

	LDVL = 1;
	LDVR = 1;
	memory->allocate(evec_left, LDVL);
	memory->allocate(evec_right, LDVR);

	JOBZ = 'N';
	UPLO = 'U';

	JOBVL = 'N';
	JOBVR = 'N';

	memory->allocate(amat, ns * ns);


	k = 0;
	for(i = 0; i < ns; ++i){
		for (j = 0; j < ns; ++j){
			amat[k++] = mat_in[i][j];
		}
	}

	//	zheev_(&JOBZ, &UPLO, &ns, amat, &ns, eval_out, WORK, &LWORK, RWORK, &INFO);

	zgeev_(&JOBVL, &JOBVR, &ns, amat, &ns, eval_out, evec_left, &LDVL, evec_right, &LDVR, WORK, &LWORK, RWORK, &INFO);

	for (i = 0; i < ns; ++i) {
		eval_out[i] = std::sqrt(eval_out[i]);
	}


}

void Interpolation::r2q(double *xk_in, unsigned int nk, unsigned int nx, unsigned int ny, unsigned int nz, unsigned int ns, std::complex<double> ***mat_r, std::complex<double> **mat_k)
{
	unsigned int ix, iy, iz;
	unsigned int icell;
	unsigned int i, j;
	double phase;
	std::complex<double> im(0.0, 1.0);


	for (i = 0; i < ns; ++i) {
		for (j = 0; j < ns; ++j) {
			mat_k[i][j] = std::complex<double>(0.0, 0.0);
		}
	}


	for (ix = 0; ix < nx; ++ix) {
		for (iy = 0; iy < ny; ++iy) {
			for (iz = 0; iz < nz; ++iz) {
				icell = iz + nz * iy + nz * ny * ix;

				phase = 2.0 * pi * (static_cast<double>(ix) * xk_in[0] + static_cast<double>(iy) * xk_in[1] + static_cast<double>(iz) * xk_in[2]);

				for (i = 0; i < ns; ++i) {
					for (j = 0; j < ns; ++j) {
						mat_k[i][j] += mat_r[i][j][icell] * exp(im*phase);
					}
				}
			}
		}
	}
}


void Interpolation::setup_damping()
{
	unsigned int i, j;
	unsigned int ik, is;
	int knum, knum0;

	double Tmin = system->Tmin;
	double Tmax = system->Tmax;
	double dT = system->dT;

	unsigned int nk = kpoint->nk;
	unsigned int ns = dynamical->neval;
	unsigned int NT = static_cast<unsigned int>((Tmax - Tmin) / dT);

	memory->allocate(damp, NT, nk, ns);

	for (ik = 0; ik < kpoint->nk_reduced; ++ik) {

		knum0 = kpoint->k_reduced[ik][0];

		for (is = 0; is < ns; ++is) {

			for (i = 0; i < NT; ++i) {
				damp[i][knum0][is] = conductivity->tau[ik * ns + is][i];
			}
		}

		for (i = 0; i < NT; ++i) {
			for (j = 1; j < kpoint->nk_equiv[ik]; ++j) {
				knum = kpoint->k_reduced[ik][j];

				for (is = 0; is < ns; ++is) {
					damp[i][knum][is] = damp[i][knum0][is];
				}
			}
		}
	}

	for (i = 0; i < NT; ++i) {
		for (ik = 0; ik < nk; ++ik) {
			for (is = 0; is < ns; ++is) {
				if (std::pow(damp[i][ik][is], 2) < eps) {
					//		if (ik == 0 && (is == 0 || is == 1 || is == 2)) {
					damp[i][ik][is] = 0.0;
				} else {
					damp[i][ik][is] = 1.0e-12 / (damp[i][ik][is] * time_ry);
				}
			}
		}
	}
}

void Interpolation::prepare_self_energy_extend()
{
	unsigned int i, j;
	unsigned int ik, is;
	int knum, knum0;

	unsigned int ns = dynamical->neval;

	memory->allocate(self_energy_extend, ntemp_ref, nk_ref, ns);

	for (ik = 0; ik < nk_reduced_ref; ++ik) {

		knum0 = k_reduced_ref[ik][0];

		for (is = 0; is < ns; ++is) {

			for (i = 0; i < ntemp_ref; ++i) {
				self_energy_extend[i][knum0][is] = self_energy[i][ik * ns + is];
			}
		}

		for (i = 0; i < ntemp_ref; ++i) {
			for (j = 1; j < nk_equiv_ref[ik]; ++j) {
				knum = k_reduced_ref[ik][j];

				for (is = 0; is < ns; ++is) {
					self_energy_extend[i][knum][is] = self_energy_extend[i][knum0][is];
				}
			}
		}
	}
}



void Interpolation::create_matrix_k(unsigned int NT, std::complex<double> ***mat_k)
{
	std::complex<double> *polarization_matrix, *mat_tmp;
	std::complex<double> *eigval_matrix, *dmat;
	std::complex<double> alpha, *beta;

	std::complex<double> **dmat_tmp;
	std::complex<double> im(0.0, 1.0);

	double xk_tmp[3];

	unsigned int ik, is, js;
	int ns = dynamical->neval;

	unsigned int ns2 = ns * ns;
	unsigned int m;


	alpha = std::complex<double>(1.0, 0.0);

	char TRANSA[] = "N";
	char TRANSB[] = "C";

	memory->allocate(polarization_matrix, ns2);
	memory->allocate(mat_tmp, ns2);
	memory->allocate(eigval_matrix, ns2);
	memory->allocate(beta, ns);
	memory->allocate(dmat, ns2);


	for (is = 0; is < ns; ++is) beta[is] = std::complex<double>(0.0, 0.0);

	for (ik = 0; ik < nk_ref; ++ik) {

		// create eigval matrix

		for (is = 0; is < ns2; ++is) eigval_matrix[is] = std::complex<double>(0.0, 0.0);

		m = 0;
		for (is = 0; is < ns; ++is) {
			for (js = 0; js < ns; ++js) {
				if (is == js) {
					// eigval_matrix[m] = std::pow(eigval[ik][is] - self_energy_extend[1][ik][is], 2);
					// eigval_matrix[m] = std::pow(eigval[ik][is], 2);
					   eigval_matrix[m] = std::pow(- self_energy_extend[1][ik][is], 2);

				}
				++m;
			}
		}

		// create polarization matrix

		m = 0;

		for (is = 0; is < ns; ++is) {
			for (js = 0; js < ns; ++js) {
				polarization_matrix[m++] = eigvec[ik][is][js];
			}
		}
		zgemm_(TRANSA, TRANSB, &ns, &ns, &ns, &alpha, eigval_matrix, &ns, polarization_matrix, &ns, beta, mat_tmp, &ns);
		zgemm_(TRANSA, TRANSA, &ns, &ns, &ns, &alpha, polarization_matrix, &ns, mat_tmp, &ns, beta, dmat, &ns);


		m = 0;
		for (is = 0; is < ns; ++is) {
			for (js = 0; js < ns; ++js) {
				mat_k[is][js][ik] = dmat[m];
				++m;
			}
		}
	}

	memory->deallocate(polarization_matrix);
	memory->deallocate(mat_tmp);
	memory->deallocate(eigval_matrix);
	memory->deallocate(beta);
	memory->deallocate(dmat);

}

void Interpolation::parse_self_energy()
{
	// Parse self energy of phonons from the given .result file.

	std::fstream fs_result;

	// Restart
	fs_result.open(writes->file_result.c_str(), std::ios::in | std::ios::out);
	if (!fs_result) {
		error->exit("setup_result_io", "Could not open file_result");
	}

	// Check the consistency

	std::string line_tmp, str_tmp;
	int natmin_tmp, nkd_tmp;
	int nksym_tmp;
	int is_classical, ismear;
	double epsilon_tmp, T1, T2, delta_T;		

	bool found_tag;

	found_tag = false;
	while (fs_result >> line_tmp)
	{
		if (line_tmp == "#SYSTEM") {
			found_tag = true;
			break;
		}
	}
	if (!found_tag) error->exit("setup_result_io", "Could not find #SYSTEM tag");

	fs_result >> natmin_tmp >> nkd_tmp;

	if (!(natmin_tmp == system->natmin && nkd_tmp == system->nkd)) {
		error->exit("setup_result_io", "SYSTEM information is not consistent");
	}

	found_tag = false;
	while (fs_result >> line_tmp)
	{
		if (line_tmp == "#KPOINT") {
			found_tag = true;
			break;
		}
	}
	if (!found_tag) error->exit("setup_result_io", "Could not find #KPOINT tag");

	fs_result >> nk1 >> nk2 >> nk3;
	fs_result >> nksym_ref;

	found_tag = false;
	while (fs_result >> line_tmp)
	{
		if (line_tmp == "#FCSINFO") {
			found_tag = true;
			break;
		}
	}
	if (!found_tag) error->exit("setup_result_io", "Could not find #FCSINFO tag");

	fs_result >> str_tmp;
	if (str_tmp != fcs_phonon->file_fcs) {
		error->warn("setup_result_io", "FCSINFO is not consistent");
	}

	found_tag = false;
	while (fs_result >> line_tmp)
	{
		if (line_tmp == "#TEMPERATURE") {
			found_tag = true;
			break;
		}
	}
	if (!found_tag) error->exit("setup_result_io", "Could not find #TEMPERATURE tag");

	fs_result >> T1 >> T2 >> delta_T;

	ntemp_ref = static_cast<unsigned int>((T2 - T1) / delta_T);

	// Extract phonon self-energies

	unsigned int i, j;
	unsigned int nk_tmp, ns_tmp;
	unsigned int nks_tmp;
	unsigned int multiplicity;
	double vel_dummy[3];
	std::complex<double> im(0.0, 1.0);

	unsigned int ns = dynamical->neval;
	double tau_tmp;

	memory->allocate(self_energy, ntemp_ref, nksym_ref * ns);

	for (i = 0; i < ntemp_ref; ++i) {
		for (j = 0; j < nksym_ref * ns; ++j) {
			self_energy[i][j] = std::complex<double>(0.0, 0.0);
		}
	}

	while (fs_result >> line_tmp) {

		if (line_tmp == "#TAU_EACH") {

			fs_result >> nk_tmp >> ns_tmp;
			fs_result >> multiplicity;

			nks_tmp = (nk_tmp - 1) * ns + ns_tmp -1;

			for (i = 0; i < multiplicity; ++i) {
				fs_result >> vel_dummy[0] >> vel_dummy[1] >> vel_dummy[2];
			}

			for (i = 0; i < ntemp_ref; ++i) {
				fs_result >> tau_tmp;
				if (std::abs(tau_tmp) < eps) {
					self_energy[i][nks_tmp] = 0.0;
				} else  {
					self_energy[i][nks_tmp] = im *  time_ry * 1.0e+12 / (2.0 * tau_tmp);
				}
			}
		}
	}

}

void Interpolation::prepare_dymat_for_interpolation()
{
	unsigned int ik, is, js;
	unsigned int i;
	unsigned int ns = dynamical->neval;

	memory->allocate(eigval, nk_ref, ns);
	memory->allocate(eigvec, nk_ref, ns, ns);

	for (ik = 0; ik < nk_ref; ++ik){
		if (fcs_phonon->is_fc2_ext) {
			dynamical->eval_k(xk_interpolate[ik], kvec_na_interpolate[ik], fcs_phonon->fc2_ext, eigval[ik], eigvec[ik], true);
		} else {
			dynamical->eval_k(xk_interpolate[ik], kvec_na_interpolate[ik], fcs_phonon->fc2, eigval[ik], eigvec[ik], true);
		}

		for (is = 0; is < ns; ++is){
			eigval[ik][is] = dynamical->freq(eigval[ik][is]);
		}
	}

	bool *flag_done;
	unsigned int nk_inv;
	std::complex<double> *evec_tmp;
	// 
	// 	if (kpoint->kpoint_mode == 2) {
	// 		modify_eigenvectors();
	// 	}

	double xk_minus[3];
	int iloc, jloc, kloc;

	double diff[3];
	double norm;
	memory->allocate(flag_done, nk_ref);
	memory->allocate(evec_tmp, ns);


 	for (ik = 0; ik < nk_ref; ++ik) flag_done[ik] = false;
//	for (ik = 0; ik < nk_ref; ++ik) flag_done[ik] = true;


	for (ik = 0; ik < nk_ref; ++ik){

		if (!flag_done[ik]) {

			for (i = 0; i < 3; ++i) xk_minus[i] = -xk_interpolate[ik][i];

			diff[0] = static_cast<double>(kpoint->nint(xk_minus[0]*static_cast<double>(nk1))) - xk_minus[0]*static_cast<double>(nk1);
			diff[1] = static_cast<double>(kpoint->nint(xk_minus[1]*static_cast<double>(nk2))) - xk_minus[1]*static_cast<double>(nk2);
			diff[2] = static_cast<double>(kpoint->nint(xk_minus[2]*static_cast<double>(nk3))) - xk_minus[2]*static_cast<double>(nk3);

		    norm = std::sqrt(diff[0]*diff[0] + diff[1]*diff[1] + diff[2]*diff[2]);
			if (norm > eps12) error->exit("prepare_dymat_for_interpolation", "Cannot find a k point.");

			iloc = (kpoint->nint(xk_minus[0]*static_cast<double>(nk1) + 2.0 * static_cast<double>(nk1))) % nk1;
			jloc = (kpoint->nint(xk_minus[1]*static_cast<double>(nk2) + 2.0 * static_cast<double>(nk2))) % nk2;
			kloc = (kpoint->nint(xk_minus[2]*static_cast<double>(nk3) + 2.0 * static_cast<double>(nk3))) % nk3;

			nk_inv = kloc + nk3 * jloc + nk2 * nk3 * iloc;  

			for (is = 0; is < ns; ++is){
				for (js = 0; js < ns; ++js){
					evec_tmp[js] =eigvec[ik][is][js];
				}

				for (js = 0; js < ns; ++js){
					eigvec[nk_inv][is][js] = std::conj(evec_tmp[js]);
				}
			}

			flag_done[ik] = true;
			flag_done[nk_inv] = true;
		}
	}

	memory->deallocate(flag_done);
	memory->deallocate(evec_tmp);
}

