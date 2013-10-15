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
#include <fftw3.h>
#include <complex>
#include <iomanip>
#include <fstream>

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
	nk1 = kpoint->nkx;
	nk2 = kpoint->nky;
	nk3 = kpoint->nkz;

	unsigned int nk_ref = nk1 * nk2 * nk3;

	int ns = dynamical->neval;
	std::complex<double> **mat_tmp;
	std::complex<double> *amat;
	double *eval_out;
	std::complex<double> *eval_complex;

	parse_self_energy();
	memory->allocate(xk_interpolate, nk_ref);

	memory->allocate(mat_k, ns, ns, nk);
	memory->allocate(mat_r, ns, ns, nk);
	memory->allocate(mat_tmp, ns, ns);
	memory->allocate(amat, ns * ns);
	memory->allocate(eval_out, ns);
	memory->allocate(eval_complex, ns);


	double xk_tmp[3];

	unsigned int i, j, k;
	unsigned int ik;

	// 	for (ik = 0; ik < nk; ++ik) {
	// 		for (k = 0; k < 3; ++k) {
	// 			xk_tmp[k] = kpoint->xk[ik][k];
	// 		}
	// 		dynamical->calc_analytic_k(xk_tmp, fcs_phonon->fc2, mat_tmp);
	// 
	// 		for (i = 0; i < ns; ++i) {
	// 			for (j = 0; j < ns; ++j){
	// 				mat_k[i][j][ik] = mat_tmp[i][j];
	// 			}
	// 		}
	// 	}

	memory->allocate(matrix_k, ns, ns, nk);
	if (mympi->my_rank == 0) {
		setup_damping();
		create_matrix_k(1);
	}

	fftw_plan plan;

	for (i = 0; i < ns; ++i) {
		for (j = 0; j < ns; ++j) {
			plan = fftw_plan_dft_3d(nk1, nk2, nk3, reinterpret_cast<fftw_complex*>(matrix_k[i][j]), 
				reinterpret_cast<fftw_complex*>(mat_r[i][j]), FFTW_FORWARD, FFTW_ESTIMATE);
			fftw_execute(plan);
			fftw_destroy_plan(plan);

			for (ik = 0; ik < nk; ++ik) mat_r[i][j][ik] /= static_cast<double>(nk);
		}
	}


	char JOBZ, UPLO;
	char JOBVL, JOBVR;
	int LDVL, LDVR;
	int INFO, LWORK;
	double *RWORK;
	std::complex<double> *WORK;
	std::complex<double> *evec_left, *evec_right;

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


	for (ik = 0; ik < nk; ++ik) {

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

		r2q(xk_tmp, nk, nk1, nk2, nk3, ns, mat_r, mat_tmp);

		std::cout << "INTERPOLATE" << std::endl;
		for (i = 0; i < ns; ++i) {
			for (j = 0; j < ns; ++j) {
				std::cout << std::setw(15) << mat_tmp[i][j].real() << std::setw(15) << mat_tmp[i][j].imag();
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;



		k = 0;
		for(i = 0; i < ns; ++i){
			for (j = 0; j < ns; ++j){
				amat[k++] = mat_tmp[i][j];
			}
		}


		//	zheev_(&JOBZ, &UPLO, &ns, amat, &ns, eval_out, WORK, &LWORK, RWORK, &INFO);

		zgeev_(&JOBVL, &JOBVR, &ns, amat, &ns, eval_complex, evec_left, &LDVL, evec_right, &LDVR, WORK, &LWORK, RWORK, &INFO);

		for (i = 0; i < ns; ++i) {
				eval_complex[i] = std::sqrt(eval_complex[i]);
		}


		for (i = 0; i < ns; ++i) {
			//		std::cout << std::setw(15) << std::sqrt(eval_out[i]) * 2.0 * 1.0e-12 / time_ry << std::endl;
			std::cout << std::setw(15) << eval_complex[i].real() << std::setw(15) << -eval_complex[i].imag();
			std::cout << std::setw(15) << dynamical->eval_phonon[ik][i];
		 	std::cout << std::setw(15) << damp[1][ik][i] << std::endl;
		}
		std::cout << std::endl;
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

void Interpolation::create_matrix_k(unsigned int NT)
{
	std::complex<double> *polarization_matrix, *mat_tmp;
	std::complex<double> *eigval_matrix, *dmat;
	std::complex<double> alpha, *beta;

	std::complex<double> **dmat_tmp;
	std::complex<double> im(0.0, 1.0);

	double xk_tmp[3];

	unsigned int ik, is, js;
	unsigned int nk = kpoint->nk;
	int ns = dynamical->neval;

	unsigned int ns2 = ns * ns;
	unsigned int m;

	memory->allocate(polarization_matrix, ns2);
	memory->allocate(mat_tmp, ns2);
	memory->allocate(eigval_matrix, ns2);
	memory->allocate(beta, ns);
	memory->allocate(dmat, ns2);

	memory->allocate(dmat_tmp, ns, ns);

	for (is = 0; is < ns; ++is) beta[is] = std::complex<double>(0.0, 0.0);


	alpha = std::complex<double>(1.0, 0.0);


	char TRANSA[] = "N";
	char TRANSB[] = "C";


	for (ik = 0; ik < nk; ++ik) {

		// create eigval matrix

		for (is = 0; is < ns2; ++is) eigval_matrix[is] = std::complex<double>(0.0, 0.0);

		m = 0;
		for (is = 0; is < ns; ++is) {
			for (js = 0; js < ns; ++js) {
				if (is == js) {
						eigval_matrix[m] = std::pow(dynamical->eval_phonon[ik][is] - im * damp[1][ik][is], 2);
					//	eigval_matrix[m] = 1.0;
					// eigval_matrix[m] = std::pow(dynamical->eval_phonon[ik][is], 2);

					//	eigval_matrix[m] = std::pow(dynamical->eval_phonon[ik][is], 2);
				}
				++m;
			}
		}

		// create polarization matrix

		m = 0;

		for (is = 0; is < ns; ++is) {
			for (js = 0; js < ns; ++js) {
				polarization_matrix[m++] = dynamical->evec_phonon[ik][is][js];
			}
		}
		zgemm_(TRANSA, TRANSB, &ns, &ns, &ns, &alpha, eigval_matrix, &ns, polarization_matrix, &ns, beta, mat_tmp, &ns);
		zgemm_(TRANSA, TRANSA, &ns, &ns, &ns, &alpha, polarization_matrix, &ns, mat_tmp, &ns, beta, dmat, &ns);


		m = 0;
		for (is = 0; is < ns; ++is) {
			for (js = 0; js < ns; ++js) {
				matrix_k[is][js][ik] = dmat[m];
				++m;
			}
		}

		// 		for (int k = 0; k < 3; ++k) {
		// 			xk_tmp[k] = kpoint->xk[ik][k];
		// 		}
		// 		dynamical->calc_analytic_k(xk_tmp, fcs_phonon->fc2, dmat_tmp);
		// 
		// 
		// 		std::cout << "ORIGINAL" << std::endl;
		// 		for (is = 0; is < ns; ++is) {
		// 			for (js = 0; js < ns; ++js) {
		// 				std::cout << std::setw(15) << dmat_tmp[is][js].real() << std::setw(15) << dmat_tmp[is][js].imag();
		// 			}
		// 			std::cout << std::endl;
		// 		}
		// 		std::cout << std::endl;
		// 
		// 		std::cout << "MODIFIED" << std::endl;
		// 
		// 		m = 0;
		// 		for (is = 0; is < ns; ++is) {
		// 			for (js = 0; js < ns; ++js) {
		// 				std::cout << std::setw(15) << dmat[m].real() << std::setw(15) << dmat[m].imag();
		// 				++m;
		// 			}
		// 			std::cout << std::endl;
		// 		}
		// 		std::cout << std::endl;
	}
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
	int nk_tmp[3], nksym_tmp;
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

	fs_result >> nk_tmp[0] >> nk_tmp[1] >> nk_tmp[2];

	nk1 = nk_tmp[0];
	nk2 = nk_tmp[1];
	nk3 = nk_tmp[2];

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

// 	if (!(T1 == system->Tmin && T2 == system->Tmax && delta_T == system->dT)) {
// 		error->exit("setup_result_io", "Temperature information is not consistent");
// 
}