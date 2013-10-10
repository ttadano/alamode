#include "mpi_common.h"
#include "interpolation.h"
#include "dynamical.h"
#include "kpoint.h"
#include "memory.h"
#include "fcs_phonon.h"
#include "system.h"
#include "conductivity.h"
#include "../alm_c++/constants.h"
#include <fftw3.h>
#include <complex>
#include <iomanip>


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

	unsigned int nk = nk1 * nk2 * nk3;

	int ns = dynamical->neval;
	std::complex<double> **mat_tmp;

	memory->allocate(mat_k, ns, ns, nk);
	memory->allocate(mat_r, ns, ns, nk);
	memory->allocate(mat_tmp, ns, ns);

	double xk_tmp[3];

	unsigned int i, j, k;
	unsigned int ik;


	for (ik = 0; ik < nk; ++ik) {
		for (k = 0; k < 3; ++k) {
			xk_tmp[k] = kpoint->xk[ik][k];
		}
		dynamical->calc_analytic_k(xk_tmp, fcs_phonon->fc2, mat_tmp);

		for (i = 0; i < ns; ++i) {
			for (j = 0; j < ns; ++j){
				mat_k[i][j][ik] = mat_tmp[i][j];
			}
		}
	}

	fftw_plan plan;

	for (i = 0; i < ns; ++i) {
		for (j = 0; j < ns; ++j) {
			plan = fftw_plan_dft_3d(nk1, nk2, nk3, reinterpret_cast<fftw_complex*>(mat_k[i][j]), reinterpret_cast<fftw_complex*>(mat_r[i][j]), FFTW_FORWARD, FFTW_ESTIMATE);
			fftw_execute(plan);
			fftw_destroy_plan(plan);

			for (ik = 0; ik < nk; ++ik) mat_r[i][j][ik] /= static_cast<double>(nk);
		}
	}

	// 	for (ik = 0; ik < nk; ++ik) {
	// 		for (i = 0; i < ns; ++i) {
	// 			for (j = 0; j < ns; ++j) {
	// 				std::cout << std::setw(15) << mat_r[i][j][ik].real();
	// 			}
	// 			std::cout << std::endl;
	// 		}
	// 		std::cout << std::endl;
	// 	}


	// 	for (ik = 0; ik < nk; ++ik) {
	// 
	// 		for (k = 0; k < 3; ++k) {
	// 			xk_tmp[k] = kpoint->xk[ik][k];
	// 		}
	// 
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
	// 
	// 		r2q(xk_tmp, nk, nk1, nk2, nk3, ns, mat_r, mat_tmp);
	// 
	// 		std::cout << "INTERPOLATE" << std::endl;
	// 		for (i = 0; i < ns; ++i) {
	// 			for (j = 0; j < ns; ++j) {
	// 				std::cout << std::setw(15) << mat_tmp[i][j].real() << std::setw(15) << mat_tmp[i][j].imag();
	// 			}
	// 			std::cout << std::endl;
	// 		}
	// 		std::cout << std::endl;
	// 
	// 	}


	if (mympi->my_rank == 0) {
		setup_damping();
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

				phase = static_cast<double>(ix) * xk_in[0] + static_cast<double>(iy) * xk_in[1] + static_cast<double>(iz) * xk_in[2];
				phase *= 2.0 * pi;

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
					damp[i][ik][is] = 4.0e-24 / std::pow(damp[i][ik][is] * time_ry, 2);
				}
			}
		}
	}


}