#include "mpi_common.h"
#include "isotope.h"
#include "memory.h"
#include "system.h"
#include "dynamical.h"
#include <iomanip>
#include "kpoint.h"
#include <complex>
#include "relaxation.h"
#include "../alm_c++/constants.h"

using namespace PHON_NS;

Isotope::Isotope(PHON *phon): Pointers(phon){};


Isotope::~Isotope(){
	if (include_isotope) {
		memory->deallocate(isotope_factor);
		memory->deallocate(gamma_isotope);
	}
};


void Isotope::setup_isotope_scattering()
{
	int i;
	int nkd = system->nkd;

	MPI_Bcast(&include_isotope, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD);

	if (include_isotope) {


		if (mympi->my_rank > 0) {
			memory->allocate(isotope_factor, nkd);
		}
		MPI_Bcast(&isotope_factor[0], nkd, MPI_DOUBLE, 0, MPI_COMM_WORLD);


		if (mympi->my_rank == 0) {
			std::cout << std::endl;
			std::cout << "ISOTOPE = 1: " << std::endl;
			std::cout << "Isotope scattering effects will be considered with the following scattering factors." << std::endl;

			for (i = 0; i < nkd; ++i) {
				std::cout << std::setw(5) << system->symbol_kd[i] << ":";
				std::cout << std::setw(15) << isotope_factor[i] << std::endl;
			}
			std::cout << std::endl;
		}

		memory->allocate(gamma_isotope, kpoint->nk, dynamical->neval);
	}
}

void Isotope::calc_isotope_selfenergy(int knum, int snum, double omega, double &ret)
{
	int iat, icrd;
	int ik, is;
	double prod;
	std::complex<double> dprod;
	int nk = kpoint->nk;
	int ns = dynamical->neval;
	int natmin = system->natmin;

	double omega1;
	dprod = 0.0;

	ret = 0.0;

	for (ik = 0; ik < nk; ++ik) {
		for (is = 0; is < ns; ++is) {

			prod = 0.0;

			for (iat = 0; iat < natmin; ++iat) {

				dprod = std::complex<double>(0.0, 0.0);
				for (icrd = 0; icrd < 3; ++icrd) {
					dprod += std::conj(dynamical->evec_phonon[ik][is][3 * iat + icrd]) * dynamical->evec_phonon[knum][snum][3 * iat + icrd];
				}

				prod += isotope_factor[iat] * std::norm(dprod);
			}

			omega1 = dynamical->eval_phonon[ik][is];

			ret += omega1 * relaxation->delta_lorentz(omega - omega1) * prod;
		}
	}

	ret *= pi * omega * 0.25 / static_cast<double>(nk);
}


void Isotope::calc_isotope_selfenergy_all()
{

	int nk = kpoint->nk;
	int ns = dynamical->neval;
	int i, j;

	int nks = nk * ns;

	double *gamma_tmp, *gamma_loc;
	double tmp, omega;

	int knum, snum;

	if (include_isotope) {

		if (mympi->my_rank == 0) {
			std::cout << "Calculating self-energies from isotope scatterings..." << std::endl;
		}
	
		memory->allocate(gamma_tmp, nks);
		memory->allocate(gamma_loc, nks);

		for (i = 0; i < nks; ++i) gamma_loc[i] = 0.0;

		for (i = mympi->my_rank; i < nks; i += mympi->nprocs) {
			knum = i / ns;
			snum = i % ns;
			omega = dynamical->eval_phonon[knum][snum];
			calc_isotope_selfenergy(knum, snum, omega, tmp);
			gamma_loc[i] = tmp;
		}

		MPI_Reduce(&gamma_loc[0], &gamma_tmp[0], nks, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		
		for (i = 0; i < nk; ++i) {
			for (j = 0; j < ns; ++j) {
				gamma_isotope[i][j] = gamma_tmp[ns * i + j];
			}
		}

		memory->deallocate(gamma_tmp);
		memory->deallocate(gamma_loc);

		if (mympi->my_rank == 0) {
			std::cout << "Done !" << std::endl;
		}
	}
}

