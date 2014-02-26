#include "mpi_common.h"
#include "isotope.h"
#include "memory.h"
#include "system.h"
#include "dynamical.h"
#include <iomanip>
#include "kpoint.h"
#include <complex>
#include "relaxation.h"
#include "constants.h"
#include "integration.h"
#include "error.h"

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

        memory->allocate(gamma_isotope, kpoint->nk_reduced, dynamical->neval);
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
    double epsilon = integration->epsilon;

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

            ret += omega1 * delta_lorentz(omega - omega1, epsilon) * prod;
            //		ret += relaxation->delta_lorentz(omega - omega1) * prod;
        }
    }


    ret *= pi * omega * 0.25 / static_cast<double>(nk);
    //	ret *= pi * omega * omega * 0.25 / static_cast<double>(nk);
}


void Isotope::calc_isotope_selfenergy_tetra(int knum, int snum, double omega, double &ret)
{
    int iat, icrd;
    int ik, is;
    double prod;
    std::complex<double> dprod;
    int nk = kpoint->nk;
    int ns = dynamical->neval;
    int natmin = system->natmin;

    ret = 0.0;

    double **eval;
    double **weight;

    memory->allocate(eval, ns, nk);
    memory->allocate(weight, ns, nk);

    for (ik = 0; ik < nk; ++ik) {
        for (is = 0; is < ns; ++is) {
            eval[is][ik] = dynamical->eval_phonon[ik][is];
        }
    }

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

            //			weight[is][ik] = prod;
            weight[is][ik] = prod * dynamical->eval_phonon[ik][is];
        }
    }

    for (is = 0; is < ns; ++is) {
        ret += integration->do_tetrahedron(eval[is], weight[is], omega);
    }

    //	ret *= pi * omega * omega * 0.25;
    ret *= pi * omega * 0.25;
}


void Isotope::calc_isotope_selfenergy_all()
{

    int nk = kpoint->nk;
    int ns = dynamical->neval;
    int i, j;

    int nks = kpoint->nk_reduced * ns;

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
            knum = kpoint->kpoint_irred_all[i / ns][0].knum;
            snum = i % ns;
            omega = dynamical->eval_phonon[knum][snum];
            calc_isotope_selfenergy_tetra(knum, snum, omega, tmp);
            //	calc_isotope_selfenergy(knum, snum, omega, tmp);
            gamma_loc[i] = tmp;
        }

        MPI_Reduce(&gamma_loc[0], &gamma_tmp[0], nks, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

        for (i = 0; i < kpoint->nk_reduced; ++i) {
            for (j = 0; j < ns; ++j) {
                gamma_isotope[i][j] = gamma_tmp[ns * i + j];
            }
        }

        memory->deallocate(gamma_tmp);
        memory->deallocate(gamma_loc);

        if (mympi->my_rank == 0) {
            std::cout << "Done !" << std::endl;
        }

        // 		double tmp2;
        // 
        // 		for (i = 0; i < kpoint->nk_reduced; ++i) {
        // 
        // 			for (int k = 0; k < ns; ++k) {
        // 				for (j = 0; j < kpoint->nk_equiv[i]; ++j) {
        // 					knum = kpoint->k_reduced[i][j];
        // 
        // 					omega = dynamical->eval_phonon[knum][k];
        // 					calc_isotope_selfenergy(knum, k, omega, tmp);
        // 					calc_isotope_selfenergy_tetra(knum, k, omega, tmp2);
        // 
        // 					std::cout << " i = " << std::setw(5) << i;
        // 					std::cout << " k = " << std::setw(5) << k;
        // 					std::cout << " j = " << std::setw(5) << j;
        // 					std::cout << " omega = " << std::setw(15) << omega;
        // 					std::cout << " ret1 = " << std::setw(15) << tmp;
        // 					std::cout << " ret2 = " << std::setw(15) << tmp2 << std::endl;
        // 				}
        // 			}
        // 		}
        // 		error->exit("hoge", "hoge");
    }
}

