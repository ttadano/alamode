/*
 isotope.cpp

 Copyright (c) 2014 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#include "mpi_common.h"
#include "isotope.h"
#include "constants.h"
#include "dynamical.h"
#include "error.h"
#include "integration.h"
#include "kpoint.h"
#include "memory.h"
#include "system.h"
#include <iomanip>
#include <complex>

using namespace PHON_NS;

Isotope::Isotope(PHON *phon) : Pointers(phon)
{
    set_default_variables();
};

Isotope::~Isotope()
{
    deallocate_variables();
};

void Isotope::set_default_variables()
{
    include_isotope = false;
    isotope_factor = nullptr;
    gamma_isotope = nullptr;
}

void Isotope::deallocate_variables()
{
    if (isotope_factor) {
        memory->deallocate(isotope_factor);
    }
    if (gamma_isotope) {
        memory->deallocate(gamma_isotope);
    }
}


void Isotope::setup_isotope_scattering()
{
    const int nkd = system->nkd;

    MPI_Bcast(&include_isotope, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (include_isotope) {

        if (mympi->my_rank == 0) {
            if (!isotope_factor) {
                memory->allocate(isotope_factor, nkd);
                set_isotope_factor_from_database(nkd,
                                                 system->symbol_kd,
                                                 isotope_factor);
            }
        }

        if (mympi->my_rank > 0) {
            memory->allocate(isotope_factor, nkd);
        }
        MPI_Bcast(&isotope_factor[0], nkd, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        if (mympi->my_rank == 0) {
            std::cout << " ISOTOPE >= 1: Isotope scattering effects will be considered" << std::endl;
            std::cout << "               with the following scattering factors." << std::endl;

            for (int i = 0; i < nkd; ++i) {
                std::cout << std::setw(5) << system->symbol_kd[i] << ":";
                std::cout << std::scientific << std::setw(17) << isotope_factor[i] << std::endl;
            }
            std::cout << std::endl;

            memory->allocate(gamma_isotope, kpoint->nk_irred, dynamical->neval);
        }
    }
}

void Isotope::calc_isotope_selfenergy(const int knum,
                                      const int snum,
                                      const double omega,
                                      double &ret) const
{
    // Compute phonon selfenergy of phonon (knum, snum) 
    // due to phonon-isotope scatterings.
    // Delta functions are replaced by smearing functions with width EPSILON.

    const auto nk = kpoint->nk;
    const auto ns = dynamical->neval;
    const auto natmin = system->natmin;
    const auto epsilon = integration->epsilon;

    ret = 0.0;

    for (auto ik = 0; ik < nk; ++ik) {
        for (auto is = 0; is < ns; ++is) {

            auto prod = 0.0;

            for (auto iat = 0; iat < natmin; ++iat) {

                auto dprod = std::complex<double>(0.0, 0.0);
                for (auto icrd = 0; icrd < 3; ++icrd) {
                    dprod += std::conj(dynamical->evec_phonon[ik][is][3 * iat + icrd])
                             * dynamical->evec_phonon[knum][snum][3 * iat + icrd];
                }
                prod += isotope_factor[system->kd[system->map_p2s[iat][0]]] * std::norm(dprod);
            }

            const auto omega1 = dynamical->eval_phonon[ik][is];

            if (integration->ismear == 0) {
                ret += omega1 * delta_lorentz(omega - omega1, epsilon) * prod;
            } else {
                ret += omega1 * delta_gauss(omega - omega1, epsilon) * prod;
            }
        }
    }

    ret *= pi * omega * 0.25 / static_cast<double>(nk);
}


void Isotope::calc_isotope_selfenergy_tetra(const int knum,
                                            const int snum,
                                            const double omega,
                                            double &ret) const
{
    // Compute phonon selfenergy of phonon (knum, snum) 
    // due to phonon-isotope scatterings.
    // This version employs the tetrahedron method.

    int ik, is;
    const auto nk = kpoint->nk;
    const auto ns = dynamical->neval;
    const auto natmin = system->natmin;

    ret = 0.0;

    double *eval;
    double *weight;

    memory->allocate(eval, nk);
    memory->allocate(weight, nk);

    for (is = 0; is < ns; ++is) {
        for (ik = 0; ik < nk; ++ik) {

            auto prod = 0.0;

            for (auto iat = 0; iat < natmin; ++iat) {

                auto dprod = std::complex<double>(0.0, 0.0);
                for (auto icrd = 0; icrd < 3; ++icrd) {
                    dprod += std::conj(dynamical->evec_phonon[ik][is][3 * iat + icrd])
                             * dynamical->evec_phonon[knum][snum][3 * iat + icrd];
                }
                prod += isotope_factor[system->kd[system->map_p2s[iat][0]]] * std::norm(dprod);
            }

            weight[ik] = prod * dynamical->eval_phonon[ik][is];
            eval[ik] = dynamical->eval_phonon[ik][is];
        }
        ret += integration->do_tetrahedron(eval, weight, omega);
    }

    ret *= pi * omega * 0.25;

    memory->deallocate(eval);
    memory->deallocate(weight);
}


void Isotope::calc_isotope_selfenergy_all() const
{
    int i;
    const auto ns = dynamical->neval;
    const auto nks = kpoint->nk_irred * ns;
    double tmp;
    double *gamma_tmp = nullptr;
    double *gamma_loc = nullptr;

    if (include_isotope) {

        if (mympi->my_rank == 0) {
            std::cout << " Calculating self-energies from isotope scatterings ... ";
        }

        if (mympi->my_rank == 0) {
            memory->allocate(gamma_tmp, nks);
        } else {
            memory->allocate(gamma_tmp, 1);
        }
        memory->allocate(gamma_loc, nks);

        for (i = 0; i < nks; ++i) gamma_loc[i] = 0.0;

        for (i = mympi->my_rank; i < nks; i += mympi->nprocs) {
            const auto knum = kpoint->kpoint_irred_all[i / ns][0].knum;
            const auto snum = i % ns;
            const auto omega = dynamical->eval_phonon[knum][snum];
            if (integration->ismear == -1) {
                calc_isotope_selfenergy_tetra(knum, snum, omega, tmp);
            } else {
                calc_isotope_selfenergy(knum, snum, omega, tmp);
            }
            gamma_loc[i] = tmp;
        }

        MPI_Reduce(&gamma_loc[0], &gamma_tmp[0], nks,
                   MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

        if (mympi->my_rank == 0) {
            for (i = 0; i < kpoint->nk_irred; ++i) {
                for (int j = 0; j < ns; ++j) {
                    gamma_isotope[i][j] = gamma_tmp[ns * i + j];
                }
            }
        }

        memory->deallocate(gamma_tmp);
        memory->deallocate(gamma_loc);

        if (mympi->my_rank == 0) {
            std::cout << "done!" << std::endl;
        }
    }
}

void Isotope::set_isotope_factor_from_database(const int nkd,
                                               const std::string *symbol_in,
                                               double *isofact_out)
{
    for (int i = 0; i < nkd; ++i) {
        const auto atom_number = system->get_atomic_number_by_name(symbol_in[i]);
        if (atom_number >= isotope_factors.size() || atom_number == -1) {
            error->exit("set_isotope_factor_from_database",
                        "The isotope factor for the given element doesn't exist in the database.\n"
                        "Therefore, please input ISOFACT manually.");
        }
        const auto isofact_tmp = isotope_factors[atom_number];
        if (isofact_tmp < -0.5) {
            error->exit("set_isotope_factor_from_database",
                        "One of the elements in the KD-tag is unstable. "
                        "Therefore, please input ISOFACT manually.");
        }
        isofact_out[i] = isofact_tmp;
    }
}
