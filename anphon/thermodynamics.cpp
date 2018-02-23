/*
 phonon_thermodynamics.cpp

 Copyright (c) 2014, 2015, 2016 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#include "mpi_common.h"
#include "dynamical.h"
#include "error.h"
#include "kpoint.h"
#include "thermodynamics.h"
#include "pointers.h"
#include "system.h"
#include "constants.h"

using namespace PHON_NS;

Thermodynamics::Thermodynamics(PHON *phon): Pointers(phon)
{
    T_to_Ryd = k_Boltzmann / Ryd;
}

Thermodynamics::~Thermodynamics()
{
};

void Thermodynamics::setup()
{
    MPI_Bcast(&classical, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD);
}

double Thermodynamics::Cv(const double omega,
                          const double T)
{
    if (std::abs(T) < eps) return 0.0;

    double x = omega / (T_to_Ryd * T);
    return k_Boltzmann * std::pow(x / (2.0 * sinh(0.5 * x)), 2);
}

double Thermodynamics::Cv_classical(const double omega,
                                    const double T)
{
    if (std::abs(T) < eps) return 0.0;

    return k_Boltzmann;
}

double Thermodynamics::fB(const double omega,
                          const double T)
{
    if (std::abs(T) < eps || omega < eps8) return 0.0;

    double x = omega / (T_to_Ryd * T);
    return 1.0 / (std::exp(x) - 1.0);
}

double Thermodynamics::fC(const double omega,
                          const double T)
{
    if (std::abs(T) < eps || omega < eps8) return 0.0;

    double x = omega / (T_to_Ryd * T);
    return 1.0 / x;
}

double Thermodynamics::Cv_tot(const double T)
{
    int i;
    unsigned int ik, is;
    unsigned int nk = kpoint->nk;
    unsigned int ns = dynamical->neval;
    double omega;
    double ret = 0.0;
    int N = nk * ns;

    if (classical) {
#pragma omp parallel for private(ik, is, omega), reduction(+:ret)
        for (i = 0; i < N; ++i) {
            ik = i / ns;
            is = i % ns;

            omega = dynamical->eval_phonon[ik][is];
            if (omega < 0.0) continue;

            ret += Cv_classical(omega, T);
        }
    } else {
#pragma omp parallel for private(ik, is, omega), reduction(+:ret)
        for (i = 0; i < N; ++i) {
            ik = i / ns;
            is = i % ns;

            omega = dynamical->eval_phonon[ik][is];
            if (omega < 0.0) continue;

            ret += Cv(omega, T);
        }
    }

    return ret / static_cast<double>(nk);
}

double Thermodynamics::Cv_Debye(const double T,
                                const double TD)
{
    unsigned int natmin = system->natmin;
    unsigned int i;
    double d_theta, theta, theta_max;
    unsigned int ntheta;

    double x, y;
    double ret;
    d_theta = 0.001;

    if (TD < eps) {
        error->exit("Cv_Debye", "Too small TD");
    }
    if (T < 0.0) {
        error->exit("Cv_Debye", "Negative T");
    }

    if (T < eps) {
        return 0.0;
    } else {
        x = TD / T;
        theta_max = atan(x);

        ret = 0.0;
        ntheta = static_cast<unsigned int>(theta_max / d_theta);

        for (i = 0; i < ntheta; ++i) {
            theta = static_cast<double>(i) * d_theta;
            y = tan(theta);
            if (y > eps) {
                ret += std::pow(y, 4) * std::exp(y)
                    / std::pow((std::exp(y) - 1.0) * cos(theta), 2);
            }
        }
        y = tan(theta_max);
        ret += 0.5 * std::pow(y, 4) * std::exp(y)
            / std::pow((std::exp(y) - 1.0) * cos(theta_max), 2);

        return 9.0 * static_cast<double>(natmin) * k_Boltzmann
            * ret * d_theta / std::pow(x, 3);
    }
}

void Thermodynamics::Debye_T(const double T,
                             double &TD)
{
    double TD_old;
    double diff_C;
    double fdegfree = 1.0 / static_cast<double>(3.0 * system->natmin);

    double Cv_tot_tmp = Cv_tot(T);

    if (T > eps) {

        do {
            diff_C = fdegfree * (Cv_tot_tmp - Cv_Debye(T, TD)) / k_Boltzmann;

            TD_old = TD;
            TD = TD_old - diff_C * 10.0;
        } while (std::abs(diff_C) > 1.0e-5);
    }
}

double Thermodynamics::internal_energy(const double T)
{
    int i;
    unsigned int ik, is;
    unsigned int nk = kpoint->nk;
    unsigned int ns = dynamical->neval;
    double omega;
    double ret = 0.0;

    int N = nk * ns;

    if (classical) {
#pragma omp parallel for private(ik, is, omega), reduction(+:ret)
        for (i = 0; i < N; ++i) {
            ik = i / ns;
            is = i % ns;
            omega = dynamical->eval_phonon[ik][is];

            if (omega < eps8) continue;

            ret += T_to_Ryd * T;
        }
        ret *= 2.0;
    } else {
#pragma omp parallel for private(ik, is, omega), reduction(+:ret)
        for (i = 0; i < N; ++i) {
            ik = i / ns;
            is = i % ns;
            omega = dynamical->eval_phonon[ik][is];

            if (omega < eps8) continue;

            ret += omega * coth_T(omega, T);
        }
    }

    return ret * 0.5 / static_cast<double>(nk);
}

double Thermodynamics::vibrational_entropy(const double T)
{
    int i;
    unsigned int ik, is;
    unsigned int nk = kpoint->nk;
    unsigned int ns = dynamical->neval;
    double omega, x;
    double ret = 0.0;

    int N = nk * ns;

    if (classical) {
#pragma omp parallel for private(ik, is, omega, x), reduction(+:ret)
        for (i = 0; i < N; ++i) {
            ik = i / ns;
            is = i % ns;
            omega = dynamical->eval_phonon[ik][is];

            if (omega < eps8 || std::abs(T) < eps) continue;

            x = omega / (T * T_to_Ryd);
            ret += std::log(x) - 1.0;
        }
    } else {
#pragma omp parallel for private(ik, is, omega, x), reduction(+:ret)
        for (i = 0; i < N; ++i) {
            ik = i / ns;
            is = i % ns;
            omega = dynamical->eval_phonon[ik][is];

            if (omega < eps8 || std::abs(T) < eps) continue;

            x = omega / (T * T_to_Ryd);
            ret += std::log(1.0 - std::exp(-x)) - x / (std::exp(x) - 1.0);
        }
    }


    return -k_Boltzmann * ret / static_cast<double>(nk);
}

double Thermodynamics::free_energy(const double T)
{
    int i;
    unsigned int ik, is;
    unsigned int nk = kpoint->nk;
    unsigned int ns = dynamical->neval;
    double omega, x;
    double ret = 0.0;

    int N = nk * ns;

    if (classical) {
#pragma omp parallel for private(ik, is, omega, x), reduction(+:ret)
        for (i = 0; i < N; ++i) {
            ik = i / ns;
            is = i % ns;
            omega = dynamical->eval_phonon[ik][is];

            if (omega < eps8) continue;

            x = omega / (T * T_to_Ryd);
            ret += std::log(x);
        }

        return T * T_to_Ryd * ret / static_cast<double>(nk);

    } else {
#pragma omp parallel for private(ik, is, omega, x), reduction(+:ret)
        for (i = 0; i < N; ++i) {
            ik = i / ns;
            is = i % ns;
            omega = dynamical->eval_phonon[ik][is];

            if (omega < eps8) continue;

            if (std::abs(T) < eps) {
                ret += 0.5 * omega;
            } else {
                x = omega / (T * T_to_Ryd);
                ret += 0.5 * x + std::log(1.0 - std::exp(-x));
            }
        }

        if (std::abs(T) < eps) return ret / static_cast<double>(nk);

        return T * T_to_Ryd * ret / static_cast<double>(nk);
    }
}

double Thermodynamics::disp2_avg(const double T,
                                 const unsigned int ns1,
                                 const unsigned int ns2)
{
    int i;
    double ret = 0.0;
    unsigned int ik, is;
    unsigned int nk = kpoint->nk;
    unsigned int ns = dynamical->neval;
    double omega;

    int N = nk * ns;

    if (classical) {
#pragma omp parallel for private(ik, is, omega), reduction(+:ret)
        for (i = 0; i < N; ++i) {
            ik = i / ns;
            is = i % ns;
            omega = dynamical->eval_phonon[ik][is];

            // Skip when omega is almost zero. 
            // (neglect divergent contributions from acoustic modes at gamma point)
            if (omega < eps8) continue;

            ret += real(dynamical->evec_phonon[ik][is][ns1]
                    * std::conj(dynamical->evec_phonon[ik][is][ns2]))
                * T * T_to_Ryd / (omega * omega);
        }
    } else {
#pragma omp parallel for private(ik, is, omega), reduction(+:ret)
        for (i = 0; i < N; ++i) {
            ik = i / ns;
            is = i % ns;
            omega = dynamical->eval_phonon[ik][is];

            // Skip when omega is almost zero. 
            // (neglect divergent contributions from acoustic modes at gamma point)
            if (omega < eps8) continue;

            ret += real(dynamical->evec_phonon[ik][is][ns1]
                    * std::conj(dynamical->evec_phonon[ik][is][ns2]))
                * (fB(omega, T) + 0.5) / omega;
        }
    }


    ret *= 1.0 / (static_cast<double>(nk)
        * std::sqrt(system->mass[system->map_p2s[ns1 / 3][0]]
            * system->mass[system->map_p2s[ns2 / 3][0]]));

    return ret;
}

double Thermodynamics::coth_T(const double omega,
                              const double T)
{
    // This function returns coth(hbar*omega/2*kB*T)

    // if T = 0.0 and omega > 0, coth(hbar*omega/(2*kB*T)) = 1.0
    if (T < eps) return 1.0;

    double x = omega / (T_to_Ryd * T * 2.0);
    return 1.0 + 2.0 / (std::exp(2.0 * x) - 1.0);
}
