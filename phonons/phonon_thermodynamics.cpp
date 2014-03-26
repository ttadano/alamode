/*
 phonon_thermodynamics.cpp

 Copyright (c) 2014 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#include "mpi_common.h"
#include <iostream>
#include <iomanip>
#include "dynamical.h"
#include "fcs_phonon.h"
#include "error.h"
#include "kpoint.h"
#include "memory.h"
#include "phonon_thermodynamics.h"
#include "phonon_velocity.h"
#include "pointers.h"
#include "system.h"
#include "constants.h"

using namespace PHON_NS;

Phonon_thermodynamics::Phonon_thermodynamics(PHON *phon): Pointers(phon) {
    T_to_Ryd = k_Boltzmann / Ryd;
}

Phonon_thermodynamics::~Phonon_thermodynamics(){};

double Phonon_thermodynamics::Cv(const double omega, const double T)
{
    double x;

    if (std::abs(T) < eps) {
        return 0.0;
    } else {
        x = omega / (T_to_Ryd * T);
        return k_Boltzmann * std::pow(x/(2.0 * sinh(0.5*x)), 2);
    }
}

double Phonon_thermodynamics::Cv_classical(const double omega, const double T)
{
    double x;

    if (std::abs(T) < eps) {
        return 0.0;
    } else {
        x = omega / (T_to_Ryd * T);
        return k_Boltzmann * std::pow(x, 2) * std::exp(-x);
    }
}

double Phonon_thermodynamics::fB(const double omega, const double T)
{
    double x;

    if (std::abs(T) < eps){
        return 0.0;
    } else {
        x = omega / (T_to_Ryd * T);
        return 1.0 / (std::exp(x) - 1.0);
    }
}

double Phonon_thermodynamics::fC(const double omega, const double T)
{
    double x;

    if (std::abs(T) < eps) {
        return 0.0;
    } else {
        x = omega / (T_to_Ryd * T);
        return std::exp(-x);
    }
}

void Phonon_thermodynamics::test_fB(const double T)
{

    unsigned int i, j;
    unsigned int nk = kpoint->nk;
    unsigned int ns = dynamical->neval;

    double omega;

    for (i = 0; i < nk; ++i){
        for (j = 0; j < ns; ++j){
            omega = dynamical->eval_phonon[i][j];
            std::cout << "omega = " << omega / (T_to_Ryd * T) << " ,fB = " << fB(omega, T) << std::endl;
        }
    }
}

double Phonon_thermodynamics::Cv_tot(const double T)
{
    unsigned int ik, is;
    unsigned int nk = kpoint->nk;
    unsigned int ns = dynamical->neval;
    double omega;

    double ret = 0.0;

    for (ik = 0; ik < nk; ++ik){
        for (is = 0; is < ns; ++is){
            omega = dynamical->eval_phonon[ik][is];
            if (omega < 0.0) continue;
            ret += Cv(omega, T);
        }
    }
    return ret / static_cast<double>(nk);
}

double Phonon_thermodynamics::Cv_Debye(const double T, const double TD)
{
    unsigned int natmin = system->natmin;
    unsigned int i;
    double d_theta, theta, theta_max;
    unsigned int ntheta;

    double x, y;
    double ret;
    d_theta = 0.001;

    if (TD < eps)  {
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
        ntheta = static_cast<unsigned int>(theta_max/d_theta);

        for (i = 0; i < ntheta; ++i){
            theta = static_cast<double>(i) * d_theta;
            y = tan(theta);
            if (y > eps){
                ret += std::pow(y, 4)*std::exp(y) / std::pow((std::exp(y) - 1.0) * cos(theta), 2);
            }
        }
        y = tan(theta_max);
        ret += 0.5 * std::pow(y, 4)*std::exp(y) / std::pow((std::exp(y) - 1.0) * cos(theta_max), 2);

        return 9.0 * static_cast<double>(natmin) * k_Boltzmann * ret * d_theta/ std::pow(x, 3);
    }
}

void Phonon_thermodynamics::Debye_T(const double T, double &TD)
{    
    double TD_old;
    double diff_C;
    double fdegfree = 1.0 / static_cast<double>(3.0 * system->natmin);

    if (T > eps) {

        do {   
            //       std::cout << "T = " << T << " , TD = " << TD << std::endl;
            diff_C = fdegfree * (Cv_tot(T) - Cv_Debye(T, TD)) / k_Boltzmann;

            TD_old = TD;
            TD = TD_old - diff_C * 10.0;
        } while(std::abs(diff_C) > 1.0e-5);
    }
}

double Phonon_thermodynamics::Internal_Energy(const double T)
{
    unsigned int ik, is;
    unsigned int nk = kpoint->nk;
    unsigned int ns = dynamical->neval;
    double omega;

    double ret = 0.0;

    for (ik = 0; ik < nk; ++ik){
        for (is = 0; is < ns; ++is){
            omega = dynamical->eval_phonon[ik][is];
            if (omega <= 0.0) {
                // exactly zero is anomalous
                continue;
            } else {
                ret += omega * coth_T(omega, T);
            }
        }
    }
    return ret / static_cast<double>(nk);
}

double Phonon_thermodynamics::disp2_avg(const double T, const unsigned int ns1, const unsigned int ns2)
{
    double ret = 0.0;
    unsigned int ik, is;
    unsigned int nk = kpoint->nk;
    unsigned int ns = dynamical->neval;
    double omega;

    for (ik = 0; ik < nk; ++ik) {
        for (is = 0; is < ns; ++is) {

            omega = dynamical->eval_phonon[ik][is];

            // Skip when omega is almost zero. 
            // (neglect divergent contributions from acoustic modes at gamma point)
            if (omega < eps8) {
                //				std::cout << "ik = " << ik << " is = " << is << " omega = " << omega << std::endl;
                continue;
            }
            ret += real(dynamical->evec_phonon[ik][is][ns1] * std::conj(dynamical->evec_phonon[ik][is][ns2])) * (fB(omega, T) + 0.5) / omega;

        }
    }

    ret *= 1.0 / (static_cast<double>(nk) * std::sqrt(system->mass[system->map_p2s[ns1/3][0]] * system->mass[system->map_p2s[ns2/3][0]]));

    // ret *= 2.0 * electron_mass / time_ry * Bohr_in_Angstrom * Bohr_in_Angstrom;
    // ret *= h_planck / (2.0 * pi); // Convert to SI unit 
    // Note that hbar (Dirac's constant) is equal to 2*Me*a0**2/time_ry

    return ret;
}

double Phonon_thermodynamics::coth_T(const double omega, const double T)
{
    if (T < eps) {
        // if T = 0.0 and omega > 0, coth(hbar*omega/(2*kB*T)) = 1.0
        return 1.0;
    } else {
        double x = omega / (T_to_Ryd * T);
        return 1.0 + 2.0 / (std::exp(2.0 * x) - 1.0);
    }
}

