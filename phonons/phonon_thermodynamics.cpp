#include "phonon_thermodynamics.h"
#include "pointers.h"
#include "../alm_c++/constants.h"
#include "kpoint.h"
#include "dynamical.h"
#include "phonon_velocity.h"
#include <iostream>

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

void Phonon_thermodynamics::test_fB(const double T)
{

    unsigned int i, j;
    unsigned int nk = kpoint->nk;
    unsigned int ns = dynamical->neval;

    double omega;

    for (i = 0; i < nk; ++i){
        for (j = 0; j < ns; ++j){
            omega = phonon_velocity->freq(dynamical->eval_phonon[i][j]);
            std::cout << "omega = " << omega / (T_to_Ryd * T) << " ,fB = " << fB(omega, T) << std::endl;
        }
    }

}
