#include "phonon_thermodynamics.h"
#include "pointers.h"
#include "../alm_c++/constants.h"

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