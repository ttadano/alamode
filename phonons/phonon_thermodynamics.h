#pragma once

#include "pointers.h"

namespace PHON_NS {
    class Phonon_thermodynamics: protected Pointers {
    public:
        Phonon_thermodynamics(class PHON *);
        ~Phonon_thermodynamics();

        double T_to_Ryd;

        double Cv(const double, const double);
        double fB(const double, const double);
    };
}