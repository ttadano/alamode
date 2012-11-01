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
        double fC(const double, const double);

        void test_fB(const double);
        double Cv_tot(const double);
        double Internal_Energy(const double);
        double Cv_Debye(const double, const double);
        double Cv_classical(const double, const double);
        void Debye_T(const double, double &);

    public:
        double coth_T(const double, const double);

    };
}