/*
 phonon_thermodynamics.h

 Copyright (c) 2014 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#pragma once

#include "pointers.h"

namespace PHON_NS
{
    class Thermodynamics: protected Pointers
    {
    public:
        Thermodynamics(class PHON *);
        ~Thermodynamics();

        double T_to_Ryd;

        double Cv(const double, const double);
        double fB(const double, const double);
        double fC(const double, const double);

        double Cv_tot(const double);
        double internal_energy(const double);
        double vibrational_entropy(const double);
        double free_energy(const double);
        double Cv_Debye(const double, const double);
        double Cv_classical(const double, const double);
        void Debye_T(const double, double &);
        double disp2_avg(const double, const unsigned int, const unsigned int);

    public:
        double coth_T(const double, const double);
    };
}
