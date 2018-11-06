/*
 phonon_thermodynamics.h

 Copyright (c) 2014 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#pragma once

#include "pointers.h"
#include <complex>

namespace PHON_NS
{
    class Thermodynamics : protected Pointers
    {
    public:
        Thermodynamics(class PHON *);

        ~Thermodynamics();

        double T_to_Ryd;
        bool classical;
        bool calc_FE_bubble;
        double *FE_bubble;

        void setup();

        double Cv(double,
                  double) const;

        double fB(double,
                  double) const;

        double fC(double,
                  double) const;

        double Cv_tot(double) const;

        double internal_energy(double) const;

        double vibrational_entropy(double) const;

        double free_energy(double) const;

        double Cv_classical(double,
                            double) const;

        double disp2_avg(double,
                         unsigned int,
                         unsigned int) const;

        double coth_T(double,
                      double) const;

        void compute_free_energy_bubble();


        void compute_FE_bubble(double **,
                               std::complex<double> ***,
                               double *) const;

        double compute_FE_bubble_SCPH(double,
                                      double **,
                                      std::complex<double> ***) const;
    };
}
