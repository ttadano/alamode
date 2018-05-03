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
                  double);

        double fB(double,
                  double);

        double fC(double,
                  double);

        double Cv_tot(double);

        double internal_energy(double);

        double vibrational_entropy(double);

        double free_energy(double);

        double Cv_Debye(double,
                        double);

        double Cv_classical(double,
                            double);

        void Debye_T(double,
                     double &);

        double disp2_avg(double,
                         unsigned int,
                         unsigned int);

        double coth_T(double,
                      double);

        void compute_free_energy_bubble();
        void compute_free_energy_bubble_SCPH(double ***,
                                             std::complex<double> ****);

        void compute_FE_bubble(double **,
                               std::complex<double> ***,
                               double *);

        void compute_FE_bubble_SCPH(double ***,
                                    std::complex<double> ****,
                                    double *);
    };
}
