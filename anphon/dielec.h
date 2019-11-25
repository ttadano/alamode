/*
 dielec.h

 Copyright (c) 2019 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#pragma once

#include "pointers.h"
#include <complex>

namespace PHON_NS
{
    class Dielec : protected Pointers
    {
    public:
        Dielec(class PHON *);
        ~Dielec();

        void init();
        void run_dielec_calculation();
        double* get_omega_grid(unsigned int &nomega) const;
        double*** get_dielectric_func() const;
        void compute_dielectric_function(const unsigned int nomega_in,
                                         double *omega_grid_in,
                                         double *eval_in, 
                                         std::complex<double> **evec_in,
                                         double ***dielec_out);
        int calc_dielectric_constant;

    private:

        void set_default_variables();
        void deallocate_variables();
        double *omega_grid;
        double ***dielec;
        unsigned int nomega;
        double emax, emin, delta_e;
    };
}
