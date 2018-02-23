/*
 isotope.h

 Copyright (c) 2014 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#include "pointers.h"


namespace PHON_NS
{
    class Isotope : protected Pointers
    {
    public:

        Isotope(class PHON *);
        ~Isotope();

        int include_isotope;
        double *isotope_factor;
        double **gamma_isotope;

        void setup_isotope_scattering();
        void calc_isotope_selfenergy_all();

    private:
        void set_default_variables();
        void deallocate_variables();

        void calc_isotope_selfenergy(int,
                                     int,
                                     double,
                                     double &);
        void calc_isotope_selfenergy_tetra(int,
                                           int,
                                           double,
                                           double &);
    };
}
