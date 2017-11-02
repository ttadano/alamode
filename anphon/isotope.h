/*
 isotope.h

 Copyright (c) 2014 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#include "mpi_common.h"
#include "pointers.h"


namespace PHON_NS
{
    class Isotope: protected Pointers
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
        void calc_isotope_selfenergy(const int, const int,
                                     const double, double &);
        void calc_isotope_selfenergy_tetra(const int, const int,
                                           const double, double &);
    };
}
