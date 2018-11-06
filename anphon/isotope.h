/*
 isotope.h

 Copyright (c) 2014 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#include "pointers.h"
#include <vector>

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
                                     double &) const;

        void calc_isotope_selfenergy_tetra(int,
                                           int,
                                           double,
                                           double &) const;

        void set_isotope_factor_from_database(const int,
                                              const std::string *,
                                              double *);

        std::vector<double> isotope_factors{
            0.00011460742533846188, 8.141017510851699e-08, 0.001458820047621947, 0.0, 0.0013539153245004177,
            7.387230454537788e-05, 1.8376657914224594e-05, 3.3588025065260784e-05, 0.0, 0.0008278312160539562, 0.0,
            0.0007398827105109379, 0.0, 0.00020070043953754622, 0.0, 0.00016512811739588583, 0.0005827012790438193,
            3.4803742597617525e-05, 0.0001640001019257936, 0.00029756377957496945, 0.0, 0.00028645556480110174,
            9.548320927761995e-07, 0.0001328704171052577, 0.0, 8.244411853712163e-05, 0.0, 0.00043070442089852476,
            0.00021093312802220393, 0.0005931533235650079, 0.00019712731536263687, 0.000586966208421238, 0.0,
            0.0004627901461335774, 0.00015627716697599108, 0.00024880742737704476, 0.00010969638987469167,
            6.099700015374489e-05, 0.0, 0.00034262903024295327, 0.0, 0.0005961083777486782, -1, 0.0004066661504740231,
            0.0, 0.0003094784411091063, 8.579847704787673e-05, 0.0002716036180261363, 1.245588189674909e-05,
            0.0003340852797872776, 6.607553852631361e-05, 0.0002839333030058624, 0.0, 0.0002675566535085368, 0.0,
            6.237013178021676e-05, 4.5917491111023726e-08, 2.2495590932891925e-05, 0.0, 0.00023237187990100274, -1,
            0.000334685954430736, 4.3279441126609935e-05, 0.0001276749037273739, 0.0, 5.198070714285335e-05, 0.0,
            7.23248017569606e-05, 0.0, 8.556028002982757e-05, 8.300794202558322e-07, 5.253850496170609e-05,
            3.6715121208243084e-09, 6.966807117351981e-05, 2.7084982818795603e-05, 7.452354225251159e-05,
            2.5378700157091918e-05, 3.4285145177491544e-05, 0.0, 6.525193204276608e-05, 1.9964351041965618e-05,
            1.9437780365209887e-05, 0.0, -1, -1, -1, -1, -1, -1, 0.0, 0.0, 1.1564592331193284e-06, -1, -1
        }; // Precalculated isotope factors. For unstable elements, the value is set to -1.
    };
}
