/*
 system.h

 Copyright (c) 2014 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#pragma once

#include "pointers.h"
#include <string>
#include <vector>

namespace ALM_NS
{
    class AtomType
    {
    public:
        int element;
        double magmom;

        bool operator<(const AtomType &a) const
        {
            if (this->element < a.element) {
                return true;
            } else if (this->element == a.element) {
                return this->magmom < a.magmom;
            } else {
                return false;
            }
        }
    };


    class System: protected Pointers
    {
    public:
        System(class ALM *);
        ~System();
        void init();
        void recips(double [3][3], double [3][3]);
        void frac2cart(double **);
        void load_reference_system();
        void load_reference_system_xml(std::string, const int, double *);

        int nat, nkd;
        int ndata, nstart, nend, nskip;
        int *kd;
        double lavec[3][3], rlavec[3][3];
        double **xcoord; // fractional coordinate
        double **x_cartesian;
        double **magmom;
        int noncollinear;
        std::string *kdname;

        unsigned int nclassatom;

        std::vector<unsigned int> *atomlist_class;
        bool lspin;
        double cell_volume;


    private:
        double volume(double [3], double [3], double [3]);
        void setup_atomic_class(int *);
    };
}
