/*
 system.h

 Copyright (c) 2014, 2015, 2016 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#pragma once

#include "pointers.h"
#include <string>
#include <vector>
#include <boost/property_tree/ptree.hpp>

namespace PHON_NS
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
            }
            if (this->element == a.element) {
                return this->magmom < a.magmom;
            }
            return false;
        }
    };

    class System: protected Pointers
    {
    public:
        System(class PHON *);
        ~System();

        void setup();

        double lavec_s[3][3], rlavec_s[3][3];
        double lavec_p[3][3], rlavec_p[3][3];
        double lavec_s_anharm[3][3], rlavec_s_anharm[3][3];
        double **xr_p, **xr_s, **xc;
        double **xr_s_anharm;
        double **magmom;
        double volume_p;

        unsigned int nat, natmin, ntran;
        unsigned int nat_anharm, ntran_anharm;
        unsigned int *kd, nkd;
        unsigned int *kd_anharm;

        unsigned int nclassatom;
        std::vector<unsigned int> *atomlist_class;

        unsigned int **map_p2s, **map_p2s_anharm;

        class Maps
        {
        public:
            unsigned int atom_num;
            unsigned int tran_num;
        };

        Maps *map_s2p, *map_s2p_anharm;

        std::string *symbol_kd;
        double *mass_kd, *mass, *mass_anharm;

        double Tmin, Tmax, dT;
        double volume(double [3], double [3], double [3]);

        bool lspin, trevsym_mag;
        int noncollinear;

    private:

        void load_system_info_from_XML();
        void recips(double [3][3], double [3][3]);
        void setup_atomic_class(unsigned int, unsigned int *, double **);
        void check_consistency_primitive_lattice();
    };
}
