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
#include <boost/property_tree/ptree.hpp>

namespace PHON_NS {
    class System: protected Pointers {
    public:
        System(class PHON *);
        ~System();

        void setup();

        double lavec_s[3][3], rlavec_s[3][3];
        double lavec_p[3][3], rlavec_p[3][3];
        double **xr_p, **xr_s, **xc;
        double volume_p;

        unsigned int nat, natmin, ntran;
        unsigned int *kd, nkd;
        unsigned int cell_dimension[3];

        unsigned int nclassatom;
        std::vector<unsigned int> *atomlist_class;

        unsigned int **map_p2s;
        class Maps {
        public:
            unsigned int atom_num;
            unsigned int tran_num;
        };
        Maps *map_s2p;

        std::string *symbol_kd;
        double *mass_kd, *mass;

        double Tmin, Tmax, dT;
        double volume(double [3], double [3], double [3]);

    private:

        unsigned int *kd_prim;
        void load_system_info();
        void load_system_info_xml();
      //  std::string get_value_from_xml(boost::property_tree::ptree, std::string);
        void recips(double [3][3], double [3][3]);
        void setup_atomic_class(unsigned int, unsigned int *);
    };
}
