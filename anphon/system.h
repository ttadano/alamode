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

    class System : protected Pointers
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
        unsigned int **map_p2s_anharm_orig;

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

        double volume(double [3],
                      double [3],
                      double [3]);

        bool lspin, trevsym_mag;
        int noncollinear;

        int get_atomic_number_by_name(const std::string);

    private:
        void set_default_variables();
        void deallocate_variables();

        void set_mass_elem_from_database(const int,
                                         const std::string *,
                                         double *);

        void load_system_info_from_XML();

        void recips(double [3][3],
                    double [3][3]);

        void setup_atomic_class(unsigned int,
                                unsigned int *,
                                double **);

        void check_consistency_primitive_lattice();

        std::vector<std::string> element_names{
            "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar",
            "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr",
            "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe",
            "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf",
            "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
            "Pa", "U", "Np", "Pu"
        };
        std::vector<double> atomic_masses{
            1.007941, 4.002602, 6.940037, 9.012183, 10.811028, 12.010736, 14.006703, 15.999405, 18.998403, 20.180046,
            22.989769, 24.305052, 26.981539, 28.085499, 30.973762, 32.064787, 35.452938, 39.947799, 39.098301,
            40.078023, 44.955908, 47.866745, 50.941465, 51.996132, 54.938044, 55.845144, 58.933194, 58.693347,
            63.546040, 65.377783, 69.723066, 72.627550, 74.921595, 78.959389, 79.903528, 83.798000, 85.467664,
            87.616644, 88.905840, 91.223642, 92.906373, 95.959789, -1, 101.064940, 102.905498, 106.415328,
            107.868150, 112.411558, 114.818087, 118.710113, 121.759784, 127.603126, 126.904472, 131.292761, 132.905452,
            137.326892, 138.905469, 140.115731, 140.907658, 144.241596, -1, 150.366356, 151.964378, 157.252131,
            158.925355, 162.499473, 164.930329, 167.259083, 168.934218, 173.054150, 174.966815, 178.484979, 180.947876,
            183.841778, 186.206705, 190.224860, 192.216052, 195.084457, 196.966569, 200.599167, 204.383413, 207.216908,
            208.980399, -1, -1, -1, -1, -1, -1, 232.038056, 231.035884, 238.028910, -1, -1
        }; // For unstable elements, the atomic mass is set to -1
    };
}
