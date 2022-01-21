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

namespace PHON_NS {
class AtomType {
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

class System : protected Pointers {
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

    class Maps {
    public:
        unsigned int atom_num;
        unsigned int tran_num;
    };

    Maps *map_s2p, *map_s2p_anharm;

    std::string *symbol_kd;
    double *mass_kd, *mass, *mass_anharm;

    double Tmin, Tmax, dT;

    double volume(const double [3],
                  const double [3],
                  const double [3]) const;

    bool lspin, trevsym_mag;
    int noncollinear;

    int get_atomic_number_by_name(const std::string &);

private:
    void set_default_variables();

    void deallocate_variables();

    void set_mass_elem_from_database(const int,
                                     const std::string *,
                                     double *);

    void load_system_info_from_XML();

    void recips(double [3][3],
                double [3][3]) const;

    void setup_atomic_class(unsigned int,
                            const unsigned int *,
                            double **);

    void check_consistency_primitive_lattice() const;

    std::vector<std::string> element_names{
            "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar",
            "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br",
            "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I",
            "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu",
            "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac",
            "Th", "Pa", "U", "Np", "Pu"
    };
    std::vector<double> atomic_masses{
            1.00794, 4.002602, 6.941, 9.0121831, 10.811, 12.0107, 14.0067, 15.9994, 18.998403163, 20.1797, 22.98976928,
            24.305, 26.9815384, 28.0855, 30.973761998, 32.065, 35.453, 39.948, 39.0983, 40.078, 44.955908, 47.867,
            50.9415, 51.9961, 54.938043, 55.845, 58.933194, 58.6934, 63.546, 65.38, 69.723, 72.63, 74.921595, 78.971,
            79.904, 83.798, 85.4678, 87.62, 88.90584, 91.224, 92.90637, 95.95, -1, 101.07, 102.90549, 106.42, 107.8682,
            112.414, 114.818, 118.71, 121.76, 127.6, 126.90447, 131.293, 132.90545196, 137.327, 138.90547, 140.116,
            140.90766, 144.242, -1, 150.36, 151.964, 157.25, 158.925354, 162.5, 164.930328, 167.259, 168.934218,
            173.045, 174.9668, 178.486, 180.94788, 183.84, 186.207, 190.23, 192.217, 195.084, 196.96657, 200.592,
            204.3833, 207.2, 208.9804, -1, -1, -1, -1, -1, -1, 232.0377, 231.03588, 238.02891, -1, -1
    }; // They are standard atomic weight recommended by CIAAW.
    // Some recent changes are not considered because of the presence of uncertainty interval.
    // For unstable elements, the atomic mass is set to -1.
};
}
