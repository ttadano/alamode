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
#include <Eigen/Core>

namespace PHON_NS {


class Cell {
public:
    Eigen::Matrix3d lattice_vector;
    Eigen::Matrix3d reciprocal_lattice_vector;
    double volume;
    size_t number_of_atoms;
    size_t number_of_elems;
    std::vector<int> kind;
    Eigen::MatrixXd x_fractional;
    Eigen::MatrixXd x_cartesian;
    int has_entry{0};
};

class Spin {
public:
    int lspin;
    int time_reversal_symm;
    int noncollinear;
    std::vector<std::vector<double>> magmom;
};

class Maps {
public:
    unsigned int atom_num;
    unsigned int tran_num;
};

class MappingTable {
public:
    std::vector<Maps> to_true_primitive;
    std::vector<std::vector<unsigned int>> from_true_primitive;
};

class System : protected Pointers {
public:
    System(class PHON *);

    ~System();

    void setup();

    const Cell &get_cell(const std::string celltype, const std::string filetype="base") const;

    const Spin &get_spin(const std::string celltype) const;

    const MappingTable &get_mapping_table(const std::string celltype, const std::string filetype="base") const;

    Eigen::Matrix3d lavec_s, rlavec_s;
    Eigen::Matrix3d lavec_p, rlavec_p;
    Eigen::Matrix3d lavec_s_anharm, rlavec_s_anharm;
    Eigen::MatrixXd xr_p, xr_s, xc;
    Eigen::MatrixXd xr_s_anharm;

    int load_primitive_from_file;
    double volume_p;

    unsigned int nat, natmin, ntran;
    unsigned int nat_anharm, ntran_anharm;
    unsigned int *kd, nkd;
    unsigned int *kd_anharm;

    unsigned int **map_p2s, **map_p2s_anharm;
    unsigned int **map_p2s_anharm_orig;

    Maps *map_s2p, *map_s2p_anharm;

    std::string *symbol_kd;
    double *mass_kd, *mass, *mass_anharm;

    double Tmin, Tmax, dT;

    double volume(const double [3],
                  const double [3],
                  const double [3]) const;

    int get_atomic_number_by_name(const std::string &);

private:

    enum LatticeType {
        Direct, Reciprocal
    };

    Cell supercell_base, supercell_fc2, supercell_fc3, supercell_fc4;
    Cell primcell_base, primcell_fc2, primcell_fc3, primcell_fc4;
    Spin spin_super_base, spin_prim_base;
    MappingTable map_scell_base, map_pcell_base;
    MappingTable map_scell_fc2, map_pcell_fc2;
    MappingTable map_scell_fc3, map_pcell_fc3;
    MappingTable map_scell_fc4, map_pcell_fc4;

    void set_default_variables();

    void deallocate_variables();

    void set_mass_elem_from_database(const int,
                                     const std::string *,
                                     double *);

    void load_system_info_from_XML();

    void load_system_info_from_file();

    void get_structure_and_mapping_table_xml(const std::string &filename,
                                             Cell &scell_out,
                                             Cell &pcell_out,
                                             Spin &spin_super_out,
                                             Spin &spin_prim_out,
                                             MappingTable &map_super_out,
                                             MappingTable &map_prim_out,
                                             std::vector<std::string> &elements) const;

    void get_structure_and_mapping_table_h5(const std::string &filename,
                                            Cell &scell_out,
                                            Cell &pcell_out,
                                            Spin &spin_super_out,
                                            Spin &spin_prim_out,
                                            MappingTable &map_super_out,
                                            MappingTable &map_prim_out,
                                            std::vector<std::string> &elements) const;

    void update_primitive_lattice();

    void recips(const Eigen::Matrix3d &mat_in,
                Eigen::Matrix3d &rmat_out) const;

    void check_consistency_primitive_lattice() const;


    double volume(const Eigen::Matrix3d &mat_in,
                  const LatticeType latttype_in) const;

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
