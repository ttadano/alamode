/*
 system.h

 Copyright (c) 2014 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory
 or http://opensource.org/licenses/mit-license.php for information.
*/

#pragma once

#include <string>
#include <vector>
#include "timer.h"

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
            }
            if (this->element == a.element) {
                return this->magmom < a.magmom;
            }
            return false;
        }
    };

    class Cell
    {
    public:
        double lattice_vector[3][3];
        double reciprocal_lattice_vector[3][3];
        double volume;
        size_t number_of_atoms;
        size_t number_of_elems;
        std::vector<int> kind;
        std::vector<std::vector<double>> x_fractional;
        std::vector<std::vector<double>> x_cartesian;
    };

    class Spin
    {
    public:
        bool lspin;
        int time_reversal_symm;
        int noncollinear;
        std::vector<std::vector<double>> magmom;
    };

    class System
    {
    public:
        System();
        ~System();
        void init(const int,
                  Timer *);
        void frac2cart(double **) const;

        void set_supercell(const double [3][3],
                           const size_t,
                           const int *,
                           const double [][3]);
        void set_kdname(const std::string *);
        void set_periodicity(const int [3]);
        void set_spin_variables(const size_t nat,
                                const bool,
                                const int,
                                const int,
                                const double (*)[3]);
        void set_str_magmom(std::string);

        const Cell& get_supercell() const;
        double*** get_x_image() const;
        int* get_exist_image() const;
        std::string* get_kdname() const;
        int* get_periodicity() const;
        const Spin& get_spin() const;
        const std::string& get_str_magmom() const;
        const std::vector<std::vector<unsigned int>>& get_atomtype_group() const;

    private:
        // Variables for geometric structure
        Cell supercell;
        std::string *kdname;
        int *is_periodic; // is_periodic[3];
        double ***x_image;
        int *exist_image;

        // Variables for spins
        Spin spin;
        std::string str_magmom;

        // concatenate atomic kind and magmom (only for collinear case)
        std::vector<std::vector<unsigned int>> atomtype_group;

        enum LatticeType { Direct, Reciprocal };

        void set_reciprocal_latt(const double [3][3],
                                 double [3][3]) const;
        void set_default_variables();
        void deallocate_variables();

        double volume(const double [3][3],
                      LatticeType) const;
        void set_atomtype_group();

        void generate_coordinate_of_periodic_images();
        void print_structure_stdout(const Cell &);
        void print_magmom_stdout() const;
    };
}
