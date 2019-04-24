/*
 alm.h

 Copyright (c) 2014, 2015, 2016 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory
 or http://opensource.org/licenses/mit-license.php for information.
*/

#pragma once

#include <string>
#include "system.h"
#include "cluster.h"
#include "fcs.h"
#include "symmetry.h"
#include "optimize.h"
#include "constraint.h"
#include "files.h"
#include "patterndisp.h"
#include "timer.h"

namespace ALM_NS
{
    class ALM
    {
    public:
        ALM();
        ~ALM();

        class Cluster *cluster{};
        class Fcs *fcs{};
        class Symmetry *symmetry{};
        class Optimize *optimize{};
        class Constraint *constraint{};
        class Files *files{};
        class Displace *displace{};
        class Timer *timer{};

        void set_run_mode(std::string run_mode_in);
        std::string get_run_mode() const;
        void set_verbosity(int verbosity_in);
        int get_verbosity() const;
        void set_output_filename_prefix(std::string prefix) const;
        void set_print_symmetry(int printsymmetry) const;
        void set_print_hessian(bool print_hessian) const;
        void set_symmetry_tolerance(double tolerance) const;
        void set_displacement_param(bool trim_dispsign_for_evenfunc) const;
        void set_displacement_basis(std::string str_disp_basis) const;
        void set_periodicity(const int is_periodic[3]) const;
        void set_cell(size_t nat,
                      const double lavec[3][3],
                      const double xcoord[][3],
                      const int kind[],
                      const std::string kdname[]) const;
        void set_magnetic_params(const size_t nat,
                                 const double (*magmom)[3],
                                 const bool lspin,
                                 const int noncollinear,
                                 const int trev_sym_mag,
                                 const std::string str_magmom) const;
        void set_displacement_and_force(const double *u_in,
                                        const double *f_in,
                                        int nat,
                                        int ndata_used) const;
        void set_constraint_type(int constraint_flag) const;
        void set_rotation_axis(std::string rotation_axis) const;
        void set_sparse_mode(int sparse_mode) const;
        //void set_fitting_filenames(std::string dfile,
        //                           std::string ffile) const;
        void define(const int maxorder,
                    const size_t nkd,
                    const int *nbody_include,
                    const double *cutoff_radii) const;
        //int get_ndata_used() const;
        size_t get_nrows_sensing_matrix() const;
        Cell get_supercell() const;
        std::string* get_kdname() const;
        Spin get_spin() const;
        void set_str_magmom(std::string);
        std::string get_str_magmom() const;
        double*** get_x_image() const;
        int* get_periodicity() const;

        const std::vector<std::vector<int>>& get_atom_mapping_by_pure_translations() const;
        int get_maxorder() const;
        int* get_nbody_include() const;

        size_t get_number_of_displacement_patterns(const int fc_order) const; // harmonic=1, ...
        void get_number_of_displaced_atoms(int *numbers,
                                           int fc_order) const; // harmonic=1, ...
        int get_displacement_patterns(int *atom_indices,
                                      double *disp_patterns,
                                      int fc_order) const;          // harmonic=1, ...
        size_t get_number_of_fc_elements(const int fc_order) const; // harmonic=2, ...
        size_t get_number_of_irred_fc_elements(const int fc_order); // harmonic=2, ...

        void get_fc_origin(double *fc_values,
                           int *elem_indices,
                           // (len(fc_value), fc_order) is flatten.
                           int fc_order) const; // harmonic=2, ...


        void get_fc_irreducible(double *fc_values,
                                int *elem_indices,
                                // (len(fc_value), fc_order) is flatten.
                                int fc_order); // harmonic=2, ...


        void get_fc_all(double *fc_values,
                        int *elem_indices,
                        // (len(fc_value), fc_order) is flatten.
                        int fc_order) const; // harmonic=2, ...

        void set_fc(double *fc_in) const;

        void get_matrix_elements(double *amat,
                                 double *bvec) const;
        void generate_force_constant();
        int run_optimize();
        void run_suggest() const;
        void run();

    private:
        class System *system{};

        std::string run_mode;
        int verbosity;

        bool structure_initialized;
        bool ready_to_fit;
        std::ofstream *ofs_alm;
        std::streambuf *coutbuf;
        void create();
        void initialize_structure();
        void initialize_interaction();
    };
}
