/*
 input_setter.h

 Copyright (c) 2014, 2015, 2016 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory
 or http://opensource.org/licenses/mit-license.php for information.
*/

#pragma once

#include "alm.h"
#include <string>

namespace ALM_NS {
    class InputSetter {
    public:
        InputSetter();

        ~InputSetter();

        void set_cell_parameter(const double a,
                                const double lavec_in[3][3]);

        void set_atomic_positions(const size_t nat_in,
                                  const int *kd_in,
                                  const double (*xcoord_in)[3]);

        void set_geometric_structure(ALM *alm);

        void set_interaction_vars(const int maxorder_in,
                                  const int *nbody_include_in);

        void set_cutoff_radii(const int maxorder_in,
                              const size_t nkd_in,
                              const std::vector<double> &cutoff_radii_in);

        void define(ALM *alm) const;

        void set_general_vars(ALM *alm,
                              std::string prefix,
                              std::string mode,
                              int verbosity,
                              std::string str_disp_basis,
                              std::string str_magmom,
                              size_t nat_in,
                              size_t nkd_in,
                              int printsymmetry,
                              const int is_periodic_in[3],
                              bool trim_dispsign_for_evenfunc,
                              bool lspin_in,
                              bool print_hessian,
                              int noncollinear_in,
                              int trevsym_in,
                              const std::string *kdname_in,
                              const double *const *magmom_in,
                              double tolerance,
                              double tolerance_constraint,
                              const std::string basis_force_constant,
                              const int nmaxsave);

        void set_optimize_vars(ALM *alm,
                               const std::vector<std::vector<double>> &u_train_in,
                               const std::vector<std::vector<double>> &f_train_in,
                               const std::vector<std::vector<double>> &u_validation_in,
                               const std::vector<std::vector<double>> &f_validation_in,
                               const OptimizerControl &optcontrol_in) const;

        void set_file_vars(ALM *alm,
                           const DispForceFile &datfile_train_in,
                           const DispForceFile &datfile_validation_in) const;

        void set_constraint_vars(ALM *alm,
                                 int constraint_flag,
                                 std::string rotation_axis,
                                 std::string fc2_file,
                                 std::string fc3_file,
                                 bool fix_harmonic,
                                 bool fix_cubic) const;

        void set_geometric_structure(ALM *alm) const;

    private:
        size_t nat, nkd;
        int *kd;
        double lavec[3][3];
        double (*xcoord)[3]; // fractional coordinate
        std::string *kdname;
        int is_periodic[3];

        bool lspin;
        double (*magmom)[3];
        int noncollinear;
        int trevsym;
        std::string str_magmom;

        int maxorder;
        int *nbody_include;
        double *cutoff_radii;
    };
}
