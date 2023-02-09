/*
 input_setter.cpp

 Copyright (c) 2014, 2015, 2016 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory
 or http://opensource.org/licenses/mit-license.php for information.
*/

#include <string>
#include "input_setter.h"
#include "memory.h"
#include "files.h"
#include "symmetry.h"
#include "optimize.h"
#include "constraint.h"
#include "patterndisp.h"
#include "alm.h"
#include "error.h"

using namespace ALM_NS;

InputSetter::InputSetter()
{
    nat_base = 0;
    nkd = 0;
    maxorder = 0;
    //kd_base = nullptr;
    //kdname = nullptr;

//    for (auto i = 0; i < 3; ++i) {
//        for (auto j = 0; j < 3; j++) {
//            lavec_base[i][j] = 0.0;
//        }
//    }
    //xcoord_base = nullptr;
    is_periodic[0] = 1;
    is_periodic[1] = 1;
    is_periodic[2] = 1;

    lspin = false;
    //magmom_base = nullptr;
    noncollinear = 0;
    trevsym = 1;
    str_magmom = "";

    nbody_include = nullptr;
    cutoff_radii = nullptr;
}

InputSetter::~InputSetter()
{
//    if (kdname) {
//        deallocate(kdname);
//    }
//    if (xcoord_base) {
//        deallocate(xcoord_base);
//    }
//    if (kd_base) {
//        deallocate(kd_base);
//    }
//    if (magmom_base) {
//        deallocate(magmom_base);
//    }
    if (nbody_include) {
        deallocate(nbody_include);
    }
    if (cutoff_radii) {
        deallocate(cutoff_radii);
    }
}

//void InputSetter::set_cell_parameter(const double a,
//                                     const double lavec_in[3][3])
//{
//    for (auto i = 0; i < 3; ++i) {
//        for (auto j = 0; j < 3; ++j) {
//            lavec_base[i][j] = a * lavec_in[i][j];
//        }
//    }
//}

void InputSetter::set_cell_parameter(const Eigen::Matrix3d &lavec_in)
{
    lavec_base_mat = lavec_in;
}

void InputSetter::set_interaction_vars(const int maxorder_in,
                                       const int *nbody_include_in)
{
    maxorder = maxorder_in;
    if (nbody_include) {
        deallocate(nbody_include);
    }
    allocate(nbody_include, maxorder);
    for (auto i = 0; i < maxorder; i++) {
        nbody_include[i] = nbody_include_in[i];
    }
}

void InputSetter::set_cutoff_radii(const int maxorder_in,
                                   const size_t nkd_in,
                                   const std::vector<double> &cutoff_radii_in)
{
    if (cutoff_radii_in.size() != (nkd_in * nkd_in * maxorder_in)) {
        exit("set_cutoff_radii",
             "Incorrect size of the input array cutoff_radii_in");
    }
    if (cutoff_radii) {
        deallocate(cutoff_radii);
    }
    allocate(cutoff_radii, maxorder_in * nkd_in * nkd_in);
    auto counter = 0;
    for (auto i = 0; i < maxorder_in; i++) {
        for (size_t j = 0; j < nkd_in; j++) {
            for (size_t k = 0; k < nkd_in; k++) {
                cutoff_radii[counter] = cutoff_radii_in[counter];
                ++counter;
            }
        }
    }
}

void InputSetter::set_general_vars(ALM *alm,
                                   const std::string &prefix,
                                   const std::string &mode,
                                   const int verbosity,
                                   const std::string &str_disp_basis,
                                   const int printsymmetry,
                                   const int is_periodic_in[3],
                                   const bool trim_dispsign_for_evenfunc,
                                   const int print_hessian,
                                   const int print_fcs_alamode,
                                   const int print_fc3_shengbte,
                                   const int print_fc2_qefc,
                                   const double tolerance,
                                   const double tolerance_constraint,
                                   const std::string &basis_force_constant,
                                   const int nmaxsave,
                                   const double fc_zero_threshold)
{
    size_t i;

    alm->set_output_filename_prefix(prefix);
    alm->set_verbosity(verbosity);
    alm->set_print_symmetry(printsymmetry);
    alm->set_symmetry_tolerance(tolerance);

    for (i = 0; i < 3; i++) {
        is_periodic[i] = is_periodic_in[i];
    }

    alm->set_fcs_save_flag("hessian", print_hessian);
    alm->set_fcs_save_flag("alamode", print_fcs_alamode);
    alm->set_fcs_save_flag("shengbte", print_fc3_shengbte);
    alm->set_fcs_save_flag("qefc", print_fc2_qefc);
    alm->set_fc_zero_threshold(fc_zero_threshold);
    alm->set_tolerance_constraint(tolerance_constraint);
    alm->set_forceconstant_basis(basis_force_constant);
    alm->set_nmaxsave(nmaxsave);

    if (mode == "suggest") {
        alm->set_displacement_basis(str_disp_basis);
        alm->set_displacement_param(trim_dispsign_for_evenfunc);
    }
}

void InputSetter::define(ALM *alm) const
{
    alm->define(maxorder,
                nkd,
                nbody_include,
                cutoff_radii);
}


void InputSetter::set_optimize_vars(ALM *alm,
                                    const std::vector<std::vector<double>> &u_train_in,
                                    const std::vector<std::vector<double>> &f_train_in,
                                    const std::vector<std::vector<double>> &u_validation_in,
                                    const std::vector<std::vector<double>> &f_validation_in,
                                    const OptimizerControl &optcontrol_in) const
{
    alm->set_u_train(u_train_in);
    alm->set_f_train(f_train_in);
    alm->set_validation_data(u_validation_in, f_validation_in);
    alm->set_optimizer_control(optcontrol_in);
}

void InputSetter::set_file_vars(ALM *alm,
                                const DispForceFile &datfile_train_in,
                                const DispForceFile &datfile_validation_in) const
{
    alm->set_datfile_train(datfile_train_in);
    alm->set_datfile_validation(datfile_validation_in);
}

void InputSetter::set_constraint_vars(ALM *alm,
                                      const int constraint_flag,
                                      const std::string &rotation_axis,
                                      const std::string &fc2_file,
                                      const std::string &fc3_file,
                                      const bool fix_harmonic,
                                      const bool fix_cubic) const
{
    alm->set_constraint_mode(constraint_flag);
    alm->set_rotation_axis(rotation_axis);
    alm->set_fc_file(2, fc2_file);
    alm->set_fc_fix(2, fix_harmonic);
    alm->set_fc_file(3, fc3_file);
    alm->set_fc_fix(3, fix_cubic);
    const auto use_algebraic_constraint = constraint_flag / 10;
    alm->set_algebraic_constraint(use_algebraic_constraint);
}


//void InputSetter::set_atomic_positions(const size_t nat_in,
//                                       const int *kd_in,
//                                       const double (*xcoord_in)[3])
//{
//    if (kd_base) {
//        deallocate(kd_base);
//    }
//    if (xcoord_base) {
//        deallocate(xcoord_base);
//    }
//    allocate(xcoord_base, nat_in);
//    allocate(kd_base, nat_in);
//
//    for (size_t i = 0; i < nat_in; ++i) {
//        kd_base[i] = kd_in[i];
//        for (auto j = 0; j < 3; ++j) {
//            xcoord_base[i][j] = xcoord_in[i][j];
//        }
//    }
//}

void InputSetter::set_atomic_positions(const Eigen::MatrixXd &positions_in)
{
    nat_base = positions_in.col(0).size();
    xcoord_base_mat = positions_in;
}

void InputSetter::set_element_types(const std::vector<int> &kd_in,
                                    const std::vector<std::string> &kdnames_in)
{
    kd_base_vec = kd_in;
    kdnames_vec = kdnames_in;
}


void InputSetter::set_transformation_matrices(const Eigen::Matrix3d &transmat_super_in,
                                              const Eigen::Matrix3d &transmat_prim_in)
{
    transmat_super = transmat_super_in;
    transmat_prim = transmat_prim_in;
}

void InputSetter::set_magnetic_vars(const int lspin_in,
                                    const Eigen::MatrixXd &magmom_in,
                                    const int noncollinear_in,
                                    const int time_reversal_symm_in)
{
    lspin = lspin_in;
    magmom_base_mat = magmom_in;
    noncollinear = noncollinear_in;
    trevsym = time_reversal_symm_in;
}


void InputSetter::set_geometric_structure(ALM *alm)
{
    double (*xcoord_base)[3]; // fractional coordinate
    double (*magmom_base)[3];
    int *kd_base;
    double lavec_base[3][3];

    allocate(xcoord_base, nat_base);
    allocate(magmom_base, nat_base);
    allocate(kd_base, nat_base);

    for (auto i = 0; i < 3; ++i) {
        for (auto j = 0; j < 3; ++j) {
            lavec_base[i][j] = lavec_base_mat(i, j);
        }
    }

    for (auto i = 0; i < nat_base; ++i) {
        for (auto j = 0; j < 3; ++j) {
            xcoord_base[i][j] = xcoord_base_mat(i, j);
            magmom_base[i][j] = magmom_base_mat(i, j);
        }
        kd_base[i] = kd_base_vec[i];
    }

    alm->set_cell(nat_base, lavec_base, xcoord_base, kd_base, kdnames_vec.data());
    alm->set_periodicity(is_periodic);

    deallocate(xcoord_base);
    deallocate(kd_base);

    double transmat_super_tmp[3][3], transmat_prim_tmp[3][3];

    for (auto i = 0; i < 3; ++i) {
        for (auto j = 0 ; j <3 ; ++j) {
            transmat_super_tmp[i][j] = transmat_super(i,j);
            transmat_prim_tmp[i][j] = transmat_prim(i,j);
        }
    }
    alm->set_transformation_matrices(transmat_super_tmp,
                                     transmat_prim_tmp);

    alm->set_magnetic_params(nat_base, magmom_base, lspin, noncollinear, trevsym, str_magmom);
    deallocate(magmom_base);
}
