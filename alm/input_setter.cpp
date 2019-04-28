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
    nat = 0;
    nkd = 0;
    maxorder = 0;
    kd = nullptr;
    kdname = nullptr;

    for (auto i = 0; i < 3; ++i) {
        for (auto j = 0; j < 3; j++) {
            lavec[i][j] = 0.0;
        }
    }
    xcoord = nullptr;
    is_periodic[0] = 1;
    is_periodic[1] = 1;
    is_periodic[2] = 1;

    lspin = false;
    magmom = nullptr;
    noncollinear = 0;
    trevsym = 1;
    str_magmom = "";

    nbody_include = nullptr;
    cutoff_radii = nullptr;
}

InputSetter::~InputSetter()
{
    if (kdname) {
        deallocate(kdname);
    }
    if (xcoord) {
        deallocate(xcoord);
    }
    if (kd) {
        deallocate(kd);
    }
    if (magmom) {
        deallocate(magmom);
    }
    if (nbody_include) {
        deallocate(nbody_include);
    }
    if (cutoff_radii) {
        deallocate(cutoff_radii);
    }
}

void InputSetter::set_cell_parameter(const double a,
                                     const double lavec_in[3][3])
{
    for (auto i = 0; i < 3; ++i) {
        for (auto j = 0; j < 3; ++j) {
            lavec[i][j] = a * lavec_in[i][j];
        }
    }
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
                                   const std::string prefix,
                                   const std::string mode,
                                   const int verbosity,
                                   const std::string str_disp_basis,
                                   const std::string str_magmom,
                                   const size_t nat_in,
                                   const size_t nkd_in,
                                   const int printsymmetry,
                                   const int is_periodic_in[3],
                                   const bool trim_dispsign_for_evenfunc,
                                   const bool lspin_in,
                                   const bool print_hessian,
                                   const int noncollinear_in,
                                   const int trevsym_in,
                                   const std::string *kdname_in,
                                   const double *const *magmom_in,
                                   const double tolerance,
                                   const double tolerance_constraint)
{
    size_t i;

    alm->files->set_prefix(prefix);
    alm->set_run_mode(mode);
    alm->set_verbosity(verbosity);
    nat = nat_in;
    nkd = nkd_in;
    alm->symmetry->set_print_symmetry(printsymmetry);
    alm->symmetry->set_tolerance(tolerance);

    if (kdname) {
        deallocate(kdname);
    }
    allocate(kdname, nkd);
    for (i = 0; i < nkd; ++i) {
        kdname[i] = kdname_in[i];
    }

    if (magmom) {
        deallocate(magmom);
    }
    allocate(magmom, nat);

    for (i = 0; i < nat; i ++) {
        for (auto j = 0; j < 3; j++) {
            magmom[i][j] = magmom_in[i][j];
        }
    }
    lspin = lspin_in;
    noncollinear = noncollinear_in;
    trevsym = trevsym_in;

    for (i = 0; i < 3; i++) {
        is_periodic[i] = is_periodic_in[i];
    }

    alm->files->print_hessian = print_hessian;
    alm->constraint->set_tolerance_constraint(tolerance_constraint);

    if (mode == "suggest") {
        alm->displace->set_disp_basis(str_disp_basis);
        alm->displace->set_trim_dispsign_for_evenfunc(trim_dispsign_for_evenfunc);
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
    alm->optimize->set_training_data(u_train_in, f_train_in);
    alm->optimize->set_validation_data(u_validation_in, f_validation_in);
    alm->optimize->set_optimizer_control(optcontrol_in);
}

void InputSetter::set_file_vars(ALM *alm,
                                const DispForceFile &datfile_train_in,
                                const DispForceFile &datfile_validation_in) const
{
    alm->files->set_datfile_train(datfile_train_in);
    alm->files->set_datfile_validation(datfile_validation_in);
}

void InputSetter::set_constraint_vars(ALM *alm,
                                      const int constraint_flag,
                                      const std::string rotation_axis,
                                      const std::string fc2_file,
                                      const std::string fc3_file,
                                      const bool fix_harmonic,
                                      const bool fix_cubic) const
{
    alm->constraint->set_constraint_mode(constraint_flag);
    alm->constraint->set_rotation_axis(rotation_axis);
    alm->constraint->set_fc_file(2, fc2_file);
    alm->constraint->set_fix_harmonic(fix_harmonic);
    alm->constraint->set_fc_file(3, fc3_file);
    alm->constraint->set_fix_cubic(fix_cubic);
}


void InputSetter::set_atomic_positions(const size_t nat_in,
                                       const int *kd_in,
                                       const double (*xcoord_in)[3])
{
    if (kd) {
        deallocate(kd);
    }
    if (xcoord) {
        deallocate(xcoord);
    }
    allocate(xcoord, nat_in);
    allocate(kd, nat_in);

    for (size_t i = 0; i < nat_in; ++i) {
        kd[i] = kd_in[i];
        for (auto j = 0; j < 3; ++j) {
            xcoord[i][j] = xcoord_in[i][j];
        }
    }
}

void InputSetter::set_geometric_structure(ALM *alm)
{
    alm->set_cell(nat, lavec, xcoord, kd, kdname);
    alm->set_periodicity(is_periodic);
    alm->set_magnetic_params(nat, magmom, lspin, noncollinear, trevsym, str_magmom);
}
