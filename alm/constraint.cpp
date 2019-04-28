/*
 constraint.cpp

 Copyright (c) 2014, 2015, 2016 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory
 or http://opensource.org/licenses/mit-license.php for information.
*/

#include "constraint.h"
#include "combination.h"
#include "constants.h"
#include "error.h"
#include "fcs.h"
#include "cluster.h"
#include "mathfunctions.h"
#include "memory.h"
#include "rref.h"
#include "symmetry.h"
#include "system.h"
#include "timer.h"
#include "xml_parser.h"
#include <iostream>
#include <iomanip>
#include <boost/bimap.hpp>
#include <algorithm>
#include <unordered_set>
#include <map>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>

using namespace ALM_NS;

Constraint::Constraint()
{
    set_default_variables();
}

Constraint::~Constraint()
{
    deallocate_variables();
}

void Constraint::set_default_variables()
{
    constraint_mode = 11;
    rotation_axis = "";
    fix_harmonic = false;
    fix_cubic = false;
    constraint_algebraic = 1;
    fc2_file = "";
    fc3_file = "";
    exist_constraint = false;
    extra_constraint_from_symmetry = false;
    const_mat = nullptr;
    const_rhs = nullptr;
    const_symmetry = nullptr;
    const_fix = nullptr;
    const_relate = nullptr;
    const_relate_rotation = nullptr;
    index_bimap = nullptr;
    number_of_constraints = 0;
    tolerance_constraint = eps8;
}

void Constraint::deallocate_variables()
{
    if (const_symmetry) {
        deallocate(const_symmetry);
    }
    if (const_fix) {
        deallocate(const_fix);
    }
    if (const_relate) {
        deallocate(const_relate);
    }
    if (const_relate_rotation) {
        deallocate(const_relate_rotation);
    }
    if (index_bimap) {
        deallocate(index_bimap);
    }
    if (const_mat) {
        deallocate(const_mat);
    }
    if (const_rhs) {
        deallocate(const_rhs);
    }
}

void Constraint::setup(const System *system,
                       const Fcs *fcs,
                       const Cluster *cluster,
                       const Symmetry *symmetry,
                       const std::string alm_mode,
                       const int verbosity,
                       Timer *timer)
{
    timer->start_clock("constraint");

    if (verbosity > 0) {
        std::cout << " CONSTRAINT" << std::endl;
        std::cout << " ==========" << std::endl << std::endl;
    }

    constraint_algebraic = constraint_mode / 10;
    constraint_mode = constraint_mode % 10;
    const auto maxorder = cluster->get_maxorder();

    if (alm_mode == "lasso") {
        if (constraint_mode > 1) {
            warn("Constraint::setup()", "Sorry, only ICONST = 11 is supported when MODE = lasso.");
            constraint_mode = 1;
        }
        constraint_algebraic = 1;
    }

    switch (constraint_mode) {
    case 0: // do nothing
        impose_inv_T = false;
        impose_inv_R = false;
        if (verbosity > 0) {
            std::cout << "  ICONST = 0: Constraint for translational/rotational invariance" << std::endl;
            std::cout << "              will NOT be considered." << std::endl;
        }
        break;
    case 1:
        impose_inv_T = true;
        impose_inv_R = false;
        if (verbosity > 0) {
            std::cout << "  ICONST = 1: Constraints for translational invariance" << std::endl;
            std::cout << "              will be considered." << std::endl;
        }
        break;
    case 2:
        impose_inv_T = true;
        impose_inv_R = true;
        exclude_last_R = true;
        if (verbosity > 0) {
            std::cout << "  ICONST = 2: Constraints for translational and rotational invariance" << std::endl;
            std::cout << "              will be considered. Axis of rotation is " << rotation_axis << std::endl;
            std::cout << "              Rotational invariance of the maximum order will be neglected" << std::endl;
        }
        break;
    case 3:
        impose_inv_T = true;
        impose_inv_R = true;
        exclude_last_R = false;
        if (verbosity > 0) {
            std::cout << "  ICONST = 3: Constraints for translational and rotational invariance" << std::endl;
            std::cout << "              will be considered. Axis of rotation is " << rotation_axis << std::endl;
        }
        break;
    default:
        exit("Constraint::setup", "invalid constraint_mode", constraint_mode);
        break;
    }

    if (verbosity > 0) std::cout << std::endl;

    if (const_symmetry) {
        deallocate(const_symmetry);
    }
    if (const_fix) {
        deallocate(const_fix);
    }
    if (const_relate) {
        deallocate(const_relate);
    }
    if (index_bimap) {
        deallocate(index_bimap);
    }

    allocate(const_fix, maxorder);
    allocate(const_relate, maxorder);
    allocate(index_bimap, maxorder);
    allocate(const_symmetry, maxorder);
    allocate(const_translation, maxorder);
    allocate(const_rotation_self, maxorder);
    allocate(const_rotation_cross, maxorder);
    allocate(const_self, maxorder);

    for (auto order = 0; order < maxorder; ++order) {
        const_translation[order].clear();
        const_rotation_self[order].clear();
        const_rotation_cross[order].clear();
        const_self[order].clear();
        const_fix[order].clear();
    }

    if (fix_harmonic) {

        if (verbosity > 0) {
            std::cout << "  FC2XML is given : Harmonic force constants will be " << std::endl;
            std::cout << "                    fixed to the values given in " << fc2_file << std::endl;
            std::cout << std::endl;
        }

        fix_forceconstants_to_file(0,
                                   symmetry,
                                   fcs,
                                   fc2_file,
                                   const_fix[0]);
    }

    fix_cubic = fix_cubic & (cluster->get_maxorder() > 1);
    if (fix_cubic) {

        if (verbosity > 0) {
            std::cout << "  FC3XML is given : Cubic force constants will be " << std::endl;
            std::cout << "                    fixed to the values given in " << fc3_file << std::endl;
            std::cout << std::endl;
        }

        fix_forceconstants_to_file(1,
                                   symmetry,
                                   fcs,
                                   fc3_file,
                                   const_fix[1]);
    }

    generate_symmetry_constraint_in_cartesian(system->get_supercell().number_of_atoms,
                                              symmetry,
                                              cluster,
                                              fcs,
                                              verbosity);


    extra_constraint_from_symmetry = false;

    for (auto order = 0; order < cluster->get_maxorder(); ++order) {
        if (!const_symmetry[order].empty()) extra_constraint_from_symmetry = true;
    }

    if (impose_inv_T) {
        // const_translation is updated.
        generate_translational_constraint(system->get_supercell(),
                                          symmetry,
                                          cluster,
                                          fcs,
                                          verbosity);
    }

    if (impose_inv_R) {
        generate_rotational_constraint(system,
                                       symmetry,
                                       cluster,
                                       fcs,
                                       verbosity,
                                       tolerance_constraint);
    }

    // Merge intra-order constrants and do reduction

    for (auto order = 0; order < maxorder; ++order) {

        const auto nparam = fcs->get_nequiv()[order].size();

        const_self[order].reserve(
            const_translation[order].size()
            + const_rotation_self[order].size()
            + const_symmetry[order].size());

        // The order of const_symmetry and const_translation
        // should not be changed since the rref_sparse is sensitive to
        // the numerical accuracy of the input matrix.
        const_self[order].insert(const_self[order].end(),
                                 const_symmetry[order].begin(),
                                 const_symmetry[order].end());

        const_self[order].insert(const_self[order].end(),
                                 const_translation[order].begin(),
                                 const_translation[order].end());

        const_self[order].insert(const_self[order].end(),
                                 const_rotation_self[order].begin(),
                                 const_rotation_self[order].end());

        rref_sparse(nparam, const_self[order], tolerance_constraint);
    }

    get_mapping_constraint(maxorder,
                           fcs->get_nequiv(),
                           const_self,
                           const_fix,
                           const_relate,
                           index_bimap);

    if (!constraint_algebraic) {

        size_t Pmax = 0;
        size_t nparams = 0;
        for (auto order = 0; order < maxorder; ++order) {
            Pmax += const_self[order].size()
                + const_rotation_cross[order].size();
        }
        if (fix_harmonic) {
            Pmax -= const_self[0].size();
            Pmax += fcs->get_nequiv()[0].size();
        }
        if (fix_cubic) {
            Pmax -= const_self[1].size();
            Pmax += fcs->get_nequiv()[1].size();
        }
        for (auto order = 0; order < maxorder; ++order) {
            nparams += fcs->get_nequiv()[order].size();
        }

        if (const_mat) {
            deallocate(const_mat);
        }
        allocate(const_mat, Pmax, nparams);

        if (const_rhs) {
            deallocate(const_rhs);
        }
        allocate(const_rhs, Pmax);

        // const_mat and const_rhs are updated.
        number_of_constraints = calc_constraint_matrix(maxorder,
                                                       fcs->get_nequiv(),
                                                       nparams);
    }

    exist_constraint
        = impose_inv_T
        || fix_harmonic
        || fix_cubic
        || extra_constraint_from_symmetry;

    if (verbosity > 0) {
        if (exist_constraint) {

            int order;

            if (impose_inv_T || impose_inv_R) {
                std::cout << "  Number of constraints [T-inv, R-inv (self), R-inv (cross)]:" << std::endl;
                for (order = 0; order < maxorder; ++order) {
                    std::cout << "   " << std::setw(8) << cluster->get_ordername(order);
                    std::cout << std::setw(5) << const_translation[order].size();
                    std::cout << std::setw(5) << const_rotation_self[order].size();
                    std::cout << std::setw(5) << const_rotation_cross[order].size();
                    std::cout << std::endl;
                }
                std::cout << std::endl;
            }

            if (extra_constraint_from_symmetry) {
                std::cout << "  There are constraints from crystal symmetry." << std::endl;
                std::cout << "  The number of such constraints for each order:" << std::endl;
                for (order = 0; order < maxorder; ++order) {
                    std::cout << "   " << std::setw(8) << cluster->get_ordername(order);
                    std::cout << std::setw(5) << const_symmetry[order].size();
                    std::cout << std::endl;
                }
                std::cout << std::endl;
            }

            if (extra_constraint_from_symmetry) {
                std::cout << "  Constraints of T-inv, R-inv (self), and those from crystal symmetry are merged."
                    << std::endl;
            } else {
                std::cout << "  Constraints of T-inv and R-inv (self) are merged." << std::endl;
            }
            std::cout << "  If there are redundant constraints, they are removed in this process." << std::endl;
            std::cout << std::endl;
            std::cout << "  Number of inequivalent constraints (self, cross) : " << std::endl;

            for (order = 0; order < maxorder; ++order) {
                std::cout << "   " << std::setw(8) << cluster->get_ordername(order);
                std::cout << std::setw(5) << const_self[order].size();
                std::cout << std::setw(5) << const_rotation_cross[order].size();
                std::cout << std::endl;
            }
            std::cout << std::endl;

            if (constraint_algebraic) {

                std::cout << "  ICONST >= 10 : Constraints will be considered algebraically."
                    << std::endl << std::endl;

                if (impose_inv_R) {
                    std::cout << "  WARNING : Inter-order constraints for rotational invariance will be neglected."
                        << std::endl;
                }

                for (order = 0; order < maxorder; ++order) {
                    std::cout << "  Number of free" << std::setw(9) << cluster->get_ordername(order)
                        << " FCs : " << index_bimap[order].size() << std::endl;
                }
                std::cout << std::endl;

            } else {

                std::cout << "  Total number of constraints = " << number_of_constraints << std::endl << std::endl;

            }
        }
        timer->print_elapsed();
        std::cout << " -------------------------------------------------------------------" << std::endl;
        std::cout << std::endl;
    }


    deallocate(const_translation);
    const_translation = nullptr;
    deallocate(const_rotation_self);
    const_rotation_self = nullptr;
    deallocate(const_rotation_cross);
    const_rotation_cross = nullptr;
    deallocate(const_self);
    const_self = nullptr;

    timer->stop_clock("constraint");
}

size_t Constraint::calc_constraint_matrix(const int maxorder,
                                          const std::vector<size_t> *nequiv,
                                          const size_t nparams) const
{
    size_t i, j;
    int order;
    double *arr_tmp;
    std::vector<ConstraintClass> const_total;

    const_total.clear();
    allocate(arr_tmp, nparams);

    size_t nshift = 0;

    for (order = 0; order < maxorder; ++order) {
        const auto nelems = nequiv[order].size();

        if (const_fix[order].empty()) {
            for (auto &p : const_self[order]) {
                for (i = 0; i < nparams; ++i) arr_tmp[i] = 0.0;
                for (const auto &it : p) {
                    arr_tmp[nshift + it.first] = it.second;
                }
                const_total.emplace_back(nparams, arr_tmp);
            }
        }
        nshift += nelems;
    }


    const auto nconst1 = const_total.size();

    // Inter-order constraints
    size_t nshift2 = 0;
    for (order = 0; order < maxorder; ++order) {
        if (order > 0) {
            if (const_fix[order - 1].empty() && const_fix[order].empty()) {
                for (auto &p : const_rotation_cross[order]) {
                    for (i = 0; i < nparams; ++i) arr_tmp[i] = 0.0;
                    for (const auto &it : p) {
                        arr_tmp[nshift2 + it.first] = it.second;
                    }
                    const_total.emplace_back(nparams, arr_tmp);
                }
            }

            nshift2 += nequiv[order - 1].size();
        }
    }
    deallocate(arr_tmp);

    if (nconst1 != const_total.size())
        remove_redundant_rows(nparams, const_total, tolerance_constraint);

    auto nconst = const_total.size();

    if (fix_harmonic) nconst += nequiv[0].size();
    if (fix_cubic) nconst += nequiv[1].size();

    for (i = 0; i < nconst; ++i) {
        for (j = 0; j < nparams; ++j) {
            const_mat[i][j] = 0.0;
        }
        const_rhs[i] = 0.0;
    }

    size_t irow = 0;
    size_t icol = 0;
    size_t ishift = 0;

    if (fix_harmonic) {

        for (const auto &p : const_fix[0]) {
            i = p.p_index_target;
            const_mat[i][i] = 1.0;
            const_rhs[i] = p.val_to_fix;
        }

        irow += const_fix[0].size();
        icol += const_fix[0].size();
        ishift += const_fix[0].size();
    }

    if (fix_cubic && maxorder > 1) {

        const auto ishift2 = nequiv[0].size();

        for (const auto &p : const_fix[1]) {
            i = p.p_index_target;
            const_mat[i + ishift][i + ishift2] = 1.0;
            const_rhs[i + ishift] = p.val_to_fix;
        }

        irow += const_fix[1].size();
        icol += const_fix[1].size();
    }

    for (auto &p : const_total) {
        for (i = 0; i < nparams; ++i) {
            const_mat[irow][i] = p.w_const[i];
        }
        ++irow;
    }
    const_total.clear();

    return nconst;
}


void Constraint::get_mapping_constraint(const int nmax,
                                        const std::vector<size_t> *nequiv,
                                        const ConstraintSparseForm *const_in,
                                        std::vector<ConstraintTypeFix> *const_fix_out,
                                        std::vector<ConstraintTypeRelate> *const_relate_out,
                                        boost::bimap<size_t, size_t> *index_bimap_out) const
{
    // If const_fix_out[order] is not empty as input, it assumes that fix_forceconstant[order] is true.
    // In this case, const_fix_out[order] is not updated.

    int order;
    size_t i;

    for (order = 0; order < nmax; ++order) {

        if (const_fix_out[order].empty()) {

            size_t p_index_target;
            std::vector<double> alpha_tmp;
            std::vector<size_t> p_index_tmp;

            for (auto p = const_in[order].rbegin(); p != const_in[order].rend(); ++p) {

                alpha_tmp.clear();
                p_index_tmp.clear();
                auto counter = 0;
                for (const auto &p2 : (*p)) {
                    if (counter == 0) {
                        p_index_target = p2.first;
                    } else {
                        alpha_tmp.push_back(p2.second);
                        p_index_tmp.push_back(p2.first);
                    }
                    ++counter;
                }

                if (!alpha_tmp.empty()) {
                    const_relate_out[order].emplace_back(p_index_target,
                                                         alpha_tmp, p_index_tmp);
                } else {
                    const_fix_out[order].emplace_back(p_index_target, 0.0);
                }
            }
        }
    }

    std::vector<int> *has_constraint;
    allocate(has_constraint, nmax);
    size_t nparam;
    for (order = 0; order < nmax; ++order) {

        nparam = nequiv[order].size();
        has_constraint[order].resize(nparam, 0);

        for (i = 0; i < const_fix_out[order].size(); ++i) {
            has_constraint[order][const_fix_out[order][i].p_index_target] = 1;
        }

        for (i = 0; i < const_relate_out[order].size(); ++i) {
            has_constraint[order][const_relate_out[order][i].p_index_target] = 2;
        }
    }

    size_t icount;

    for (order = 0; order < nmax; ++order) {
        nparam = nequiv[order].size();

        icount = 0;
        for (i = 0; i < nparam; ++i) {
            if (has_constraint[order][i] == 0) {
                index_bimap_out[order].insert(
                    boost::bimap<size_t, size_t>::value_type(icount, i));
                ++icount;
            }
        }
    }

    deallocate(has_constraint);
}

int Constraint::get_constraint_mode() const
{
    return constraint_mode;
}

void Constraint::set_constraint_mode(const int constraint_mode_in)
{
    constraint_mode = constraint_mode_in;
}

size_t Constraint::get_number_of_constraints() const
{
    return number_of_constraints;
}

std::string Constraint::get_fc_file(const int order) const
{
    switch (order) {
    case 2:
        return fc2_file;
    case 3:
        return fc3_file;
    default:
        return "";
    }
}

void Constraint::set_fc_file(const int order,
                             const std::string fc_file)
{
    switch (order) {
    case 2:
        fc2_file = fc_file;
        break;
    case 3:
        fc3_file = fc_file;
        break;
    default:
        break;
    }
}

bool Constraint::get_fix_harmonic() const
{
    return fix_harmonic;
}

void Constraint::set_fix_harmonic(const bool fix_harmonic_in)
{
    fix_harmonic = fix_harmonic_in;
}

bool Constraint::get_fix_cubic() const
{
    return fix_cubic;
}

void Constraint::set_fix_cubic(const bool fix_cubic_in)
{
    fix_cubic = fix_cubic_in;
}

int Constraint::get_constraint_algebraic() const
{
    return constraint_algebraic;
}

double** Constraint::get_const_mat() const
{
    return const_mat;
}

double* Constraint::get_const_rhs() const
{
    return const_rhs;
}

double Constraint::get_tolerance_constraint() const
{
    return tolerance_constraint;
}

void Constraint::set_tolerance_constraint(const double tol)
{
    tolerance_constraint = tol;
}

bool Constraint::get_exist_constraint() const
{
    return exist_constraint;
}

bool Constraint::get_extra_constraint_from_symmetry() const
{
    return extra_constraint_from_symmetry;
}

std::string Constraint::get_rotation_axis() const
{
    return rotation_axis;
}

void Constraint::set_rotation_axis(const std::string rotation_axis_in)
{
    rotation_axis = rotation_axis_in;
}

const ConstraintSparseForm& Constraint::get_const_symmetry(const int order) const
{
    return const_symmetry[order];
}

const std::vector<ConstraintTypeFix>& Constraint::get_const_fix(const int order) const
{
    return const_fix[order];
}

void Constraint::set_const_fix_val_to_fix(const int order,
                                          const size_t idx,
                                          const double val)
{
    const_fix[order][idx].val_to_fix = val;
}

const std::vector<ConstraintTypeRelate>& Constraint::get_const_relate(const int order) const
{
    return const_relate[order];
}

const boost::bimap<size_t, size_t>& Constraint::get_index_bimap(const int order) const
{
    return index_bimap[order];
}

void Constraint::generate_symmetry_constraint_in_cartesian(const size_t nat,
                                                           const Symmetry *symmetry,
                                                           const Cluster *cluster,
                                                           const Fcs *fcs,
                                                           const int verbosity) const
{
    // Create constraint matrices arising from the crystal symmetry.

    const auto maxorder = cluster->get_maxorder();
    auto has_constraint_from_symm = false;
    std::vector<std::vector<double>> const_tmp;

    for (auto isym = 0; isym < symmetry->get_nsym(); ++isym) {
        if (!symmetry->get_SymmData()[isym].compatible_with_cartesian) {
            has_constraint_from_symm = true;
            break;
        }
    }

    has_constraint_from_symm = has_constraint_from_symm & (verbosity > 0);

    if (has_constraint_from_symm) {
        std::cout << "  Generating constraints from crystal symmetry ..." << std::endl;
    }

    for (auto order = 0; order < maxorder; ++order) {
        if (has_constraint_from_symm) {
            std::cout << "   " << std::setw(8) << cluster->get_ordername(order) << " ...";
        }

        fcs->get_constraint_symmetry(nat,
                                     symmetry,
                                     order,
                                     "Cartesian",
                                     fcs->get_fc_table()[order],
                                     fcs->get_nequiv()[order].size(),
                                     tolerance_constraint,
                                     const_symmetry[order], true);

        if (has_constraint_from_symm) {
            std::cout << " done." << std::endl;
        }
    }

    if (has_constraint_from_symm) {
        std::cout << "  Finished !" << std::endl << std::endl;
    }
}


void Constraint::generate_translational_constraint(const Cell &supercell,
                                                   const Symmetry *symmetry,
                                                   const Cluster *cluster,
                                                   const Fcs *fcs,
                                                   const int verbosity) const
{
    // Create constraint matrix for the translational invariance (aka acoustic sum rule).

    if (verbosity > 0) {
        std::cout << "  Generating constraints for translational invariance ..." << std::endl;
    }

    for (auto order = 0; order < cluster->get_maxorder(); ++order) {

        if (verbosity > 0)
            std::cout << "   " << std::setw(8) << cluster->get_ordername(order) << " ...";

        const auto nparams = fcs->get_nequiv()[order].size();

        if (nparams == 0) {
            if (verbosity > 0) std::cout << "  No parameters! Skipped." << std::endl;
            continue;
        }

        get_constraint_translation(supercell,
                                   symmetry,
                                   cluster,
                                   fcs,
                                   order,
                                   fcs->get_fc_table()[order],
                                   fcs->get_nequiv()[order].size(),
                                   const_translation[order], true);

        if (verbosity > 0) std::cout << " done." << std::endl;
    }

    if (verbosity > 0) std::cout << "  Finished !" << std::endl << std::endl;
}


void Constraint::get_constraint_translation(const Cell &supercell,
                                            const Symmetry *symmetry,
                                            const Cluster *cluster,
                                            const Fcs *fcs,
                                            const int order,
                                            const std::vector<FcProperty> &fc_table,
                                            const size_t nparams,
                                            ConstraintSparseForm &const_out,
                                            const bool do_rref) const
{
    // Generate equality constraint for the acoustic sum rule.

    int i, j;
    int iat, jat, icrd, jcrd;
    int idata;
    int loc_nonzero;

    int *ind;
    int *intarr, *intarr_copy;
    int **xyzcomponent;

    int ixyz;
    const auto natmin = symmetry->get_nat_prim();
    const auto nat = supercell.number_of_atoms;

    unsigned int isize;

    std::vector<int> data;
    std::unordered_set<FcProperty> list_found;
    std::unordered_set<FcProperty>::iterator iter_found;
    std::vector<std::vector<int>> data_vec;
    std::vector<FcProperty> list_vec;
    std::vector<int> const_now;

    typedef std::vector<ConstraintIntegerElement> ConstEntry;
    std::vector<ConstEntry> constraint_all;

    ConstEntry const_tmp;

    if (order < 0) return;

    if (nparams == 0) return;

    allocate(ind, order + 2);

    // Create force constant table for search

    list_found.clear();

    for (const auto &p : fc_table) {
        for (i = 0; i < order + 2; ++i) {
            ind[i] = p.elems[i];
        }
        if (list_found.find(FcProperty(order + 2, p.sign,
                                       ind, p.mother)) != list_found.end()) {
            exit("get_constraint_translation", "Duplicate interaction list found");
        }
        list_found.insert(FcProperty(order + 2, p.sign,
                                     ind, p.mother));
    }

    deallocate(ind);

    // Generate xyz component for each order

    const auto nxyz = static_cast<int>(std::pow(static_cast<double>(3), order + 1));
    allocate(xyzcomponent, nxyz, order + 1);
    fcs->get_xyzcomponent(order + 1, xyzcomponent);

    allocate(intarr, order + 2);
    allocate(intarr_copy, order + 2);

    const_now.resize(nparams);

    for (i = 0; i < natmin; ++i) {

        iat = symmetry->get_map_p2s()[i][0];

        // Generate atom pairs for each order

        if (order == 0) {


            for (icrd = 0; icrd < 3; ++icrd) {

                intarr[0] = 3 * iat + icrd;

                for (jcrd = 0; jcrd < 3; ++jcrd) {

                    // Reset the temporary array for another constraint
                    for (j = 0; j < nparams; ++j) const_now[j] = 0;

                    for (jat = 0; jat < 3 * nat; jat += 3) {
                        intarr[1] = jat + jcrd;

                        iter_found = list_found.find(FcProperty(order + 2, 1.0,
                                                                intarr, 1));

                        //  If found an IFC
                        if (iter_found != list_found.end()) {
                            // Round the coefficient to integer
                            const_now[(*iter_found).mother] += nint((*iter_found).sign);
                        }

                    }
                    // Add to the constraint list
                    if (!is_allzero(const_now, loc_nonzero)) {
                        if (const_now[loc_nonzero] < 0) {
                            for (j = 0; j < nparams; ++j) const_now[j] *= -1;
                        }
                        const_tmp.clear();
                        for (j = 0; j < nparams; ++j) {
                            if (std::abs(const_now[j]) > 0) {
                                const_tmp.emplace_back(j, const_now[j]);
                            }
                        }
                        constraint_all.emplace_back(const_tmp);
                    }
                }
            }

        } else {

            // Anharmonic cases

            auto intlist(cluster->get_interaction_pair(order, i));
            std::sort(intlist.begin(), intlist.end());

            data_vec.clear();
            // Generate data_vec that contains possible interacting clusters.
            // Each cluster contains (order + 1) atoms, and the last atom index
            // will be treated seperately below.
            CombinationWithRepetition<int> g2(intlist.begin(), intlist.end(), order);
            do {
                data = g2.now();

                intarr[0] = iat;

                for (isize = 0; isize < data.size(); ++isize) {
                    intarr[isize + 1] = data[isize];
                }

                if (cluster->satisfy_nbody_rule(order + 1, intarr, order)) {
                    if (cluster->is_incutoff(order + 1, intarr, order, supercell.kind)) {
                        // Add to list if the atoms interact with each other.
                        data_vec.push_back(data);
                    }
                }

            } while (g2.next());

            const auto ndata = data_vec.size();

            // Use openmp for acceleration if possible
#ifdef _OPENMP
#pragma omp parallel
#endif
            {
                int *intarr_omp, *intarr_copy_omp;

                allocate(intarr_omp, order + 2);
                allocate(intarr_copy_omp, order + 2);

                std::vector<int> data_omp;
                std::vector<int> const_now_omp;

                ConstEntry const_tmp_omp;
                std::vector<ConstEntry> constraint_list_omp;

                const_now_omp.resize(nparams);
#ifdef _OPENMP
#pragma omp for private(isize, ixyz, jcrd, j, jat, iter_found, loc_nonzero), schedule(guided), nowait
#endif
                for (idata = 0; idata < ndata; ++idata) {

                    data_omp = data_vec[idata];

                    intarr_omp[0] = iat;
                    for (isize = 0; isize < data_omp.size(); ++isize) {
                        intarr_omp[isize + 1] = data_omp[isize];
                    }

                    // Loop for xyz component
                    for (ixyz = 0; ixyz < nxyz; ++ixyz) {
                        // Loop for the xyz index of the last atom
                        for (jcrd = 0; jcrd < 3; ++jcrd) {

                            // Reset the temporary array for another constraint
                            for (j = 0; j < nparams; ++j) const_now_omp[j] = 0;

                            // Loop for the last atom index
                            for (jat = 0; jat < 3 * nat; jat += 3) {
                                intarr_omp[order + 1] = jat / 3;

                                if (cluster->satisfy_nbody_rule(order + 2, intarr_omp, order)) {
                                    for (j = 0; j < order + 1; ++j) {
                                        intarr_copy_omp[j] = 3 * intarr_omp[j] + xyzcomponent[ixyz][j];
                                    }
                                    intarr_copy_omp[order + 1] = jat + jcrd;

                                    sort_tail(order + 2, intarr_copy_omp);

                                    iter_found = list_found.find(FcProperty(order + 2, 1.0,
                                                                            intarr_copy_omp, 1));
                                    if (iter_found != list_found.end()) {
                                        const_now_omp[(*iter_found).mother] += nint((*iter_found).sign);
                                    }

                                }
                            } // close loop jat

                            // Add the constraint to the private array
                            if (!is_allzero(const_now_omp, loc_nonzero)) {
                                if (const_now_omp[loc_nonzero] < 0) {
                                    for (j = 0; j < nparams; ++j) const_now_omp[j] *= -1;
                                }

                                const_tmp_omp.clear();
                                for (j = 0; j < nparams; ++j) {
                                    if (std::abs(const_now_omp[j]) > 0) {
                                        const_tmp_omp.emplace_back(j, const_now_omp[j]);
                                    }
                                }
                                if (const_tmp_omp.empty()) {
                                    std::cout << "This cannot happen" << std::endl;
                                }
                                constraint_list_omp.emplace_back(const_tmp_omp);
                            }
                        }
                    }

                } // close idata (openmp main loop)

                deallocate(intarr_omp);
                deallocate(intarr_copy_omp);

                // Merge vectors
#pragma omp critical
                {
                    for (const auto &it : constraint_list_omp) {
                        constraint_all.emplace_back(it);
                    }
                }
                constraint_list_omp.clear();
            } // close openmp

            intlist.clear();
        } // close if
    }     // close loop i

    deallocate(xyzcomponent);
    deallocate(intarr);
    deallocate(intarr_copy);

    std::sort(constraint_all.begin(), constraint_all.end());
    constraint_all.erase(std::unique(constraint_all.begin(),
                                     constraint_all.end()),
                         constraint_all.end());

    typedef std::map<size_t, double> ConstDoubleEntry;
    ConstDoubleEntry const_tmp2;
    auto division_factor = 1.0;
    int counter;
    const_out.clear();

    for (const auto &it : constraint_all) {
        const_tmp2.clear();
        counter = 0;
        for (const auto &it2 : it) {
            if (counter == 0) {
                division_factor = 1.0 / it2.val;
            }
            const_tmp2[it2.col] = it2.val * division_factor;
            ++counter;
        }
        const_out.emplace_back(const_tmp2);
    }
    constraint_all.clear();
    if (do_rref) rref_sparse(nparams, const_out, eps8);
}

void Constraint::generate_rotational_constraint(const System *system,
                                                const Symmetry *symmetry,
                                                const Cluster *cluster,
                                                const Fcs *fcs,
                                                const int verbosity,
                                                const double tolerance)
{
    // Create constraints for the rotational invariance

    if (verbosity > 0)
        std::cout << "  Generating constraints for rotational invariance ..." << std::endl;

#ifdef _DEBUG
    std::ofstream ofs_constraint;
    ofs_constraint.open("CONSTRAINT", std::ios::out);
#endif

    int i, j;
    int iat, jat;
    int icrd, jcrd;
    int order;
    const auto maxorder = cluster->get_maxorder();
    const auto natmin = symmetry->get_nat_prim();
    int mu, nu;
    int ixyz, nxyz{0}, nxyz2;
    int mu_lambda, lambda;
    int levi_factor;

    int *ind;
    int **xyzcomponent = nullptr;
    int **xyzcomponent2 = nullptr;
    size_t *nparams, nparam_sub;
    int *interaction_index, *interaction_atom;
    int *interaction_tmp;
    int loc_nonzero;

    std::vector<double> arr_constraint;
    std::vector<double> arr_constraint_self;
    std::vector<double> arr_constraint_lower;

    bool valid_rotation_axis[3][3];

    double vec_for_rot[3];

    std::vector<int> interaction_list;

    std::unordered_set<FcProperty> list_found;
    std::unordered_set<FcProperty> list_found_last;
    std::unordered_set<FcProperty>::iterator iter_found;

    CombinationWithRepetition<int> g;

    std::vector<int> atom_tmp;
    std::vector<std::vector<int>> cell_dummy;
    std::set<InteractionCluster>::iterator iter_cluster;

    typedef std::vector<ConstraintDoubleElement> ConstEntry;
    ConstEntry const_tmp;
    std::vector<ConstEntry> *const_self_vec, *const_cross_vec;

    allocate(const_self_vec, maxorder);
    allocate(const_cross_vec, maxorder);

    setup_rotation_axis(valid_rotation_axis);

    allocate(ind, maxorder + 1);
    allocate(nparams, maxorder);

    for (order = 0; order < maxorder; ++order) {

        nparams[order] = fcs->get_nequiv()[order].size();

        if (order == 0) {
            if (verbosity > 0) {
                std::cout << "   Constraints between " << std::setw(8)
                    << "1st-order IFCs (which are zero) and "
                    << std::setw(8) << cluster->get_ordername(order) << " ...";
            }

            nparam_sub = nparams[order];
        } else {
            if (verbosity > 0) {
                std::cout << "   Constraints between " << std::setw(8)
                    << cluster->get_ordername(order - 1) << " and "
                    << std::setw(8) << cluster->get_ordername(order) << " ...";
            }

            nparam_sub = nparams[order] + nparams[order - 1];
        }
        arr_constraint.resize(nparam_sub);
        arr_constraint_self.resize(nparams[order]);

        if (order > 0) arr_constraint_lower.resize(nparams[order - 1]);

        allocate(interaction_atom, order + 2);
        allocate(interaction_index, order + 2);
        allocate(interaction_tmp, order + 2);
        const_self_vec[order].clear();
        const_cross_vec[order].clear();

        if (order > 0) {
            list_found_last = list_found;
            nxyz = static_cast<int>(pow(static_cast<double>(3), order));
            allocate(xyzcomponent, nxyz, order);
            fcs->get_xyzcomponent(order, xyzcomponent);
        }

        list_found.clear();

        for (auto p = fcs->get_fc_table()[order].begin(); p != fcs->get_fc_table()[order].end(); ++p) {
            for (i = 0; i < order + 2; ++i) {
                ind[i] = (*p).elems[i];
            }
            list_found.insert(FcProperty(order + 2, (*p).sign,
                                         ind, (*p).mother));
        }

        for (i = 0; i < natmin; ++i) {

            iat = symmetry->get_map_p2s()[i][0];

            interaction_atom[0] = iat;

            if (order == 0) {

                auto interaction_list_now(cluster->get_interaction_pair(order, i));
                std::sort(interaction_list_now.begin(), interaction_list_now.end());

                // Special treatment for harmonic force constants

                for (icrd = 0; icrd < 3; ++icrd) {

                    interaction_index[0] = 3 * iat + icrd;

                    for (mu = 0; mu < 3; ++mu) {

                        for (nu = 0; nu < 3; ++nu) {

                            if (!valid_rotation_axis[mu][nu]) continue;

                            // Clear history

                            for (j = 0; j < nparam_sub; ++j) arr_constraint[j] = 0.0;

                            for (auto &iter_list : interaction_list_now) {

                                jat = iter_list;
                                interaction_index[1] = 3 * jat + mu;
                                iter_found = list_found.find(FcProperty(order + 2, 1.0,
                                                                        interaction_index, 1));

                                atom_tmp.clear();
                                atom_tmp.push_back(jat);
                                cell_dummy.clear();
                                iter_cluster = cluster->get_interaction_cluster(order, i).find(
                                    InteractionCluster(atom_tmp, cell_dummy));

                                if (iter_cluster == cluster->get_interaction_cluster(order, i).end()) {
                                    exit("generate_rotational_constraint",
                                         "cluster not found ...");
                                } else {
                                    for (j = 0; j < 3; ++j) vec_for_rot[j] = 0.0;

                                    const auto nsize_equiv = (*iter_cluster).cell.size();

                                    for (j = 0; j < nsize_equiv; ++j) {
                                        for (auto k = 0; k < 3; ++k) {
                                            vec_for_rot[k]
                                                += system->get_x_image()[(*iter_cluster).cell[j][0]][jat][k];
                                        }
                                    }

                                    for (j = 0; j < 3; ++j) {
                                        vec_for_rot[j] /= static_cast<double>(nsize_equiv);
                                    }
                                }


                                if (iter_found != list_found.end()) {
                                    arr_constraint[(*iter_found).mother] += (*iter_found).sign * vec_for_rot[nu];
                                }

                                // Exchange mu <--> nu and repeat again.
                                // Note that the sign is inverted (+ --> -) in the summation

                                interaction_index[1] = 3 * jat + nu;
                                iter_found = list_found.find(FcProperty(order + 2, 1.0,
                                                                        interaction_index, 1));
                                if (iter_found != list_found.end()) {
                                    arr_constraint[(*iter_found).mother]
                                        -= (*iter_found).sign * vec_for_rot[mu];
                                }
                            }

                            if (!is_allzero(arr_constraint, tolerance, loc_nonzero)) {
                                // Add to constraint list
                                if (arr_constraint[loc_nonzero] < 0.0) {
                                    for (j = 0; j < nparam_sub; ++j) arr_constraint[j] *= -1.0;
                                }
                                const_tmp.clear();
                                for (j = 0; j < nparam_sub; ++j) {
                                    if (std::abs(arr_constraint[j]) >= tolerance) {
                                        const_tmp.emplace_back(j, arr_constraint[j]);
                                    }
                                }
                                const_self_vec[order].emplace_back(const_tmp);
                            }

                        } // nu
                    }     // mu
                }
            } else {

                // Constraint between different orders

                auto interaction_list_now(cluster->get_interaction_pair(order, i));
                auto interaction_list_old(cluster->get_interaction_pair(order - 1, i));
                std::sort(interaction_list_now.begin(), interaction_list_now.end());
                std::sort(interaction_list_old.begin(), interaction_list_old.end());

                for (icrd = 0; icrd < 3; ++icrd) {

                    interaction_index[0] = 3 * iat + icrd;

                    const CombinationWithRepetition<int> g_now(interaction_list_now.begin(),
                                                               interaction_list_now.end(), order);
                    const CombinationWithRepetition<int> g_old(interaction_list_old.begin(),
                                                               interaction_list_old.end(), order);

                    // m    -th order --> (m-1)-th order
                    // (m-1)-th order -->     m-th order
                    // 2-different directions to find all constraints

                    for (unsigned int direction = 0; direction < 2; ++direction) {

                        if (direction == 0) {
                            g = g_now;
                            interaction_list = interaction_list_now;
                        } else {
                            g = g_old;
                            interaction_list = interaction_list_old;
                        }

                        // Loop for the interacting pairs

                        do {
                            auto data = g.now();

                            for (size_t idata = 0; idata < data.size(); ++idata) {
                                interaction_atom[idata + 1] = data[idata];
                            }

                            for (ixyz = 0; ixyz < nxyz; ++ixyz) {

                                for (j = 0; j < order; ++j)
                                    interaction_index[j + 1] = 3 * interaction_atom[j + 1] + xyzcomponent[ixyz][j];

                                for (mu = 0; mu < 3; ++mu) {

                                    for (nu = 0; nu < 3; ++nu) {

                                        if (!valid_rotation_axis[mu][nu]) continue;

                                        // Search for a new constraint below

                                        for (j = 0; j < nparam_sub; ++j) arr_constraint[j] = 0.0;

                                        // Loop for m_{N+1}, a_{N+1}
                                        for (auto &iter_list : interaction_list) {
                                            jat = iter_list;

                                            interaction_atom[order + 1] = jat;
                                            if (!cluster->is_incutoff(order + 2,
                                                                      interaction_atom,
                                                                      order,
                                                                      system->get_supercell().kind))
                                                continue;

                                            atom_tmp.clear();

                                            for (j = 1; j < order + 2; ++j) {
                                                atom_tmp.push_back(interaction_atom[j]);
                                            }
                                            std::sort(atom_tmp.begin(), atom_tmp.end());

                                            iter_cluster = cluster->get_interaction_cluster(order, i).find(
                                                InteractionCluster(atom_tmp,
                                                                   cell_dummy));
                                            if (iter_cluster != cluster->get_interaction_cluster(order, i).end()) {

                                                int iloc = -1;

                                                for (j = 0; j < atom_tmp.size(); ++j) {
                                                    if (atom_tmp[j] == jat) {
                                                        iloc = j;
                                                        break;
                                                    }
                                                }

                                                if (iloc == -1) {
                                                    exit("generate_rotational_constraint", "This cannot happen.");
                                                }

                                                for (j = 0; j < 3; ++j) vec_for_rot[j] = 0.0;

                                                const auto nsize_equiv = (*iter_cluster).cell.size();

                                                for (j = 0; j < nsize_equiv; ++j) {
                                                    for (auto k = 0; k < 3; ++k) {
                                                        vec_for_rot[k] += system->get_x_image()[(*iter_cluster).cell[j][
                                                            iloc]][jat][k];
                                                    }
                                                }

                                                for (j = 0; j < 3; ++j) {
                                                    vec_for_rot[j] /= static_cast<double>(nsize_equiv);
                                                }
                                            }


                                            // mu, nu

                                            interaction_index[order + 1] = 3 * jat + mu;
                                            for (j = 0; j < order + 2; ++j) interaction_tmp[j] = interaction_index[j];

                                            sort_tail(order + 2, interaction_tmp);

                                            iter_found = list_found.
                                                find(FcProperty(order + 2, 1.0, interaction_tmp, 1));
                                            if (iter_found != list_found.end()) {
                                                arr_constraint[nparams[order - 1] + (*iter_found).mother]
                                                    += (*iter_found).sign * vec_for_rot[nu];
                                            }

                                            // Exchange mu <--> nu and repeat again.

                                            interaction_index[order + 1] = 3 * jat + nu;
                                            for (j = 0; j < order + 2; ++j) interaction_tmp[j] = interaction_index[j];

                                            sort_tail(order + 2, interaction_tmp);

                                            iter_found = list_found.
                                                find(FcProperty(order + 2, 1.0, interaction_tmp, 1));
                                            if (iter_found != list_found.end()) {
                                                arr_constraint[nparams[order - 1] + (*iter_found).mother]
                                                    -= (*iter_found).sign * vec_for_rot[mu];
                                            }
                                        }

                                        for (lambda = 0; lambda < order + 1; ++lambda) {

                                            mu_lambda = interaction_index[lambda] % 3;

                                            for (jcrd = 0; jcrd < 3; ++jcrd) {

                                                for (j = 0; j < order + 1; ++j)
                                                    interaction_tmp[j] = interaction_index[j
                                                    ];

                                                interaction_tmp[lambda] = 3 * interaction_atom[lambda] + jcrd;

                                                levi_factor = 0;

                                                for (j = 0; j < 3; ++j) {
                                                    levi_factor += levi_civita(j, mu, nu) * levi_civita(
                                                        j, mu_lambda, jcrd);
                                                }

                                                if (levi_factor == 0) continue;

                                                sort_tail(order + 1, interaction_tmp);

                                                iter_found = list_found_last.find(FcProperty(order + 1, 1.0,
                                                                                             interaction_tmp, 1));
                                                if (iter_found != list_found_last.end()) {
                                                    arr_constraint[(*iter_found).mother]
                                                        += (*iter_found).sign * static_cast<double>(levi_factor);
                                                }
                                            }
                                        }

                                        if (!is_allzero(arr_constraint, tolerance, loc_nonzero)) {

                                            // A Candidate for another constraint found !
                                            // Add to the appropriate set

                                            if (arr_constraint[loc_nonzero] < 0.0) {
                                                for (j = 0; j < nparam_sub; ++j) arr_constraint[j] *= -1.0;
                                            }
                                            for (j = 0; j < nparams[order]; ++j) {
                                                arr_constraint_self[j] = arr_constraint[j + nparams[order - 1]];
                                            }
                                            for (j = 0; j < nparams[order - 1]; ++j) {
                                                arr_constraint_lower[j] = arr_constraint[j];
                                            }

                                            const_tmp.clear();

                                            if (is_allzero(arr_constraint_self, tolerance, loc_nonzero)) {
                                                // If all elements of the "order"th order is zero,
                                                // the constraint is intraorder of the "order-1"th order.
                                                for (j = 0; j < nparams[order - 1]; ++j) {
                                                    if (std::abs(arr_constraint_lower[j]) >= tolerance) {
                                                        const_tmp.emplace_back(j, arr_constraint_lower[j]);
                                                    }
                                                }
                                                const_self_vec[order - 1].emplace_back(const_tmp);

                                            } else if (is_allzero(arr_constraint_lower, tolerance, loc_nonzero)) {
                                                // If all elements of the "order-1"th order is zero,
                                                // the constraint is intraorder of the "order"th order.
                                                for (j = 0; j < nparams[order]; ++j) {
                                                    if (std::abs(arr_constraint_self[j]) >= tolerance) {
                                                        const_tmp.emplace_back(j, arr_constraint_self[j]);
                                                    }
                                                }
                                                const_self_vec[order].emplace_back(const_tmp);

                                            } else {
                                                // If nonzero elements exist in both of the "order-1" and "order",
                                                // the constraint is intrerorder.

                                                for (j = 0; j < nparam_sub; ++j) {
                                                    if (std::abs(arr_constraint[j]) >= tolerance) {
                                                        const_tmp.emplace_back(j, arr_constraint[j]);
                                                    }
                                                }
                                                const_cross_vec[order].emplace_back(const_tmp);
                                            }
                                        }

                                    } // nu
                                }     // mu

                            } // ixyz

                        } while (g.next());

                    } // direction
                }     // icrd
            }

            // Additional constraint for the last order.
            // All IFCs over maxorder-th order are neglected.

            if (order == maxorder - 1 && !exclude_last_R) {

                auto interaction_list_now(cluster->get_interaction_pair(order, i));
                std::sort(interaction_list_now.begin(), interaction_list_now.end());

                nxyz2 = static_cast<int>(pow(static_cast<double>(3), order + 1));
                allocate(xyzcomponent2, nxyz2, order + 1);
                fcs->get_xyzcomponent(order + 1, xyzcomponent2);

                for (icrd = 0; icrd < 3; ++icrd) {

                    interaction_index[0] = 3 * interaction_atom[0] + icrd;

                    CombinationWithRepetition<int> g_now(interaction_list_now.begin(),
                                                         interaction_list_now.end(), order + 1);
                    do {

                        auto data = g_now.now();

                        for (auto idata = 0; idata < data.size(); ++idata)
                            interaction_atom[idata + 1] = data[idata];

                        for (ixyz = 0; ixyz < nxyz2; ++ixyz) {

                            for (j = 0; j < order + 1; ++j)
                                interaction_index[j + 1] = 3 * interaction_atom[j + 1] + xyzcomponent2[ixyz][j];

                            for (mu = 0; mu < 3; ++mu) {

                                for (nu = 0; nu < 3; ++nu) {

                                    if (!valid_rotation_axis[mu][nu]) continue;

                                    for (j = 0; j < nparams[order]; ++j) arr_constraint_self[j] = 0.0;

                                    for (lambda = 0; lambda < order + 2; ++lambda) {

                                        mu_lambda = interaction_index[lambda] % 3;

                                        for (jcrd = 0; jcrd < 3; ++jcrd) {

                                            for (j = 0; j < order + 2; ++j)
                                                interaction_tmp[j] = interaction_index[j];

                                            interaction_tmp[lambda] = 3 * interaction_atom[lambda] + jcrd;

                                            levi_factor = 0;
                                            for (j = 0; j < 3; ++j) {
                                                levi_factor += levi_civita(j, mu, nu) * levi_civita(j, mu_lambda, jcrd);
                                            }

                                            if (levi_factor == 0) continue;

                                            sort_tail(order + 2, interaction_tmp);

                                            iter_found = list_found.find(FcProperty(order + 2, 1.0,
                                                                                    interaction_tmp, 1));
                                            if (iter_found != list_found.end()) {
                                                arr_constraint_self[(*iter_found).mother]
                                                    += (*iter_found).sign * static_cast<double>(levi_factor);
                                            }
                                        } // jcrd
                                    }     // lambda

                                    if (!is_allzero(arr_constraint_self, tolerance, loc_nonzero)) {
                                        if (arr_constraint_self[loc_nonzero] < 0.0) {
                                            for (j = 0; j < nparams[order]; ++j) arr_constraint_self[j] *= -1.0;
                                        }
                                        const_tmp.clear();
                                        for (j = 0; j < nparams[order]; ++j) {
                                            if (std::abs(arr_constraint_self[j]) >= tolerance) {
                                                const_tmp.emplace_back(j, arr_constraint_self[j]);
                                            }
                                        }
                                        const_self_vec[order].emplace_back(const_tmp);
                                    }

                                } // nu
                            }     // mu

                        } // ixyz


                    } while (g_now.next());

                } // icrd

                deallocate(xyzcomponent2);
            }
        } // iat

        if (verbosity > 0) std::cout << " done." << std::endl;

        if (order > 0) {
            deallocate(xyzcomponent);
        }
        deallocate(interaction_tmp);
        deallocate(interaction_index);
        deallocate(interaction_atom);
    } // order

    int counter;
    std::map<size_t, double> const_copy;
    auto division_factor = 1.0;

    for (order = 0; order < maxorder; ++order) {
        // Sort & unique
        std::sort(const_self_vec[order].begin(), const_self_vec[order].end());
        const_self_vec[order].erase(std::unique(const_self_vec[order].begin(),
                                                const_self_vec[order].end()),
                                    const_self_vec[order].end());
        std::sort(const_cross_vec[order].begin(), const_cross_vec[order].end());
        const_cross_vec[order].erase(std::unique(const_cross_vec[order].begin(),
                                                 const_cross_vec[order].end()),
                                     const_cross_vec[order].end());

        // Copy to the return variable
        for (const auto &it : const_self_vec[order]) {
            const_copy.clear();
            counter = 0;
            for (const auto &it2 : it) {
                if (counter == 0) {
                    division_factor = 1.0 / it2.val;
                }
                const_copy[it2.col] = it2.val * division_factor;
                ++counter;
            }
            const_rotation_self[order].emplace_back(const_copy);
        }
        const_self_vec[order].clear();

        for (const auto &it : const_cross_vec[order]) {
            const_copy.clear();
            counter = 0;
            for (const auto &it2 : it) {
                if (counter == 0) {
                    division_factor = 1.0 / it2.val;
                }
                const_copy[it2.col] = it2.val * division_factor;
                ++counter;
            }
            const_rotation_cross[order].emplace_back(const_copy);
        }
        const_cross_vec[order].clear();

        //  Perform rref
        rref_sparse(nparams[order],
                    const_rotation_self[order],
                    eps6);

        if (order > 0) {
            rref_sparse(nparams[order - 1] + nparams[order],
                        const_rotation_cross[order],
                        eps6);
        }

    }

    if (verbosity > 0) std::cout << "  Finished !" << std::endl << std::endl;

    deallocate(ind);
    deallocate(nparams);
    deallocate(const_self_vec);
    deallocate(const_cross_vec);
}


int Constraint::levi_civita(const int i,
                            const int j,
                            const int k) const
{
    return (j - i) * (k - i) * (k - j) / 2;
}


void Constraint::setup_rotation_axis(bool flag[3][3])
{
    for (auto mu = 0; mu < 3; ++mu) {
        for (auto nu = 0; nu < 3; ++nu) {
            if (mu == nu) {
                flag[mu][nu] = false;
            } else {
                flag[mu][nu] = true;
            }
        }
    }
    std::sort(rotation_axis.begin(), rotation_axis.end());

    if (rotation_axis == "x") {
        flag[0][1] = false;
        flag[1][0] = false;
        flag[0][2] = false;
        flag[2][0] = false;
    } else if (rotation_axis == "y") {
        flag[0][1] = false;
        flag[1][0] = false;
        flag[1][2] = false;
        flag[2][1] = false;
    } else if (rotation_axis == "z") {
        flag[0][2] = false;
        flag[2][0] = false;
        flag[1][2] = false;
        flag[2][1] = false;
    } else if (rotation_axis == "xy") {
        flag[0][1] = false;
        flag[1][0] = false;
    } else if (rotation_axis == "yz") {
        flag[1][2] = false;
        flag[2][1] = false;
    } else if (rotation_axis == "xz") {
        flag[0][2] = false;
        flag[2][0] = false;
    } else if (rotation_axis == "xyz") {
        // do nothing
    } else {
        warn("setup_rotation_axis",
             "Invalid rotation_axis. Default value(xyz) will be used.");
    }
}


void Constraint::fix_forceconstants_to_file(const int order,
                                            const Symmetry *symmetry,
                                            const Fcs *fcs,
                                            const std::string file_to_fix,
                                            std::vector<ConstraintTypeFix> &const_out) const
{
    using namespace boost::property_tree;
    ptree pt;

    try {
        read_xml(file_to_fix, pt);
    }
    catch (std::exception &e) {
        if (order == 0) {
            auto str_error = "Cannot open file FC2XML ( " + file_to_fix + " )";
        } else if (order == 1) {
            auto str_error = "Cannot open file FC3XML ( " + file_to_fix + " )";
        }
        exit("fix_forceconstants_to_file", "Failed to open ", file_to_fix.c_str());
    }

    const auto nat_ref = boost::lexical_cast<size_t>(
        get_value_from_xml(pt, "Data.Structure.NumberOfAtoms"));
    const auto ntran_ref = boost::lexical_cast<size_t>(
        get_value_from_xml(pt, "Data.Symmetry.NumberOfTranslations"));
    const auto natmin_ref = nat_ref / ntran_ref;

    if (natmin_ref != symmetry->get_nat_prim()) {
        exit("fix_forceconstants_to_file",
             "The number of atoms in the primitive cell is not consistent.");
    }

    const auto nfcs = fcs->get_nequiv()[order].size();

    if (order == 0) {
        const auto nfcs_ref = boost::lexical_cast<size_t>(
            get_value_from_xml(pt, "Data.ForceConstants.HarmonicUnique.NFC2"));

        if (nfcs_ref != nfcs) {
            exit("load_reference_system_xml",
                 "The number of harmonic force constants is not consistent.");
        }
    } else if (order == 1) {
        const auto nfcs_ref = boost::lexical_cast<size_t>(
            get_value_from_xml(pt, "Data.ForceConstants.CubicUnique.NFC3"));

        if (nfcs_ref != nfcs) {
            exit("load_reference_system_xml",
                 "The number of cubic force constants is not consistent.");
        }
    }

    int **intpair_ref;
    double *fcs_ref;

    allocate(fcs_ref, nfcs);
    allocate(intpair_ref, nfcs, 3);

    int counter = 0;

    if (order == 0) {
        BOOST_FOREACH(const ptree::value_type & child_, pt.get_child("Data.ForceConstants.HarmonicUnique")) {
            if (child_.first == "FC2") {
                const auto &child = child_.second;
                const auto str_intpair = child.get<std::string>("<xmlattr>.pairs");
                const auto str_multiplicity = child.get<std::string>("<xmlattr>.multiplicity");

                std::istringstream is(str_intpair);
                is >> intpair_ref[counter][0] >> intpair_ref[counter][1];
                fcs_ref[counter] = boost::lexical_cast<double>(child.data());
                ++counter;
            }
        }
    } else if (order == 1) {
        BOOST_FOREACH(const ptree::value_type & child_, pt.get_child("Data.ForceConstants.CubicUnique")) {
            if (child_.first == "FC3") {
                const auto &child = child_.second;
                const auto str_intpair = child.get<std::string>("<xmlattr>.pairs");
                const auto str_multiplicity = child.get<std::string>("<xmlattr>.multiplicity");

                std::istringstream is(str_intpair);
                is >> intpair_ref[counter][0] >> intpair_ref[counter][1] >> intpair_ref[counter][2];
                fcs_ref[counter] = boost::lexical_cast<double>(child.data());
                ++counter;
            }
        }
    }

    const auto nterms = order + 2;

    std::unordered_set<FcProperty> list_found;
    std::unordered_set<FcProperty>::iterator iter_found;

    list_found.clear();

    for (auto &list_tmp : fcs->get_fc_table()[order]) {
        list_found.insert(FcProperty(list_tmp));
    }

    for (auto i = 0; i < nfcs; ++i) {
        iter_found = list_found.find(FcProperty(nterms, 1.0,
                                                intpair_ref[i], 1));
        if (iter_found == list_found.end()) {
            exit("fix_forceconstants_to_file",
                 "Cannot find equivalent force constant, number: ",
                 i + 1);
        }
        const_out.emplace_back(ConstraintTypeFix((*iter_found).mother, fcs_ref[i]));
    }
    deallocate(intpair_ref);
    deallocate(fcs_ref);

    list_found.clear();
}


bool Constraint::is_allzero(const int n,
                            const double *arr,
                            const int nshift) const
{
    for (auto i = nshift; i < n; ++i) {
        if (std::abs(arr[i]) > eps10) {
            return false;
        }
    }
    return true;
}

bool Constraint::is_allzero(const std::vector<int> &vec,
                            int &loc) const
{
    loc = -1;
    for (auto i = 0; i < vec.size(); ++i) {
        if (std::abs(vec[i]) > 0) {
            loc = i;
            return false;
        }
    }
    return true;
}

bool Constraint::is_allzero(const std::vector<double> &vec,
                            const double tol,
                            int &loc,
                            const int nshift) const
{
    loc = -1;
    const auto n = vec.size();
    for (auto i = nshift; i < n; ++i) {
        if (std::abs(vec[i]) > tol) {
            loc = i;
            return false;
        }
    }
    return true;
}


void Constraint::remove_redundant_rows(const size_t n,
                                       std::vector<ConstraintClass> &Constraint_vec,
                                       const double tolerance) const
{
    size_t i, j;

    auto nparam = n;
    const auto nconst = Constraint_vec.size();
    double *arr_tmp;
    double **mat_tmp;

    size_t nrank;

    if (nconst > 0) {

        allocate(mat_tmp, nconst, nparam);

        i = 0;

        for (auto &p : Constraint_vec) {
            for (j = 0; j < nparam; ++j) {
                mat_tmp[i][j] = p.w_const[j];
            }
            ++i;
        }

        rref(nconst, nparam, mat_tmp, nrank, tolerance);

        allocate(arr_tmp, nparam);

        Constraint_vec.clear();

        for (i = 0; i < nrank; ++i) {
            for (j = 0; j < nparam; ++j) arr_tmp[j] = 0.0;
            long iloc = -1;
            for (j = 0; j < nparam; ++j) {
                if (std::abs(mat_tmp[i][j]) < tolerance) {
                    arr_tmp[j] = 0.0;
                } else {
                    arr_tmp[j] = mat_tmp[i][j];
                }
                if (std::abs(arr_tmp[j]) >= tolerance) {
                    iloc = j;
                }
            }

            if (iloc != -1) {
                Constraint_vec.emplace_back(nparam, arr_tmp);
            }
        }

        deallocate(mat_tmp);
        deallocate(arr_tmp);

    }
}

void Constraint::print_constraint(const ConstraintSparseForm &const_in) const
{
    const auto nconst = const_in.size();
    auto counter = 0;
    std::cout << std::endl;
    std::cout << "TOTAL CONST SIZE :" << std::setw(6) << nconst << std::endl;
    for (const auto &it : const_in) {
        std::cout << "CONST : " << std::setw(5) << counter + 1 << std::endl;
        for (const auto &it2 : it) {
            std::cout << std::setw(5) << it2.first + 1;
            std::cout << std::setw(15) << it2.second << std::endl;
        }
        std::cout << std::endl;
        ++counter;
    }
}
