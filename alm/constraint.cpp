/*
 constraint.cpp

 Copyright (c) 2014, 2015, 2016 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#include <iomanip>
#include "constraint.h"
#include "interaction.h"
#include "memory.h"
#include "fcs.h"
#include "symmetry.h"
#include "system.h"
#include "timer.h"
#include "combination.h"
#include "constants.h"
#include "error.h"
#include <boost/bimap.hpp>
#include "mathfunctions.h"

#ifdef _USE_EIGEN
#include <Eigen/Dense>
#endif

using namespace ALM_NS;

Constraint::Constraint(ALM *alm) : Pointers(alm)
{
}

Constraint::~Constraint()
{
    if (exist_constraint && alm->mode == "fitting") {

        memory->deallocate(const_symmetry);

        if (constraint_algebraic) {
            memory->deallocate(const_fix);
            memory->deallocate(const_relate);
            memory->deallocate(index_bimap);
        } else {
            memory->deallocate(const_mat);
            memory->deallocate(const_rhs);
        }
    }
}

void Constraint::setup()
{
    std::cout << " CONSTRAINT" << std::endl;
    std::cout << " ==========" << std::endl << std::endl;

    constraint_algebraic = constraint_mode / 10;
    constraint_mode = constraint_mode % 10;

    switch (constraint_mode) {
    case 0: // do nothing
        impose_inv_T = false;
        impose_inv_R = false;
        std::cout << "  ICONST = 0: Constraint for translational/rotational invariance" << std::endl;
        std::cout << "              will NOT be considered." << std::endl;
        break;
    case 1:
        impose_inv_T = true;
        impose_inv_R = false;
        std::cout << "  ICONST = 1: Constraints for translational invariance" << std::endl;
        std::cout << "              will be considered." << std::endl;
        break;
    case 2:
        impose_inv_T = true;
        impose_inv_R = true;
        exclude_last_R = true;
        std::cout << "  ICONST = 2: Constraints for translational and rotational invariance" << std::endl;
        std::cout << "              will be considered. Axis of rotation is " << rotation_axis << std::endl;
        std::cout << "              Rotational invariance of the maximum order will be neglected" << std::endl;
        break;
    case 3:
        impose_inv_T = true;
        impose_inv_R = true;
        exclude_last_R = false;
        std::cout << "  ICONST = 3: Constraints for translational and rotational invariance" << std::endl;
        std::cout << "              will be considered. Axis of rotation is " << rotation_axis << std::endl;
        break;
    default:
        error->exit("Constraint::setup", "invalid constraint_mode", constraint_mode);
        break;
    }

    std::cout << std::endl;

    if (fix_harmonic) {
        std::cout << "  FC2XML is given : Harmonic force constants will be " << std::endl;
        std::cout << "                    fixed to the values given in " << fc2_file << std::endl;
        std::cout << std::endl;
    }

    fix_cubic = fix_cubic & (interaction->maxorder > 1);
    if (fix_cubic) {
        std::cout << "  FC3XML is given : Cubic force constants will be " << std::endl;
        std::cout << "                    fixed to the values given in " << fc3_file << std::endl;
        std::cout << std::endl;
    }

    extra_constraint_from_symmetry = false;
    memory->allocate(const_symmetry, interaction->maxorder);
    constraint_from_symmetry(const_symmetry);
    for (int order = 0; order < interaction->maxorder; ++order) {
        if (const_symmetry[order].size() > 0) extra_constraint_from_symmetry = true;
    }

    if (impose_inv_T || fix_harmonic || fix_cubic || extra_constraint_from_symmetry) {
        exist_constraint = true;
    } else {
        exist_constraint = false;
    }

    if (exist_constraint) {

        int i;
        int maxorder = interaction->maxorder;
        int Pmax, order;
        int N = 0;

        for (i = 0; i < maxorder; ++i) {
            N += fcs->ndup[i].size();
        }

        memory->allocate(const_translation, maxorder);
        memory->allocate(const_rotation_self, maxorder);
        memory->allocate(const_rotation_cross, maxorder);

        for (order = 0; order < maxorder; ++order) {
            const_translation[order].clear();
            const_rotation_self[order].clear();
            const_rotation_cross[order].clear();
        }

        if (impose_inv_T) {
            translational_invariance();
        }
        if (impose_inv_R) {
            rotational_invariance();
        }

        if (impose_inv_T || impose_inv_R) {
            std::cout << "  Number of constraints [T-inv, R-inv (self), R-inv (cross)]:" << std::endl;
            for (order = 0; order < maxorder; ++order) {
                std::cout << "   " << std::setw(8) << interaction->str_order[order];
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
                std::cout << "   " << std::setw(8) << interaction->str_order[order];
                std::cout << std::setw(5) << const_symmetry[order].size();
                std::cout << std::endl;
            }
            std::cout << std::endl;
        }

        memory->allocate(const_self, maxorder);
        for (order = 0; order < maxorder; ++order) const_self[order].clear();

        int nparam;
        double *arr_tmp;

        std::vector<ConstraintClass>::const_iterator it_const;

        // Merge translational and rotational invariance excluding order-crossing constraints

        for (order = 0; order < maxorder; ++order) {

            nparam = fcs->ndup[order].size();
            memory->allocate(arr_tmp, nparam);

            for (it_const = const_translation[order].begin();
                 it_const != const_translation[order].end(); ++it_const) {
                for (i = 0; i < nparam; ++i) arr_tmp[i] = (*it_const).w_const[i];
                const_self[order].push_back(ConstraintClass(nparam, arr_tmp));
            }

            if (const_rotation_self[order].size() > 0) {
                for (it_const = const_rotation_self[order].begin();
                     it_const != const_rotation_self[order].end(); ++it_const) {
                    for (i = 0; i < nparam; ++i) arr_tmp[i] = (*it_const).w_const[i];
                    const_self[order].push_back(ConstraintClass(nparam, arr_tmp));
                }
                remove_redundant_rows(nparam, const_self[order], eps8);
            }

            if (const_symmetry[order].size() > 0) {
                for (it_const = const_symmetry[order].begin();
                     it_const != const_symmetry[order].end(); ++it_const) {
                    for (i = 0; i < nparam; ++i) arr_tmp[i] = (*it_const).w_const[i];
                    const_self[order].push_back(ConstraintClass(nparam, arr_tmp));
                }
                remove_redundant_rows(nparam, const_self[order], eps8);
            }

            memory->deallocate(arr_tmp);
            const_translation[order].clear();
            const_rotation_self[order].clear();
        }

        if (extra_constraint_from_symmetry) {
            std::cout << "  Constraints of T-inv, R-inv (self), and those from crystal symmetry are merged." << std::endl;
        } else {
            std::cout << "  Constraints of T-inv and R-inv (self) are merged." << std::endl;
        }
        std::cout << "  If there are redundant constraints, they are removed in this process." << std::endl;
        std::cout << std::endl;
        std::cout << "  Number of inequivalent constraints (self, cross) : " << std::endl;

        for (order = 0; order < maxorder; ++order) {
            std::cout << "   " << std::setw(8) << interaction->str_order[order];
            std::cout << std::setw(5) << const_self[order].size();
            std::cout << std::setw(5) << const_rotation_cross[order].size();
            std::cout << std::endl;
        }
        std::cout << std::endl;

        if (constraint_algebraic) {

            std::cout << "  ICONST >= 10 : Constraints will be considered algebraically."
                << std::endl << std::endl;

            if (impose_inv_R) {
                std::cout << "  WARNING : Inter-order constraints for rotational invariance will be neglected." << std::endl;
            }

            memory->allocate(const_fix, maxorder);
            memory->allocate(const_relate, maxorder);
            memory->allocate(index_bimap, maxorder);

            get_mapping_constraint(maxorder, const_self, const_fix,
                                   const_relate, index_bimap, false);

            for (order = 0; order < maxorder; ++order) {
                std::cout << "  Number of free" << std::setw(9) << interaction->str_order[order]
                    << " FCs : " << index_bimap[order].size() << std::endl;
            }
            std::cout << std::endl;

        } else {

            Pmax = 0;
            for (order = 0; order < maxorder; ++order) {
                Pmax += const_self[order].size() + const_rotation_cross[order].size();
            }
            if (fix_harmonic) {
                Pmax -= const_self[0].size();
                Pmax += fcs->ndup[0].size();
            }
            if (fix_cubic) {
                Pmax -= const_self[1].size();
                Pmax += fcs->ndup[1].size();
            }
            memory->allocate(const_mat, Pmax, N);
            memory->allocate(const_rhs, Pmax);

            calc_constraint_matrix(N, P);
            std::cout << "  Total number of constraints = " << P << std::endl << std::endl;

        }

        memory->deallocate(const_translation);
        memory->deallocate(const_rotation_self);
        memory->deallocate(const_rotation_cross);
        memory->deallocate(const_self);

        timer->print_elapsed();
        std::cout << " -------------------------------------------------------------------" << std::endl;
        std::cout << std::endl;
    }
}

void Constraint::calc_constraint_matrix(const int N, int &P)
{
    int i, j;
    int maxorder = interaction->maxorder;
    int order;
    int nconst1;
    int icol, irow;
    double *arr_tmp;
    std::vector<ConstraintClass> const_total;

    const_total.clear();
    memory->allocate(arr_tmp, N);

    int nshift = 0;

    // Intra-order constraints
    for (order = 0; order < maxorder; ++order) {
        int nparam = fcs->ndup[order].size();

        if ((order == 0 && !fix_harmonic) || (order == 1 && !fix_cubic) || order > 1) {
            for (i = 0; i < N; ++i) arr_tmp[i] = 0.0;

            for (std::vector<ConstraintClass>::iterator p = const_self[order].begin();
                 p != const_self[order].end(); ++p) {
                ConstraintClass const_now = *p;
                for (i = 0; i < nparam; ++i) {
                    arr_tmp[nshift + i] = const_now.w_const[i];
                }
                const_total.push_back(ConstraintClass(N, arr_tmp));
            }
        }
        nshift += nparam;
    }
    nconst1 = const_total.size();

    // order-crossing constraints
    int nshift2 = 0;
    for (order = 0; order < maxorder; ++order) {
        if (order > 0) {
            int nparam2 = fcs->ndup[order - 1].size() + fcs->ndup[order].size();
            for (i = 0; i < N; ++i) arr_tmp[i] = 0.0;
            for (std::vector<ConstraintClass>::iterator p = const_rotation_cross[order].begin();
                 p != const_rotation_cross[order].end(); ++p) {
                ConstraintClass const_now = *p;
                for (i = 0; i < nparam2; ++i) {
                    arr_tmp[nshift2 + i] = const_now.w_const[i];
                }
                const_total.push_back(ConstraintClass(N, arr_tmp));
            }
            nshift2 += fcs->ndup[order - 1].size();
        }
    }
    memory->deallocate(arr_tmp);

    if (nconst1 != const_total.size())
        remove_redundant_rows(N, const_total, eps8);

    P = const_total.size();

    if (fix_harmonic) {
        std::cout << "  Harmonic force constants will be fixed to the values " << std::endl;
        std::cout << "  of the reference " << fc2_file << std::endl;
        std::cout << "  Constraint for HARMONIC IFCs will be updated accordingly."
            << std::endl << std::endl;
        P += fcs->ndup[0].size();
    }

    if (fix_cubic) {
        std::cout << "  Cubic force constants will be fixed to the values " << std::endl;
        std::cout << "  of the reference " << fc3_file << std::endl;
        std::cout << "  Constraint for ANHARM3 IFCs will be updated accordingly."
            << std::endl << std::endl;
        P += fcs->ndup[1].size();
    }

    for (i = 0; i < P; ++i) {
        for (j = 0; j < N; ++j) {
            const_mat[i][j] = 0.0;
        }
        const_rhs[i] = 0.0;
    }

    irow = 0;
    icol = 0;

    if (fix_harmonic) {
        double *const_rhs_tmp;
        int nfcs_tmp = fcs->ndup[0].size();
        memory->allocate(const_rhs_tmp, nfcs_tmp);
        system->load_reference_system_xml(fc2_file, 0, const_rhs_tmp);

        for (i = 0; i < nfcs_tmp; ++i) {
            const_mat[i][i] = 1.0;
            const_rhs[i] = const_rhs_tmp[i];
        }

        irow += nfcs_tmp;
        icol += nfcs_tmp;
        memory->deallocate(const_rhs_tmp);
    }

    if (fix_cubic) {
        int ishift = fcs->ndup[0].size();
        double *const_rhs_tmp;
        int nfcs_tmp = fcs->ndup[1].size();
        memory->allocate(const_rhs_tmp, nfcs_tmp);
        system->load_reference_system_xml(fc3_file, 1, const_rhs_tmp);

        for (i = 0; i < nfcs_tmp; ++i) {
            const_mat[i + ishift][i + ishift] = 1.0;
            const_rhs[i + ishift] = const_rhs_tmp[i];
        }

        irow += nfcs_tmp;
        icol += nfcs_tmp;
        memory->deallocate(const_rhs_tmp);
    }

    for (std::vector<ConstraintClass>::iterator p = const_total.begin();
         p != const_total.end(); ++p) {
        for (i = 0; i < N; ++i) {
            const_mat[irow][i] = (*p).w_const[i];
        }
        ++irow;
    }
    const_total.clear();
}


void Constraint::get_mapping_constraint(const int nmax,
                                        std::vector<ConstraintClass> *const_in,
                                        std::vector<ConstraintTypeFix> *const_fix_out,
                                        std::vector<ConstraintTypeRelate> *const_relate_out,
                                        boost::bimap<int, int> *index_bimap_out,
                                        const bool is_suggest_mode)
{
    int order;
    unsigned int i;
    double const_tol = eps8;

    bool *fix_forceconstant;
    std::string *file_forceconstant;

    memory->allocate(fix_forceconstant, nmax);
    memory->allocate(file_forceconstant, nmax);

    for (i = 0; i < nmax; ++i) {
        fix_forceconstant[i] = false;
        file_forceconstant[i] = "";
    }

    if (!is_suggest_mode) {
        fix_forceconstant[0] = fix_harmonic;
        fix_forceconstant[1] = fix_cubic;
        if (fix_forceconstant[0]) file_forceconstant[0] = fc2_file;
        if (fix_forceconstant[1]) file_forceconstant[1] = fc3_file;
    }

    int nparam;
    for (order = 0; order < nmax; ++order) {

        nparam = fcs->ndup[order].size();

        if (fix_forceconstant[order]) {

            double *const_rhs_tmp;
            memory->allocate(const_rhs_tmp, nparam);
            system->load_reference_system_xml(file_forceconstant[order],
                                              order, const_rhs_tmp);

            for (i = 0; i < nparam; ++i) {
                const_fix_out[order].push_back(ConstraintTypeFix(i, const_rhs_tmp[i]));
            }
            memory->deallocate(const_rhs_tmp);

        } else {

            int p_index_target;
            std::vector<double> alpha_tmp;
            std::vector<unsigned int> p_index_tmp;

            for (std::vector<ConstraintClass>::reverse_iterator p = const_in[order].rbegin();
                 p != const_in[order].rend(); ++p) {
                p_index_target = -1;
                for (i = 0; i < nparam; ++i) {
                    if (std::abs((*p).w_const[i]) > const_tol) {
                        p_index_target = i;
                        break;
                    }
                }

                if (p_index_target == -1) {
                    error->exit("get_mapping_constraint",
                                "No finite entry found in the constraint.");
                }

                alpha_tmp.clear();
                p_index_tmp.clear();

                for (i = p_index_target + 1; i < nparam; ++i) {
                    if (std::abs((*p).w_const[i]) > const_tol) {
                        alpha_tmp.push_back((*p).w_const[i]);
                        p_index_tmp.push_back(i);
                    }
                }

                if (alpha_tmp.size() > 0) {
                    const_relate_out[order].push_back(
                        ConstraintTypeRelate(p_index_target,
                                             alpha_tmp, p_index_tmp));
                } else {
                    const_fix_out[order].push_back(
                        ConstraintTypeFix(p_index_target, 0.0));
                }
            }
        }
    }


    std::vector<int> *has_constraint;
    memory->allocate(has_constraint, nmax);

    for (order = 0; order < nmax; ++order) {

        nparam = fcs->ndup[order].size();

        for (i = 0; i < nparam; ++i) {
            has_constraint[order].push_back(0);
        }

        for (i = 0; i < const_fix_out[order].size(); ++i) {
            has_constraint[order][const_fix_out[order][i].p_index_target] = 1;
        }

        for (i = 0; i < const_relate_out[order].size(); ++i) {
            has_constraint[order][const_relate_out[order][i].p_index_target] = 2;
        }
    }

    int icount;

    for (order = 0; order < nmax; ++order) {
        nparam = fcs->ndup[order].size();

        icount = 0;
        for (i = 0; i < nparam; ++i) {

            if (!has_constraint[order][i]) {
                index_bimap_out[order].insert(
                    boost::bimap<int, int>::value_type(icount, i));
                ++icount;
            }
        }
    }

    memory->deallocate(has_constraint);
}


void Constraint::constraint_from_symmetry(std::vector<ConstraintClass> *const_out)
{
    // Create constraint matrices arising from the crystal symmetry.

    int i;
    unsigned int isym;
    int ixyz, nxyz;
    int order;
    int maxorder = interaction->maxorder;

    int *index_tmp;
    int **xyzcomponent;
    int nparams;

    double *arr_constraint;
    bool has_constraint_from_symm = false;
    std::set<FcProperty> list_found;
    std::vector<std::vector<double>> const_mat;

    for (isym = 0; isym < symmetry->nsym; ++isym) {
        if (symmetry->sym_available[isym]) continue;
        has_constraint_from_symm = true;
    }

    for (order = 0; order < maxorder; ++order) const_out[order].clear();

    if (has_constraint_from_symm) {
        std::cout << "  Generating constraints from crystal symmetry ..." << std::endl;
    }

    memory->allocate(index_tmp, maxorder + 1);

    const_mat.clear();

    for (order = 0; order < maxorder; ++order) {

        nparams = fcs->ndup[order].size();

        if (has_constraint_from_symm) {
            std::cout << "   " << std::setw(8) << interaction->str_order[order] << " ...";
            if (nparams == 0) {
                std::cout << "  No parameters! Skipped." << std::endl;
                continue;
            }
        }

        // Generate temporary list of parameters
        list_found.clear();
        for (auto p = fcs->fc_set[order].begin(); p != fcs->fc_set[order].end(); ++p) {
            for (i = 0; i < order + 2; ++i) index_tmp[i] = (*p).elems[i];
            list_found.insert(FcProperty(order + 2, (*p).coef,
                                         index_tmp, (*p).mother));
        }

        nxyz = static_cast<int>(std::pow(static_cast<double>(3), order + 2));
        memory->allocate(xyzcomponent, nxyz, order + 2);
        fcs->get_xyzcomponent(order + 2, xyzcomponent);

        int nfcs = fcs->fc_set[order].size();

#ifdef _OPENMP
#pragma omp parallel
#endif
        {
            int j;
            int i_prim;
            int loc_nonzero;
            int *ind;
            int *atm_index, *atm_index_symm;
            int *xyz_index;
            double c_tmp;

            std::set<FcProperty>::iterator iter_found;
            std::vector<double> const_now_omp;
            std::vector<std::vector<double>> const_omp;

            memory->allocate(ind, order + 2);
            memory->allocate(atm_index, order + 2);
            memory->allocate(atm_index_symm, order + 2);
            memory->allocate(xyz_index, order + 2);

            const_omp.clear();
            const_now_omp.resize(nparams);

#ifdef _OPENMP
#pragma omp for private(i, isym, ixyz)
#endif
            for (int ii = 0; ii < nfcs; ++ii) {
                FcProperty list_tmp = fcs->fc_set[order][ii];

                for (i = 0; i < order + 2; ++i) {
                    atm_index[i] = list_tmp.elems[i] / 3;
                    xyz_index[i] = list_tmp.elems[i] % 3;
                }

                for (isym = 0; isym < symmetry->nsym; ++isym) {

                    if (symmetry->sym_available[isym]) continue;

                    for (i = 0; i < order + 2; ++i)
                        atm_index_symm[i] = symmetry->map_sym[atm_index[i]][isym];
                    if (!fcs->is_inprim(order + 2, atm_index_symm)) continue;

                    for (i = 0; i < nparams; ++i) const_now_omp[i] = 0.0;

                    const_now_omp[list_tmp.mother] = -list_tmp.coef;

                    for (ixyz = 0; ixyz < nxyz; ++ixyz) {
                        for (i = 0; i < order + 2; ++i)
                            ind[i] = 3 * atm_index_symm[i] + xyzcomponent[ixyz][i];

                        i_prim = fcs->min_inprim(order + 2, ind);
                        std::swap(ind[0], ind[i_prim]);
                        fcs->sort_tail(order + 2, ind);

                        iter_found = list_found.find(FcProperty(order + 2, 1.0, ind, 1));
                        if (iter_found != list_found.end()) {
                            c_tmp = fcs->coef_sym(order + 2, isym, xyz_index, xyzcomponent[ixyz]);
                            const_now_omp[(*iter_found).mother] += (*iter_found).coef * c_tmp;
                        }
                    }

                    if (!is_allzero(const_now_omp, loc_nonzero)) {
                        if (const_now_omp[loc_nonzero] < 0.0) {
                            for (j = 0; j < nparams; ++j) const_now_omp[j] *= -1.0;
                        }
                        const_omp.push_back(const_now_omp);
                    }

                } // close isym loop


                // sort-->uniq the array
                std::sort(const_omp.begin(), const_omp.end());
                const_omp.erase(std::unique(const_omp.begin(), const_omp.end()),
                                const_omp.end());

                // Merge vectors
#pragma omp critical
                {
                    for (std::vector<std::vector<double>>::iterator it = const_omp.begin();
                         it != const_omp.end(); ++it) {
                        const_mat.push_back(*it);
                    }
                }
                const_omp.clear();

            } // close ii loop

            memory->deallocate(ind);
            memory->deallocate(atm_index);
            memory->deallocate(atm_index_symm);
            memory->deallocate(xyz_index);

        } // close openmp region

        memory->allocate(arr_constraint, nparams);
        for (std::vector<std::vector<double>>::reverse_iterator it = const_mat.rbegin();
             it != const_mat.rend(); ++it) {
            for (i = 0; i < (*it).size(); ++i) {
                arr_constraint[i] = (*it)[i];
            }
            const_out[order].push_back(ConstraintClass(nparams,
                                                       arr_constraint));
        }
        const_mat.clear();

        memory->deallocate(xyzcomponent);
        memory->deallocate(arr_constraint);

        remove_redundant_rows(nparams, const_out[order], eps8);

        if (has_constraint_from_symm) {
            std::cout << " done." << std::endl;
        }
    } // close loop order

    memory->deallocate(index_tmp);

    if (has_constraint_from_symm) {
        std::cout << "  Finished !" << std::endl << std::endl;
    }
}

void Constraint::translational_invariance()
{
    // Create constraint matrix for the translational invariance (aka acoustic sum rule).

    int i, j;
    int iat, jat, icrd, jcrd;
    int idata;
    int order;
    int loc_nonzero;
    int maxorder = interaction->maxorder;

    int *ind;
    int *intarr, *intarr_copy;
    int **xyzcomponent;

    int ixyz, nxyz;
    int natmin = symmetry->natmin;
    int nat = system->nat;
    int nparams;

    unsigned int isize;
    double *arr_constraint;

    std::vector<int> intlist, data;
    std::set<FcProperty> list_found;
    std::set<FcProperty>::iterator iter_found;
    std::vector<std::vector<int>> data_vec;
    std::vector<FcProperty> list_vec;
    std::vector<FcProperty>::iterator iter_vec;
    std::vector<int> const_now;
    std::vector<std::vector<int>> const_mat;

    std::cout << "  Generating constraints for translational invariance ..." << std::endl;

    memory->allocate(ind, maxorder + 1);

    for (order = 0; order < maxorder; ++order) {

        std::cout << "   " << std::setw(8) << interaction->str_order[order] << " ...";

        const_mat.clear();

        nparams = fcs->ndup[order].size();

        if (nparams == 0) {
            std::cout << "  No parameters! Skipped." << std::endl;
            continue;
        }

        // Make interaction list

        list_found.clear();

        for (std::vector<FcProperty>::iterator p = fcs->fc_set[order].begin();
             p != fcs->fc_set[order].end(); ++p) {
            for (i = 0; i < order + 2; ++i) {
                ind[i] = (*p).elems[i];
            }
            if (list_found.find(FcProperty(order + 2, (*p).coef,
                                           ind, (*p).mother)) != list_found.end()) {
                error->exit("translational invariance", "Duplicate interaction list found");
            }
            list_found.insert(FcProperty(order + 2, (*p).coef,
                                         ind, (*p).mother));
        }

        // Generate xyz component for each order

        nxyz = static_cast<int>(std::pow(static_cast<double>(3), order + 1));
        memory->allocate(xyzcomponent, nxyz, order + 1);
        fcs->get_xyzcomponent(order + 1, xyzcomponent);

        memory->allocate(arr_constraint, nparams);
        memory->allocate(intarr, order + 2);
        memory->allocate(intarr_copy, order + 2);

        const_now.resize(nparams);

        for (i = 0; i < natmin; ++i) {

            iat = symmetry->map_p2s[i][0];

            // Generate atom pairs for each order

            if (order == 0) {
                for (icrd = 0; icrd < 3; ++icrd) {

                    intarr[0] = 3 * iat + icrd;

                    for (jcrd = 0; jcrd < 3; ++jcrd) {

                        // Reset the temporary array for another constraint
                        for (j = 0; j < nparams; ++j) const_now[j] = 0;

                        for (jat = 0; jat < 3 * nat; jat += 3) {
                            intarr[1] = jat + jcrd;

                            iter_found = list_found.find(
                                FcProperty(order + 2, 1.0, intarr, 1));

                            //  If found a IFC
                            if (iter_found != list_found.end()) {
                                // Round the coefficient to integer
                                const_now[(*iter_found).mother] += nint((*iter_found).coef);
                            }

                        }
                        // Add to the constraint list
                        if (!is_allzero(const_now, loc_nonzero)) {
                            if (const_now[loc_nonzero] < 0) {
                                for (j = 0; j < nparams; ++j) const_now[j] *= -1;
                            }
                            const_mat.push_back(const_now);
                        }
                    }
                }

            } else {

                // Anharmonic cases

                for (j = 0; j < interaction->interaction_pair[order][i].size(); ++j) {
                    intlist.push_back(interaction->interaction_pair[order][i][j]);
                }
                std::sort(intlist.begin(), intlist.end());

                data_vec.clear();
                // Generate data_vec that contains possible interaction clusters.
                // Each cluster contains (order + 1) atoms, and the last atom index
                // will be treated seperately below.
                CombinationWithRepetition<int> g2(intlist.begin(), intlist.end(), order);
                do {
                    data = g2.now();

                    intarr[0] = iat;

                    for (isize = 0; isize < data.size(); ++isize) {
                        intarr[isize + 1] = data[isize];
                    }

                    if (interaction->nbody(order + 1, intarr)
                        <= interaction->nbody_include[order]) {
                        if (interaction->is_incutoff(order + 1, intarr, order)) {
                            // Add to list if the atoms interact with each other.
                            data_vec.push_back(data);
                        }
                    }

                } while (g2.next());

                int ndata = data_vec.size();

                // Use openmp for acceleration if possible
#ifdef _OPENMP
#pragma omp parallel
#endif
                {
                    int *intarr_omp, *intarr_copy_omp;

                    memory->allocate(intarr_omp, order + 2);
                    memory->allocate(intarr_copy_omp, order + 2);

                    std::vector<std::vector<int>> const_omp;
                    std::vector<int> data_omp;
                    std::vector<int> const_now_omp;

                    const_omp.clear();
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

                                    if (interaction->nbody(order + 2, intarr_omp)
                                        <= interaction->nbody_include[order]) {
                                        for (j = 0; j < order + 1; ++j) {
                                            intarr_copy_omp[j] = 3 * intarr_omp[j] + xyzcomponent[ixyz][j];
                                        }
                                        intarr_copy_omp[order + 1] = jat + jcrd;

                                        fcs->sort_tail(order + 2, intarr_copy_omp);

                                        iter_found = list_found.find(FcProperty(order + 2, 1.0,
                                                                                intarr_copy_omp, 1));
                                        if (iter_found != list_found.end()) {
                                            const_now_omp[(*iter_found).mother] += nint((*iter_found).coef);
                                        }

                                    }
                                } // close loop jat

                                // Add the constraint to the private array
                                if (!is_allzero(const_now_omp, loc_nonzero)) {
                                    if (const_now_omp[loc_nonzero] < 0) {
                                        for (j = 0; j < nparams; ++j) const_now_omp[j] *= -1;
                                    }
                                    const_omp.push_back(const_now_omp);
                                }
                            }
                        }
                        // sort-->uniq the array
                        std::sort(const_omp.begin(), const_omp.end());
                        const_omp.erase(std::unique(const_omp.begin(), const_omp.end()),
                                        const_omp.end());

                        // Merge vectors
#pragma omp critical
                        {
                            for (std::vector<std::vector<int>>::iterator it = const_omp.begin();
                                 it != const_omp.end(); ++it) {
                                const_mat.push_back(*it);
                            }
                        }
                        const_omp.clear();

                    }// close idata (openmp main loop)

                    memory->deallocate(intarr_omp);
                    memory->deallocate(intarr_copy_omp);

                } // close openmp

                intlist.clear();
            } // close if

            // sort--> uniq the array (to save memory consumption)
            std::sort(const_mat.begin(), const_mat.end());
            const_mat.erase(std::unique(const_mat.begin(), const_mat.end()),
                            const_mat.end());
            //            timer->print_elapsed();
        } // close loop i

        memory->deallocate(xyzcomponent);
        memory->deallocate(intarr);
        memory->deallocate(intarr_copy);
        // Copy to constraint class 

        const_translation[order].clear();
        for (std::vector<std::vector<int>>::reverse_iterator it = const_mat.rbegin();
             it != const_mat.rend(); ++it) {
            for (i = 0; i < (*it).size(); ++i) {
                arr_constraint[i] = static_cast<double>((*it)[i]);
            }
            const_translation[order].push_back(ConstraintClass(nparams,
                                                               arr_constraint));
        }
        const_mat.clear();
        memory->deallocate(arr_constraint);

        remove_redundant_rows(nparams, const_translation[order], eps8);

        std::cout << " done." << std::endl;
    } // close loop order
    memory->deallocate(ind);

    std::cout << "  Finished !" << std::endl << std::endl;
}

void Constraint::rotational_invariance()
{
    // Create constraints for the rotational invariance

    std::cout << "  Generating constraints for rotational invariance ..." << std::endl;

#ifdef _DEBUG
    std::ofstream ofs_constraint;
    ofs_constraint.open("CONSTRAINT", std::ios::out);
#endif

    int i, j;
    int iat, jat;
    int icrd, jcrd;
    int order;
    int maxorder = interaction->maxorder;
    int natmin = symmetry->natmin;
    int mu, nu;
    int ixyz, nxyz, nxyz2;
    int mu_lambda, lambda;
    int levi_factor;

    int *ind;
    int **xyzcomponent, **xyzcomponent2;
    int *nparams, nparam_sub;
    int *interaction_index, *interaction_atom;
    int *interaction_tmp;

    double *arr_constraint;
    double *arr_constraint_self;

    bool valid_rotation_axis[3][3];

    double vec_for_rot[3];

    std::vector<int> interaction_list, interaction_list_old, interaction_list_now;

    std::set<FcProperty> list_found;
    std::set<FcProperty> list_found_last;
    std::set<FcProperty>::iterator iter_found;

    CombinationWithRepetition<int> g;

    std::vector<int> atom_tmp;
    std::vector<std::vector<int>> cell_dummy;
    std::set<MinimumDistanceCluster>::iterator iter_cluster;

    setup_rotation_axis(valid_rotation_axis);

    memory->allocate(ind, maxorder + 1);
    memory->allocate(nparams, maxorder);

    for (order = 0; order < maxorder; ++order) {

        nparams[order] = fcs->ndup[order].size();

        if (order == 0) {
            std::cout << "   Constraints between " << std::setw(8)
                << "1st-order IFCs (which are zero) and "
                << std::setw(8) << interaction->str_order[order] << " ...";
            nparam_sub = nparams[order];
        } else {
            std::cout << "   Constraints between " << std::setw(8)
                << interaction->str_order[order - 1] << " and "
                << std::setw(8) << interaction->str_order[order] << " ...";
            nparam_sub = nparams[order] + nparams[order - 1];
        }

        memory->allocate(arr_constraint, nparam_sub);
        memory->allocate(arr_constraint_self, nparams[order]);
        memory->allocate(interaction_atom, order + 2);
        memory->allocate(interaction_index, order + 2);
        memory->allocate(interaction_tmp, order + 2);

        if (order > 0) {
            list_found_last = list_found;
            nxyz = static_cast<int>(pow(static_cast<double>(3), order));
            memory->allocate(xyzcomponent, nxyz, order);
            fcs->get_xyzcomponent(order, xyzcomponent);
        }

        list_found.clear();

        for (std::vector<FcProperty>::iterator p = fcs->fc_set[order].begin();
             p != fcs->fc_set[order].end(); ++p) {
            for (i = 0; i < order + 2; ++i) {
                ind[i] = (*p).elems[i];
            }
            list_found.insert(FcProperty(order + 2, (*p).coef,
                                         ind, (*p).mother));
        }

        for (i = 0; i < natmin; ++i) {

            iat = symmetry->map_p2s[i][0];

            interaction_atom[0] = iat;

            if (order == 0) {

                interaction_list_now.clear();
                for (j = 0; j < interaction->interaction_pair[order][i].size(); ++j) {
                    interaction_list_now.push_back(interaction->interaction_pair[order][i][j]);
                }
                std::sort(interaction_list_now.begin(), interaction_list_now.end());

                // Special treatment for harmonic force constants

                for (icrd = 0; icrd < 3; ++icrd) {

                    interaction_index[0] = 3 * iat + icrd;

                    for (mu = 0; mu < 3; ++mu) {

                        for (nu = 0; nu < 3; ++nu) {

                            if (!valid_rotation_axis[mu][nu]) continue;

                            // Clear history

                            for (j = 0; j < nparam_sub; ++j) arr_constraint[j] = 0.0;

                            for (std::vector<int>::iterator iter_list = interaction_list_now.begin();
                                 iter_list != interaction_list_now.end(); ++iter_list) {

                                jat = *iter_list;
                                interaction_index[1] = 3 * jat + mu;
                                iter_found = list_found.find(FcProperty(order + 2, 1.0,
                                                                        interaction_index, 1));

                                atom_tmp.clear();
                                atom_tmp.push_back(jat);
                                cell_dummy.clear();
                                iter_cluster = interaction->mindist_cluster[order][i].find(
                                    MinimumDistanceCluster(atom_tmp, cell_dummy));

                                if (iter_cluster == interaction->mindist_cluster[order][i].end()) {
                                    error->exit("rotational_invariance",
                                                "interaction not found ...");
                                } else {
                                    for (j = 0; j < 3; ++j) vec_for_rot[j] = 0.0;

                                    int nsize_equiv = (*iter_cluster).cell.size();

                                    for (j = 0; j < nsize_equiv; ++j) {
                                        for (int k = 0; k < 3; ++k) {
                                            vec_for_rot[k]
                                                += interaction->x_image[(*iter_cluster).cell[j][0]][jat][k];
                                        }
                                    }

                                    for (j = 0; j < 3; ++j) {
                                        vec_for_rot[j] /= static_cast<double>(nsize_equiv);
                                    }
                                }


                                if (iter_found != list_found.end()) {
                                    arr_constraint[(*iter_found).mother]
                                        += (*iter_found).coef * vec_for_rot[nu];
                                }

                                // Exchange mu <--> nu and repeat again. 
                                // Note that the sign is inverted (+ --> -) in the summation

                                interaction_index[1] = 3 * jat + nu;
                                iter_found = list_found.find(FcProperty(order + 2, 1.0,
                                                                        interaction_index, 1));
                                if (iter_found != list_found.end()) {
                                    arr_constraint[(*iter_found).mother]
                                        -= (*iter_found).coef * vec_for_rot[mu];
                                }
                            }

                            if (!is_allzero(nparam_sub, arr_constraint)) {
                                // Add to constraint list
                                const_rotation_self[order].push_back(
                                    ConstraintClass(nparam_sub, arr_constraint));
                            }

                        } // nu
                    } // mu
                }
            } else {

                // Constraint between different orders

                interaction_list_old.clear();
                interaction_list_now.clear();

                for (j = 0; j < interaction->interaction_pair[order][i].size(); ++j) {
                    interaction_list_now.push_back(interaction->interaction_pair[order][i][j]);
                }
                for (j = 0; j < interaction->interaction_pair[order - 1][i].size(); ++j) {
                    interaction_list_old.push_back(interaction->interaction_pair[order - 1][i][j]);
                }

                std::sort(interaction_list_now.begin(), interaction_list_now.end());
                std::sort(interaction_list_old.begin(), interaction_list_old.end());

                for (icrd = 0; icrd < 3; ++icrd) {

                    interaction_index[0] = 3 * iat + icrd;

                    CombinationWithRepetition<int> g_now(interaction_list_now.begin(),
                                                         interaction_list_now.end(), order);
                    CombinationWithRepetition<int> g_old(interaction_list_old.begin(),
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
                            std::vector<int> data = g.now();

                            for (int idata = 0; idata < data.size(); ++idata)
                                interaction_atom[idata + 1] = data[idata];

                            for (ixyz = 0; ixyz < nxyz; ++ixyz) {

                                for (j = 0; j < order; ++j)
                                    interaction_index[j + 1] = 3 * interaction_atom[j + 1] + xyzcomponent[ixyz][j];

                                for (mu = 0; mu < 3; ++mu) {

                                    for (nu = 0; nu < 3; ++nu) {

                                        if (!valid_rotation_axis[mu][nu]) continue;

                                        // Search for a new constraint below

                                        for (j = 0; j < nparam_sub; ++j) arr_constraint[j] = 0.0;

                                        // Loop for m_{N+1}, a_{N+1}
                                        for (std::vector<int>::iterator iter_list = interaction_list.begin();
                                             iter_list != interaction_list.end(); ++iter_list) {
                                            jat = *iter_list;

                                            interaction_atom[order + 1] = jat;
                                            if (!interaction->is_incutoff(order + 2, interaction_atom, order)) continue;

                                            atom_tmp.clear();

                                            for (j = 1; j < order + 2; ++j) {
                                                atom_tmp.push_back(interaction_atom[j]);
                                            }
                                            std::sort(atom_tmp.begin(), atom_tmp.end());

                                            iter_cluster = interaction->mindist_cluster[order][i].find(
                                                MinimumDistanceCluster(atom_tmp, cell_dummy));
                                            if (iter_cluster != interaction->mindist_cluster[order][i].end()) {

                                                int iloc = -1;

                                                for (j = 0; j < atom_tmp.size(); ++j) {
                                                    if (atom_tmp[j] == jat) {
                                                        iloc = j;
                                                        break;
                                                    }
                                                }

                                                if (iloc == -1) {
                                                    error->exit("rotational_invariance", "This cannot happen.");
                                                }

                                                for (j = 0; j < 3; ++j) vec_for_rot[j] = 0.0;

                                                int nsize_equiv = (*iter_cluster).cell.size();

                                                for (j = 0; j < nsize_equiv; ++j) {
                                                    for (int k = 0; k < 3; ++k) {
                                                        vec_for_rot[k] += interaction->x_image[(*iter_cluster).cell[j][iloc]][jat][k];
                                                    }
                                                }

                                                for (j = 0; j < 3; ++j) {
                                                    vec_for_rot[j] /= static_cast<double>(nsize_equiv);
                                                }
                                            }


                                            // mu, nu

                                            interaction_index[order + 1] = 3 * jat + mu;
                                            for (j = 0; j < order + 2; ++j) interaction_tmp[j] = interaction_index[j];

                                            fcs->sort_tail(order + 2, interaction_tmp);

                                            iter_found = list_found.find(FcProperty(order + 2, 1.0, interaction_tmp, 1));
                                            if (iter_found != list_found.end()) {
                                                arr_constraint[nparams[order - 1] + (*iter_found).mother]
                                                    += (*iter_found).coef * vec_for_rot[nu];
                                            }

                                            // Exchange mu <--> nu and repeat again.

                                            interaction_index[order + 1] = 3 * jat + nu;
                                            for (j = 0; j < order + 2; ++j) interaction_tmp[j] = interaction_index[j];

                                            fcs->sort_tail(order + 2, interaction_tmp);

                                            iter_found = list_found.find(FcProperty(order + 2, 1.0, interaction_tmp, 1));
                                            if (iter_found != list_found.end()) {
                                                arr_constraint[nparams[order - 1] + (*iter_found).mother]
                                                    -= (*iter_found).coef * vec_for_rot[mu];
                                            }
                                        }

                                        for (lambda = 0; lambda < order + 1; ++lambda) {

                                            mu_lambda = interaction_index[lambda] % 3;

                                            for (jcrd = 0; jcrd < 3; ++jcrd) {

                                                for (j = 0; j < order + 1; ++j) interaction_tmp[j] = interaction_index[j];

                                                interaction_tmp[lambda] = 3 * interaction_atom[lambda] + jcrd;

                                                levi_factor = 0;

                                                for (j = 0; j < 3; ++j) {
                                                    levi_factor += levi_civita(j, mu, nu) * levi_civita(j, mu_lambda, jcrd);
                                                }

                                                if (levi_factor == 0) continue;

                                                fcs->sort_tail(order + 1, interaction_tmp);

                                                iter_found = list_found_last.find(FcProperty(order + 1, 1.0,
                                                                                             interaction_tmp, 1));
                                                if (iter_found != list_found_last.end()) {
                                                    arr_constraint[(*iter_found).mother]
                                                        += (*iter_found).coef * static_cast<double>(levi_factor);
                                                }
                                            }
                                        }

                                        if (!is_allzero(nparam_sub, arr_constraint)) {

                                            // A Candidate for another constraint found !
                                            // Add to the appropriate set

                                            if (is_allzero(nparam_sub, arr_constraint, nparams[order - 1])) {
                                                const_rotation_self[order - 1].push_back(
                                                    ConstraintClass(nparams[order - 1], arr_constraint));
                                            } else if (is_allzero(nparams[order - 1], arr_constraint)) {
                                                const_rotation_self[order].push_back(
                                                    ConstraintClass(nparam_sub, arr_constraint, nparams[order - 1]));
                                            } else {
                                                const_rotation_cross[order].push_back(
                                                    ConstraintClass(nparam_sub, arr_constraint));
                                            }
                                        }

                                    } // nu
                                } // mu

                            } // ixyz

                        } while (g.next());

                    } // direction
                } // icrd
            }

            // Additional constraint for the last order.
            // All IFCs over maxorder-th order are neglected.

            if (order == maxorder - 1 && !exclude_last_R) {

                nxyz2 = static_cast<int>(pow(static_cast<double>(3), order + 1));
                memory->allocate(xyzcomponent2, nxyz2, order + 1);
                fcs->get_xyzcomponent(order + 1, xyzcomponent2);

                for (icrd = 0; icrd < 3; ++icrd) {

                    interaction_index[0] = 3 * interaction_atom[0] + icrd;

                    CombinationWithRepetition<int> g_now(interaction_list_now.begin(),
                                                         interaction_list_now.end(), order + 1);
                    do {

                        std::vector<int> data = g_now.now();

                        for (int idata = 0; idata < data.size(); ++idata)
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

                                            fcs->sort_tail(order + 2, interaction_tmp);

                                            iter_found = list_found.find(FcProperty(order + 2, 1.0,
                                                                                    interaction_tmp, 1));
                                            if (iter_found != list_found.end()) {
                                                arr_constraint_self[(*iter_found).mother]
                                                    += (*iter_found).coef * static_cast<double>(levi_factor);
                                            }
                                        } // jcrd
                                    } // lambda

                                    if (!is_allzero(nparams[order], arr_constraint_self)) {
                                        const_rotation_self[order].push_back(
                                            ConstraintClass(nparams[order], arr_constraint_self));
                                    }

                                } // nu
                            } // mu

                        } // ixyz


                    } while (g_now.next());

                } // icrd

                memory->deallocate(xyzcomponent2);
            }
        } // iat

        std::cout << " done." << std::endl;

        if (order > 0) {
            memory->deallocate(xyzcomponent);
        }
        memory->deallocate(arr_constraint);
        memory->deallocate(arr_constraint_self);
        memory->deallocate(interaction_tmp);
        memory->deallocate(interaction_index);
        memory->deallocate(interaction_atom);
    } // order

    for (order = 0; order < maxorder; ++order) {
        remove_redundant_rows(nparam_sub, const_rotation_cross[order], eps6);
        remove_redundant_rows(nparams[order], const_rotation_self[order], eps6);
    }

    std::cout << "  Finished !" << std::endl << std::endl;

    memory->deallocate(ind);
    memory->deallocate(nparams);
}

void Constraint::remove_redundant_rows(const int n,
                                       std::vector<ConstraintClass> &Constraint_vec,
                                       const double tolerance)
{
#ifdef _USE_EIGEN_DISABLED

    // This function doesn't make the reduced row echelon form of the constraint matrix.
    // It just returns the image of the matrix, though they are similar.

    using namespace Eigen;
    int i;
    int nrow = n;
    int ncol = Constraint_vec.size();
    double *arr_tmp;
    double **mat;

    if (ncol > 0) {
        memory->allocate(arr_tmp, nrow);
        MatrixXd mat_tmp(nrow, ncol);

        int icol = 0;

        for (std::vector<ConstraintClass>::iterator p = Constraint_vec.begin(); p != Constraint_vec.end(); ++p) {
            ConstraintClass const_now = *p;
            for (i = 0; i < nrow; ++i) {
                mat_tmp(i, icol) = const_now.w_const[i];
            }
            ++icol;
        }

        FullPivLU<MatrixXd> lu_decomp(mat_tmp);
        lu_decomp.setThreshold(tolerance);
        int nrank = lu_decomp.rank();
        MatrixXd c_reduced = lu_decomp.image(mat_tmp);

        memory->allocate(mat, nrank, nrow);

        i = 0;
        for (icol = 0; icol < nrank; ++icol) {
            for (int irow = 0; irow < nrow; ++irow) {
                mat[icol][irow] = c_reduced(irow, icol);
            }
        }
    //      std::cout << "nrank = " << nrank << std::flush << std::endl;
        rref(nrank, nrow, mat, nrank, tolerance);

        Constraint_vec.clear();

        for (i = 0; i < nrank; ++i) {
            for (int j = 0; j < i; ++j) arr_tmp[j] = 0.0;

            for (int j = i; j < nrow; ++j) {
                arr_tmp[j] = mat[i][j];
            }
            Constraint_vec.push_back(ConstraintClass(nrow, arr_tmp));
        }

        memory->deallocate(mat);

    // for (icol = 0; icol < nrank; ++icol) {
    //     for (int irow = 0; irow < nrow; ++irow) {
    //         arr_tmp[irow] = c_reduced(irow, icol);
    //     }

    //     Constraint_Set.insert(ConstraintClass(nrow, arr_tmp));
    // }

        memory->deallocate(arr_tmp);
    }
#else 
    int i, j;

    int nparam = n;
    int nconst = Constraint_vec.size();
    double *arr_tmp;
    double **mat_tmp;

    int nrank;

    if (nconst > 0) {

        memory->allocate(mat_tmp, nconst, nparam);

        i = 0;

        for (std::vector<ConstraintClass>::iterator p = Constraint_vec.begin();
             p != Constraint_vec.end(); ++p) {
            for (j = 0; j < nparam; ++j) {
                mat_tmp[i][j] = (*p).w_const[j];
            }
            ++i;
        }

        rref(nconst, nparam, mat_tmp, nrank, tolerance);

        /*
        // 		// Transpose matrix A 

        memory->allocate(arr_tmp, nconst * nparam);

        k = 0;

        for (j = 0; j < nparam; ++j) {
        for (i = 0; i < nconst; ++i) {
        arr_tmp[k++] = mat_tmp[i][j];
        }
        }

        // Perform LU decomposition

        int nmin = std::min<int>(nconst, nparam);
        memory->allocate(ipiv, nmin);

        dgetrf_(&nconst, &nparam, arr_tmp, &nconst, ipiv, &INFO);

        k = 0;

        for (j = 0; j < nparam; ++j) {
        for (i = 0; i < nconst; ++i) {
        mat_tmp[i][j] = arr_tmp[k++];
        }
        }

        memory->deallocate(arr_tmp);
        memory->deallocate(ipiv);
        */

        memory->allocate(arr_tmp, nparam);

        Constraint_vec.clear();

        for (i = 0; i < nrank; ++i) {
            for (j = 0; j < i; ++j) arr_tmp[j] = 0.0;

            for (j = i; j < nparam; ++j) {
                arr_tmp[j] = mat_tmp[i][j];
            }
            Constraint_vec.push_back(ConstraintClass(nparam, arr_tmp));
        }

        memory->deallocate(mat_tmp);
        memory->deallocate(arr_tmp);
    }

#endif
}

int Constraint::levi_civita(const int i, const int j, const int k)
{
    int epsilon = (j - i) * (k - i) * (k - j) / 2;
    return epsilon;
}

bool Constraint::is_allzero(const int n, const double *arr, const int nshift)
{
    for (int i = nshift; i < n; ++i) {
        if (std::abs(arr[i]) > eps10) {
            return false;
        }
    }
    return true;
}

bool Constraint::is_allzero(const std::vector<int> vec, int &loc)
{
    loc = -1;
    for (int i = 0; i < vec.size(); ++i) {
        if (std::abs(vec[i]) > 0) {
            loc = i;
            return false;
        }
    }
    return true;
}

bool Constraint::is_allzero(const std::vector<double> vec, int &loc)
{
    loc = -1;
    for (int i = 0; i < vec.size(); ++i) {
        if (std::abs(vec[i]) > eps) {
            loc = i;
            return false;
        }
    }
    return true;
}

void Constraint::setup_rotation_axis(bool flag[3][3])
{
    unsigned int mu, nu;

    for (mu = 0; mu < 3; ++mu) {
        for (nu = 0; nu < 3; ++nu) {
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
        error->warn("setup_rotation_axis",
                    "Invalid rotation_axis. Default value(xyz) will be used.");
    }
}


void Constraint::rref(int nrows,
                      int ncols,
                      double **mat,
                      int &nrank,
                      double tolerance)
{
    // Return the reduced row echelon form (rref) of matrix mat.
    // In addition, rank of the matrix is estimated.

    int irow, icol, jrow, jcol;
    int pivot;
    double tmp, *arr;

    memory->allocate(arr, ncols);

    nrank = 0;

    icol = 0;

    for (irow = 0; irow < nrows; ++irow) {

        pivot = irow;

        while (std::abs(mat[pivot][icol]) < tolerance) {
            ++pivot;

            if (pivot == nrows) {
                pivot = irow;
                ++icol;

                if (icol == ncols) break;
            }
        }

        if (icol == ncols) break;

        if (std::abs(mat[pivot][icol]) > tolerance) ++nrank;

        if (pivot != irow) {
#pragma omp parallel for private(tmp)
            for (jcol = icol; jcol < ncols; ++jcol) {
                tmp = mat[pivot][jcol];
                mat[pivot][jcol] = mat[irow][jcol];
                mat[irow][jcol] = tmp;
            }
        }

        tmp = mat[irow][icol];
        tmp = 1.0 / tmp;
#pragma omp parallel for
        for (jcol = icol; jcol < ncols; ++jcol) {
            mat[irow][jcol] *= tmp;
        }

        for (jrow = 0; jrow < nrows; ++jrow) {
            if (jrow == irow) continue;

            tmp = mat[jrow][icol];
#pragma omp parallel for
            for (jcol = icol; jcol < ncols; ++jcol) {
                mat[jrow][jcol] -= tmp * mat[irow][jcol];
            }
        }
    }

    memory->deallocate(arr);
}
