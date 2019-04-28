/*
 constraint.h

 Copyright (c) 2014, 2015, 2016 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory
 or http://opensource.org/licenses/mit-license.php for information.
*/

#pragma once

#include <boost/bimap.hpp>
#include <utility>
#include <vector>
#include <string>
#include <iomanip>
#include <map>
#include "constants.h"
#include "fcs.h"
#include "cluster.h"
#include "system.h"
#include "timer.h"

namespace ALM_NS
{
    class ConstraintClass
    {
    public:
        std::vector<double> w_const;

        ConstraintClass() = default;

        ConstraintClass(const ConstraintClass &a) = default;

        ConstraintClass(std::vector<double> vec) : w_const(std::move(vec)) { }

        ConstraintClass(const int n,
                        const double *arr,
                        const int nshift = 0)
        {
            for (auto i = nshift; i < n; ++i) {
                w_const.push_back(arr[i]);
            }
        }

        bool operator<(const ConstraintClass &a) const
        {
            return std::lexicographical_compare(w_const.begin(), w_const.end(),
                                                a.w_const.begin(), a.w_const.end());
        }
    };

    class ConstraintTypeFix
    {
    public:
        size_t p_index_target;
        double val_to_fix;

        ConstraintTypeFix(const size_t index_in,
                          const double val_in) :
            p_index_target(index_in), val_to_fix(val_in) { }
    };

    class ConstraintTypeRelate
    {
    public:
        size_t p_index_target;
        std::vector<double> alpha;
        std::vector<size_t> p_index_orig;

        ConstraintTypeRelate(const size_t index_in,
                             std::vector<double> alpha_in,
                             std::vector<size_t> p_index_in) :
            p_index_target(index_in), alpha(std::move(alpha_in)), p_index_orig(std::move(p_index_in)) { }
    };

    inline bool equal_within_eps12(const std::vector<double> &a,
                                   const std::vector<double> &b)
    {
        const auto n = a.size();
        const auto m = b.size();
        if (n != m) return false;
        for (size_t i = 0; i < n; ++i) {
            if (std::abs(a[i] - b[i]) > eps12) return false;
        }
        return true;
    }

    class ConstraintIntegerElement
    {
        // For sparse representation
    public:
        size_t col;
        int val;

        ConstraintIntegerElement(const size_t col_in,
                                 const int val_in) :
            col(col_in), val(val_in) {}
    };

    // Operator for sort
    inline bool operator<(const std::vector<ConstraintIntegerElement> &obj1,
                          const std::vector<ConstraintIntegerElement> &obj2)
    {
        const auto len1 = obj1.size();
        const auto len2 = obj2.size();
        const auto min = (std::min)(len1, len2);

        for (size_t i = 0; i < min; ++i) {
            if (obj1[i].col < obj2[i].col) {
                return true;
            }
            if (obj1[i].col > obj2[i].col) {
                return false;
            }
            if (obj1[i].val < obj2[i].val) {
                return true;
            }
            if (obj1[i].val > obj2[i].val) {
                return false;
            }
        }
        return false;
    }

    // Operator for unique
    inline bool operator==(const std::vector<ConstraintIntegerElement> &obj1,
                           const std::vector<ConstraintIntegerElement> &obj2)
    {
        const auto len1 = obj1.size();
        const auto len2 = obj2.size();
        if (len1 != len2) return false;

        for (size_t i = 0; i < len1; ++i) {
            if (obj1[i].col != obj2[i].col || obj1[i].val != obj2[i].val) {
                return false;
            }
        }
        return true;
    }

    class ConstraintDoubleElement
    {
        // For sparse representation
    public:
        size_t col;
        double val;

        ConstraintDoubleElement(const size_t col_in,
                                const double val_in) :
            col(col_in), val(val_in) {}
    };

    // Operator for sort
    inline bool operator<(const std::vector<ConstraintDoubleElement> &obj1,
                          const std::vector<ConstraintDoubleElement> &obj2)
    {
        const auto len1 = obj1.size();
        const auto len2 = obj2.size();
        const auto min = (std::min)(len1, len2);

        for (size_t i = 0; i < min; ++i) {
            if (obj1[i].col < obj2[i].col) {
                return true;
            }
            if (obj1[i].col > obj2[i].col) {
                return false;
            }
            if (obj1[i].val < obj2[i].val) {
                return true;
            }
            if (obj1[i].val > obj2[i].val) {
                return false;
            }
        }
        return false;
    }

    // Operator for unique
    inline bool operator==(const std::vector<ConstraintDoubleElement> &obj1,
                           const std::vector<ConstraintDoubleElement> &obj2)
    {
        const auto len1 = obj1.size();
        const auto len2 = obj2.size();
        if (len1 != len2) return false;

        for (size_t i = 0; i < len1; ++i) {
            if (obj1[i].col != obj2[i].col || (std::abs(obj1[i].val - obj2[i].val) > 1.0e-10)) {
                return false;
            }
        }
        return true;
    }

    inline bool operator<(const std::map<size_t, double> &obj1,
                          const std::map<size_t, double> &obj2)
    {
        return obj1.begin()->first < obj2.begin()->first;
    }

    class Constraint
    {
    public:
        Constraint();
        ~Constraint();

        void setup(const System *system,
                   const Fcs *fcs,
                   const Cluster *cluster,
                   const Symmetry *symmetry,
                   const std::string alm_mode,
                   const int verbosity,
                   Timer *timer);

        void get_mapping_constraint(const int nmax,
                                    const std::vector<size_t> *nequiv,
                                    const ConstraintSparseForm *const_in,
                                    std::vector<ConstraintTypeFix> *const_fix_out,
                                    std::vector<ConstraintTypeRelate> *const_relate_out,
                                    boost::bimap<size_t, size_t> *index_bimap_out) const;

        int get_constraint_mode() const;
        void set_constraint_mode(const int);
        size_t get_number_of_constraints() const;
        std::string get_fc_file(const int) const;
        void set_fc_file(const int,
                         const std::string);
        bool get_fix_harmonic() const;
        void set_fix_harmonic(const bool);
        bool get_fix_cubic() const;
        void set_fix_cubic(const bool);
        int get_constraint_algebraic() const;

        double** get_const_mat() const;
        double* get_const_rhs() const;

        double get_tolerance_constraint() const;
        void set_tolerance_constraint(const double);

        bool get_exist_constraint() const;
        bool get_extra_constraint_from_symmetry() const;

        std::string get_rotation_axis() const;
        void set_rotation_axis(const std::string);

        const ConstraintSparseForm& get_const_symmetry(const int) const;
        const std::vector<ConstraintTypeFix>& get_const_fix(const int) const;
        void set_const_fix_val_to_fix(const int order,
                                      const size_t idx,
                                      const double val);
        const std::vector<ConstraintTypeRelate>& get_const_relate(const int) const;
        const boost::bimap<size_t, size_t>& get_index_bimap(const int) const;

    private:

        int constraint_mode;
        size_t number_of_constraints;
        std::string fc2_file, fc3_file;
        bool fix_harmonic, fix_cubic;
        int constraint_algebraic;

        double **const_mat;
        double *const_rhs;

        double tolerance_constraint;

        bool exist_constraint;
        bool extra_constraint_from_symmetry;

        std::string rotation_axis;
        ConstraintSparseForm *const_symmetry;
        std::vector<ConstraintTypeFix> *const_fix;
        std::vector<ConstraintTypeRelate> *const_relate;
        std::vector<ConstraintTypeRelate> *const_relate_rotation;
        boost::bimap<size_t, size_t> *index_bimap;

        bool impose_inv_T, impose_inv_R, exclude_last_R;

        ConstraintSparseForm *const_translation;
        ConstraintSparseForm *const_rotation_self;
        ConstraintSparseForm *const_rotation_cross;
        ConstraintSparseForm *const_self;

        void set_default_variables();
        void deallocate_variables();

        int levi_civita(const int,
                        const int,
                        const int) const;

        void generate_rotational_constraint(const System *,
                                            const Symmetry *,
                                            const Cluster *,
                                            const Fcs *,
                                            const int,
                                            const double);

        // const_mat and const_rhs are updated.
        size_t calc_constraint_matrix(const int maxorder,
                                      const std::vector<size_t> *nequiv,
                                      const size_t nparams) const;

        void print_constraint(const ConstraintSparseForm &) const;

        void setup_rotation_axis(bool [3][3]);
        bool is_allzero(const int,
                        const double *,
                        const int nshift = 0) const;
        bool is_allzero(const std::vector<int> &,
                        int &) const;
        bool is_allzero(const std::vector<double> &,
                        const double,
                        int &,
                        const int nshift = 0) const;


        void remove_redundant_rows(const size_t n,
                                   std::vector<ConstraintClass> &Constraint_vec,
                                   const double tolerance = eps12) const;

        // const_symmetry is updated.
        void generate_symmetry_constraint_in_cartesian(const size_t nat,
                                                       const Symmetry *symmetry,
                                                       const Cluster *cluster,
                                                       const Fcs *fcs,
                                                       const int verbosity) const;

        void get_constraint_translation(const Cell &supercell,
                                        const Symmetry *symmetry,
                                        const Cluster *cluster,
                                        const Fcs *fcs,
                                        const int order,
                                        const std::vector<FcProperty> &fc_table,
                                        const size_t nparams,
                                        ConstraintSparseForm &const_out,
                                        const bool do_rref = false) const;

        // const_translation is updated.
        void generate_translational_constraint(const Cell &,
                                               const Symmetry *,
                                               const Cluster *,
                                               const Fcs *,
                                               const int) const;

        void fix_forceconstants_to_file(const int,
                                        const Symmetry *,
                                        const Fcs *,
                                        const std::string,
                                        std::vector<ConstraintTypeFix> &) const;
    };
}
