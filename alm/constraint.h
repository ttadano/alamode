/*
 constraint.h

 Copyright (c) 2014, 2015, 2016 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#pragma once

#include <iostream>
#include <vector>
#include <set>
#include <string>
#include "pointers.h"
#include "constants.h"
#include "interaction.h"
#include "symmetry.h"
#include "fcs.h"
#include <boost/bimap.hpp>

namespace ALM_NS
{
    class ConstraintClass
    {
    public:
        std::vector<double> w_const;

        ConstraintClass();

        ConstraintClass(const ConstraintClass &a)
        {
            for (std::vector<double>::const_iterator p = a.w_const.begin();
                 p != a.w_const.end(); ++p) {
                w_const.push_back(*p);
            }
        }

        ConstraintClass(const int n, const double *arr, const int nshift = 0)
        {
            for (int i = nshift; i < n; ++i) {
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
        unsigned int p_index_target;
        double val_to_fix;

        ConstraintTypeFix(const unsigned int index_in, const double val_in)
        {
            p_index_target = index_in;
            val_to_fix = val_in;
        }
    };

    class ConstraintTypeRelate
    {
    public:
        unsigned int p_index_target;
        std::vector<double> alpha;
        std::vector<unsigned int> p_index_orig;

        ConstraintTypeRelate(const unsigned int index_in,
                             const std::vector<double> alpha_in,
                             const std::vector<unsigned int> p_index_in)
        {
            p_index_target = index_in;
            std::copy(alpha_in.begin(), alpha_in.end(), std::back_inserter(alpha));
            std::copy(p_index_in.begin(), p_index_in.end(), std::back_inserter(p_index_orig));
        }
    };

    inline bool equal_within_eps12(const std::vector<double> &a, const std::vector<double> &b)
    {
        int n = a.size();
        int m = b.size();
        if (n != m) return false;
        double res = 0.0;
        for (int i = 0; i < n; ++i) {
            if (std::abs(a[i] - b[i]) > eps12) return false;
        }
        //        if (std::sqrt(res)>eps12) return false;
        return true;
    }

    class Constraint: protected Pointers
    {
    public:
        Constraint(class ALM *);
        ~Constraint();

        void setup();

        int constraint_mode;
        int P;
        std::string fc2_file, fc3_file;
        bool fix_harmonic, fix_cubic;
        int constraint_algebraic;

        double **const_mat;
        double *const_rhs;
        double tolerance_constraint;

        bool exist_constraint;
        bool extra_constraint_from_symmetry;
        std::string rotation_axis;
        std::vector<ConstraintClass> *const_symmetry;

        std::vector<ConstraintTypeFix> *const_fix;
        std::vector<ConstraintTypeRelate> *const_relate;
        boost::bimap<int, int> *index_bimap;

        //void constraint_from_symmetry(std::vector<ConstraintClass> *);
        void get_symmetry_constraint(const int, const std::set<IntList>,
                                     const std::vector<SymmetryOperation>,
                                     const std::string,
                                     const std::vector<FcProperty>,
                                     const std::vector<int>,
                                     std::vector<ConstraintClass> &);

        void get_mapping_constraint(const int, std::vector<int> *,
                                    std::vector<ConstraintClass> *,
                                    std::vector<ConstraintTypeFix> *,
                                    std::vector<ConstraintTypeRelate> *,
                                    boost::bimap<int, int> *, const bool);

    private:

        bool impose_inv_T, impose_inv_R, exclude_last_R;

        std::vector<ConstraintClass> *const_translation;
        std::vector<ConstraintClass> *const_rotation_self;
        std::vector<ConstraintClass> *const_rotation_cross;

        std::vector<ConstraintClass> *const_self;

        int levi_civita(const int, const int, const int);

        void translational_invariance();
        void rotational_invariance();
        void calc_constraint_matrix(const int, int &);

        void setup_rotation_axis(bool [3][3]);
        bool is_allzero(const int, const double *, const int nshift = 0);
        bool is_allzero(const std::vector<int>, int &);
        bool is_allzero(const std::vector<double>, const double, int &);

        void remove_redundant_rows(const int, std::vector<ConstraintClass> &,
                                   const double tolerance = eps12);

        void rref(int, int, double **, int &, double tolerance = eps12);
        void rref(std::vector<std::vector<double>> &, const double tolerance = eps12);

        void generate_symmetry_constraint_in_cartesian(std::vector<ConstraintClass> *);
    };

    extern "C"
    {
        void dgetrf_(int *m, int *n, double *a, int *lda, int *ipiv, int *info);
    }
}
