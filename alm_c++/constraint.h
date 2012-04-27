#pragma once

#include <iostream>
#include <vector>
#include <set>
#include <string>
#include "pointers.h"

namespace ALM_NS
{
    class ConstraintClass {
    public:
        std::vector<double> w_const;

        ConstraintClass();
        ConstraintClass(const ConstraintClass &a){
            for(std::vector<double>::const_iterator p = a.w_const.begin(); p != a.w_const.end(); ++p){
                w_const.push_back(*p);
            }
        }
        ConstraintClass(const int n, const double *arr, const int nshift=0){
            for(int i = nshift; i < n; ++i){
                w_const.push_back(arr[i]);
            }
        }
    };

    inline bool operator<(const ConstraintClass a, const ConstraintClass b){
        return std::lexicographical_compare(a.w_const.begin(), a.w_const.end(), b.w_const.begin(), b.w_const.end());
    }

    class Constraint: protected Pointers {
    public:
        Constraint(class ALM *);
        ~Constraint();

        void setup();

        int constraint_mode;
        int P;
        std::string fc2_file;

        double **const_mat;
        double *const_rhs;

    private:

        int levi_civita(const int, const int, const int);

        void translational_invariance();
        void rotational_invariance();
        void calc_constraint_matrix(const int, int &);
        
        std::set<ConstraintClass> *const_translation;
        std::set<ConstraintClass> *const_rotation_self;
        std::set<ConstraintClass> *const_rotation_cross;

        bool is_allzero(const int, const double *, const int nshift = 0);
        void remove_redundant_rows(const int, std::set<ConstraintClass> &);
    };

}
