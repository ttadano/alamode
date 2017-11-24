/*
 fcs.h

 Copyright (c) 2014--2017 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#pragma once

#include "pointers.h"
#include <vector>
#include <set>
#include <algorithm>
#include "symmetry.h"
#include "interaction.h"

namespace ALM_NS
{
    class FcProperty
    {
    public:
        std::vector<int> elems; // flattened index of (iatom, icoordinate) in the supercell
        double sign; // factor (+1 or -1) to convert the mother FC to the child
        int mother;

        FcProperty();

        FcProperty(const FcProperty &obj)
        {
            sign = obj.sign;
            mother = obj.mother;
            for (auto it = obj.elems.begin(); it != obj.elems.end(); ++it) {
                elems.push_back(*it);
            }
        };

        FcProperty(const int n, const double c, const int *arr, const int m)
        {
            sign = c;
            mother = m;
            for (int i = 0; i < n; ++i) {
                elems.push_back(arr[i]);
            }
        }

        bool operator<(const FcProperty &a) const
        {
            return std::lexicographical_compare(elems.begin(), elems.end(),
                                                a.elems.begin(), a.elems.end());
        }
    };

    class ForceConstantTable
    {
    public:
        double fc_value;
        int multiplicity;
        std::vector<FcProperty> fclist;
        ForceConstantTable();
    };

    class Fcs: protected Pointers
    {
    public:
        Fcs(class ALM *);
        ~Fcs();

        void init();

        std::vector<int> *nequiv;
        std::vector<FcProperty> *fc_table;
        std::vector<FcProperty> *fc_zeros;

        std::string easyvizint(const int);
        void get_xyzcomponent(int, int **);
        void sort_tail(const int, int *);

        bool is_inprim(const int, const int *);
        bool is_inprim(const int);
        int min_inprim(const int, const int *);
        double coef_sym(const int, const int, const int *, const int *);
        double coef_sym(const int, double **, const int *, const int *);

        void generate_force_constant_table(const int,
                                           const std::set<IntList>,
                                           const std::vector<SymmetryOperation>,
                                           std::string,
                                           std::vector<FcProperty> &,
                                           std::vector<int> &,
                                           std::vector<FcProperty> &,
                                           const bool);

    private:

        bool is_ascending(const int, const int *);
    };
}
