/*
 fcs.h

 Copyright (c) 2014--2017 Terumasa Tadano


 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory
 or http://opensource.org/licenses/mit-license.php for information.
*/

#pragma once

#include <vector>
#include <set>
#include <map>
#include <algorithm>
#include "cluster.h"
#include "symmetry.h"
#include "timer.h"

using ConstraintSparseForm = std::vector<std::map<size_t, double>>;

namespace ALM_NS
{
    class FcProperty
    {
    public:
        std::vector<int> elems; // flattened index of (iatom, icoordinate) in the supercell
        double sign;            // factor (+1 or -1) to convert the mother FC to the child
        size_t mother;          // index of the reducible force constants

        FcProperty();

        FcProperty(const FcProperty &obj) = default;

        FcProperty(const int n,
                   const double c,
                   const int *arr,
                   const size_t m)
        {
            sign = c;
            mother = m;
            for (auto i = 0; i < n; ++i) {
                elems.push_back(arr[i]);
            }
        }

        bool operator<(const FcProperty &a) const
        {
            return std::lexicographical_compare(elems.begin(), elems.end(),
                                                a.elems.begin(), a.elems.end());
        }

        bool operator==(const FcProperty &a) const
        {
            auto n = elems.size();
            auto n_ = a.elems.size();
            if (n != n_) return false;
            for (size_t i = 0; i < n; ++i) {
                if (elems[i] != a.elems[i]) return false;
            }
            return true;
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

    class Fcs
    {
    public:
        Fcs();
        ~Fcs();

        void init(const Cluster *cluster,
                  const Symmetry *symmetry,
                  const size_t number_of_atoms,
                  const int verbosity,
                  Timer *timer);

        void get_xyzcomponent(int,
                              int **) const;
        void generate_force_constant_table(const int,
                                           const size_t nat,
                                           const std::set<IntList> &,
                                           const Symmetry *,
                                           const std::string,
                                           std::vector<FcProperty> &,
                                           std::vector<size_t> &,
                                           std::vector<FcProperty> &,
                                           const bool) const;

        void get_constraint_symmetry(const size_t nat,
                                     const Symmetry *symmetry,
                                     const int order,
                                     const std::string basis,
                                     const std::vector<FcProperty> &fc_table_in,
                                     const size_t nparams,
                                     const double tolerance,
                                     ConstraintSparseForm &const_out,
                                     const bool do_rref = false) const;

        std::vector<size_t>* get_nequiv() const;
        std::vector<FcProperty>* get_fc_table() const;

    private:
        std::vector<size_t> *nequiv;       // stores duplicate number of irreducible force constants
        std::vector<FcProperty> *fc_table; // all force constants
        std::vector<FcProperty> *fc_zeros; // zero force constants (due to space group symm.)

        bool store_zeros;
        void set_default_variables();
        void deallocate_variables();
        bool is_ascending(int,
                          const int *) const;
        bool is_inprim(const int n,
                       const int *arr,
                       const size_t natmin,
                       const std::vector<std::vector<int>> &map_p2s) const;
        bool is_inprim(const int n,
                       const size_t natmin,
                       const std::vector<std::vector<int>> &map_p2s) const;
        bool is_allzero(const std::vector<double> &,
                        double,
                        int &) const;
        void get_available_symmop(const size_t nat,
                                  const Symmetry *symmetry,
                                  const std::string basis,
                                  int &nsym_avail,
                                  int **mapping_symm,
                                  double ***rotation,
                                  const bool use_compatible) const;
        int get_minimum_index_in_primitive(const int n,
                                           const int *arr,
                                           const size_t nat,
                                           const size_t natmin,
                                           const std::vector<std::vector<int>> &map_p2s) const;
        double coef_sym(const int,
                        const double *const *,
                        const int *,
                        const int *) const;
    };
}

// Define a hash function for FcProperty class
// Use boost::hash_combine
namespace std
{
    template <>
    struct hash<ALM_NS::FcProperty>
    {
        std::size_t operator ()(ALM_NS::FcProperty const &obj) const
        {
            hash<int> hasher;
            size_t seed = 0;
            for (auto i : obj.elems) {
                seed ^= hasher(i) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
            }
            return seed;
        }
    };
}
