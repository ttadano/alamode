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
#include <unordered_map>
#include <algorithm>
#include "cluster.h"
#include "symmetry.h"
#include "timer.h"
#include <Eigen/Core>

// By default, we use unordered_map for better performance of rref_sparse (in rref.cpp).
#ifdef _USE_MAP_FOR_CONSTRAINT
using MapConstraintElement = std::map<size_t, double>;
#else
using MapConstraintElement = std::unordered_map<size_t, double>;
#endif
using ConstraintSparseForm = std::vector<MapConstraintElement>;

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

        static bool compare_atom_index(const FcProperty &a,
                                       const FcProperty &b)
        {
            const auto n1 = a.elems.size();
            const auto n2 = b.elems.size();
            if (n1 != n2) return n1 < n2;

            std::vector<int> atom_array1(n1), atom_array2(n2);

            for (auto i = 0; i < n1; ++i) {
                atom_array1[i] = a.elems[i] / 3;
                atom_array2[i] = b.elems[i] / 3;
            }

            return std::lexicographical_compare(atom_array1.begin(), atom_array1.end(),
                                                atom_array2.begin(), atom_array2.end());
        }
    };

    class ForceConstantTable
    {
    public:
        double fc_value;
        std::vector<int> atoms, coords, flattenarray;
        bool is_ascending_order; // true if the elements except the first element is sorted in ascending order.
        ForceConstantTable();

        ForceConstantTable(const ForceConstantTable &obj) = default;

        ForceConstantTable(const int nelems,
                           const double fc_in,
                           const int *atoms_in,
                           const int *coords_in)
        {
            atoms.resize(nelems);
            coords.resize(nelems);
            flattenarray.resize(nelems);
            fc_value = fc_in;
            for (auto i = 0; i < nelems; ++i) {
                atoms[i] = atoms_in[i];
                coords[i] = coords_in[i];
                flattenarray[i] = 3 * atoms_in[i] + coords_in[i];
            }
            is_ascending_order = true;

            for (auto i = 1; i < nelems - 1; ++i) {
                if (flattenarray[i] > flattenarray[i + 1]) {
                    is_ascending_order = false;
                    break;
                }
            }
        }

        bool operator<(const ForceConstantTable &a) const
        {
            return std::lexicographical_compare(flattenarray.begin(), flattenarray.end(),
                                                a.flattenarray.begin(), a.flattenarray.end());
        }
    };

    class Fcs
    {
    public:
        Fcs();
        ~Fcs();

        void init(const Cluster *cluster,
                  const Symmetry *symmetry,
                  const Cell &supercell,
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

        void get_constraint_symmetry_in_integer(const size_t nat,
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
        std::vector<ForceConstantTable>* get_fc_cart() const;
        std::vector<size_t> get_nfc_cart(const int permutation) const;

        void set_forceconstant_basis(const std::string preferred_basis_in);
        std::string get_forceconstant_basis() const;
        Eigen::Matrix3d get_basis_conversion_matrix() const;

        void set_forceconstant_cartesian(const int maxorder,
                                         double *param_in);

    private:
        std::vector<size_t> *nequiv;       // stores duplicate number of irreducible force constants
        std::vector<FcProperty> *fc_table; // all force constants in preferred_basis
        std::vector<FcProperty> *fc_zeros; // zero force constants (due to space group symm.)

        std::vector<ForceConstantTable> *fc_cart; // all force constants in Cartesian coordinate
        std::vector<size_t> nfc_cart_permu; // Number of nonzero elements with permutation
        std::vector<size_t> nfc_cart_nopermu; // Number of nonzero elements without permutation

        std::string preferred_basis; // "Cartesian" or "Lattice"
        Eigen::Matrix3d basis_conversion_matrix;

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
        bool is_allzero(const std::vector<int> &,
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
                        const double * const *,
                        const int *,
                        const int *) const;

        void set_basis_conversion_matrix(const Cell &supercell);
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
