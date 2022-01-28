/*
 cluster.h

 Copyright (c) 2014 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory
 or http://opensource.org/licenses/mit-license.php for information.
*/

#pragma once

#include <string>
#include <utility>
#include <vector>
#include <set>
#include <iterator>
#include <algorithm>
#include "constants.h"
#include "system.h"
#include "symmetry.h"
#include "timer.h"

#include "mathfunctions.h"


namespace ALM_NS {
class IntList {
public:
    std::vector<int> iarray;

    IntList() = default;

    ~IntList() = default;

    IntList(const IntList &a) = default;

    IntList(const int n,
            const int *arr)
    {
        for (auto i = 0; i < n; ++i) {
            iarray.push_back(arr[i]);
        }
    }

    bool operator<(const IntList &a) const
    {
        return std::lexicographical_compare(iarray.begin(), iarray.end(),
                                            a.iarray.begin(), a.iarray.end());
    }

    bool operator==(const IntList &a) const
    {
        const auto n = iarray.size();
        const auto n_ = a.iarray.size();
        if (n != n_) return false;
        for (size_t i = 0; i < n; ++i) {
            if (iarray[i] != a.iarray[i]) return false;
        }
        return true;
    }
};

class RelativeVectors {
public:
    int order;
    std::vector<std::vector<double>> relvecs_cartesian;
    std::vector<std::vector<double>> relvecs_fractional;

    RelativeVectors() = default;

    RelativeVectors(int order_in)
    {
        relvecs_cartesian.clear();
        order = order_in;
    }

    void clear_relvecs()
    {
        int i;
        for (i = 0; i < order + 2; i++) {
            relvecs_fractional[i].clear();
            relvecs_cartesian[i].clear();
        }
    }

    void make_fractional_from_cartesian(const double reciprocal_lat[3][3])
    {

        int i, j, xyz, xyz2;
        std::vector<double> vectmp;

        relvecs_fractional.clear();

        for (i = 0; i < order + 1; i++) {

            vectmp.clear();
            for (xyz = 0; xyz < 3; xyz++) {
                vectmp.push_back(0.0);
                for (xyz2 = 0; xyz2 < 3; xyz2++) {
                    vectmp[xyz] += reciprocal_lat[xyz][xyz2] / (2.0 * pi) * relvecs_cartesian[i][xyz2];
                }
            }

            relvecs_fractional.push_back(vectmp);
        }
    }

    bool is_equal(RelativeVectors &another, double threshold)
    {
        std::vector<int> is_checked(order + 1);
        int is_found = 0;
        int is_equal_flg;

        int i, j, xyz;
        for (i = 0; i < order + 1; i++) {
            is_found = 0;
            // find corresponding relative vector
            for (j = 0; j < order + 1; j++) {
                if (is_checked[j] == 1) {
                    continue;
                }
                is_equal_flg = 1;
                for (xyz = 0; xyz < 3; xyz++) {
                    if (std::fabs(relvecs_cartesian[i][xyz] - another.relvecs_cartesian[j][xyz]) > threshold) {
                        is_equal_flg = 0;
                        break;
                    }
                }
                if (is_equal_flg == 1) {
                    is_found = 1;
                    is_checked[j] = 1;
                    break;
                }
            }
            // if not found
            if (is_found == 0) {
                return false;
            }
        }
        return true;
    }
};

class DistInfo {
public:
    int cell;
    double dist;
    double relvec[3]{};

    DistInfo() = default;

    ~DistInfo() = default;

    DistInfo(const DistInfo &obj) = default;

    DistInfo(const int n,
             const double d,
             const double x[3])
    {
        cell = n;
        dist = d;
        for (auto i = 0; i < 3; ++i) relvec[i] = x[i];
    }

    bool operator<(const DistInfo &a) const
    {
        return dist < a.dist;
    }
};

class DistList
    // This class is used only in print_neighborlist. Can be replaced by a more generalic function.
{
public:
    size_t atom;
    double dist;

    DistList() = default;

    DistList(const size_t atom_,
             const double dist_) : atom(atom_), dist(dist_) {};

    bool operator<(const DistList &a) const
    {
        if (std::abs(dist - a.dist) > eps8) {
            return dist < a.dist;
        }
        return atom < a.atom;
    }
};

class MinDistList {
public:
    std::vector<int> cell;
    std::vector<double> dist;

    MinDistList() = default;

    MinDistList(std::vector<int> cell_in,
                std::vector<double> dist_in)
            : cell(std::move(cell_in)), dist(std::move(dist_in)) {};

    static bool compare_sum_distance(const MinDistList &a,
                                     const MinDistList &b)
    {
        double dist_a = 0;
        double dist_b = 0;

        for (auto i: a.dist) {
            dist_a += i;
        }
        for (auto i: b.dist) {
            dist_b += i;
        }
        return dist_a < dist_b;
    }

    static bool compare_max_distance(const MinDistList &a,
                                     const MinDistList &b)
    {
        // This function works properly when dvec_a.size() > 0 and dvec_b.size() > 0
        std::vector<double> dvec_a, dvec_b;
        std::copy(a.dist.begin(), a.dist.end(), std::back_inserter(dvec_a));
        std::copy(b.dist.begin(), b.dist.end(), std::back_inserter(dvec_b));
        const auto max_dist_a = *std::max_element(dvec_a.begin(), dvec_a.end());
        const auto max_dist_b = *std::max_element(dvec_b.begin(), dvec_b.end());

        return max_dist_a < max_dist_b;
    }
};

class InteractionCluster {
public:
    std::vector<int> atom;
    std::vector<std::vector<int>> cell;
    double distmax;

    InteractionCluster() = default;

    InteractionCluster(std::vector<int> atom_in,
                       std::vector<std::vector<int>> cell_in,
                       const double dist_in)
            : atom(std::move(atom_in)), cell(std::move(cell_in)), distmax(dist_in) {};

    InteractionCluster(std::vector<int> atom_in,
                       std::vector<std::vector<int>> cell_in)
            : atom(std::move(atom_in)), cell(std::move(cell_in)), distmax(0.0) {};


    bool operator<(const InteractionCluster &a) const
    {
        return lexicographical_compare(atom.begin(), atom.end(),
                                       a.atom.begin(), a.atom.end());
    }
};

class Cluster {
public:
    Cluster();

    ~Cluster();

    void init(const System *system,
              const Symmetry *symmetry,
              const int mirror_image_conv,
              const int verbosity,
              Timer *timer);

    bool satisfy_nbody_rule(const int nelem,
                            const int *arr,
                            const int order) const;

    bool is_incutoff(const int n,
                     const int *atomnumlist,
                     const size_t order,
                     const std::vector<int> &kd) const;

    void define(const int maxorder_in,
                const size_t nkd,
                const int *nbody_include_in,
                const double *cutoff_radii_in);

    int get_maxorder() const;

    int *get_nbody_include() const;

    std::string get_ordername(const unsigned int order) const;

    const std::set<IntList> &get_cluster_list(const unsigned int order) const;

    const std::vector<int> &get_interaction_pair(const unsigned int order,
                                                 const size_t atom_index) const;

    const std::set<InteractionCluster> &get_interaction_cluster(const unsigned int order,
                                                                const size_t atom_index) const;

private:

    int maxorder;
    int *nbody_include;
    double ***cutoff_radii;
    std::set<IntList> *cluster_list;
    std::vector<int> **interaction_pair; // List of atoms inside the cutoff radius for each order
    std::set<InteractionCluster> **interaction_cluster;

    std::vector<DistInfo> **distall;       // Distance of all pairs (i,j) under the PBC
    std::vector<DistInfo> **mindist_pairs; // All pairs (i,j) with the minimum distance
    // Interacting many-body clusters with mirrow image information

    void set_default_variables();

    void deallocate_variables();

    // can be made const function, but mindist_pairs is modified
    // in this function.
    void get_pairs_of_minimum_distance(const size_t nat,
                                       const double *const *const *xc_in,
                                       const int *exist) const;

    void generate_interaction_information_by_cutoff(const size_t nat,
                                                    const size_t natmin,
                                                    const std::vector<int> &kd,
                                                    const std::vector<std::vector<int>> &map_p2s,
                                                    const double *const *rc,
                                                    std::vector<int> *interaction_list) const;

    void set_interaction_by_cutoff(const size_t nat,
                                   const std::vector<int> &kd,
                                   const size_t nat_prim,
                                   const std::vector<std::vector<int>> &map_p2s) const;

    void print_neighborlist(const size_t,
                            const size_t,
                            const std::vector<std::vector<int>> &,
                            const std::vector<int> &,
                            const std::string *) const;

    void print_interaction_information(const size_t natmin,
                                       const std::vector<std::vector<int>> &map_p2s,
                                       const std::vector<int> &kd,
                                       const std::string *kdname,
                                       const std::vector<int> *const *interaction_list) const;

    double distance(const double *,
                    const double *) const;

    int nbody(const int,
              const int *) const;

    void calc_interaction_clusters(const size_t natmin,
                                   const std::vector<int> &kd,
                                   const std::vector<std::vector<int>> &map_p2s,
                                   const double *const *const *x_image,
                                   const int *exist,
                                   const int mirror_image_conv) const;

    void set_interaction_cluster(const int order,
                                 const size_t natmin,
                                 const std::vector<int> &kd,
                                 const std::vector<std::vector<int>> &map_p2s,
                                 const std::vector<int> *interaction_pair_in,
                                 const double *const *const *x_image,
                                 const int *exist,
                                 const int mirror_image_conv,
                                 std::set<InteractionCluster> *interaction_cluster_out) const;

    void cell_combination(const std::vector<std::vector<int>> &,
                          const size_t,
                          const std::vector<int> &,
                          std::vector<std::vector<int>> &) const;

    void generate_pairs(const size_t natmin,
                        const std::vector<std::vector<int>> &map_p2s,
                        std::set<IntList> *pair_out) const;

    void check_permutation_symmetry(const System *system,
                                    const Symmetry *symmetry,
                                    int order);

    void make_symnum_tran_to_prim(const System *system,
                                  const Symmetry *symmetry,
                                  std::vector<int> &symnum_tran_to_prim);

    bool is_inprim(const int iat, // atom index in supercell
                   const size_t natmin,
                   const std::vector<std::vector<int>> &map_p2s) const;
};
}

namespace std {
template<>
struct hash<ALM_NS::IntList> {
    std::size_t operator()(ALM_NS::IntList const &obj) const
    {
        hash<int> hasher;
        size_t seed = 0;
        for (auto i: obj.iarray) {
            seed ^= hasher(i) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }
};
}
