/*
 interaction.h

 Copyright (c) 2014 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#pragma once

#include <string>
#include <vector>
#include <set>
#include <iterator>
#include <algorithm>
#include "pointers.h"
#include "constants.h"


namespace ALM_NS
{
    class IntList
    {
    public:
        std::vector<int> iarray;

        IntList()
        {
            iarray.clear();
        }

        IntList(const IntList &a)
        {
            for (std::vector<int>::const_iterator p = a.iarray.begin(); p != a.iarray.end(); ++p) {
                iarray.push_back(*p);
            }
        }

        IntList(const int n, const int *arr)
        {
            for (int i = 0; i < n; ++i) {
                iarray.push_back(arr[i]);
            }
        }

        bool operator<(const IntList &a) const
        {
            return std::lexicographical_compare(iarray.begin(), iarray.end(),
                                                a.iarray.begin(), a.iarray.end());
        }
    };

    class InteractionCluster
    {
    public:
        std::vector<double> x;

        InteractionCluster();

        InteractionCluster(const double *arr)
        {
            for (int i = 0; i < 3; ++i) {
                x.push_back(arr[i]);
            }
        }

        bool operator<(const InteractionCluster &a) const
        {
            return std::lexicographical_compare(x.begin(), x.end(),
                                                a.x.begin(), a.x.end());
        }
    };


    class DistInfo
    {
    public:
        int cell;
        double dist;
        double relvec[3];

        DistInfo();

        DistInfo(const int n, const double d, const double x[3])
        {
            cell = n;
            dist = d;
            for (int i = 0; i < 3; ++i) relvec[i] = x[i];
        }

        DistInfo(const DistInfo &obj)
        {
            cell = obj.cell;
            dist = obj.dist;
            for (int i = 0; i < 3; ++i) relvec[i] = obj.relvec[i];
        }

        bool operator<(const DistInfo &a) const
        {
            return dist < a.dist;
        }
    };

    class DistList
    {
    public:
        int atom;
        double dist;

        DistList();

        DistList(int n, double dist_tmp)
        {
            atom = n;
            dist = dist_tmp;
        }

        bool operator<(const DistList &a) const
        {
            if (std::abs(dist - a.dist) > eps8) {
                return dist < a.dist;
            } else {
                return atom < a.atom;
            }
        }
    };

    class MinDistList
    {
    public:
        std::vector<int> cell;
        std::vector<double> dist;

        MinDistList();

        MinDistList(const std::vector<int> cell_in,
                    const std::vector<double> dist_in)
        {
            for (std::vector<int>::const_iterator it = cell_in.begin(); it != cell_in.end(); ++it) {
                cell.push_back(*it);
            }
            for (std::vector<double>::const_iterator it = dist_in.begin(); it != dist_in.end(); ++it) {
                dist.push_back(*it);
            }
        }

        static bool compare_sum_distance(const MinDistList &a, const MinDistList &b)
        {
            double dist_a = 0;
            double dist_b = 0;

            for (int i = 0; i < a.dist.size(); ++i) {
                dist_a += a.dist[i];
            }
            for (int i = 0; i < b.dist.size(); ++i) {
                dist_b += b.dist[i];
            }
            return dist_a < dist_b;
        }

        static bool compare_max_distance(const MinDistList &a, const MinDistList &b)
        {
            // This function works properly when dvec_a.size() > 0 and dvec_b.size() > 0
            std::vector<double> dvec_a, dvec_b;
            std::copy(a.dist.begin(), a.dist.end(), std::back_inserter(dvec_a));
            std::copy(b.dist.begin(), b.dist.end(), std::back_inserter(dvec_b));
            double max_dist_a = *std::max_element(dvec_a.begin(), dvec_a.end());
            double max_dist_b = *std::max_element(dvec_b.begin(), dvec_b.end());

            return max_dist_a < max_dist_b;
        }
    };

    class MinimumDistanceCluster
    {
    public:
        std::vector<int> atom;
        std::vector<std::vector<int>> cell;
        double distmax;

        MinimumDistanceCluster();

        MinimumDistanceCluster(const std::vector<int> atom_in,
                               const std::vector<std::vector<int>> cell_in,
                               const double dist_in)
        {
            for (int i = 0; i < atom_in.size(); ++i) {
                atom.push_back(atom_in[i]);
            }
            for (int i = 0; i < cell_in.size(); ++i) {
                cell.push_back(cell_in[i]);
            }
            distmax = dist_in;
        }

        MinimumDistanceCluster(const std::vector<int> atom_in,
                               const std::vector<std::vector<int>> cell_in)
        {
            for (int i = 0; i < atom_in.size(); ++i) {
                atom.push_back(atom_in[i]);
            }
            for (int i = 0; i < cell_in.size(); ++i) {
                cell.push_back(cell_in[i]);
            }
        }

        bool operator<(const MinimumDistanceCluster &a) const
        {
            return lexicographical_compare(atom.begin(), atom.end(),
                                           a.atom.begin(), a.atom.end());
        }
    };

    class Interaction: protected Pointers
    {
    public:
        Interaction(class ALM *);
        ~Interaction();

        int is_periodic[3];
        int nneib;
        int maxorder;

        int *nbody_include;

        double ***rcs;
        double ***x_image;
        int *exist_image;

        std::string *str_order;
        std::vector<DistInfo> **distall;
        std::vector<DistInfo> **mindist_pairs;
        std::set<IntList> *pairs;
        std::vector<int> **interaction_pair;
        std::set<MinimumDistanceCluster> **mindist_cluster;

        void init();
        double distance(double *, double *);
        int nbody(const int, const int *);
        bool is_incutoff(const int, int *, const int);
        bool is_incutoff2(const int, int *, const int);

        template <typename T>
        void insort(int n, T *arr)
        {
            int i, j;
            T tmp;

            for (i = 1; i < n; ++i) {
                tmp = arr[i];
                for (j = i - 1; j >= 0 && arr[j] > tmp; --j) {
                    arr[j + 1] = arr[j];
                }
                arr[j + 1] = tmp;
            }
        }

    private:
        void generate_coordinate_of_periodic_images(const unsigned int, double **,
                                                    const int [3], double ***, int *);

        void get_pairs_of_minimum_distance(int, double ***, int *,
                                           std::vector<DistInfo> **,
                                           std::vector<DistInfo> **);

        void print_neighborlist(std::vector<DistInfo> **);
        void search_interactions(std::vector<int> **, std::set<IntList> *);
        void set_ordername();

        void calc_mindist_clusters(std::vector<int> **,
                                   std::vector<DistInfo> **,
                                   std::vector<DistInfo> **,
                                   int *, std::set<MinimumDistanceCluster> **);

        void calc_mindist_clusters2(std::vector<int> **,
                                    std::vector<DistInfo> **,
                                    std::vector<DistInfo> **,
                                    int *, std::set<MinimumDistanceCluster> **);

        void cell_combination(std::vector<std::vector<int>>,
                              int, std::vector<int>,
                              std::vector<std::vector<int>> &);

        void generate_pairs(std::set<IntList> *, std::set<MinimumDistanceCluster> **);
    };
}
