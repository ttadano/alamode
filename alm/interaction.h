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
#include "listcomparison.h"
#include "pointers.h"
#include "constants.h"


namespace ALM_NS {
    class InteractionCluster {
    public:
        std::vector<double> x;

        InteractionCluster();
        InteractionCluster(const double *arr){
            for (int i = 0; i < 3; ++i){
                x.push_back(arr[i]);
            }
        }
    };

    inline bool operator<(const InteractionCluster a, const InteractionCluster b){
        return std::lexicographical_compare(a.x.begin(), a.x.end(), b.x.begin(), b.x.end());
    }

    class DistInfo {
    public:
        int cell;
        double dist;
        double relvec[3];

        DistInfo();
        DistInfo(const int n, const double d, const double x[3]) {
            cell = n;
            dist = d;
            for (int i = 0; i < 3; ++i) relvec[i] = x[i];
        }

        DistInfo(const DistInfo &obj) {
            cell = obj.cell;
            dist = obj.dist;
            for (int i = 0; i < 3; ++i) relvec[i] = obj.relvec[i];
        }
    };

    inline bool operator<(const DistInfo a, const DistInfo b) {
        return a.dist < b.dist;
    }

    class DistList {
    public:
        int atom;
        double dist;

        DistList();
        DistList(int n, double dist_tmp) {
            atom = n;
            dist = dist_tmp;
        }
    };
    inline bool operator<(const DistList a, const DistList b) {
        if (std::abs(a.dist - b.dist) > eps8) {
            return a.dist < b.dist;
        } else {
            return a.atom < b.atom;
        }
    }

    class MinDistList {
    public:
        std::vector<int> cell;
        std::vector<double> dist;

        MinDistList();
        MinDistList(const std::vector<int> cell_in, const std::vector<double> dist_in) {
            for (std::vector<int>::const_iterator it = cell_in.begin(); it != cell_in.end(); ++it) {
                cell.push_back(*it);
            }
            for (std::vector<double>::const_iterator it = dist_in.begin(); it != dist_in.end(); ++it) {
                dist.push_back(*it);
            }
        }

        static bool compare_sum_distance(const MinDistList &a, const MinDistList &b) {
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
    };

    class MinimumDistanceCluster {
    public:
        std::vector<int> atom;
        std::vector<std::vector<int> > cell;

        MinimumDistanceCluster();

        MinimumDistanceCluster(const std::vector<int> atom_in, const std::vector<std::vector<int> > cell_in) {
            for (int i = 0; i < atom_in.size(); ++i) {
                atom.push_back(atom_in[i]);
            }
            for (int i = 0; i < cell_in.size(); ++i) {
                cell.push_back(cell_in[i]);
            }
        }

        friend bool operator<(const MinimumDistanceCluster &a, const MinimumDistanceCluster &b) {
            return lexicographical_compare(a.atom.begin(), a.atom.end(), b.atom.begin(), b.atom.end());
        }
    };

    class Interaction: protected Pointers {
    public:
        Interaction(class ALM *);
        ~Interaction();

        int is_periodic[3];
        int nneib;
        int maxorder;

        int *nbody_include;

        double ***rcs;
        double ***xcrd;

        std::string *str_order;

        std::vector<DistInfo> **mindist_pairs;
        std::set<IntList> *pairs;
        std::vector<int> **interaction_pair;
        std::set<MinimumDistanceCluster> **mindist_cluster;

        void init();
        double distance(double *, double *);
        bool is_incutoff(const int, int *);

        template <typename T>
        void insort(int n, T *arr)
        {
            int i, j;
            T tmp;

            for (i = 1; i < n; ++i){
                tmp = arr[i];
                for (j = i - 1; j >= 0 && arr[j] > tmp; --j){
                    arr[j + 1] = arr[j];
                }
                arr[j + 1] = tmp;
            }
        }

    private:
        int nsize[3];
        void get_pairs_of_minimum_distance(int, double **, std::vector<DistInfo> **);
        void print_neighborlist(std::vector<DistInfo> **);
        void search_interactions(std::vector<int> **, std::set<IntList> *);
        void set_ordername();
        void calc_mindist_clusters(std::vector<int> **, std::vector<DistInfo> **, std::set<MinimumDistanceCluster> **);
        int nbody(const int, const int *);
        void cell_combination(std::vector<std::vector<int> >, int, std::vector<int>, std::vector<std::vector<int> > &);

    };
}
