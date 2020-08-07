/*
kpoint.h

Copyright (c) 2014 Terumasa Tadano

This file is distributed under the terms of the MIT license.
Please see the file 'LICENCE.txt' in the root directory
or http://opensource.org/licenses/mit-license.php for information.
*/

#pragma once

#include "pointers.h"
#include <string>
#include <vector>
#include <map>

namespace PHON_NS {
    class KpointList {
    public:
        std::vector<double> kval;
        unsigned int knum;

        KpointList() {};

        KpointList(const KpointList &obj) : kval(obj.kval), knum(obj.knum) {};

        KpointList(const unsigned int knum_in,
                   const std::vector<double> &vec)
                : kval(vec), knum(knum_in) {};
    };

    class KpointInp {
    public:
        std::vector<std::string> kpelem;

        KpointInp() {};

        KpointInp(const std::vector<std::string> &obj) : kpelem(obj) {};
    };

    class KpointPlaneGeometry {
    public:
        double xk_origin[3];
        double xk_edges[2][3];
        int npoints[2];

        KpointPlaneGeometry() {};

        KpointPlaneGeometry(const double *xk_origin_in,
                            const double *xk_edge1_in,
                            const double *xk_edge2_in,
                            const int *num)
        {
            for (int i = 0; i < 3; ++i) {
                xk_origin[i] = xk_origin_in[i];
                xk_edges[0][i] = xk_edge1_in[i];
                xk_edges[1][i] = xk_edge2_in[i];
            }
            for (int i = 0; i < 2; ++i) {
                npoints[i] = num[i];
            }
        }
    };

    class KpointPlane {
    public:
        double k[3];
        int n[2];

        KpointPlane() {};

        KpointPlane(const double *xk_in,
                    const int *n_in)
        {
            for (int i = 0; i < 3; ++i) k[i] = xk_in[i];
            for (int i = 0; i < 2; ++i) n[i] = n_in[i];
        }
    };

    class KpointPlaneTriangle {
    public:
        int index;
        int knum[3];

        KpointPlaneTriangle() {};

        KpointPlaneTriangle(int index_in,
                            const int *nk_in)
        {
            index = index_in;

            for (int i = 0; i < 3; ++i) {
                knum[i] = nk_in[i];
            }
        }
    };

    class KsList {
    public:
        std::vector<int> ks;
        int symnum;

        KsList();

        KsList(const KsList &a) : ks(a.ks), symnum(a.symnum) {};

        KsList(const int n,
               int *ks_in,
               const int sym)
        {
            for (int i = 0; i < n; ++i) {
                ks.push_back(ks_in[i]);
            }
            symnum = sym;
        }

        bool operator<(const KsList &obj) const
        {
            return std::lexicographical_compare(ks.begin(), ks.end(),
                                                obj.ks.begin(), obj.ks.end());
        }
    };

    class KsListGroup {
    public:
        std::vector<KsList> group;

        KsListGroup();

        KsListGroup(const std::vector<KsList> &a) : group(a) {};
    };

    class Kpoint : protected Pointers {
    public:
        Kpoint(class PHON *);

        ~Kpoint();

        void kpoint_setups(std::string);

        int kpoint_mode;
        unsigned int nkx, nky, nkz;
        unsigned int nk;
        unsigned int *knum_minus;

        double **xk;
        double *kaxis;
        double **kvec_na;

        std::vector<KpointInp> kpInp;
        std::vector<double> weight_k; // Weight of each kpoint used for integration
        std::vector<std::vector<KpointList>> kpoint_irred_all;

        unsigned int nplanes;
        std::vector<KpointPlane> *kp_planes;
        std::vector<KpointPlaneGeometry> kp_plane_geometry;
        std::vector<KpointPlaneTriangle> *kp_planes_tri;
        unsigned int nk_irred; // Number of irreducible k points
        std::map<int, int> kmap_to_irreducible;
        std::vector<int> *small_group_of_k;


        int get_knum(double,
                     double,
                     double) const;

        int get_knum(const double [3],
                     const unsigned int [3]) const;

        void generate_irreducible_kmap(int *,
                                       unsigned int &,
                                       std::vector<int> &,
                                       unsigned int,
                                       unsigned int,
                                       unsigned int,
                                       double **,
                                       int,
                                       int ***) const;

        void gen_kmesh(bool,
                       const unsigned int [3],
                       double **,
                       std::vector<std::vector<KpointList>> &) const;

        void get_small_group_k(const double *,
                               std::vector<int> &,
                               double [3][3]) const;

        int knum_sym(int,
                     int) const;

        void get_commensurate_kpoints(const double [3][3],
                                      const double [3][3],
                                      std::vector<std::vector<double>> &) const;

        void get_unique_triplet_k(const int,
                                  const bool,
                                  const bool,
                                  std::vector<KsListGroup> &,
                                  const int sign = -1);

    private:
        void set_default_variables();

        void deallocate_variables();

        void setup_kpoint_given(const std::vector<KpointInp> &,
                                unsigned int &,
                                double **&,
                                double **&) const;

        void setup_kpoint_band(const std::vector<KpointInp> &,
                               unsigned int &,
                               double **&,
                               double **&,
                               double *&) const;

        void setup_kpoint_mesh(const std::vector<KpointInp> &,
                               unsigned int &,
                               unsigned int &,
                               unsigned int &,
                               unsigned int &,
                               double **&,
                               double **&,
                               bool,
                               std::vector<std::vector<KpointList>> &) const;

        void setup_kpoint_plane(const std::vector<KpointInp> &,
                                unsigned int &,
                                std::vector<KpointPlane> *&);

        void reduce_kpoints(unsigned int,
                            double **,
                            const unsigned int [3],
                            std::vector<std::vector<KpointList>> &) const;

        void gen_nkminus(unsigned int,
                         unsigned int *,
                         double **) const;

        void gen_kpoints_plane(const std::vector<KpointInp> &,
                               std::vector<KpointPlane> *,
                               std::vector<KpointPlaneTriangle> *);

        bool in_first_BZ(const double *) const;

        void mpi_broadcast_kpoint_vector(std::vector<std::vector<KpointList>> &) const;

        void mpi_broadcast_kplane_vector(unsigned int,
                                         std::vector<KpointPlane> *&) const;

        void calc_small_groups_k_irred(std::vector<int> *);

        std::vector<int> get_small_group_of_k(int) const;
    };
}
