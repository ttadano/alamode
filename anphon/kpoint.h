/*
kpoint.h

Copyright (c) 2014 Terumasa Tadano

This file is distributed under the terms of the MIT license.
Please see the file 'LICENCE.txt' in the root directory
or http://opensource.org/licenses/mit-license.php for information.
*/

#pragma once

#include "pointers.h"
#include <iostream>
#include <string>
#include <vector>
#include <set>
#include <map>
#include "constants.h"

namespace PHON_NS
{
    class KpointList
    {
    public:
        std::vector<double> kval;
        unsigned int knum;

        KpointList()
        {
        };

        KpointList(const KpointList &obj)
        {
            knum = obj.knum;
            for (auto it = obj.kval.cbegin(); it != obj.kval.cend(); ++it) {
                kval.push_back(*it);
            }
        }

        KpointList(const unsigned int knum_in, const std::vector<double> vec)
        {
            knum = knum_in;
            for (auto it = vec.cbegin(); it != vec.cend(); ++it) {
                kval.push_back(*it);
            }
        }
    };

    class KpointInp
    {
    public:
        std::vector<std::string> kpelem;

        KpointInp()
        {
        };

        KpointInp(const std::vector<std::string> &obj)
        {
            for (auto it = obj.cbegin(); it != obj.cend(); ++it) {
                kpelem.push_back(*it);
            }
        }
    };

    class KpointPlaneGeometry
    {
    public:
        double xk_origin[3];
        double xk_edges[2][3];
        int npoints[2];

        KpointPlaneGeometry()
        {
        };

        KpointPlaneGeometry(double *xk_origin_in,
                            double *xk_edge1_in,
                            double *xk_edge2_in,
                            int *num)
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

    class KpointPlane
    {
    public:
        double k[3];
        int n[2];

        KpointPlane()
        {
        };

        KpointPlane(double *xk_in, int *n_in)
        {
            for (int i = 0; i < 3; ++i) k[i] = xk_in[i];
            for (int i = 0; i < 2; ++i) n[i] = n_in[i];
        }
    };

    class KpointPlaneTriangle
    {
    public:
        int index;
        int knum[3];

        KpointPlaneTriangle()
        {
        };

        KpointPlaneTriangle(int index_in, int *nk_in)
        {
            index = index_in;

            for (int i = 0; i < 3; ++i) {
                knum[i] = nk_in[i];
            }
        }
    };

    class Kpoint : protected Pointers
    {
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
        std::vector<double> weight_k;
        std::vector<std::vector<KpointList>> kpoint_irred_all;

        unsigned int nplanes;
        std::vector<KpointPlane> *kp_planes;
        std::vector<KpointPlaneGeometry> kp_plane_geometry;
        std::vector<KpointPlaneTriangle> *kp_planes_tri;
        unsigned int nk_reduced;
        std::map<int, int> kmap_to_irreducible;
        std::vector<int> *small_group_of_k;


        int get_knum(const double, const double, const double);
        int get_knum(const double [3], const unsigned int [3]);

        void generate_irreducible_kmap(int *, unsigned int &,
                                       std::vector<int> &,
                                       const unsigned int, const unsigned int, const unsigned int,
                                       double **, const int, int ***);

        void gen_kmesh(const bool, const unsigned int [3],
                       double **,
                       std::vector<std::vector<KpointList>> &);

        void get_small_group_k(double *, std::vector<int> &, double [3][3]);
        int knum_sym(const int, const int);


    private:
        void setup_kpoint_given(std::vector<KpointInp> &,
                                unsigned int &,
                                double **&, double **&);

        void setup_kpoint_band(std::vector<KpointInp> &,
                               unsigned int &,
                               double **&, double **&, double *&);

        void setup_kpoint_mesh(std::vector<KpointInp> &,
                               unsigned int &, unsigned int &, unsigned int &, unsigned int &,
                               double **&, double **&,
                               const bool,
                               std::vector<std::vector<KpointList>> &);

        void setup_kpoint_plane(std::vector<KpointInp> &,
                                unsigned int &,
                                std::vector<KpointPlane> *&);

        void reduce_kpoints(const unsigned int, double **,
                            const unsigned int [3],
                            std::vector<std::vector<KpointList>> &);

        void gen_nkminus(const unsigned int, unsigned int *, double **);

        void gen_kpoints_plane(std::vector<KpointInp>,
                               std::vector<KpointPlane> *,
                               std::vector<KpointPlaneTriangle> *);

        bool in_first_BZ(double *);

        void mpi_broadcast_kpoint_vector(std::vector<std::vector<KpointList>> &);
        void mpi_broadcast_kplane_vector(const unsigned int, std::vector<KpointPlane> *&);
        void calc_small_groups_k_irred(std::vector<int> *);
        std::vector<int> get_small_group_of_k(const int);
    };
}
