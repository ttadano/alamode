/*
anharmonic_core.h

Copyright (c) 2014 Terumasa Tadano

This file is distributed under the terms of the MIT license.
Please see the file 'LICENCE.txt' in the root directory
or http://opensource.org/licenses/mit-license.php for information.
*/

#pragma once

#include "pointers.h"
#include <complex>
#include <vector>
#include "fcs_phonon.h"

namespace PHON_NS
{
    class KsList
    {
    public:
        std::vector<int> ks;
        int symnum;

        KsList();

        KsList(const KsList &a) : ks(a.ks), symnum(a.symnum) {};

        KsList(const int n, int *ks_in, const int sym)
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

    class KsListGroup
    {
    public:
        std::vector<KsList> group;

        KsListGroup();

        KsListGroup(const std::vector<KsList> &a) : group(a) {};
    };

    class KsListMode
    {
    public:
        double xk[3];
        int nmode;

        KsListMode();

        KsListMode(double xk_in[3], const int n)
        {
            for (int i = 0; i < 3; ++i) xk[i] = xk_in[i];
            nmode = n;
        }
    };

    class KpointListWithCoordinate
    {
    public:
        double xk[3];
        double x, y;
        int plane;
        int selection_type;

        KpointListWithCoordinate();

        KpointListWithCoordinate(const std::vector<double> &a,
                                 const double x_in,
                                 const double y_in,
                                 const int plane_in,
                                 const int selection_type_in)
        {
            for (int i = 0; i < 3; ++i) xk[i] = a[i];
            x = x_in;
            y = y_in;
            plane = plane_in;
            selection_type = selection_type_in;
        }
    };


    class AnharmonicCore : protected Pointers
    {
    public:
        AnharmonicCore(class PHON *);
        ~AnharmonicCore();

        void setup();
        void calc_damping_smearing(const unsigned int, double *, const double,
                                   const unsigned int, const unsigned int,
                                   double *);

        void calc_damping_tetrahedron(const unsigned int, double *, const double,
                                      const unsigned int, const unsigned int,
                                      double *);

        int quartic_mode;
        bool use_tuned_ver;
        bool use_triplet_symmetry;
        bool **is_imaginary;

        std::complex<double> V3(const unsigned int [3]);
        std::complex<double> V4(const unsigned int [4]);


        std::complex<double> V3_mode(int, double *, double *,
                                     int, int, double **,
                                     std::complex<double> ***);

        void prepare_relative_vector(const std::vector<FcsArrayWithCell> &,
                                     const unsigned int, double ***);

        void prepare_group_of_force_constants(const std::vector<FcsArrayWithCell> &,
                                              const unsigned int, int &,
                                              std::vector<double> *&);

        void detect_imaginary_branches(double **);


        void get_unique_triplet_k(const int,
                                  const bool,
                                  const bool,
                                  std::vector<KsListGroup> &);

        void calc_self3omega_tetrahedron(const double, double **,
                                         std::complex<double> ***,
                                         const unsigned int, const unsigned int,
                                         const unsigned int, double *, double *);


    private:
        void set_default_variables();
        void deallocate_variables();
        std::complex<double> im;

        double ***vec_for_v3, *invmass_for_v3;
        double ***vec_for_v4, *invmass_for_v4;
        int **evec_index;
        int **evec_index4;

        bool sym_permutation;

        void setup_cubic();
        void setup_quartic();
        void store_exponential_for_acceleration(const int nk_in[3], int &,
                                                std::complex<double> *,
                                                std::complex<double> ***);
        int ngroup;
        int ngroup2;
        std::vector<double> *fcs_group;
        std::vector<double> *fcs_group2;
        std::complex<double> *exp_phase, ***exp_phase3;

        int nk_grid[3];
        int nk_represent;
        unsigned int tune_type;
        double dnk[3];
    };
}
