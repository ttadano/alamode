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
    class KsListMode
    {
    public:
        double xk[3];
        int nmode;

        KsListMode();

        KsListMode(double xk_in[3],
                   const int n)
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

    class RelativeVector
    {
    public:
        double vecs[3][3];

        RelativeVector();

        // Constructor for cubic term
        RelativeVector(const double vec1[3],
                       const double vec2[3])
        {
            for (int i = 0; i < 3; ++i) {
                vecs[0][i] = vec1[i];
                vecs[1][i] = vec2[i];
                vecs[2][i] = 0.0;
            }
        }

        // Constructor for quartic term
        RelativeVector(const double vec1[3],
                       const double vec2[3],
                       const double vec3[3])
        {
            for (int i = 0; i < 3; ++i) {
                vecs[0][i] = vec1[i];
                vecs[1][i] = vec2[i];
                vecs[2][i] = vec3[i];
            }
        }
    };


    class AnharmonicCore : protected Pointers
    {
    public:
        AnharmonicCore(class PHON *);

        ~AnharmonicCore();

        void setup();

        void calc_damping_smearing(unsigned int,
                                   double *,
                                   double,
                                   unsigned int,
                                   unsigned int,
                                   double *);

        void calc_damping_tetrahedron(unsigned int,
                                      double *,
                                      double,
                                      unsigned int,
                                      unsigned int,
                                      double *);

        int quartic_mode;
        bool use_tuned_ver;
        bool use_triplet_symmetry;

        std::complex<double> V3(const unsigned int [3]);
        std::complex<double> V4(const unsigned int [4]);

        std::complex<double> V3(const unsigned int [3],
                                double **,
                                std::complex<double> ***);

        std::complex<double> V4(const unsigned int [4],
                                double **,
                                std::complex<double> ***);

        std::complex<double> V3_mode(int,
                                     double *,
                                     double *,
                                     int,
                                     int,
                                     double **,
                                     std::complex<double> ***);

        void prepare_relative_vector(const std::vector<FcsArrayWithCell> &,
                                     unsigned int,
                                     double ***);

        void prepare_relative_vector(const std::vector<FcsArrayWithCell> &,
                                     unsigned int,
                                     int,
                                     std::vector<double> *,
                                     std::vector<RelativeVector> *&);

        void prepare_group_of_force_constants(const std::vector<FcsArrayWithCell> &,
                                              unsigned int,
                                              int &,
                                              std::vector<double> *&);


        void calc_self3omega_tetrahedron(double,
                                         double **,
                                         std::complex<double> ***,
                                         unsigned int,
                                         unsigned int,
                                         unsigned int,
                                         double *,
                                         double *);


    private:
        void set_default_variables();
        void deallocate_variables();

        std::complex<double> im;

        double *invmass_v3;
        double *invmass_v4;
        int **evec_index_v3;
        int **evec_index_v4;
        int ngroup_v3;
        int ngroup_v4;
        std::vector<double> *fcs_group_v3;
        std::vector<double> *fcs_group_v4;
        std::complex<double> *exp_phase, ***exp_phase3;
        std::complex<double> *phi3_reciprocal, *phi4_reciprocal;
        std::vector<RelativeVector> *relvec_v3, *relvec_v4;

        int nk_grid[3];
        int nk_represent;
        unsigned int tune_type;
        double dnk[3];

        bool sym_permutation;

        int kindex_phi3_stored[2] = {-1, -1};
        int kindex_phi4_stored[3] = {-1, -1, -1};

        void setup_cubic();
        void setup_quartic();

        void store_exponential_for_acceleration(const int nk_in[3],
                                                int &,
                                                std::complex<double> *,
                                                std::complex<double> ***);

        void calc_phi3_reciprocal(unsigned int,
                                  unsigned int,
                                  std::complex<double> *);

        void calc_phi4_reciprocal(unsigned int,
                                  unsigned int,
                                  unsigned int,
                                  std::complex<double> *);
    };
}
