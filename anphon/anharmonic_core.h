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


    class AnharmonicCore : protected Pointers
    {
    public:
        AnharmonicCore(class PHON *);

        ~AnharmonicCore();

        void setup();

        void calc_damping_smearing(const unsigned int,
                                   double *,
                                   const double,
                                   const unsigned int,
                                   const unsigned int,
                                   double *);

        void calc_damping_tetrahedron(const unsigned int,
                                      double *,
                                      const double,
                                      const unsigned int,
                                      const unsigned int,
                                      double *);

        void calc_damping_tetrahedron2(const unsigned int,
                                       double *,
                                       const double,
                                       const unsigned int,
                                       const unsigned int,
                                       double *);

        int quartic_mode;
        bool use_tuned_ver;
        bool use_triplet_symmetry;

        std::complex<double> V3(const unsigned int [3]);

        std::complex<double> V3(const unsigned int [3],
                                std::complex<double> *);


        std::complex<double> V4(const unsigned int [4]);


        std::complex<double> V3_mode(int,
                                     double *,
                                     double *,
                                     int,
                                     int,
                                     double **,
                                     std::complex<double> ***);

        void prepare_relative_vector(const std::vector<FcsArrayWithCell> &,
                                     const unsigned int,
                                     double ***);

        void prepare_group_of_force_constants(const std::vector<FcsArrayWithCell> &,
                                              const unsigned int,
                                              int &,
                                              std::vector<double> *&);


        void calc_self3omega_tetrahedron(const double,
                                         double **,
                                         std::complex<double> ***,
                                         const unsigned int,
                                         const unsigned int,
                                         const unsigned int,
                                         double *,
                                         double *);


    private:
        void set_default_variables();
        void deallocate_variables();

        std::complex<double> im;

        double ***relvec_v3, *invmass_v3;
        double ***relvec_v4, *invmass_v4;
        int **evec_index_v3;
        int **evec_index_v4;

        bool sym_permutation;

        void setup_cubic();
        void setup_quartic();

        void store_exponential_for_acceleration(const int nk_in[3],
                                                int &,
                                                std::complex<double> *,
                                                std::complex<double> ***);

        void phase_V3(const unsigned int,
                      const unsigned int,
                      std::complex<double> *);

        int ngroup_v3;
        int ngroup_v4;
        std::vector<double> *fcs_group_v3;
        std::vector<double> *fcs_group_v4;
        std::complex<double> *exp_phase, ***exp_phase3;

        int nk_grid[3];
        int nk_represent;
        unsigned int tune_type;
        double dnk[3];
    };
}
