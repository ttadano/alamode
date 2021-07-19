/*
 integration.h

 Copyright (c) 2014 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#pragma once

#include "pointers.h"
#include "constants.h"
#include "kpoint.h"
#include <vector>

namespace PHON_NS {
    struct tetra_pair {
        double e;
        double f;
    };

    inline bool operator<(const tetra_pair &a,
                          const tetra_pair &b)
    {
        return a.e < b.e;
    }

    struct TetraWithKnum {
        double e;
        int knum;
    };

    inline bool operator<(const TetraWithKnum &a,
                          const TetraWithKnum &b)
    {
        return a.e < b.e;
    }

    class Integration : protected Pointers {
    public:
        Integration(class PHON *);

        ~Integration();

        bool use_tetrahedron;
        int ismear;
        int ismear_4ph;
        double epsilon;

        void setup_integration();

        double do_tetrahedron(const double *,
                              const double *,
                              double);

        double dos_integration(double *,
                               double);

        void calc_weight_tetrahedron(int,
                                     const int *,
                                     double *,
                                     const double *,
                                     double);

        void calc_weight_smearing(const std::vector<std::vector<KpointList>> &,
                                  double *,
                                  double *,
                                  double,
                                  int) const;

        void calc_weight_smearing(int,
                                  int,
                                  const int *,
                                  double *,
                                  double *,
                                  double,
                                  int) const;

        // overload for 3ph or 4ph
        void adaptive_smearing( int,int,int,int,
                                double *);

        void adaptive_smearing( int,int,int,int,
                                int,int,double *);

    private:
        void set_default_variables();

        void deallocate_variables();

        unsigned int ntetra;
        int **tetras;
        // for adaptive smearing
        double ***vel;
        double **dq;

        void prepare_tetrahedron(int,
                                 int,
                                 int) const;

        void prepare_adaptivesmearing();

        inline double fij(double,
                          double,
                          double) const;

        inline double volume(const int *) const;

        std::vector<tetra_pair> tetra_data;

        inline double refold(double) const;

        void insertion_sort(double *,
                            int *,
                            int) const;
    };

    inline double delta_lorentz(const double omega,
                                const double epsilon)
    {
        return inverse_pi * epsilon / (omega * omega + epsilon * epsilon);
    }

    inline double delta_gauss(const double omega,
                              const double epsilon)
    {
        return std::exp(-omega * omega / (epsilon * epsilon)) / (epsilon * std::sqrt(pi));
    }
}
