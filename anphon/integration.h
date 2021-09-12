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
#include "memory.h"
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

    class TetraNodes {
     public:
        TetraNodes() {
            nk1 = 0; nk2 = 0; nk3 = 0;
            ntetra = 0;
            tetras = nullptr;
        };
        TetraNodes(unsigned int nk1_in,
                   unsigned int nk2_in,
                   unsigned int nk3_in) {
            nk1 = nk1_in; nk2 = nk2_in; nk3 = nk3_in;
            ntetra = 6 * nk1 * nk2 * nk3;
            allocate(tetras, ntetra, 4);
        };
        ~TetraNodes() {
            if (tetras) deallocate(tetras);
        }

        void setup();
        unsigned int get_ntetra() const;
        unsigned int **get_tetras() const;

     private:
        unsigned int nk1, nk2, nk3;
        unsigned int ntetra;
        unsigned int **tetras;
    };

    class Integration : protected Pointers {
    public:
        Integration(class PHON *);

        ~Integration();

        int ismear; // ismear = -1: tetrahedron, ismear = 0: gaussian
        double epsilon;

        void setup_integration();

        double do_tetrahedron(const double *,
                              const double *,
                              double);

        void calc_weight_tetrahedron(const unsigned int,
                                     const unsigned int *,
                                     double *,
                                     const double *,
                                     double);

    void calc_weight_tetrahedron(const unsigned int nk_irreducible,
                                 const unsigned int *map_to_irreducible_k,
                                 const double *energy,
                                 const double e_ref,
                                 const unsigned int ntetra,
                                 const unsigned int * const * tetras,
                                 double *weight) const;

       void calc_weight_smearing(const unsigned int nk,
                                 const unsigned int nk_irreducible,
                                 const unsigned int *map_to_irreducible_k,
                                 const double *energy,
                                 const double e_ref,
                                 const int smearing_method,
                                 double *weight) const;

    private:
        void set_default_variables();

        void deallocate_variables();

        unsigned int ntetra;
        int **tetras;

        void prepare_tetrahedron(int,
                                 int,
                                 int) const;

        inline double fij(double,
                          double,
                          double) const;

       // inline double volume(const int *) const;

        std::vector<tetra_pair> tetra_data;

       // inline double refold(double) const;

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
