/*
 phonon_dos.h

 Copyright (c) 2014 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#pragma once

#include "pointers.h"
#include "kpoint.h"
#include <vector>
#include <complex>

namespace PHON_NS
{
    class Dos : protected Pointers
    {
    public:
        Dos(class PHON *);
        ~Dos();

        void setup();
        void calc_dos_all();

        bool flag_dos;
        bool compute_dos;
        bool projected_dos, two_phonon_dos;
        int scattering_phase_space;

        int n_energy;
        double emin, emax, delta_e;
        double *energy_dos;
        double *dos_phonon;
        double **pdos_phonon;
        double ***dos2_phonon;
        double total_sps3, ***sps3_mode;
        double ****sps3_with_bose;


        void calc_dos_scph(double ***, double **);

    private:
        void set_default_variables();
        void deallocate_variables();

        unsigned int nk_irreducible;
        int *kmap_irreducible;
        std::vector<int> k_irreducible;

        void calc_dos(unsigned int, int *,
                      double **, unsigned int, double *,
                      double *, unsigned int, int,
                      std::vector<std::vector<KpointList>> &);

        void calc_atom_projected_dos(unsigned int, double **,
                                     unsigned int, double *, double **,
                                     unsigned int, unsigned int, int,
                                     std::complex<double> ***);

        void calc_two_phonon_dos(unsigned int, double *, double ***, int,
                                 const std::vector<std::vector<KpointList>> &);

        void calc_total_scattering_phase_space(double **, int,
                                               const std::vector<std::vector<KpointList>> &,
                                               double ***, double &);

        void calc_scattering_phase_space_with_Bose(double **, int,
                                                   const std::vector<std::vector<KpointList>> &,
                                                   double ****);

        void calc_scattering_phase_space_with_Bose_mode(unsigned int,
                                                        unsigned int,
                                                        unsigned int,
                                                        double, double **, double *,
                                                        unsigned int *, int, double **);
    };
}
