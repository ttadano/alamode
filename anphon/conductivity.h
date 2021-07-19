/*
 conductivity.h

 Copyright (c) 2014, 2015, 2016 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#pragma once

#include "pointers.h"
#include <vector>
#include <set>
#include <complex>

namespace PHON_NS {
    class Conductivity : protected Pointers {
    public:
        Conductivity(class PHON *);

        ~Conductivity();

        void setup_kappa();

        void prepare_restart();

        void calc_anharmonic_imagself();

        void compute_kappa();

        int calc_kappa_spec;
        unsigned int ntemp;
        double **damping3;
        double **damping4;
        double ***kappa;
        double ***kappa_spec;
        double ***kappa_coherent;
        double *Temperature;
        int calc_coherent;
        int fph_rta;

    private:
        void set_default_variables();

        void deallocate_variables();

        double ***vel;
        std::complex<double> ****velmat;
        unsigned int nk, ns;
        int nshift_restart, nshift_restart4;
        std::vector<int> vks_l, vks_done, vks_done4;
        std::set<int> vks_job, vks_job4;
        std::string file_coherent_elems;

        void calc_anharmonic_imagself3();
        void calc_anharmonic_imagself4();
        void calc_anharmonic_combined();

        void write_result_gamma(unsigned int,
                                unsigned int,
                                double ***,
                                double **) const;

        
        void write_result_gamma(unsigned int,
                                unsigned int,
                                double ***,
                                double **,
                                int) const;

        void write_result_gamma34(unsigned int,
                                unsigned int,
                                double ***,
                                double **,
                                double **) const;

        void average_self_energy_at_degenerate_point(int,
                                                     int,
                                                     double **) const;

        void compute_frequency_resolved_kappa(int,
                                              double ****,
                                              int);

        void compute_kappa_intraband(double ***kappa_intra, double **lifetime);

        void compute_kappa_coherent(double ***kappa_coherent, double **gamma_total) const;
    };
}
