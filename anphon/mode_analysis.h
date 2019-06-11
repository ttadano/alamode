/*
 mode_analysis.h

Copyright (c) 2018 Terumasa Tadano

This file is distributed under the terms of the MIT license.
Please see the file 'LICENCE.txt' in the root directory
or http://opensource.org/licenses/mit-license.php for information.
*/

#pragma once

#include "pointers.h"
#include <complex>
#include <vector>
#include <string>
#include "anharmonic_core.h"
#include "kpoint.h"

namespace PHON_NS
{
    class ModeAnalysis : protected Pointers
    {
    public:
        ModeAnalysis(class PHON *);

        ~ModeAnalysis();

        void run_mode_analysis();

        void setup_mode_analysis();

        bool ks_analyze_mode;
        bool calc_imagpart;
        bool calc_realpart;
        bool calc_fstate_omega;
        bool calc_fstate_k;
        int print_V3;
        bool spectral_func;

        std::string ks_input;
        std::vector<unsigned int> kslist;

    private:
        void set_default_variables();
        void deallocate_variables() const;

        std::vector<KsListMode> kslist_fstate_k;

        void calc_frequency_resolved_final_state(const unsigned int,
                                                 double *,
                                                 const double,
                                                 const unsigned int,
                                                 const double *,
                                                 const unsigned int,
                                                 const unsigned int,
                                                 double ***) const;

        void calc_frequency_resolved_final_state_tetrahedron(const unsigned int,
                                                             double *,
                                                             const double,
                                                             const unsigned int,
                                                             const double *,
                                                             const unsigned int,
                                                             const unsigned int,
                                                             double ***) const;

        void print_momentum_resolved_final_state(const unsigned int,
                                                 double *,
                                                 const double);

        void print_frequency_resolved_final_state(const unsigned int,
                                                  double *);

        void print_V3_elements();

        void print_Phi3_elements();

        void calc_V3norm2(const unsigned int,
                          const unsigned int,
                          double **) const;

        void calc_Phi3(unsigned int,
                       unsigned int,
                       const std::vector<KsListGroup> &,
                       std::complex<double> **) const;

        void print_selfenergy(const int,
                              double *);

        void print_spectral_function(const int,
                                     double *);
    };
}
