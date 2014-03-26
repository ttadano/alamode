/*
 gruneisen.h

 Copyright (c) 2014 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#pragma once

#include "pointers.h"
#include <string>
#include <complex>
#include <vector>
#include "fcs_phonon.h"

namespace PHON_NS {

    class FcsClassGru {
    public:
        std::vector<int> elems;
        double fcs_val;

        FcsClassGru();
        FcsClassGru(const int n, int *arr, double fcs) {
            for (int i = 0; i < n; ++i) {
                elems.push_back(arr[i]);
            }
            fcs_val = fcs;
        }

    };

    inline bool operator<(const FcsClassGru &a, const FcsClassGru &b) {
        //	return std::lexicographical_compare(a.elems.begin(), a.elems.end() - 1, b.elems.begin(), b.elems.end() - 1);
        return std::lexicographical_compare(a.elems.begin(), a.elems.end(), b.elems.begin(), b.elems.end());
    }

    class Gruneisen: protected Pointers {
    public:
        Gruneisen(class PHON *);
        ~Gruneisen();

        double delta_a;
        bool print_gruneisen;
        bool print_newfcs;
        void setup();
        std::complex<double> **gruneisen;
        void calc_gruneisen();
        void calc_gruneisen2();
        void calc_gruneisen3();
        void finish_gruneisen();
        void write_newinfo_all();


    private:
        double ****fc2_plus, ****fc2_minus;
        double ****dfc2;
        std::vector<FcsClassExtent> fc2_plus_ext, fc2_minus_ext;
        std::vector<FcsClassGru> fc3_plus, fc3_minus;
        void prepare_delta_fc2();
        void prepare_newfc2();
        void prepare_newfc3();
        void write_newinfo(std::ifstream &, std::ofstream &, const double, double ****, std::vector<FcsClassExtent>, std::vector<FcsClassGru>);
        void calc_dfc2_reciprocal(std::complex<double> **, double *);
        void calc_pressure();

    };
}
