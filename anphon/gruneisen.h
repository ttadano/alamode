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

namespace PHON_NS
{
    //     class FcsClassGru {
    //     public:
    //         std::vector<int> elems;
    //         double fcs_val;
    // 
    //         FcsClassGru();
    //         FcsClassGru(const int n, int *arr, double fcs) {
    //             for (int i = 0; i < n; ++i) {
    //                 elems.push_back(arr[i]);
    //             }
    //             fcs_val = fcs;
    //         }
    //     };
    // 
    //     inline bool operator<(const FcsClassGru &a, const FcsClassGru &b) {
    //         return std::lexicographical_compare(a.elems.begin(), a.elems.end(), b.elems.begin(), b.elems.end());
    //     }

    class FcsAlignedForGruneisen
    {
    public:
        std::vector<AtomCellSuper> pairs;
        double fcs_val;

        FcsAlignedForGruneisen()
        {
        };

        FcsAlignedForGruneisen(const double fcs_in, const std::vector<AtomCellSuper> pairs_in)
        {
            fcs_val = fcs_in;

            for (auto it = pairs_in.cbegin(); it != pairs_in.cend(); ++it) {
                pairs.push_back(*it);
            }
        }
    };

    inline bool operator<(const FcsAlignedForGruneisen &a, const FcsAlignedForGruneisen &b)
    {
        std::vector<unsigned int> array_a, array_b;
        array_a.clear();
        array_b.clear();
        int len = a.pairs.size();
        for (int i = 0; i < len - 1; ++i) {
            array_a.push_back(a.pairs[i].index);
            //            array_a.push_back(a.pairs[i].index/3);
            array_a.push_back(a.pairs[i].tran);
            //            array_a.push_back(a.pairs[i].cell_s);
            //            array_a.push_back(a.pairs[i].index%3);
            array_b.push_back(b.pairs[i].index);
            //            array_b.push_back(b.pairs[i].index/3);
            array_b.push_back(b.pairs[i].tran);
            //            array_b.push_back(b.pairs[i].cell_s);
            //            array_b.push_back(b.pairs[i].index%3);
        }
        for (int i = 0; i < len - 1; ++i) {
            array_a.push_back(a.pairs[i].cell_s);
            array_b.push_back(b.pairs[i].cell_s);
        }

        array_a.push_back(a.pairs[len - 1].index);
        array_a.push_back(a.pairs[len - 1].tran);
        array_b.push_back(b.pairs[len - 1].index);
        array_b.push_back(b.pairs[len - 1].tran);
        return std::lexicographical_compare(array_a.begin(), array_a.end(),
                                            array_b.begin(), array_b.end());
    }

    class Gruneisen: protected Pointers
    {
    public:
        Gruneisen(class PHON *);
        ~Gruneisen();

        double delta_a;
        bool print_gruneisen;
        bool print_newfcs;
        void setup();
        std::complex<double> **gruneisen;
        void calc_gruneisen();
        void finish_gruneisen();
        void write_new_fcsxml_all();

    private:
        double **xshift_s;
        std::vector<FcsArrayWithCell> delta_fc2, delta_fc3;
        void prepare_delta_fcs(const std::vector<FcsArrayWithCell>,
                               std::vector<FcsArrayWithCell> &);
        void calc_dfc2_reciprocal(std::complex<double> **, double *);
        void write_new_fcsxml(const std::string, const double);
        std::string double2string(const double);
        //  double calc_stress_energy2(const std::vector<FcsArrayWithCell>);
        void calc_stress_energy3(const std::vector<FcsArrayWithCell>, double ****);
        void print_stress_energy();
    };
}
