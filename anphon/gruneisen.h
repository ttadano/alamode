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
    class FcsAlignedForGruneisen
    {
    public:
        std::vector<AtomCellSuper> pairs;
        double fcs_val;

        FcsAlignedForGruneisen() { };

        FcsAlignedForGruneisen(const double fcs_in,
                               const std::vector<AtomCellSuper> &pairs_in)
            : pairs(pairs_in), fcs_val(fcs_in) {};

        bool operator<(const FcsAlignedForGruneisen &obj) const
        {
            std::vector<unsigned int> array_a, array_b;
            array_a.clear();
            array_b.clear();
            int len = pairs.size();
            for (int i = 0; i < len - 1; ++i) {
                array_a.push_back(pairs[i].index);
                array_a.push_back(pairs[i].tran);
                array_b.push_back(obj.pairs[i].index);
                array_b.push_back(obj.pairs[i].tran);
            }
            for (int i = 0; i < len - 1; ++i) {
                array_a.push_back(pairs[i].cell_s);
                array_b.push_back(obj.pairs[i].cell_s);
            }

            array_a.push_back(pairs[len - 1].index);
            array_a.push_back(pairs[len - 1].tran);
            array_b.push_back(obj.pairs[len - 1].index);
            array_b.push_back(obj.pairs[len - 1].tran);
            return std::lexicographical_compare(array_a.begin(), array_a.end(),
                                                array_b.begin(), array_b.end());
        }
    };

    //inline bool operator<(const FcsAlignedForGruneisen &a, const FcsAlignedForGruneisen &b)
    //{
    //    std::vector<unsigned int> array_a, array_b;
    //    array_a.clear();
    //    array_b.clear();
    //    int len = a.pairs.size();
    //    for (int i = 0; i < len - 1; ++i) {
    //        array_a.push_back(a.pairs[i].index);
    //        //            array_a.push_back(a.pairs[i].index/3);
    //        array_a.push_back(a.pairs[i].tran);
    //        //            array_a.push_back(a.pairs[i].cell_s);
    //        //            array_a.push_back(a.pairs[i].index%3);
    //        array_b.push_back(b.pairs[i].index);
    //        //            array_b.push_back(b.pairs[i].index/3);
    //        array_b.push_back(b.pairs[i].tran);
    //        //            array_b.push_back(b.pairs[i].cell_s);
    //        //            array_b.push_back(b.pairs[i].index%3);
    //    }
    //    for (int i = 0; i < len - 1; ++i) {
    //        array_a.push_back(a.pairs[i].cell_s);
    //        array_b.push_back(b.pairs[i].cell_s);
    //    }

    //    array_a.push_back(a.pairs[len - 1].index);
    //    array_a.push_back(a.pairs[len - 1].tran);
    //    array_b.push_back(b.pairs[len - 1].index);
    //    array_b.push_back(b.pairs[len - 1].tran);
    //    return std::lexicographical_compare(array_a.begin(), array_a.end(),
    //                                        array_b.begin(), array_b.end());
    //}

    class Gruneisen : protected Pointers
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
        void write_new_fcsxml_all();

    private:
        void set_default_variables();
        void deallocate_variables();

        double **xshift_s;
        std::vector<FcsArrayWithCell> delta_fc2, delta_fc3;

        void prepare_delta_fcs(const std::vector<FcsArrayWithCell> &,
                               std::vector<FcsArrayWithCell> &);

        void calc_dfc2_reciprocal(std::complex<double> **,
                                  const double *);

        void write_new_fcsxml(std::string,
                              double);

        std::string double2string(double);

        //  double calc_stress_energy2(const std::vector<FcsArrayWithCell>);
        void calc_stress_energy3(std::vector<FcsArrayWithCell>,
                                 double ****);

        void print_stress_energy();
    };
}
