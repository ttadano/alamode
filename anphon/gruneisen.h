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
class FcsAlignedForGruneisen : public FcsArrayWithCell {
public:
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

struct sort_by_heading_indices {
    inline bool operator()(const FcsArrayWithCell &a, FcsArrayWithCell &b)
    {
        std::vector<int> array_a, array_b;
        array_a.clear();
        array_b.clear();
        int len = a.pairs.size();
        for (int i = 0; i < len - 1; ++i) {
            array_a.push_back(a.pairs[i].index);
            array_b.push_back(b.pairs[i].index);
        }
        // The components of relvec should be integers,
        // so let's convert their types into int for sort.
        for (int i = 0; i < len - 2; ++i) {
            for (auto j = 0; j < 3; ++j) {
                array_a.push_back(nint(a.relvecs[i][j]));
                array_b.push_back(nint(b.relvecs[i][j]));
            }
        }
        // Register the last index
        array_a.push_back(a.pairs[len - 1].index);
        array_b.push_back(b.pairs[len - 1].index);
        for (auto j = 0; j < 3; ++j) {
            array_a.push_back(nint(a.relvecs[len - 2][j]));
            array_b.push_back(nint(b.relvecs[len - 2][j]));
        }

        return std::lexicographical_compare(array_a.begin(), array_a.end(),
                                            array_b.begin(), array_b.end());
    }
};


class Gruneisen : protected Pointers {
public:
    Gruneisen(class PHON *);

    ~Gruneisen();

    double delta_a;
    bool print_gruneisen;
    bool print_newfcs;

    void setup();

    std::complex<double> **gruneisen_bs;
    std::complex<double> **gruneisen_dos;

    void calc_gruneisen();

    void write_new_fcsxml_all();

private:
    void set_default_variables();

    void deallocate_variables();

    double **xshift_s;
    std::vector<FcsArrayWithCell> delta_fc2, delta_fc3;

    void prepare_delta_fcs(const std::vector<FcsArrayWithCell> &,
                           std::vector<FcsArrayWithCell> &,
                           const int) const;

    void calc_dfc2_reciprocal(std::complex<double> **,
                               const double *);

    void write_new_fcsxml(const std::string &,
                          double);

    std::string double2string(double) const;

    //  double calc_stress_energy2(const std::vector<FcsArrayWithCell>);
    void calc_stress_energy3(std::vector<FcsArrayWithCell>,
                             double ****);

    void print_stress_energy();
};
}
