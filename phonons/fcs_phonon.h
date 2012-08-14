#pragma once

#include "pointers.h"
#include <string>
#include <vector>
#include <set>

namespace PHON_NS {

    class FcsClass {
    public:
        std::vector<int> elems;
        double fcs_val;

        FcsClass();
        FcsClass(const FcsClass &obj){
            fcs_val = obj.fcs_val;
            for (std::vector<int>::const_iterator it = obj.elems.begin(); it != obj.elems.end(); ++it){
            elems.push_back(*it);
            }
        }

        FcsClass(const int n, const double val, const int *arr)
        {
            fcs_val = val;
            for (int i = 0; i < n; ++i){
            elems.push_back(arr[i]);
            }
        }

    };

    inline bool operator<(const FcsClass a, const FcsClass b){
        return std::lexicographical_compare(a.elems.begin(), a.elems.end(), b.elems.begin(), b.elems.end());
    }

    class Fcs_phonon: protected Pointers {
    public:
        Fcs_phonon(class PHON *);
        ~Fcs_phonon();

        void setup();

        std::string file_fcs;
        double ****fc2;
        std::vector<FcsClass> *force_constant;

    private:
        void load_fc2();
        void load_fcN();
        unsigned int coordinate_index(const char);
    };
}