#pragma once

#include "pointers.h"
#include <string>
#include <vector>
#include <set>

namespace PHON_NS {

    struct Triplet {
        unsigned int atom, cell, xyz;
    };

    class FcsClass {
    public:
        std::vector<Triplet> elems;
        double fcs_val;

        FcsClass(){};
        FcsClass(const FcsClass &obj){
            fcs_val = obj.fcs_val;
            for (std::vector<Triplet>::const_iterator it = obj.elems.begin(); it != obj.elems.end(); ++it){
                elems.push_back(*it);
            }
        }

        FcsClass(const unsigned int n, const double val, const Triplet *arr)
        {
            fcs_val = val;
            for (unsigned int i = 0; i < n; ++i){
                elems.push_back(arr[i]);
            }
        }

        FcsClass(const double val, const std::vector<Triplet> vec)
        {
            fcs_val = val;
            for (std::vector<Triplet>::const_iterator it = vec.begin(); it != vec.end(); ++it){
                elems.push_back(*it);
            }
        }

    };

    class Fcs_phonon: protected Pointers {
    public:
        Fcs_phonon(class PHON *);
        ~Fcs_phonon();

        void setup(std::string);
        unsigned int maxorder;
        std::string file_fcs;
        double ****fc2;
        std::vector<FcsClass> *force_constant;

    private:
        void load_fc2();
        void load_fcs();
        unsigned int coordinate_index(const char);
    };
}
