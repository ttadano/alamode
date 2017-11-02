/*
 fcs_phonon.h

 Copyright (c) 2014 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#pragma once

#include "pointers.h"
#include <string>
#include <vector>
#include <set>

namespace PHON_NS
{
    struct Triplet
    {
        unsigned int atom, cell, xyz;
    };

    class FcsClass
    {
    public:
        std::vector<Triplet> elems;
        double fcs_val;

        FcsClass()
        {
        };

        FcsClass(const FcsClass &obj)
        {
            fcs_val = obj.fcs_val;
            for (auto it = obj.elems.cbegin(); it != obj.elems.cend(); ++it) {
                elems.push_back(*it);
            }
        }

        FcsClass(const unsigned int n, const double val, const Triplet *arr)
        {
            fcs_val = val;
            for (unsigned int i = 0; i < n; ++i) {
                elems.push_back(arr[i]);
            }
        }

        FcsClass(const double val, const std::vector<Triplet> vec)
        {
            fcs_val = val;
            for (auto it = vec.cbegin(); it != vec.cend(); ++it) {
                elems.push_back(*it);
            }
        }
    };

    inline bool operator<(const FcsClass &a, const FcsClass &b)
    {
        std::vector<int> a_tmp, b_tmp;
        a_tmp.clear();
        b_tmp.clear();
        for (int i = 0; i < a.elems.size(); ++i) {
            a_tmp.push_back(3 * a.elems[i].atom + a.elems[i].xyz);
            b_tmp.push_back(3 * b.elems[i].atom + b.elems[i].xyz);
        }
        return lexicographical_compare(a_tmp.begin(), a_tmp.end(), b_tmp.begin(), b_tmp.end());
    }

    class FcsClassExtent
    {
    public:
        unsigned int atm1, atm2;
        unsigned int xyz1, xyz2;
        unsigned int cell_s;
        double fcs_val;

        FcsClassExtent()
        {
        };

        FcsClassExtent(const FcsClassExtent &obj)
        {
            atm1 = obj.atm1;
            atm2 = obj.atm2;
            xyz1 = obj.xyz1;
            xyz2 = obj.xyz2;
            cell_s = obj.cell_s;
            fcs_val = obj.fcs_val;
        }

        bool operator==(const FcsClassExtent &a) const
        {
            return (this->atm1 == a.atm1) & (this->atm2 == a.atm2)
                & (this->xyz1 == a.xyz1) & (this->xyz2 == a.xyz2)
                & (this->cell_s == a.cell_s);
        }
    };

    struct AtomCellSuper
    {
        unsigned int index;
        unsigned int tran;
        unsigned int cell_s;
    };

    inline bool operator<(const AtomCellSuper &a, const AtomCellSuper &b)
    {
        return a.index < b.index;
    }

    class FcsArrayWithCell
    {
    public:
        std::vector<AtomCellSuper> pairs;
        double fcs_val;

        FcsArrayWithCell()
        {
        };

        FcsArrayWithCell(const double fcs_in, const std::vector<AtomCellSuper> pairs_in)
        {
            fcs_val = fcs_in;

            for (auto it = pairs_in.cbegin(); it != pairs_in.cend(); ++it) {
                pairs.push_back(*it);
            }
        }
    };

    inline bool operator<(const FcsArrayWithCell &a, const FcsArrayWithCell &b)
    {
        std::vector<unsigned int> index_a, index_b;
        index_a.clear();
        index_b.clear();
        for (int i = 0; i < a.pairs.size(); ++i) {
            index_a.push_back(a.pairs[i].index);
            index_b.push_back(b.pairs[i].index);
        }
        for (int i = 0; i < a.pairs.size(); ++i) {
            index_a.push_back(a.pairs[i].tran);
            index_a.push_back(a.pairs[i].cell_s);
            index_b.push_back(b.pairs[i].tran);
            index_b.push_back(b.pairs[i].cell_s);
        }
        return lexicographical_compare(index_a.begin(), index_a.end(), index_b.begin(), index_b.end());
    }


    class Fcs_phonon: protected Pointers
    {
    public:
        Fcs_phonon(class PHON *);
        ~Fcs_phonon();

        void setup(std::string);
        unsigned int maxorder;
        std::string file_fcs, file_fc2;

        std::vector<FcsClass> *force_constant;
        std::vector<FcsArrayWithCell> *force_constant_with_cell;
        std::vector<FcsClassExtent> fc2_ext;

        bool update_fc2;

    private:
        bool require_cubic;
        bool require_quartic;

        void load_fc2_xml();
        void load_fcs_xml();

        void examine_translational_invariance(const int, const unsigned int, const unsigned int,
                                              double *,
                                              std::vector<FcsClassExtent> &,
                                              std::vector<FcsArrayWithCell> *);


        void MPI_Bcast_fc_class(const unsigned int);
        void MPI_Bcast_fcs_array(const unsigned int);
        void MPI_Bcast_fc2_ext();
    };
}
