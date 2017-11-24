/*
 patterndisp.h

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

namespace ALM_NS
{
    class DispAtomSet
    {
    public:
        std::vector<int> atomset;

        DispAtomSet();

        DispAtomSet(std::vector<int> vec)
        {
            for (auto it = vec.begin(); it != vec.end(); ++it) {
                atomset.push_back((*it));
            }
        }
    };

    class DirectionVec
    {
    public:
        double direction[3];

        DirectionVec();

        DirectionVec(double vec_in[3])
        {
            for (int i = 0; i < 3; ++i) direction[i] = vec_in[i];
        }
    };

    class DispDirectionHarmonic
    {
    public:
        int atom;
        std::vector<DirectionVec> directionlist;

        DispDirectionHarmonic();

        DispDirectionHarmonic(int n, std::vector<DirectionVec> list_in)
        {
            atom = n;
            for (auto it = list_in.begin(); it != list_in.end(); ++it) {
                directionlist.push_back(*it);
            }
        }
    };

    class AtomWithDirection
    {
    public:
        std::vector<int> atoms;
        std::vector<double> directions;

        AtomWithDirection();

        AtomWithDirection(std::vector<int> a, std::vector<double> b)
        {
            for (unsigned int i = 0; i < a.size(); ++i) {
                atoms.push_back(a[i]);
            }
            for (unsigned int i = 0; i < b.size(); ++i) {
                directions.push_back(b[i]);
            }
        }

        AtomWithDirection(int n, int *a, double **b)
        {
            for (int i = 0; i < n; ++i) {
                atoms.push_back(a[i]);
                for (int j = 0; j < 3; ++j) {
                    directions.push_back(b[i][j]);
                }
            }
        }
    };


    inline bool operator<(const DispAtomSet &a, const DispAtomSet &b)
    {
        return std::lexicographical_compare(a.atomset.begin(), a.atomset.end(),
                                            b.atomset.begin(), b.atomset.end());
    }


    class IndexWithSign
    {
    public:
        int ind, sign;

        IndexWithSign();

        IndexWithSign(const int ind_in, const int sign_in)
        {
            ind = ind_in;
            sign = sign_in;
        }
    };

    inline bool operator<(const IndexWithSign &a, const IndexWithSign &b)
    {
        return a.ind < b.ind;
    }

    class Displace: protected Pointers
    {
    public:
        Displace(class ALM *);
        ~Displace();

        bool trim_dispsign_for_evenfunc;

        std::string disp_basis;
        std::vector<AtomWithDirection> *pattern_all;
        void gen_displacement_pattern();

    private:
        std::vector<DispDirectionHarmonic> disp_harm, disp_harm_best;
        void generate_pattern_all(const int,
                                  std::vector<AtomWithDirection> *,
                                  std::set<DispAtomSet> *,
                                  const std::string);

        void generate_signvecs(const int,
                               std::vector<std::vector<int>> &,
                               std::vector<int>);

        void find_unique_sign_pairs(const int,
                                    std::vector<std::vector<int>>,
                                    std::vector<int>,
                                    std::vector<std::vector<int>> &,
                                    const std::string);
    };
}
