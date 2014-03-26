/*
 listcomparison.h

 Copyright (c) 2014 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#pragma once

#include <vector>
#include <iostream>
#include <iomanip>
#include <algorithm>

namespace ALM_NS
{
    class IntList{
    public:
        std::vector<int> iarray;

        IntList() {
            iarray.clear();
        }

        IntList(const IntList &a){
            for(std::vector<int>::const_iterator p = a.iarray.begin(); p != a.iarray.end(); ++p){
                iarray.push_back(*p);
            }
        }

        IntList(const int n, const int *arr){
            for (int i = 0; i < n; ++i){
                iarray.push_back(arr[i]);
            }
        }
    };

    // inline options below are necessary for successful compilation. (why?)

    inline bool operator<(const IntList a, const IntList b){
        return std::lexicographical_compare(a.iarray.begin(), a.iarray.end(), b.iarray.begin(), b.iarray.end());
    };

    inline std::ostream &operator<<(std::ostream &s, const IntList &o){

        for (unsigned int i = 0; i < o.iarray.size(); ++i){
            s << std::setw(5) << o.iarray[i];
        }
        s << std::endl;
        return s;
    }
}
