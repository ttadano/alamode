/*
 combination.h

 Copyright (c) 2014 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#pragma once

#include<vector>
#include<set>
#include "../external/combination.hpp"

namespace ALM_NS {
template<class TYPE>
class CombinationWithRepetition {
    std::vector<TYPE> vec;
    unsigned int ndim;


public:
    CombinationWithRepetition() {};

    template<class InputIter>
    CombinationWithRepetition(InputIter begin,
                              InputIter end,
                              const unsigned int n)
    {
        ndim = n;

        // remove redundunt elements
        std::set<TYPE> set_tmp;
        for (InputIter iter = begin; iter != end; ++iter) set_tmp.insert(*iter);

        vec.clear();

        for (auto iter = set_tmp.begin(); iter != set_tmp.end(); ++iter) {
            for (unsigned int i = 0; i < ndim; i++) {
                vec.push_back(*iter);
            }
        }
    }

    bool next()
    {
        return boost::next_combination(vec.begin(), vec.begin() + ndim, vec.end());
    }

    std::vector<TYPE> now() const
    {
        return std::vector<TYPE>(vec.begin(), vec.begin() + ndim);
    }

    unsigned int size() const
    {
        const unsigned int n = vec.size() / ndim;
        const auto r = ndim;
        return factorial(n + r - 1, n - 1) / factorial(r);
    }

private:
    unsigned int factorial(const unsigned int max,
                           const unsigned int min = 1) const
    {
        unsigned int result = 1;
        for (auto i = min + 1; i <= max; ++i) result *= i;
        return result;
    }
};
}
