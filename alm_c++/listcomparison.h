#ifndef ALM_LISTCOMPARISON_HEADER
#define ALM_LISTCOMPARISON_HEADER

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

        IntList(const int n, const int *arr){
            for (int i = 0; i < n; ++i){
                iarray.push_back(arr[i]);
            }
        }
    };

    // 以下の演算子は inlineオプションをつけないと何故か多重定義判定されてしまいコンパイルできない．
    // 
    // inline bool operator<(const IntList a, const IntList b){
    //     return a.iarray < b.iarray;
    // };

    inline bool operator<(const IntList a, const IntList b){
                return std::lexicographical_compare(a.iarray.begin(), a.iarray.end(), b.iarray.begin(), b.iarray.end());
    };

    /*inline bool operator<(const std::vector<int> a, const std::vector<int> b){
    std::vector<int>::const_iterator aiter = a.begin(), biter = b.begin();
    do {
    if (*aiter > *biter) {
    return false;
    ++aiter;
    ++biter;
    }
    } while(aiter != a.end() && biter != b.end());

    aiter = a.begin();
    biter = b.begin();
    do {
    if (*aiter != *biter){
    return true;
    ++aiter;
    ++biter;
    }
    } while(aiter != a.end() && biter != b.end());
    return false;
    }*/

    inline std::ostream &operator<<(std::ostream &s, const IntList &o){

        for (unsigned int i = 0; i < o.iarray.size(); ++i){
            s << std::setw(5) << o.iarray[i];
        }
        s << std::endl;
        return s;
    }
}

#endif