#ifndef ALM_LIST_COMPARISON
#define ALM_LIST_COMPARISON

#include <vector>
#include <iostream>
#include <iomanip>

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

    bool operator<(const IntList a, const IntList b){
        return a.iarray < b.iarray;
    };

    bool operator<(std::vector<int> a, std::vector<int> b){
        std::vector<int>::iterator aiter = a.begin(), biter = b.begin();
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
    }

    std::ostream &operator<<(std::ostream &s, const IntList &o){

        for (int i = 0; i < o.iarray.size(); ++i){
            s << std::setw(5) << o.iarray[i];
        }
        s << std::endl;
        return s;
    }
}

#endif