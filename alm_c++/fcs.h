#ifndef ALM_FCS_HEADER
#define ALM_FCS_HEADER

#include "pointers.h"
#include <vector>
#include <set>
#include "listcomparison.h"
#include <algorithm>

namespace ALM_NS {

    class FcProperty {
    public:
        std::vector<int> elems;
        double coef;
        int mother;

        FcProperty();
        FcProperty(const FcProperty &obj){
            coef = obj.coef;
            mother = obj.mother;
            for (std::vector<int>::const_iterator it = obj.elems.begin(); it != obj.elems.end(); ++it){
                elems.push_back(*it);
            }
        };

        FcProperty(const int n, const double c, const int *arr, const int m)
        {
            coef = c;
            mother = m;
            for (int i = 0; i < n; ++i){
                elems.push_back(arr[i]);
            }
        }    
    };


    inline bool operator<(const FcProperty a, const FcProperty b){
        return std::lexicographical_compare(a.elems.begin(), a.elems.end(), b.elems.begin(), b.elems.end());
    }

    class Fcs: protected Pointers{
    public:
        Fcs(class ALM *);
        ~Fcs();

        void init();

        int *nzero;

        std::set<IntList> *pairs;
        std::vector<int> *ndup;
        std::vector<FcProperty> *fc_set;

        std::string easyvizint(const int);
        void get_xyzcomponent(int, int **);
        void sort_tail(const int, int *); 

    private:
        int *nints;
        void read_pairs(int);
        void generate_fclists(int);
        int min_inprim(const int, const int *);
        bool is_inprim(const int, const int *);
        bool is_inprim(const int);
        double coef_sym(const int, const int, const int *, const int *);

        bool is_ascending(const int, const int *);
        


    };

}

#endif