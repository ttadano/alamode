#ifndef ALM_FCS_HEADER
#define ALM_FCS_HEADER

#include "pointers.h"
#include <vector>
#include <set>
#include "listcomparison.h"

namespace ALM_NS {

    class FcSet{
    public:
        class FcList {
        public:
            std::vector<int> elems;
            double coef;

            FcList();
            FcList(const int n, const double c, const int *arr)
            {
                coef = c;
                for (int i = 0; i < n; ++i){
                    elems.push_back(arr[i]);
                }
            }    
        };

        std::vector<FcList> *fcset;
    };


    class Fcs: protected Pointers{
    public:
        Fcs(class ALM *);
        ~Fcs();

        void init();

        //        std::vector<Pairs> *pairs;
        std::set<IntList> *pairs;
    private:
        int *nints;
        void read_pairs(int);
        void generate_fclists(int);
        int min_inprim(const int, const int *);
        bool is_inprim(const int, const int *);
        double coef_sym(const int, const int, const int *, const int *);
        void get_xyzcomponent(int, int **);
        bool is_ascending(const int, const int *);
       };

}

#endif