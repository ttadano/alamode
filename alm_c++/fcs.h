#ifndef ALM_FCS_HEADER
#define ALM_FCS_HEADER

#include "pointers.h"
#include <vector>
#include <set>
#include "listcomparison.h"

namespace ALM_NS {

    class FcProperty {
    public:
        std::vector<int> elems;
        double coef;

        FcProperty();
        FcProperty(const int n, const double c, const int *arr)
        {
            coef = c;
            for (int i = 0; i < n; ++i){
                elems.push_back(arr[i]);
            }
        }    
    };

    //class FcSet {
    //public:
    //    std::vector<FcProperty> fc_property;

    //    FcSet();
    //    FcSet(std::vector<FcProperty> fc_tmp){
    //        for (std::vector<FcProperty>::iterator iter = fc_tmp.begin(); iter != fc_tmp.end(); ++iter){
    //            fc_property.push_back(*iter);
    //        }
    //    }
    //};


    class Fcs: protected Pointers{
    public:
        Fcs(class ALM *);
        ~Fcs();

        void init();

        int *nzero;

        std::set<IntList> *pairs;
        std::vector<int> *ndup;
        std::vector<FcProperty> *fc_set;

    private:
        int *nints;
        void read_pairs(int);
        void generate_fclists(int);
        int min_inprim(const int, const int *);
        bool is_inprim(const int, const int *);
        double coef_sym(const int, const int, const int *, const int *);
        void get_xyzcomponent(int, int **);
        bool is_ascending(const int, const int *);
        void sort_tail(const int, int *); 

        std::string easyvizint(const int);
    };

}

#endif