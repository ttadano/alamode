#ifndef ALM_FITTING_HEADER
#define ALM_FITTING_HEADER

#include "pointers.h"
#include <vector>
#include <set>

namespace ALM_NS {

    class Constraint {
    public:
        std::vector<double> w_const;

        Constraint();
        Constraint(const int n, const double *arr){
            for(int i = 0; i < n; ++i){
            w_const.push_back(arr[i]);
            }
        
        }
    
    };

   inline bool operator<(const Constraint a, const Constraint b){
        return std::lexicographical_compare(a.w_const.begin(), a.w_const.end(), b.w_const.begin(), b.w_const.end());
   }

   /* inline bool operator<(const Constraint a, const Constraint b){
        std::vector<double>::const_iterator aiter = a.w_const.begin(), biter = b.w_const.begin();
        do {
            if(*aiter > *biter + 1.0e-15){
            return false;
            ++aiter;
            ++biter;
            }
        } while(aiter != a.w_const.end() && biter != b.w_const.end());

        aiter = a.w_const.begin();
        biter = b.w_const.begin();

        do {
            if(std::abs(*aiter - *biter) > 1.0e-15){
                return true;
                ++aiter;
                ++biter;
            }
        } while(aiter != a.w_const.end() && biter != b.w_const.end());
        return false;
    }*/

    class Fitting: protected Pointers {
    public:
        Fitting(class ALM *);
        ~Fitting();

        void fitmain();

        int constraint;
   
    private:

        int inprim_index(const int);
        void wrtfcs(const double *);
        void fit_without_constraints(int, int);
        void fit_with_constraints(int, int);
        void calc_matrix_elements(const int, const int, const int, 
            const int, const int, const int, const int);
        double gamma(const int, const int *);
        int factorial(const int);

        void translational_invariance();
        void rotational_invariance();
        bool is_allzero(const int, const double *);

        double **amat;
        double *fsum;
        std::set<Constraint> *const_translation;
    };

    extern "C" void dgelss_(int *m, int *n, int *nrhs, double *a, int *lda,	
        double *b, int *ldb, double *s, double *rcond, int *rank,
        double *work,	int *lwork, int *info);

}
#endif