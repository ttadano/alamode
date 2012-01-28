#ifndef ALM_FITTING_HEADER
#define ALM_FITTING_HEADER

#include "pointers.h"

namespace ALM_NS {
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
        void fit_with_constraints();
        void calc_matrix_elements(const int, const int, const int, 
            const int, const int, const int, const int);
        double gamma(const int, const int *);
        int factorial(const int);

        double **amat;
        double *fsum;
    };

    extern "C" void dgelss_(int *m, int *n, int *nrhs, double *a, int *lda,	
        double *b, int *ldb, double *s, double *rcond, int *rank,
        double *work,	int *lwork, int *info);

}
#endif