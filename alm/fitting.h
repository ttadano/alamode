#pragma once

#include "pointers.h"
#include <vector>
#include <set>
#include <string>
#ifdef _VSL
#include "mkl_vsl.h"
#endif

namespace ALM_NS {

    class Fitting: protected Pointers {
    public:
        Fitting(class ALM *);
        ~Fitting();

        void fitmain();
        // int rank(const int, const int, double **);
        int rank(int, int, double *);
        int rank2(const int, const int, double **);

        // int getRankEigen(const int, const int,const  int);

        double *params;
        double *fc2_ref;
        unsigned int nboot;

        double **u, **f;

#ifdef _VSL
        VSLStreamStatePtr stream;
        int brng;
#endif
        unsigned int seed;

    private:

        int inprim_index(const int);
        void wrtfcs(const double *);
        void data_multiplier(const int, const int, const int, const int, const int, int &, const int);
        void fit_without_constraints(int, int, int, double **, double *);
        void fit_with_constraints(int, int, int, int, double **, double *, double **, double *);
        void fit_consecutively(int, int, const int, const int,
            const int, const int, double **, double *, double **, double *);
        void calc_matrix_elements(const int, const int, const int, 
            const int, const int, const int, const int, double **, double *);
        void fit_bootstrap(int, int, int, int, int, double **, double *, double **, double *);
        double gamma(const int, const int *);
        int factorial(const int);       
    };

    extern "C" {

        void dgelss_(int *m, int *n, int *nrhs, double *a, int *lda,	
            double *b, int *ldb, double *s, double *rcond, int *rank,
            double *work, int *lwork, int *info);

        void dgglse_(int *m, int *n, int *p, double *a, int *lda,
            double *b, int *ldb, double *c, double *d, double *x,
            double *work, int *lwork, int *info);

        void dgesdd_(const char *jobz, int *m, int *n, double *a, int *lda,
            double *s, double *u, int *ldu, double *vt, int *ldvt, double *work,
            int *lwork, int *iwork, int *info);

        void dgeqrf_(int *m, int *n, double *a, double *tau,
            double *work, int *lwork, int *info);
    }

}
