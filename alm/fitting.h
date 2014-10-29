/*
 fitting.h

 Copyright (c) 2014 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

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

        double *params;
        unsigned int nboot;
        unsigned int seed;

#ifdef _VSL
        VSLStreamStatePtr stream;
        int brng;
#endif

    private:

        int inprim_index(const int);
        void wrtfcs(const double *);
        void data_multiplier(const int, const int, const int, const int, const int, int &, const int, double **&, double **&);
        void fit_without_constraints(int, int, int, double **, double *);
        void fit_with_constraints(int, int, int, int, double **, double *, double **, double *);
        void fit_consecutively(int, int, const int, const int,
            const int, const int, double **, double *, double **, double *);
        void calc_matrix_elements(const int, const int, const int, 
            const int, const int, const int, const int, double **, double **, double **, double *);
        void fit_bootstrap(int, int, int, int, int, double **, double *, double **, double *);
        double gamma(const int, const int *);
        int factorial(const int);
        int rankSVD(const int, const int, double *, const double);
        int rankQRD(const int, const int, double *, const double);
        int rankSVD2(const int, const int, double **, const double);
#ifdef _USE_EIGEN
        int getRankEigen(const int, const int, double **);
#endif
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

        void dgeqp3_(int *m, int *n, double *a, int *lda, int *jpvt,
            double *tau, double *work, int *lwork, int *info);
    }

}
