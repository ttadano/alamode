#pragma once

#include "pointers.h"
#include "fcs_phonon.h"
#include <vector>
#include <complex>
#include <string>

namespace PHON_NS {
    class Dynamical: protected Pointers {
    public:
        Dynamical(class PHON *);
        ~Dynamical();

        void diagonalize_dynamical_all();
        void finish_dynamical();

        unsigned int neval;
        bool eigenvectors, nonanalytic;
        bool print_eigenvectors;

        std::string file_born;
        double na_sigma;

        double **eval_phonon;
        std::complex<double> ***evec_phonon;

        void setup_dynamical(std::string);

        void eval_k(double *, double *, double ****, double *, std::complex<double> **, bool);
        void eval_k(double *, double *, std::vector<FcsClassExtent>, double *, std::complex<double> **, bool);


        double fold(double);
        double freq(const double);

    private:
        //        void calc_analytic();

        void load_born();
        void calc_analytic_k(double *, double ****, std::complex<double> **);
        void calc_analytic_k(double *, std::vector<FcsClassExtent>, std::complex<double> **);
        void calc_nonanalytic_k(double *, double *, double **);
        void modify_eigenvectors();
        void modify_eigenvectors_sym();

        double **xshift_s;
        char UPLO;
        std::complex<double> ***dymat;
        double dielec[3][3];
        double ***borncharge;
    };

    extern "C" {
        void zheev_(const char *jobz, const char *uplo, int *n,	std::complex<double> *a, int *lda, 
            double *w, std::complex<double> *work, int *lwork, double *rwork, int *info);
        void zgemm_(const char *transa, const char *transb, int *m, int *n, int *k, 
            std::complex<double> *alpha, std::complex<double> *a, int *lda, std::complex<double> *b, int *ldb, 
            std::complex<double> *beta, std::complex<double> *c, int *ldc);
    }
}
