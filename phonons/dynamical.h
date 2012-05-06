#pragma once

#include "pointers.h"
#include <complex>
#include <string>

namespace PHON_NS {
    class Dynamical: protected Pointers {
    public:
        Dynamical(class PHON *);
        ~Dynamical();

        void calc_dynamical_matrix();
        void diagonalize_dynamical();

        unsigned int neval;
        bool eigenvectors, nonanalytic;

        std::string file_born;
        double na_sigma;

        double **eval_phonon;
        std::complex<double> ***dymat;

        void eval_k(double *, double *);
        void setup_dynamical();

    private:
        void calc_analytic();
        void calc_analytic_k(std::complex<double> **, double *);
        void calc_nonanalytic();
        double fold(double);
    };

    extern "C" {
    void zheev_(const char *jobz, const char *uplo, int *n,	std::complex<double> *a, int *lda, 
        double *w, std::complex<double> *work, int *lwork, double *rwork, int *info);
    }
}