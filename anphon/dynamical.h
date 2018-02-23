/*
 dynamical.h

 Copyright (c) 2014 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#pragma once

#include "pointers.h"
#include "fcs_phonon.h"
#include <vector>
#include <complex>
#include <string>

namespace PHON_NS
{
    class DistWithCell
    {
    public:
        int cell;
        double dist;

        DistWithCell();

        DistWithCell(const int n, const double d) : cell(n), dist(d) {};
    };

    inline bool operator<(const DistWithCell a,
                          const DistWithCell b)
    {
        return a.dist < b.dist;
    }

    class Dynamical : protected Pointers
    {
    public:
        Dynamical(class PHON *);
        ~Dynamical();

        unsigned int neval;
        bool eigenvectors;
        bool print_eigenvectors;
        unsigned int symmetrize_borncharge;
        unsigned int nonanalytic;
        bool participation_ratio;
        unsigned int band_connection;

        std::string file_born;
        double na_sigma;

        double **eval_phonon;
        int **index_bconnect;
        std::complex<double> ***evec_phonon;
        double dielec[3][3];
        double ***borncharge;

        void diagonalize_dynamical_all();

        void setup_dynamical(std::string);

        void eval_k(double *, double *,
                    std::vector<FcsClassExtent>,
                    double *, std::complex<double> **, bool);
        void modify_eigenvectors();
        void eval_k_ewald(double *, double *,
                          std::vector<FcsClassExtent>,
                          double *, std::complex<double> **, bool,
                          const int);


        double fold(const double);
        double freq(const double);

        void calc_participation_ratio_all(std::complex<double> ***,
                                          double **,
                                          double ***);

        void calc_analytic_k(double *,
                             const std::vector<FcsClassExtent> &,
                             std::complex<double> **);
        void calc_nonanalytic_k(double *, double *,
                                std::complex<double> **);
        void calc_nonanalytic_k2(double *, double *,
                                 std::complex<double> **);

        void calc_analytic_k_ewald(double *,
                                   std::vector<FcsClassExtent>,
                                   std::complex<double> **);

    private:
        void set_default_variables();
        void deallocate_variables();
        void load_born(const unsigned int);

        void prepare_mindist_list(std::vector<int> **);
        void calc_atomic_participation_ratio(std::complex<double> *, double *);
        double distance(double *, double *);
        void connect_band_by_eigen_similarity(std::complex<double> ***, int **);

        double **xshift_s;
        char UPLO;
        std::complex<double> ***dymat;
        std::vector<int> **mindist_list;
    };

    extern "C" {
    void zheev_(const char *jobz, const char *uplo, int *n, std::complex<double> *a, int *lda,
                double *w, std::complex<double> *work, int *lwork, double *rwork, int *info);
    void zgemm_(const char *transa, const char *transb, int *m, int *n, int *k,
                std::complex<double> *alpha, std::complex<double> *a, int *lda, std::complex<double> *b, int *ldb,
                std::complex<double> *beta, std::complex<double> *c, int *ldc);
    }
}
