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
#include <Eigen/Core>

namespace PHON_NS {
    class DistWithCell {
    public:
        int cell;
        double dist;

        DistWithCell();

        DistWithCell(const int n,
                     const double d) : cell(n), dist(d) {};
    };

    inline bool operator<(const DistWithCell a,
                          const DistWithCell b)
    {
        return a.dist < b.dist;
    }

    class Dynamical : protected Pointers {
    public:
        Dynamical(class PHON *);

        ~Dynamical();

        unsigned int neval{};
        bool eigenvectors{};
        bool print_eigenvectors{};
        unsigned int symmetrize_borncharge{};
        unsigned int nonanalytic{};
        bool participation_ratio{};
        unsigned int band_connection{};

        std::string file_born;
        double na_sigma{};

        double **eval_phonon{};
        int **index_bconnect{};
        std::complex<double> ***evec_phonon{};
        double dielec[3][3]{};
        double ***borncharge{};

        bool **is_imaginary{};

        void diagonalize_dynamical_all();

        void setup_dynamical();

        void setup_dielectric(const unsigned int verbosity = 1);

        void eval_k(double *,
                    double *,
                    const std::vector<FcsClassExtent> &,
                    double *,
                    std::complex<double> **,
                    bool) const;

        void modify_eigenvectors() const;

        void eval_k_ewald(double *,
                          double *,
                          const std::vector<FcsClassExtent> &,
                          double *,
                          std::complex<double> **,
                          bool) const;

        double fold(const double) const;

        double freq(const double) const;

        void calc_participation_ratio_all(std::complex<double> ***,
                                          double **,
                                          double ***) const;

        void calc_analytic_k(const double *,
                             const std::vector<FcsClassExtent> &,
                             std::complex<double> **) const;

        void calc_nonanalytic_k(double *,
                                double *,
                                std::complex<double> **) const;

        void calc_nonanalytic_k2(const double *,
                                 double *,
                                 std::complex<double> **) const;

        void calc_analytic_k_ewald(double *,
                                   std::vector <FcsClassExtent>,
                                   std::complex<double> **);

        void project_degenerate_eigenvectors(double *xk_in,
                                             const std::vector <std::vector<double>> &project_directions,
                                             std::complex<double> **evec_out) const;

        std::vector <std::vector<double>> get_projection_directions() const;

        void set_projection_directions(const std::vector <std::vector<double>> projections_in);

    private:
        void set_default_variables();

        void deallocate_variables();

        void load_born(const unsigned int flag_symmborn,
                       const unsigned int verbosity = 1);

        void prepare_mindist_list(std::vector<int> **) const;

        void calc_atomic_participation_ratio(std::complex<double> *,
                                             double *) const;

        double distance(double *,
                        double *) const;

        void connect_band_by_eigen_similarity(std::complex<double> ***,
                                              int **) const;

        void detect_imaginary_branches(double **);

        std::vector <std::vector<double>> projection_directions;

        int transform_eigenvectors(double *xk_in,
                                   std::vector<double> perturb_direction,
                                   const double dk,
                                   Eigen::MatrixXcd &evec_sub) const;

        double **xshift_s;
        char UPLO{};
        std::complex<double> ***dymat{};
        std::vector<int> **mindist_list{};
    };

    extern "C" {
    void zheev_(const char *jobz,
                const char *uplo,
                int *n,
                std::complex<double> *a,
                int *lda,
                double *w,
                std::complex<double> *work,
                int *lwork,
                double *rwork,
                int *info);

    void zgemm_(const char *transa,
                const char *transb,
                int *m,
                int *n,
                int *k,
                std::complex<double> *alpha,
                std::complex<double> *a,
                int *lda,
                std::complex<double> *b,
                int *ldb,
                std::complex<double> *beta,
                std::complex<double> *c,
                int *ldc);
    }
}
