/*
 scph.h

 Copyright (c) 2015 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#pragma once

#include "pointers.h"
#include "kpoint.h"
#include <complex>
#include <Eigen/Dense>

namespace PHON_NS
{
    class DistList
    {
    public:
        unsigned int cell_s;
        double dist;

        DistList();

        DistList(const unsigned int cell_s_,
                 const double dist_) : cell_s(cell_s_), dist(dist_) {};

        bool operator<(const DistList &obj) const
        {
            return dist < obj.dist;
        }
    };

    struct ShiftCell
    {
    public:
        int sx, sy, sz;
    };

    struct MinimumDistList
    {
    public:
        double dist;
        std::vector<ShiftCell> shift;
    };

    struct KpointSymmetry
    {
    public:
        unsigned int symmetry_op;
        unsigned int knum_irred_orig;
        unsigned int knum_orig;
    };

    class Scph : protected Pointers
    {
    public:

        Scph(class PHON *phon);

        ~Scph();

        unsigned int kmesh_scph[3];
        unsigned int kmesh_interpolate[3];
        unsigned int ialgo;

        bool restart_scph;
        bool warmstart_scph;
        bool lower_temp;
        double tolerance_scph;

        double **xk_scph, **kvec_na_scph;
        double **xk_interpolate;
        std::vector<std::vector<KpointList>> kp_irred_scph;
        std::vector<std::vector<KpointList>> kp_irred_interpolate;

        void exec_scph();
        void setup_scph();

        double mixalpha;
        unsigned int maxiter;
        bool print_self_consistent_fc2;
        bool selfenergy_offdiagonal;
        bool relax_coordinate;

    private:

        // Information of kmesh for SCPH calculation
        unsigned int nk_scph;
        unsigned int nk_interpolate;
        int *kmap_interpolate_to_scph;

        // Information for calculating the ph-ph interaction coefficients
        double ***vec_for_v3, *invmass_for_v3;
        double ***vec_for_v4, *invmass_for_v4;
        int **evec_index3;
        int **evec_index4;
        int ngroup, ngroup2;
        std::vector<double> *fcs_group;
        std::vector<double> *fcs_group2;
        std::complex<double> *exp_phase, ***exp_phase3;
        int nk_grid[3];
        int nk_represent;
        unsigned int tune_type;
        double dnk[3];

        // Information of harmonic dynamical matrix
        std::complex<double> im;
        double **omega2_harmonic;
        std::complex<double> ***evec_harmonic;
        MinimumDistList ***mindist_list_scph;

        // Local variables for handling symmetry of dynamical matrix
        std::complex<double> ****mat_transform_sym;
        std::vector<int> *small_group_at_k;
        std::vector<int> *symop_minus_at_k;
        KpointSymmetry *kpoint_map_symmetry;


        void set_default_variables();
        void deallocate_variables();

        void setup_kmesh();
        void setup_eigvecs();
        void setup_pp_interaction();
        void setup_transform_ifc();
        void setup_transform_symmetry();


        void load_scph_dymat_from_file(std::complex<double> ****);
        void store_scph_dymat_to_file(std::complex<double> ****);

        void exec_scph_main(std::complex<double> ****);
        void postprocess(std::complex<double> ****delta_dymat_scph);

        void compute_V4_elements_mpi_over_kpoint(std::complex<double> ***,
                                                 std::complex<double> ***,
                                                 bool,
                                                 bool);

        void compute_V4_elements_mpi_over_band(std::complex<double> ***,
                                               std::complex<double> ***,
                                               bool);

        void compute_V3_elements_mpi_over_kpoint(std::complex<double> ***,
                                                 std::complex<double> ***,
                                                 bool);

        void zerofill_elements_acoustic_at_gamma(double **,
                                                 std::complex<double> ***,
                                                 const int) const;

        void calc_new_dymat_with_evec(std::complex<double> ***,
                                      double **,
                                      std::complex<double> ***);

        void compute_anharmonic_frequency(std::complex<double> ***,
                                          double **,
                                          std::complex<double> ***,
                                          double,
                                          std::vector<int> *,
                                          bool &,
                                          std::complex<double> ***,
                                          bool);

        void exec_interpolation(const unsigned int [3],
                                std::complex<double> ***,
                                unsigned int,
                                double **,
                                double **,
                                double **,
                                std::complex<double> ***);

        void r2q(const double *,
                 unsigned int,
                 unsigned int,
                 unsigned int,
                 unsigned int,
                 std::complex<double> ***,
                 std::complex<double> **) const;

        void diagonalize_interpolated_matrix(std::complex<double> **,
                                             double *,
                                             std::complex<double> **,
                                             bool) const;

        void find_degeneracy(std::vector<int> *,
                             unsigned int,
                             const std::vector<std::vector<KpointList>> &,
                             double **) const;

        double distance(double *,
                        double *) const;

        void symmetrize_dynamical_matrix(unsigned int,
                                         Eigen::MatrixXcd &) const;

        void replicate_dymat_for_all_kpoints(std::complex<double> ***) const;

        void duplicate_xk_boundary(double *,
                                   std::vector<std::vector<double>> &) const;

        void write_anharmonic_correction_fc2(std::complex<double> ****,
                                             unsigned int);

        void mpi_bcast_complex(std::complex<double> ****,
                               int,
                               int,
                               int) const;


        void compute_free_energy_bubble_SCPH(const unsigned int [3],
                                             std::complex<double> ****);
    };

    extern "C" {
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
