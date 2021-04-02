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
#include "anharmonic_core.h"
#include <complex>
#include <Eigen/Dense>

namespace PHON_NS {
    class DistList {
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

    struct ShiftCell {
    public:
        int sx, sy, sz;
    };

    struct MinimumDistList {
    public:
        double dist;
        std::vector<ShiftCell> shift;
    };

    struct KpointSymmetry {
    public:
        unsigned int symmetry_op;
        unsigned int knum_irred_orig;
        unsigned int knum_orig;
    };

    class Scph : protected Pointers {
    public:

        Scph(class PHON *phon);

        ~Scph();

        unsigned int kmesh_scph[3];
        unsigned int kmesh_interpolate[3];
        unsigned int ialgo;
        unsigned int bubble;

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
        double *invmass_v3;
        double *invmass_v4;
        int **evec_index_v3;
        int **evec_index_v4;
        int ngroup_v3, ngroup_v4;
        std::vector<RelativeVector> *relvec_v3, *relvec_v4;
        std::complex<double> *phi3_reciprocal;
        int kindex_phi3_stored[2] = {-1, -1};

        std::vector<double> *fcs_group_v3;
        std::vector<double> *fcs_group_v4;
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

        void postprocess(std::complex<double> ****delta_dymat_scph,
                         std::complex<double> ****delta_dymat_scph_plus_bubble);

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
                                          bool &,
                                          std::complex<double> ***,
                                          bool,
                                          const unsigned int verbosity);

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

        void find_degeneracy(std::vector<int> *degeneracy_out,
                             unsigned int nk_in,
                             double **eval_in) const;

        static double distance(double *,
                               double *);

        void symmetrize_dynamical_matrix(unsigned int,
                                         Eigen::MatrixXcd &) const;

        void replicate_dymat_for_all_kpoints(std::complex<double> ***) const;

        static void duplicate_xk_boundary(double *,
                                          std::vector<std::vector<double>> &);

        void write_anharmonic_correction_fc2(std::complex<double> ****delta_dymat,
                                             const unsigned int NT,
                                             const int type = 0);

        static void mpi_bcast_complex(std::complex<double> ****data,
                                      const unsigned int NT,
                                      const unsigned int nk,
                                      const unsigned int ns);


        void compute_free_energy_bubble_SCPH(const unsigned int [3],
                                             std::complex<double> ****);

        void bubble_correction(std::complex<double> ****,
                               std::complex<double> ****);

        std::complex<double> V3_this(const unsigned int ks[3],
                                     double **eval,
                                     std::complex<double> ***evec);

        void calc_phi3_reciprocal_this(const unsigned int ik1,
                                             const unsigned int ik2,
                                             std::complex<double> *ret);

        std::vector<std::complex<double>> get_bubble_selfenergy(const unsigned int nk_in,
                                                                      const unsigned int ns_in,
                                                                      const unsigned int kmesh_in[3],
                                                                      double **xk_in,
                                                                      double **eval_in,
                                                                      std::complex<double> ***evec_in,
                                                                      const unsigned int knum,
                                                                      const unsigned int snum,
                                                                      const double temp_in,
                                                                      const std::vector<std::complex<double>> &omegalist);



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
