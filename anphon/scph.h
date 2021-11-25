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
#include "gruneisen.h"
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

    void exec_scph();

    void setup_scph();

    double mixalpha;
    unsigned int maxiter;
    bool print_self_consistent_fc2;
    bool selfenergy_offdiagonal;
    bool relax_coordinate;

    // variables related to structural optimization
    int relax_algo;
    int max_str_iter;
    double str_conv_tol;
    int set_init_str;
    double mixing_beta;
    double alpha_steepest_decent;

private:

    // Information of kmesh for SCPH calculation
    KpointMeshUniform *kmesh_coarse = nullptr;
    KpointMeshUniform *kmesh_dense = nullptr;
    int *kmap_interpolate_to_scph = nullptr;

    // Information for calculating the ph-ph interaction coefficients
    std::complex<double> *phi3_reciprocal, *phi4_reciprocal;

    // Phase shift
    PhaseFactorStorage *phase_factor_scph;

    // Information of harmonic dynamical matrix
    double **omega2_harmonic;
    std::complex<double> ***evec_harmonic;
    MinimumDistList ***mindist_list_scph;

    // Local variables for handling symmetry of dynamical matrix
    std::complex<double> ****mat_transform_sym;
    std::vector<int> *symop_minus_at_k;
    KpointSymmetry *kpoint_map_symmetry;

    int compute_Cv_anharmonic;

    void set_default_variables();

    void deallocate_variables();

    void setup_kmesh();

    void setup_eigvecs();

    void setup_pp_interaction();

    void setup_transform_ifc();

    void setup_transform_symmetry();

    void load_scph_dymat_from_file(std::complex<double> ****);

    // from here
    void load_scph_dymat_from_file(std::complex<double> ****,
                                       std::string);

    void store_scph_dymat_to_file(const std::complex<double> *const *const *const *dymat_in,
                                    std::string);

    void zerofill_harmonic_dymat_renormalize(std::complex<double> ****, 
                                            unsigned int);
                                                
    void exec_scph_main(std::complex<double> ****);

    void exec_scph_relax_main(std::complex<double> ****,
                                std::complex<double> ****);

    void exec_scph_relax_cell_coordinate_main(std::complex<double> ****,
                                          std::complex<double> ****);

    void read_C1_array(double *);

    void read_elastic_constants(double **, 
                                  double ***);

    void read_initial_q0(double *);

    void read_initial_strain(double **);

    void read_Tdep_initial_q0_from_u(double *, int);

    // void read_str_opt_input(int &str_opt_algo, 
    //                             int &max_str_loop, 
    //                             double &alpha_steepest_decent, 
    //                             double &mixing_beta, 
    //                             double &dq0_threashold);

    void read_cell_opt_input(double &,
                             double &);

    void calculate_u0(double *, double *);

    // void print_force(std::complex<double> *);
    void calculate_force_in_real_space(std::complex<double> *, double *);

    // void postprocess(std::complex<double> ****delta_dymat_scph,
    //                  std::complex<double> ****delta_dymat_scph_plus_bubble);
    void postprocess(std::complex<double> ****,
                        std::complex<double> ****,
                        std::complex<double> ****);
                         
    void compute_V4_elements_mpi_over_kpoint(std::complex<double> ***,
                                             std::complex<double> ***,
                                             bool,
                                             bool);
 
    void compute_V4_elements_mpi_over_band(std::complex<double> ***,
                                           std::complex<double> ***,
                                           bool);

    void compute_V3_elements_mpi_over_kpoint(std::complex<double> ***v3_out,
                                             const std::complex<double> *const *const *evec_in,
                                             const bool self_offdiag);

    void compute_V3_elements_for_given_IFCs(std::complex<double> ***v3_out,
                                            const int ngroup_v3_in,
                                            std::vector<double> *fcs_group_v3_in,
                                            std::vector<RelativeVector> *relvec_v3_in,
                                            double *invmass_v3_in,
                                            int **evec_index_v3_in,
                                            const std::complex<double> *const *const *evec_in,
                                            const bool self_offdiag);

    void zerofill_elements_acoustic_at_gamma(double **,
                                             std::complex<double> ***,
                                             const int) const;

    void calc_new_dymat_with_evec(std::complex<double> ***,
                                  double **,
                                  std::complex<double> ***);
    void compute_del_v1_strain_from_harmonic(std::complex<double> **,
                                               std::complex<double> ***);

    void compute_del_v1_strain_from_cubic(std::complex<double> **,
                                               std::complex<double> ***);

    void compute_del_v1_strain_from_quartic(std::complex<double> **,
                                               std::complex<double> ***);

    void compute_del_v2_strain_from_cubic(std::complex<double> ***,
                                          std::complex<double> ***);
    
    void compute_del_v2_strain_from_quartic(std::complex<double> ***,
                                            std::complex<double> ***);

    void compute_del_v3_strain_from_quartic(std::complex<double> ****,
                                            std::complex<double> ***);

    void compute_del_v_strain_in_real_space1(const std::vector<FcsArrayWithCell> &,
                                               std::vector<FcsArrayWithCell> &,
                                               const int,
                                               const int,
                                               const int);
    
    void compute_del_v_strain_in_real_space2(const std::vector<FcsArrayWithCell> &,
                                               std::vector<FcsArrayWithCell> &,
                                               const int,
                                               const int,
                                               const int,
                                               const int,
                                               const int);

    FcsClassExtent from_FcsArrayWithCell_to_FcsClassExtent(const FcsArrayWithCell &);

    void calculate_del_v0_strain_with_strain_displace(std::complex<double> *, 
                                               double *,
                                               double **,
                                               double ***,
                                               double **,
                                               double **,
                                               std::complex<double> **,
                                               std::complex<double> **,
                                               std::complex<double> **,
                                               std::complex<double> ***,
                                               std::complex<double> ***,
                                               std::complex<double> ****,
                                               double *);

    void calculate_eta_tensor(double **, 
                                double **);

    void renormalize_v0_from_strain(double &, 
                                    double , 
                                    double **, 
                                    double *,
                                    double **, 
                                    double ***);

    void renormalize_v1_array_from_strain(std::complex<double> *, 
                                            std::complex<double> *,
                                            std::complex<double> **, 
                                            std::complex<double> **, 
                                            std::complex<double> **, 
                                            double **);

    void renormalize_v2_array_from_strain(std::complex<double> **, 
                                            std::complex<double> ***, 
                                            std::complex<double> ***,
                                            double **);

    void renormalize_v3_array_from_strain(std::complex<double> ***, 
                                            std::complex<double> ***, 
                                            std::complex<double> ****,
                                            double **);

    void renormalize_v1_array(std::complex<double> *, 
                            std::complex<double> *, 
                            std::complex<double> **,
                            std::complex<double> ***, 
                            std::complex<double> ***,
                            double *);

    void renormalize_v2_array(std::complex<double> **, 
                            std::complex<double> **,
                            std::complex<double> ***, 
                            std::complex<double> ***,  
                            double *);

    void renormalize_v3_array(std::complex<double> ***,
                            std::complex<double> ***, 
                            std::complex<double> ***, 
                            double *);

    void renormalize_v0(double &,
                        double ,
                        std::complex<double> *,
                        std::complex<double> **,
                        std::complex<double> ***,
                        std::complex<double> ***,
                        double *);

    void compute_anharmonic_v1_array(std::complex<double> *, 
                                        std::complex<double> ***, 
                                        std::complex<double> ***, 
                                        double ** , 
                                        const double ,
                                        std::complex<double> *);

    void compute_anharmonic_del_v0_strain(std::complex<double> *, 
                                            std::complex<double> *,
                                            std::complex<double> ***,
                                            std::complex<double> ***,
                                            std::complex<double> ****,
                                            double **,
                                            double *,
                                            std::complex<double> ***, 
                                            double ** , 
                                            const double);

    void compute_anharmonic_frequency(std::complex<double> ***,
                                        double **,
                                        std::complex<double> ***,
                                        double,
                                        bool &,
                                        std::complex<double> ***,
                                        bool,
                                        std::complex<double> **,
                                        const unsigned int verbosity);

    void compute_renormalized_harmonic_frequency(double **,
                                    std::complex<double> ***,
                                    std::complex<double> **,
                                    const unsigned int );

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

    void print_distance_harmonic_IFC();

    void write_anharmonic_correction_fc2(std::complex<double> ****delta_dymat,
                                            const unsigned int NT,
                                            const int type = 0);

    static void mpi_bcast_complex(std::complex<double> ****data,
                                    const unsigned int NT,
                                    const unsigned int nk,
                                    const unsigned int ns);

    void get_derivative_central_diff(const double delta_t,
                                        const unsigned int nk,
                                        double **omega0,
                                        double **omega2,
                                        double **domega_dt);

    void compute_free_energy_bubble_SCPH(const unsigned int [3],
                                            std::complex<double> ****);

    void bubble_correction(std::complex<double> ****,
                            std::complex<double> ****);

    std::vector<std::complex<double>> get_bubble_selfenergy(const KpointMeshUniform *kmesh_in,
                                                            const unsigned int ns_in,
                                                            const double *const *eval_in,
                                                            const std::complex<double> *const *const *evec_in,
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