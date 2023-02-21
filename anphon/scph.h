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
    int relax_coordinate;

    // variables related to structural optimization
    int relax_algo;
    int max_str_iter;
    double coord_conv_tol;
    double mixbeta_coord;
    double alpha_steepest_decent;
    double cell_conv_tol;
    double mixbeta_cell;

    int set_init_str;
    int cooling_u0_index; // used if set_init_str is 3 
    double cooling_u0_thr; // used if set_init_str is 3
    double add_hess_diag;
    double stat_pressure;

    // optimization scheme used in QHA
    int qha_scheme;
    
    // options of IFC renormalization
    int renorm_3to2nd;
    int renorm_2to1st;
    int renorm_anharmto1st;

    // initial strain and displacement
    double **init_u_tensor;
    double *init_u0;
    int natmin_tmp;

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

    void load_scph_dymat_from_file(std::complex<double> ****,
                                       std::string);

    void store_scph_dymat_to_file(const std::complex<double> *const *const *const *dymat_in,
                                    std::string);

    void zerofill_harmonic_dymat_renormalize(std::complex<double> ****, 
                                            unsigned int);
                                                
    void exec_scph_main(std::complex<double> ****);

    void exec_scph_relax_cell_coordinate_main(std::complex<double> ****,
                                          std::complex<double> ****);

    void exec_QHA_relax_main(std::complex<double> ****);

    void exec_perturbative_QHA();

    void read_C1_array(double * const);

    void read_elastic_constants(double * const * const, 
                                  double * const* const* const);

    void set_initial_q0(double * const);

    void set_initial_strain(double * const * const);

    void read_cell_opt_input(double &,
                             double &,
                             double &,
                             double &);

    void calculate_u0(const double * const, double * const);

    void calculate_force_in_real_space(const std::complex<double> * const, 
                                       double * const);

    void transform_to_real_space_at_Gamma(const std::complex<double> * const, 
                                          double * const);

    void update_cell_coordinate(double *q0,
                                double *u0,
                                double **u_tensor,
                                const std::complex<double> * const v1_array_atT,
                                const double * const* const omega2_array,
                                const std::complex<double> * const del_v0_strain_atT,
                                const double * const* const C2_array,
                                const std::complex<double> * const* const* const cmat_convert,
                                const std::vector<int> &harm_optical_modes,
                                double *delta_q0,
                                double *delta_u0,
                                double *delta_u_tensor,
                                double &du0,
                                double &du_tensor);

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
                                             const std::complex<double> * const* const* const);
 
    void compute_del_v1_strain_from_cubic(std::complex<double> **,
                                          const std::complex<double> * const* const* const);

    void compute_del_v1_strain_from_quartic(std::complex<double> **,
                                            const std::complex<double> * const* const* const);

    void compute_del_v2_strain_from_cubic(std::complex<double> ***,
                                          const std::complex<double> * const* const* const);
    
    void compute_del_v2_strain_from_quartic(std::complex<double> ***,
                                            const std::complex<double> * const* const* const);

    void compute_del_v3_strain_from_quartic(std::complex<double> ****,
                                            const std::complex<double> * const* const* const);

    void calculate_del_v2_strain_from_cubic_by_finite_difference(const std::complex<double> * const* const* const,
                                                                 std::complex<double> ***);

    void read_del_v2_strain_from_cubic_in_kspace(const std::complex<double> * const* const* const,
                                                 std::complex<double> ***);
    
    void calculate_del_v2_strain_from_cubic_by_finite_difference_from_allmode(const std::complex<double> * const* const* const,
                                                                              std::complex<double> ***);

    void calculate_del_v1_strain_from_harmonic_by_finite_difference_from_allmode(std::complex<double> **,
                                                                                 const std::complex<double> * const* const* const); 
                                                                
    void make_supercell_mapping_by_symmetry_operations(int **);

    void make_inverse_translation_mapping(int **);

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
                                               double *,
                                               double );

                                            
    void calculate_del_v1_strain_with_strain_displace(std::complex<double> **, 
                                                      double **,
                                                      std::complex<double> **,
                                                      std::complex<double> **,
                                                      std::complex<double> **,
                                                      std::complex<double> ***,
                                                      std::complex<double> ***,
                                                      std::complex<double> ****,
                                                      double *);

    void calculate_C2_array_with_strain_displace(double **, 
                                                 double **,
                                                 double **,
                                                 double **,
                                                 double ***,
                                                 std::complex<double> **,
                                                 std::complex<double> **,
                                                 std::complex<double> ***,
                                                 double *);

    void calculate_eta_tensor(double **, 
                              const double * const * const);

    void renormalize_v0_from_strain(double &, 
                                    double , 
                                    double **, 
                                    double *,
                                    double **, 
                                    double ***,
                                    double **,
                                    const double);

    void renormalize_v1_array_from_strain(std::complex<double> *, 
                                          const std::complex<double> * const ,
                                          const std::complex<double> * const* const, 
                                          const std::complex<double> * const* const, 
                                          const std::complex<double> * const* const, 
                                          const double * const* const);

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
                                     std::complex<double> *, 
                                     std::complex<double> ***, 
                                     std::complex<double> ***, 
                                     double ** , 
                                     const double);

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

    int get_xyz_string(const int, std::string&);

    // QHA
    void compute_cmat(std::complex<double> ***,
                      const std::complex<double> * const* const* const);

    void calc_v1_array_vib(std::complex<double> *, 
                           std::complex<double> ***,
                           const double);

    void calc_del_v0_strain_vib(std::complex<double> *, 
                                std::complex<double> ***, 
                                double);

    void compute_del_v_strain(std::complex<double> **,
                              std::complex<double> **,
                              std::complex<double> **,
                              std::complex<double> ***,
                              std::complex<double> ***,
                              std::complex<double> ****,
                              std::complex<double> ***,
                              int );

    void set_elastic_constants(double *C1_array,
                               double **C2_array,
                               double ***C3_array);


    void set_init_structure_atT(double *q0,
                                double **u_tensor,
                                double *u0,
                                bool &converged_prev,
                                const int set_init_str,
                                const int i_temp_loop);

    void write_resfile_header(std::ofstream &fout_q0,
                              std::ofstream &fout_u0,
                              std::ofstream &fout_v0,
                              std::ofstream &fout_u_tensor);

    void write_stepresfile_header_atT(std::ofstream &fout_step_q0,
                                      std::ofstream &fout_step_u0,
                                      std::ofstream &fout_step_u_tensor,
                                      const double temp);
    
    void write_stepresfile(double *q0,
                           double **u_tensor,
                           double *u0,
                           const int i_str_loop,
                           std::ofstream &fout_step_q0,
                           std::ofstream &fout_step_u0,
                           std::ofstream &fout_step_u_tensor);

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