/*
 relaxation.h

 Copyright (c) 2022 Ryota Masuki, Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory
 or http://opensource.org/licenses/mit-license.php for information.
*/

#pragma once
#include "pointers.h"
#include <complex>
#include "kpoint.h"
#include "fcs_phonon.h"

namespace PHON_NS {
class Relaxation : protected Pointers {
public:

    Relaxation(class PHON *phon);
    ~Relaxation();

    int relax_str;

    // initial strain and displacement
    double **init_u_tensor = nullptr;
    std::vector<double> init_u0;

    // zero-th order term of the potential energy surface
    double *V0 = nullptr;

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

    int renorm_3to2nd;
    int renorm_2to1st;
    int renorm_34to1st;
    std::string strain_IFC_dir;

    void setup_relaxation();

    void compute_del_v_strain(const unsigned int nk,
                              std::complex<double> **,
                              std::complex<double> **,
                              std::complex<double> **,
                              std::complex<double> ***,
                              std::complex<double> ***,
                              std::complex<double> ****,
                              std::complex<double> ***,
                              int);

    void load_V0_from_file();
    void store_V0_to_file();

    void set_init_structure_atT(double *q0,
                                double **u_tensor,
                                double *u0,
                                bool &converged_prev,
                                int &str_diverged,
                                const int i_temp_loop,
                                double **omega2_harmonic,
                                std::complex<double> ***evec_harmonic);


    void set_elastic_constants(double *C1_array,
                               double **C2_array,
                               double ***C3_array);

    void renormalize_v0_from_umn(double &,
                                 double,
                                 double **,
                                 double *,
                                 double **,
                                 double ***,
                                 double **,
                                 const double);

    void renormalize_v1_from_umn(std::complex<double> *,
                                 const std::complex<double> *const,
                                 const std::complex<double> *const *const,
                                 const std::complex<double> *const *const,
                                 const std::complex<double> *const *const,
                                 const double *const *const);

    void renormalize_v2_from_umn(std::complex<double> **,
                                 std::complex<double> ***,
                                 std::complex<double> ***,
                                 double **);

    void renormalize_v3_from_umn(std::complex<double> ***,
                                 std::complex<double> ***,
                                 std::complex<double> ****,
                                 double **);

    void renormalize_v1_from_q0(std::complex<double> *,
                                std::complex<double> *,
                                std::complex<double> **,
                                std::complex<double> ***,
                                std::complex<double> ***,
                                double *);

    void renormalize_v2_from_q0(std::complex<double> **,
                                std::complex<double> **,
                                std::complex<double> ***,
                                std::complex<double> ***,
                                double *);

    void renormalize_v3_from_q0(std::complex<double> ***,
                                std::complex<double> ***,
                                std::complex<double> ***,
                                double *);

    void renormalize_v0_from_q0(double &,
                                double,
                                std::complex<double> *,
                                std::complex<double> **,
                                std::complex<double> ***,
                                std::complex<double> ***,
                                double *);

    void calculate_u0(const double *const q0, double *const u0,
                      double **omega2_harmonic,
                      std::complex<double> ***evec_harmonic);

    void update_cell_coordinate(double *,
                                double *,
                                double **,
                                const std::complex<double> *const,
                                const double *const *const,
                                const std::complex<double> *const,
                                const double *const *const,
                                const std::complex<double> *const *const *const,
                                const std::vector<int> &,
                                double *,
                                double *,
                                double *,
                                double &,
                                double &,
                                double **omega2_harmonic,
                                std::complex<double> ***evec_harmonic);

    void check_str_divergence(int &diverged,
                              const double *const q0,
                              const double *const u0,
                              const double *const *const u_tensor);



private:

    void set_default_variables();
    void deallocate_variables();

    void read_C1_array(double *const);

    void read_elastic_constants(double *const *const,
                                double *const *const *const);

    void set_initial_q0(double *const q0, std::complex<double> ***evec_harmonic);


    void set_initial_strain(double *const *const);


    void compute_del_v1_del_umn(std::complex<double> **,
                                const std::complex<double> *const *const *const);

    void compute_del2_v1_del_umn2(std::complex<double> **,
                                  const std::complex<double> *const *const *const);

    void compute_del3_v1_del_umn3(std::complex<double> **,
                                  const std::complex<double> *const *const *const);



    void compute_del_v2_del_umn(std::complex<double> ***,
                                const std::complex<double> *const *const *const,
                                const unsigned int nk,
                                const unsigned int nk_interpolate,
                                double **xk_in);

    void compute_del2_v2_del_umn2(std::complex<double> ***,
                                  const std::complex<double> *const *const *const,
                                  const unsigned int nk,
                                  const unsigned int nk_interpolate,
                                  double **xk_in);

    void compute_del_v3_del_umn(std::complex<double> ****,
                                const std::complex<double> *const *const *const,
                                const unsigned int nk,
                                const unsigned int nk_interpolate);

    void calculate_delv2_delumn_finite_difference(const std::complex<double> *const *const *const,
                                                  std::complex<double> ***,
                                                  const unsigned int nk,
                                                  const unsigned int nk_interpolate,
                                                  const int kmesh_interpolate[3]);

    void read_del_v2_del_umn_in_kspace(const std::complex<double> *const *const *const,
                                       std::complex<double> ***,
                                       const unsigned int nk,
                                       const unsigned int nk_interpolate);

    void calculate_delv1_delumn_finite_difference(std::complex<double> **,
                                                  const std::complex<double> *const *const *const);

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


    void make_supercell_mapping_by_symmetry_operations(int **);

    void make_inverse_translation_mapping(int **);











};
}



