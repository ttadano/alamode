/*
 phonon_velocity.h

 Copyright (c) 2014 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#pragma once

#include "pointers.h"
#include "fcs_phonon.h"
#include "kpoint.h"
#include <vector>
#include <complex>

namespace PHON_NS {
class PhononVelocity : protected Pointers {
 public:
    PhononVelocity(class PHON *);

    ~PhononVelocity();

    void setup_velocity();

    void phonon_vel_k(const double *,
                      double **) const;

    void phonon_vel_k2(const double *,
                       const double *,
                       std::complex<double> **,
                       double **) const;

    void get_phonon_group_velocity_mesh(const KpointMeshUniform &kmesh_in,
                                        const double lavec_p[3][3],
                                        const std::vector<FcsClassExtent> &fc2_ext_in,
                                        const bool irreducible_only,
                                        double ***phvel3_out) const;

    void get_phonon_group_velocity_mesh_mpi(const KpointMeshUniform &kmesh_in,
                                            const double lavec_p[3][3],
                                            const std::vector<FcsClassExtent> &fc2_ext_in,
                                            double ***phvel3_out) const;

    void calc_phonon_velmat_mesh(std::complex<double> ****velmat_out) const;

    void get_phonon_group_velocity_bandstructure(const KpointBandStructure *kpoint_bs_in,
                                                 const double lavec_p[3][3],
                                                 const double rlavec_p[3][3],
                                                 const std::vector<FcsClassExtent> &fc2_ext_in,
                                                 double **phvel_out) const;

    void velocity_matrix_analytic(const double *xk_in,
                                  const std::vector<FcsClassExtent> &fc2_in,
                                  const double *omega_in,
                                  std::complex<double> **evec_in,
                                  std::complex<double> ***velmat_out) const;

    bool print_velocity;
    std::complex<double> ***velmat;

 private:

    double **xshift_s;

    double diff(const double *,
                unsigned int,
                double) const;

    void set_default_variables();

    void deallocate_variables();

    void calc_derivative_dynmat_k(const double *,
                                  const std::vector<FcsClassExtent> &,
                                  std::complex<double> ***) const;

    void diagonalize_hermite_mat(int,
                                 std::complex<double> **,
                                 double *) const;

    bool print_velocity_xyz;
};
}
