/*
 phonon_dos.h

 Copyright (c) 2014 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#pragma once

#include "pointers.h"
#include "dynamical.h"
#include "integration.h"
#include "kpoint.h"
#include <vector>
#include <complex>

namespace PHON_NS {
class Dos : protected Pointers {
public:
    Dos(class PHON *);

    ~Dos();

    void setup();

    void calc_dos_all();

    bool flag_dos;
    bool compute_dos;
    bool projected_dos, two_phonon_dos;
    int scattering_phase_space;

    int n_energy;
    double emin, emax, delta_e;
    double *energy_dos;
    double *dos_phonon;
    double **pdos_phonon;
    double ***dos2_phonon;
    double total_sps3, ***sps3_mode;
    double ****sps3_with_bose;

    TetraNodes *tetra_nodes_dos;
    KpointMeshUniform *kmesh_dos;
    DymatEigenValue *dymat_dos;

    void calc_dos_from_given_frequency(const KpointMeshUniform *kmesh_in,
                                       const double *const *eval_in,
                                       const unsigned int ntetra_in,
                                       const unsigned int *const *tetras_in,
                                       double *dos_out) const;

    void set_dos_energy_grid();

private:
    void set_default_variables();

    void deallocate_variables();

    void calc_dos(const unsigned int nk,
                  const unsigned int nk_irreducible,
                  const unsigned int *map_k,
                  const double *const *eval,
                  const unsigned int n,
                  const double *energy,
                  const unsigned int neval,
                  const int smearing_method,
                  const unsigned int ntetra,
                  const unsigned int *const *tetras,
                  double *ret) const;

    void calc_atom_projected_dos(const unsigned int nk,
                                 double *const *eval,
                                 const unsigned int n,
                                 const double *energy,
                                 double **ret,
                                 const unsigned int neval,
                                 const unsigned int natmin,
                                 const int smearing_method,
                                 std::complex<double> ***evec) const;

    void calc_two_phonon_dos(double *const *eval,
                             const unsigned int n,
                             const double *energy,
                             const int smearing_method,
                             double ***ret) const;

    void calc_total_scattering_phase_space(double *const *eval_in,
                                           const int smearing_method,
                                           double ***ret_mode,
                                           double &ret) const;

    void calc_scattering_phase_space_with_Bose(const double *const *eval_in,
                                               const int smearing_method,
                                               double ****ret) const;

    void calc_scattering_phase_space_with_Bose_mode(const unsigned int nk,
                                                    const unsigned int ns,
                                                    const unsigned int N,
                                                    const double omega,
                                                    const double *const *eval,
                                                    const double *temperature,
                                                    const unsigned int *k_pair,
                                                    const int smearing_method,
                                                    double **ret) const;
};
}
