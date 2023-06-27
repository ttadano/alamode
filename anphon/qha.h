/*
 qha.h

 Copyright (c) 2022 Ryota Masuki, Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory
 or http://opensource.org/licenses/mit-license.php for information.
*/

#pragma once
#include "pointers.h"
#include "anharmonic_core.h"
#include "scph.h"
#include <complex>
#include "kpoint.h"

namespace PHON_NS {
class Qha : protected Pointers {
public:

    Qha(class PHON *phon);

    ~Qha();

    unsigned int kmesh_qha[3];
    unsigned int kmesh_interpolate[3];

    // optimization scheme used in QHA
    int qha_scheme;

    bool restart_qha;

    void exec_qha_optimization();


private:

    void set_default_variables();

    void deallocate_variables();

    void exec_QHA_relax_main(std::complex<double> ****,
                             std::complex<double> ****);

    void exec_perturbative_QHA(std::complex<double> ****,
                               std::complex<double> ****);

    void calc_del_v0_del_umn_vib(std::complex<double> *,
                                 std::complex<double> ***,
                                 double);


    void calculate_del_v1_del_umn_renorm(std::complex<double> **,
                                         double **,
                                         std::complex<double> **,
                                         std::complex<double> **,
                                         std::complex<double> **,
                                         std::complex<double> ***,
                                         std::complex<double> ***,
                                         std::complex<double> ****,
                                         double *);

    void calculate_C2_array_renorm(double **,
                                   double **,
                                   double **,
                                   double **,
                                   double ***,
                                   std::complex<double> **,
                                   std::complex<double> **,
                                   std::complex<double> ***,
                                   double *);

    void calculate_C2_array_ZSISA(double **,
                                  double **,
                                  std::complex<double> **,
                                  double **);

    void compute_ZSISA_stress(double **,
                              std::complex<double> *,
                              std::complex<double> ***,
                              double **,
                              std::complex<double> *,
                              std::complex<double> **,
                              std::complex<double> *,
                              std::vector<int> &);

    void compute_vZSISA_stress(std::complex<double> *,
                               double **,
                               std::complex<double> *,
                               std::complex<double> *,
                               double **);

    KpointMeshUniform *kmesh_coarse = nullptr;
    KpointMeshUniform *kmesh_dense = nullptr;

    // Information of harmonic dynamical matrix
    double **omega2_harmonic;
    std::complex<double> ***evec_harmonic;
    MinimumDistList ***mindist_list_qha;

};
}



