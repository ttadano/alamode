/*
 integration.cpp

 Copyright (c) 2014-2021 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#include "mpi_common.h"
#include "integration.h"
#include "error.h"
#include "kpoint.h"
#include "mathfunctions.h"
#include "memory.h"
#include "system.h"
#include "phonon_dos.h"
#include "phonon_velocity.h"
#include "dynamical.h"
#include "anharmonic_core.h"
#include "fcs_phonon.h"
#include <iomanip>
#include <vector>
#include <algorithm>
#include <cmath>

using namespace PHON_NS;

Integration::Integration(PHON *phon) : Pointers(phon)
{
    set_default_variables();
}

Integration::~Integration()
{
    deallocate_variables();
};

void Integration::set_default_variables()
{
    ismear = -1; // for 3ph scattering
    ismear_4ph = 1; // for 4ph scattering
    epsilon = 10.0;
    epsilon_4ph = 10.0;
    adaptive_sigma = nullptr;
    adaptive_sigma4 = nullptr;
}

void Integration::deallocate_variables()
{
    delete adaptive_sigma;
    delete adaptive_sigma4;
}

void Integration::setup_integration()
{
    MPI_Bcast(&ismear, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&ismear_4ph, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (mympi->my_rank == 0) {
        std::cout << std::endl;
        if (ismear == -1) {
            std::cout << " ISMEAR = -1: Tetrahedron method will be used." << std::endl;
        } else if (ismear == 0) {
            std::cout << " ISMEAR = 0: Lorentzian broadening with epsilon = "
                      << std::fixed << std::setprecision(2) << epsilon << " (cm^-1)" << std::endl;
        } else if (ismear == 1) {
            std::cout << " ISMEAR = 1: Gaussian broadening with epsilon = "
                      << std::fixed << std::setprecision(2) << epsilon << " (cm^-1)" << std::endl;
        } else if (ismear == 2) {
            std::cout << " ISMEAR = 2: adaptive method will be used." << std::endl;
        } else {
            exit("setup_relaxation", "Invalid ismear");
        }
        std::cout << std::endl;

        if (anharmonic_core->quartic_mode) {

            std::cout << " -- QUARTIC > 0: 4ph scattering will calculated with flag ISMEAR_4PH --" << std::endl;
            if (ismear_4ph == -1) {
                std::cout << "Tetrahedron method is not implemented, changed to adaptive smearing !" << std::endl;
                ismear_4ph = 2;
            } else if (ismear_4ph == 0) {
                std::cout << " ISMEAR_4PH = 0: Lorentzian broadening with epsilon = "
                          << std::fixed << std::setprecision(2) << epsilon << " (cm^-1)" << std::endl;
            } else if (ismear_4ph == 1) {
                std::cout << " ISMEAR_4PH = 1: Gaussian broadening with epsilon = "
                          << std::fixed << std::setprecision(2) << epsilon << " (cm^-1)" << std::endl;
            } else if (ismear_4ph == 2) {
                std::cout << " ISMEAR_4PH = 2: Adaptive smearing method will be used " << std::endl;
            } else {
                exit("setup_relaxation", "Invalid ismear_4ph");
            }

            std::cout << std::endl;
        }
    }

    prepare_adaptivesmearing();

    epsilon *= time_ry / Hz_to_kayser; // Convert epsilon to a.u.
    epsilon_4ph *= time_ry / Hz_to_kayser; // Convert epsilon to a.u.
    MPI_Bcast(&epsilon, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&epsilon_4ph, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

void Integration::prepare_adaptivesmearing()
{
    if (ismear == 2) {
        adaptive_sigma = new AdaptiveSmearingSigma(dos->kmesh_dos->nk,
                                                   dynamical->neval);
        adaptive_sigma->setup(phonon_velocity,
                              dos->kmesh_dos,
                              system->lavec_p,
                              system->rlavec_p,
                              fcs_phonon->fc2_ext);
    }
}

void TetraNodes::setup()
{
    // This menber function creates node information of the tetrahedra.

    const auto nk23 = nk2 * nk3;

    for (int i = 0; i < nk1; ++i) {
        for (int j = 0; j < nk2; ++j) {
            for (int k = 0; k < nk3; ++k) {

                const auto ii = (i + 1) % nk1;
                const auto jj = (j + 1) % nk2;
                const auto kk = (k + 1) % nk3;

                const auto n1 = k + j * nk3 + i * nk23;
                const auto n2 = k + j * nk3 + ii * nk23;
                const auto n3 = k + jj * nk3 + i * nk23;
                const auto n4 = k + jj * nk3 + ii * nk23;
                const auto n5 = kk + j * nk3 + i * nk23;
                const auto n6 = kk + j * nk3 + ii * nk23;
                const auto n7 = kk + jj * nk3 + i * nk23;
                const auto n8 = kk + jj * nk3 + ii * nk23;

                auto m = 6 * (k + j * nk3 + i * nk23);

                tetras[m][0] = n1;
                tetras[m][1] = n2;
                tetras[m][2] = n3;
                tetras[m][3] = n6;

                ++m;

                tetras[m][0] = n2;
                tetras[m][1] = n3;
                tetras[m][2] = n4;
                tetras[m][3] = n6;

                ++m;

                tetras[m][0] = n1;
                tetras[m][1] = n3;
                tetras[m][2] = n5;
                tetras[m][3] = n6;

                ++m;

                tetras[m][0] = n3;
                tetras[m][1] = n4;
                tetras[m][2] = n6;
                tetras[m][3] = n8;

                ++m;

                tetras[m][0] = n3;
                tetras[m][1] = n6;
                tetras[m][2] = n7;
                tetras[m][3] = n8;

                ++m;

                tetras[m][0] = n3;
                tetras[m][1] = n5;
                tetras[m][2] = n6;
                tetras[m][3] = n7;
            }
        }
    }
}

unsigned int TetraNodes::get_ntetra() const
{
    return this->ntetra;
}

unsigned int **TetraNodes::get_tetras() const
{
    return this->tetras;
}

double Integration::do_tetrahedron(const double *energy,
                                   const double *f,
                                   const unsigned int ntetra,
                                   const unsigned int *const *tetras,
                                   const double e_ref)
{
    /*
    This function returns the summation of the given function f_{k}
    over the k-points which have the energy "e_ref" using the tetrahedron method.

    Ret(e_ref) = \int f(k) \delta(e_ref - energy(k))

    */

    auto ret = 0.0;
    double I1, I2, I3, I4;

    const auto frac3 = 1.0 / 3.0;
    double g;

    tetra_pair pair{};

    for (unsigned int i = 0; i < ntetra; ++i) {

        tetra_data.clear();

        for (unsigned int j = 0; j < 4; ++j) {
            const auto knum = tetras[i][j];
            pair.e = energy[knum];
            pair.f = f[knum];
            tetra_data.push_back(pair);
        }

        std::sort(tetra_data.begin(), tetra_data.end());

        const auto e1 = tetra_data[0].e;
        const auto e2 = tetra_data[1].e;
        const auto e3 = tetra_data[2].e;
        const auto e4 = tetra_data[3].e;

        const auto f1 = tetra_data[0].f;
        const auto f2 = tetra_data[1].f;
        const auto f3 = tetra_data[2].f;
        const auto f4 = tetra_data[3].f;

        if (e3 <= e_ref && e_ref < e4) {
            g = 3.0 * std::pow(e4 - e_ref, 2) / ((e4 - e1) * (e4 - e2) * (e4 - e3));

            I1 = frac3 * fij(e1, e4, e_ref);
            I2 = frac3 * fij(e2, e4, e_ref);
            I3 = frac3 * fij(e3, e4, e_ref);
            I4 = frac3 * (fij(e4, e1, e_ref) + fij(e4, e2, e_ref) + fij(e4, e3, e_ref));

            ret += g * (I1 * f1 + I2 * f2 + I3 * f3 + I4 * f4);

        } else if (e2 <= e_ref && e_ref < e3) {
            g = 3.0 * (e2 - e1 + 2.0 * (e_ref - e2) - (e4 + e3 - e2 - e1)
                                                      * std::pow(e_ref - e2, 2) / ((e3 - e2) * (e4 - e2))) /
                ((e3 - e1) * (e4 - e1));

            I1 = frac3 * fij(e1, e4, e_ref) * g + fij(e1, e3, e_ref) * fij(e3, e1, e_ref) * fij(e2, e3, e_ref) / (e4 -
                                                                                                                  e1);
            I2 = frac3 * fij(e2, e3, e_ref) * g + std::pow(fij(e2, e4, e_ref), 2) * fij(e3, e2, e_ref) / (e4 - e1);
            I3 = frac3 * fij(e3, e2, e_ref) * g + std::pow(fij(e3, e1, e_ref), 2) * fij(e2, e3, e_ref) / (e4 - e1);
            I4 = frac3 * fij(e4, e1, e_ref) * g + fij(e4, e2, e_ref) * fij(e2, e4, e_ref) * fij(e3, e2, e_ref) / (e4 -
                                                                                                                  e1);

            ret += I1 * f1 + I2 * f2 + I3 * f3 + I4 * f4;

        } else if (e1 <= e_ref && e_ref < e2) {
            g = 3.0 * std::pow(e_ref - e1, 2) / ((e2 - e1) * (e3 - e1) * (e4 - e1));

            I1 = frac3 * (fij(e1, e2, e_ref) + fij(e1, e3, e_ref) + fij(e1, e4, e_ref));
            I2 = frac3 * fij(e2, e1, e_ref);
            I3 = frac3 * fij(e3, e1, e_ref);
            I4 = frac3 * fij(e4, e1, e_ref);

            ret += g * (I1 * f1 + I2 * f2 + I3 * f3 + I4 * f4);

        }
    }

    return ret / static_cast<double>(ntetra);
}

void Integration::calc_weight_tetrahedron(const unsigned int nk_irreducible,
                                          const unsigned int *map_to_irreducible_k,
                                          const double *energy,
                                          const double e_ref,
                                          const unsigned int ntetra,
                                          const unsigned int *const *tetras,
                                          double *weight) const
{
    int i;

    double g;
    double e_tmp[4];
    int sort_arg[4], kindex[4];

    for (i = 0; i < nk_irreducible; ++i) weight[i] = 0.0;

    for (i = 0; i < ntetra; ++i) {

        for (int j = 0; j < 4; ++j) {
            e_tmp[j] = energy[tetras[i][j]];
            kindex[j] = map_to_irreducible_k[tetras[i][j]];
        }

        insertion_sort(e_tmp, sort_arg, 4);
        const auto e1 = e_tmp[0];
        const auto e2 = e_tmp[1];
        const auto e3 = e_tmp[2];
        const auto e4 = e_tmp[3];

        const auto k1 = kindex[sort_arg[0]];
        const auto k2 = kindex[sort_arg[1]];
        const auto k3 = kindex[sort_arg[2]];
        const auto k4 = kindex[sort_arg[3]];

        auto I1 = 0.0;
        auto I2 = 0.0;
        auto I3 = 0.0;
        auto I4 = 0.0;

        if (e3 <= e_ref && e_ref < e4) {
            g = std::pow(e4 - e_ref, 2) / ((e4 - e1) * (e4 - e2) * (e4 - e3));

            I1 = g * fij(e1, e4, e_ref);
            I2 = g * fij(e2, e4, e_ref);
            I3 = g * fij(e3, e4, e_ref);
            I4 = g * (fij(e4, e1, e_ref) + fij(e4, e2, e_ref) + fij(e4, e3, e_ref));

        } else if (e2 <= e_ref && e_ref < e3) {
            g = (e2 - e1 + 2.0 * (e_ref - e2) - (e4 + e3 - e2 - e1)
                                                * std::pow(e_ref - e2, 2) / ((e3 - e2) * (e4 - e2))) /
                ((e3 - e1) * (e4 - e1));

            I1 = g * fij(e1, e4, e_ref) + fij(e1, e3, e_ref) * fij(e3, e1, e_ref) * fij(e2, e3, e_ref) / (e4 - e1);
            I2 = g * fij(e2, e3, e_ref) + std::pow(fij(e2, e4, e_ref), 2) * fij(e3, e2, e_ref) / (e4 - e1);
            I3 = g * fij(e3, e2, e_ref) + std::pow(fij(e3, e1, e_ref), 2) * fij(e2, e3, e_ref) / (e4 - e1);
            I4 = g * fij(e4, e1, e_ref) + fij(e4, e2, e_ref) * fij(e2, e4, e_ref) * fij(e3, e2, e_ref) / (e4 - e1);

        } else if (e1 <= e_ref && e_ref < e2) {
            g = std::pow(e_ref - e1, 2) / ((e2 - e1) * (e3 - e1) * (e4 - e1));

            I1 = g * (fij(e1, e2, e_ref) + fij(e1, e3, e_ref) + fij(e1, e4, e_ref));
            I2 = g * fij(e2, e1, e_ref);
            I3 = g * fij(e3, e1, e_ref);
            I4 = g * fij(e4, e1, e_ref);

        }
        weight[k1] += I1;
        weight[k2] += I2;
        weight[k3] += I3;
        weight[k4] += I4;
    }
    auto factor = 1.0 / static_cast<double>(ntetra);
    for (i = 0; i < nk_irreducible; ++i) weight[i] *= factor;
}

void Integration::calc_weight_smearing(const unsigned int nk,
                                       const unsigned int nk_irreducible,
                                       const unsigned int *map_to_irreducible_k,
                                       const double *energy,
                                       const double e_ref,
                                       const int smearing_method,
                                       double *weight) const
{
    int i;

    const auto epsilon_kayser = this->epsilon * Hz_to_kayser / time_ry;
    const auto invnk = 1.0 / static_cast<double>(nk);

    for (i = 0; i < nk_irreducible; ++i) weight[i] = 0.0;

    if (smearing_method == 0) {
        for (i = 0; i < nk; ++i) {
            weight[map_to_irreducible_k[i]] += delta_lorentz(e_ref - energy[i], epsilon_kayser);
        }
    } else if (smearing_method == 1) {
        for (i = 0; i < nk; ++i) {
            weight[map_to_irreducible_k[i]] += delta_gauss(e_ref - energy[i], epsilon_kayser);
        }
    }

    for (i = 0; i < nk_irreducible; ++i) weight[i] *= invnk;
}

double Integration::fij(const double ei,
                        const double ej,
                        const double e) const
{
    return (e - ej) / (ei - ej);
}

void Integration::insertion_sort(double *a,
                                 int *ind,
                                 int n) const
{
    int i;

    for (i = 0; i < n; ++i) ind[i] = i;

    for (i = 1; i < n; ++i) {
        double tmp = a[i];

        int j = i;
        while (j > 0 && tmp < a[j - 1]) {
            a[j] = a[j - 1];
            ind[j] = ind[j - 1];
            --j;
        }
        a[j] = tmp;
        ind[j] = i;
    }
}

void AdaptiveSmearingSigma::setup(const PhononVelocity *phvel_class,
                                  const KpointMeshUniform *kmesh_in,
                                  const double lavec_p_in[3][3],
                                  const double rlavec_p_in[3][3],
                                  const std::vector<FcsClassExtent> &fc2_ext_in)
{
    phvel_class->get_phonon_group_velocity_mesh(*kmesh_in,
                                                lavec_p_in,
                                                fc2_ext_in,
                                                false,
                                                vel);

    for (auto u = 0; u < 3; u++) {
        for (auto a = 0; a < 3; a++) {
            dq[u][a] = rlavec_p_in[u][a] / static_cast<double>(kmesh_in->nk_i[u]);
        }
    }
}

void AdaptiveSmearingSigma::get_sigma(const unsigned int k1,
                                      const unsigned int s1,
                                      double &sigma_out)
{
    double parts;
    double tmp;

    parts = 0;

    for (auto u = 0; u < 3; ++u) {

        tmp = 0;
        for (auto a = 0; a < 3; ++a) {
            tmp += vel[k1][s1][a] * dq[u][a];
        }

        parts += std::pow(tmp, 2);
    }

    sigma_out = std::max(2.0e-5, std::sqrt(parts / 12)); // for (w1 - w2)
}

void AdaptiveSmearingSigma::get_sigma(const unsigned int k1,
                                      const unsigned int s1,
                                      const unsigned int k2,
                                      const unsigned int s2,
                                      double sigma_out[2])
{
    double parts[2];
    double tmp[2];
    int i;

    for (i = 0; i < 2; ++i) parts[i] = 0;

    for (auto u = 0; u < 3; ++u) {

        for (i = 0; i < 2; ++i) tmp[i] = 0;

        for (auto a = 0; a < 3; ++a) {
            tmp[0] += (vel[k1][s1][a] - vel[k2][s2][a]) * dq[u][a];
            tmp[1] += (vel[k1][s1][a] + vel[k2][s2][a]) * dq[u][a];
        }

        for (i = 0; i < 2; ++i) parts[i] += std::pow(tmp[i], 2);
    }

    sigma_out[0] = std::max(2.0e-5, std::sqrt((parts[0]) / 12)); // for (w1 - w2 - w3)
    sigma_out[1] = std::max(2.0e-5, std::sqrt((parts[1]) / 12)); // for (w1 + w3 - w3)
    // 2.0e-5 ry ~ 3 cm^-1
}

void AdaptiveSmearingSigma::get_sigma(const unsigned int k2,
                                      const unsigned int s2,
                                      const unsigned int k3,
                                      const unsigned int s3,
                                      const unsigned int k4,
                                      const unsigned int s4,
                                      double sigma_out[2])
{
    double vel_diff;
    double parts[3];
    double tmp[3];
    int i;
    for (i = 0; i < 3; ++i) parts[i] = 0;

    for (auto u = 0; u < 3; ++u) {

        for (i = 0; i < 3; ++i) tmp[i] = 0;

        for (auto a = 0; a < 3; ++a) {
            tmp[0] += (vel[k2][s2][a] - vel[k4][s4][a]) * dq[u][a];
            tmp[1] += (vel[k3][s3][a] - vel[k4][s4][a]) * dq[u][a];
            tmp[2] += (vel[k2][s2][a] + vel[k4][s4][a]) * dq[u][a];
        }

        for (i = 0; i < 3; ++i) parts[i] += std::pow(tmp[i], 2);
    }

    sigma_out[0] = std::max(2.0e-5, std::sqrt((parts[0] + parts[1]) / 12));  // for delta(w1 - w2 - w3 - w4)
    sigma_out[1] = std::max(2.0e-5,
                            std::sqrt((parts[2] + parts[1]) /
                                      12));  // for delta(w1 + w2 - w3 - w4) and (w1 - w2 + w3 + w4)
}
