/*
 integration.cpp

 Copyright (c) 2014 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#include "mpi_common.h"
#include "integration.h"
#include "kpoint.h"
#include "memory.h"
#include "system.h"
#include "error.h"
#include <iomanip>
#include <vector>
#include <algorithm>
#include <cmath>
#include "mathfunctions.h"

using namespace PHON_NS;

Integration::Integration(PHON *phon): Pointers(phon)
{
}

Integration::~Integration()
{
};

void Integration::setup_integration()
{
    MPI_Bcast(&ismear, 1, MPI_INT, 0, MPI_COMM_WORLD);

    unsigned int nk = kpoint->nk;
    unsigned int nkx = kpoint->nkx;
    unsigned int nky = kpoint->nky;
    unsigned int nkz = kpoint->nkz;

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
        } else {
            error->exit("setup_relaxation", "Invalid ksum_mode");
        }
        std::cout << std::endl;
    }

    if (ismear == -1) {
        ntetra = 6 * nk;
        memory->allocate(tetras, ntetra, 4);
        prepare_tetrahedron(nkx, nky, nkz);
    }

    epsilon *= time_ry / Hz_to_kayser;
    MPI_Bcast(&epsilon, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

void Integration::finish_integration()
{
    if (ismear == -1) {
        memory->deallocate(tetras);
    }
}

void Integration::prepare_tetrahedron(const int nk1, const int nk2, const int nk3)
{
    int i, j, k, ii, jj, kk;
    int n1, n2, n3, n4, n5, n6, n7, n8;
    int m;

    int nk23 = nk2 * nk3;

    for (i = 0; i < nk1; ++i) {
        for (j = 0; j < nk2; ++j) {
            for (k = 0; k < nk3; ++k) {

                ii = (i + 1) % nk1;
                jj = (j + 1) % nk2;
                kk = (k + 1) % nk3;

                n1 = k + j * nk3 + i * nk23;
                n2 = k + j * nk3 + ii * nk23;
                n3 = k + jj * nk3 + i * nk23;
                n4 = k + jj * nk3 + ii * nk23;
                n5 = kk + j * nk3 + i * nk23;
                n6 = kk + j * nk3 + ii * nk23;
                n7 = kk + jj * nk3 + i * nk23;
                n8 = kk + jj * nk3 + ii * nk23;

                m = 6 * (k + j * nk3 + i * nk23);

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

double Integration::do_tetrahedron(double *energy, double *f, const double e_ref)
{
    /*
    This function returns the summation of the given function f_{k}
    over the k-points which have the energy "e_ref" using the tetrahedron method.

    Ret(e_ref) = \int f(k) \delta(e_ref - energy(k))

    */

    unsigned int i, j;
    int knum;

    double ret = 0.0;
    double e1, e2, e3, e4;
    double I1, I2, I3, I4;
    double f1, f2, f3, f4;

    double frac3 = 1.0 / 3.0;
    double g, vol, vol_tot;

    tetra_pair pair;

    vol_tot = 0.0;

    for (i = 0; i < ntetra; ++i) {

        tetra_data.clear();

        for (j = 0; j < 4; ++j) {
            knum = tetras[i][j];
            pair.e = energy[knum];
            pair.f = f[knum];
            tetra_data.push_back(pair);
        }

        std::sort(tetra_data.begin(), tetra_data.end());

        e1 = tetra_data[0].e;
        e2 = tetra_data[1].e;
        e3 = tetra_data[2].e;
        e4 = tetra_data[3].e;

        f1 = tetra_data[0].f;
        f2 = tetra_data[1].f;
        f3 = tetra_data[2].f;
        f4 = tetra_data[3].f;

        vol = volume(tetras[i]);
        vol_tot += vol;

        if (e3 <= e_ref && e_ref < e4) {
            g = 3.0 * std::pow(e4 - e_ref, 2) / ((e4 - e1) * (e4 - e2) * (e4 - e3));

            I1 = frac3 * fij(e1, e4, e_ref);
            I2 = frac3 * fij(e2, e4, e_ref);
            I3 = frac3 * fij(e3, e4, e_ref);
            I4 = frac3 * (fij(e4, e1, e_ref) + fij(e4, e2, e_ref) + fij(e4, e3, e_ref));

            ret += vol * g * (I1 * f1 + I2 * f2 + I3 * f3 + I4 * f4);

        } else if (e2 <= e_ref && e_ref < e3) {
            g = 3.0 * ((e2 - e1) + 2.0 * (e_ref - e2) - (e4 + e3 - e2 - e1)
                * std::pow((e_ref - e2), 2) / ((e3 - e2) * (e4 - e2))) / ((e3 - e1) * (e4 - e1));

            I1 = frac3 * fij(e1, e4, e_ref) * g + fij(e1, e3, e_ref) * fij(e3, e1, e_ref) * fij(e2, e3, e_ref) / (e4 - e1);
            I2 = frac3 * fij(e2, e3, e_ref) * g + std::pow(fij(e2, e4, e_ref), 2) * fij(e3, e2, e_ref) / (e4 - e1);
            I3 = frac3 * fij(e3, e2, e_ref) * g + std::pow(fij(e3, e1, e_ref), 2) * fij(e2, e3, e_ref) / (e4 - e1);
            I4 = frac3 * fij(e4, e1, e_ref) * g + fij(e4, e2, e_ref) * fij(e2, e4, e_ref) * fij(e3, e2, e_ref) / (e4 - e1);

            ret += vol * (I1 * f1 + I2 * f2 + I3 * f3 + I4 * f4);

        } else if (e1 <= e_ref && e_ref < e2) {
            g = 3.0 * std::pow(e_ref - e1, 2) / ((e2 - e1) * (e3 - e1) * (e4 - e1));

            I1 = frac3 * (fij(e1, e2, e_ref) + fij(e1, e3, e_ref) + fij(e1, e4, e_ref));
            I2 = frac3 * fij(e2, e1, e_ref);
            I3 = frac3 * fij(e3, e1, e_ref);
            I4 = frac3 * fij(e4, e1, e_ref);

            ret += vol * g * (I1 * f1 + I2 * f2 + I3 * f3 + I4 * f4);

        }

    }

    return ret / vol_tot;
}

double Integration::dos_integration(double *energy, const double e_ref)
{
    double dos_ret = 0.0;

    unsigned int i, j;
    std::vector<double> e_tetra;
    double e1, e2, e3, e4;
    double vol, vol_tot;

    vol_tot = 0.0;

    for (i = 0; i < ntetra; ++i) {
        e_tetra.clear();
        for (j = 0; j < 4; ++j) {
            e_tetra.push_back(energy[tetras[i][j]]);
        }
        std::sort(e_tetra.begin(), e_tetra.end());

        e1 = e_tetra[0];
        e2 = e_tetra[1];
        e3 = e_tetra[2];
        e4 = e_tetra[3];

        vol = volume(tetras[i]);
        vol_tot += vol;

        if (e3 <= e_ref && e_ref < e4) {
            dos_ret += vol * (3.0 * std::pow((e4 - e_ref), 2) / ((e4 - e1) * (e4 - e2) * (e4 - e3)));
        } else if (e2 <= e_ref && e_ref < e3) {
            dos_ret += vol * 3.0 * ((e2 - e1) + 2.0 * (e_ref - e2) - (e4 + e3 - e2 - e1)
                * std::pow((e_ref - e2), 2) / ((e3 - e2) * (e4 - e2))) / ((e3 - e1) * (e4 - e1));
        } else if (e1 <= e_ref && e_ref < e2) {
            dos_ret += vol * 3.0 * std::pow((e_ref - e1), 2) / ((e2 - e1) * (e3 - e1) * (e4 - e1));
        }
    }

    return dos_ret / vol_tot;
}


void Integration::calc_weight_tetrahedron(const int nk_irreducible,
                                          int *map_to_irreducible_k,
                                          double *weight,
                                          double *energy,
                                          const double e_ref)
{
    int i, j;
    double vol;
    double vol_tot;

    double g;
    double I1, I2, I3, I4;
    int k1, k2, k3, k4;
    double e1, e2, e3, e4;

    vol_tot = 0.0;
    // std::vector<TetraWithKnum> tetra_data;
    // TetraWithKnum pair;
    double e_tmp[4];
    int sort_arg[4], kindex[4];

    for (i = 0; i < nk_irreducible; ++i) weight[i] = 0.0;

    for (i = 0; i < ntetra; ++i) {

        for (j = 0; j < 4; ++j) {
            e_tmp[j] = energy[tetras[i][j]];
            kindex[j] = map_to_irreducible_k[tetras[i][j]];
        }

        insertion_sort(e_tmp, sort_arg, 4);
        e1 = e_tmp[0];
        e2 = e_tmp[1];
        e3 = e_tmp[2];
        e4 = e_tmp[3];

        k1 = kindex[sort_arg[0]];
        k2 = kindex[sort_arg[1]];
        k3 = kindex[sort_arg[2]];
        k4 = kindex[sort_arg[3]];


        /*
        tetra_data.clear();

        for (j = 0; j < 4; ++j) {
            pair.e = energy[tetras[i][j]];
            pair.knum = map_to_irreducible_k[tetras[i][j]];
            tetra_data.push_back(pair);
        }

        std::sort(tetra_data.begin(), tetra_data.end());

        e1 = tetra_data[0].e;
        e2 = tetra_data[1].e;
        e3 = tetra_data[2].e;
        e4 = tetra_data[3].e;

        k1 = tetra_data[0].knum;
        k2 = tetra_data[1].knum;
        k3 = tetra_data[2].knum;
        k4 = tetra_data[3].knum;
        */
        vol = volume(tetras[i]);
        vol_tot += vol;

        I1 = 0.0;
        I2 = 0.0;
        I3 = 0.0;
        I4 = 0.0;

        if (e3 <= e_ref && e_ref < e4) {
            // g = 3.0 * std::pow(e4 - e_ref, 2) / ((e4 - e1) * (e4 - e2) * (e4 - e3));
            g = std::pow(e4 - e_ref, 2) / ((e4 - e1) * (e4 - e2) * (e4 - e3));

            I1 = g * fij(e1, e4, e_ref);
            I2 = g * fij(e2, e4, e_ref);
            I3 = g * fij(e3, e4, e_ref);
            I4 = g * (fij(e4, e1, e_ref) + fij(e4, e2, e_ref) + fij(e4, e3, e_ref));

        } else if (e2 <= e_ref && e_ref < e3) {
            //  g = 3.0 * ((e2 - e1) + 2.0 * (e_ref - e2) - (e4 + e3 - e2 - e1)
            //      * std::pow((e_ref - e2), 2) / ((e3 - e2) * (e4 - e2))) / ((e3 - e1) * (e4 - e1));
            g = ((e2 - e1) + 2.0 * (e_ref - e2) - (e4 + e3 - e2 - e1)
                * std::pow((e_ref - e2), 2) / ((e3 - e2) * (e4 - e2))) / ((e3 - e1) * (e4 - e1));

            I1 = g * fij(e1, e4, e_ref) + fij(e1, e3, e_ref) * fij(e3, e1, e_ref) * fij(e2, e3, e_ref) / (e4 - e1);
            I2 = g * fij(e2, e3, e_ref) + std::pow(fij(e2, e4, e_ref), 2) * fij(e3, e2, e_ref) / (e4 - e1);
            I3 = g * fij(e3, e2, e_ref) + std::pow(fij(e3, e1, e_ref), 2) * fij(e2, e3, e_ref) / (e4 - e1);
            I4 = g * fij(e4, e1, e_ref) + fij(e4, e2, e_ref) * fij(e2, e4, e_ref) * fij(e3, e2, e_ref) / (e4 - e1);

        } else if (e1 <= e_ref && e_ref < e2) {
            //  g = 3.0 * std::pow(e_ref - e1, 2) / ((e2 - e1) * (e3 - e1) * (e4 - e1));
            g = std::pow(e_ref - e1, 2) / ((e2 - e1) * (e3 - e1) * (e4 - e1));

            I1 = g * (fij(e1, e2, e_ref) + fij(e1, e3, e_ref) + fij(e1, e4, e_ref));
            I2 = g * fij(e2, e1, e_ref);
            I3 = g * fij(e3, e1, e_ref);
            I4 = g * fij(e4, e1, e_ref);

        }
        weight[k1] += vol * I1;
        weight[k2] += vol * I2;
        weight[k3] += vol * I3;
        weight[k4] += vol * I4;
    }

    for (i = 0; i < nk_irreducible; ++i) weight[i] /= vol_tot;
}

void PHON_NS::Integration::calc_weight_smearing(const std::vector<std::vector<KpointList>> &kpinfo,
                                                double *weight,
                                                double *energy,
                                                const double e_ref,
                                                const int smearing_method)
{
    unsigned int i;
    unsigned int knum;

    double epsilon = this->epsilon * Hz_to_kayser / time_ry;

    if (smearing_method == 0) {
        for (i = 0; i < kpinfo.size(); ++i) {
            knum = kpinfo[i][0].knum;
            weight[i] = kpoint->weight_k[i] * delta_lorentz(e_ref - energy[knum], epsilon);
        }
    } else if (smearing_method == 1) {
        for (i = 0; i < kpinfo.size(); ++i) {
            knum = kpinfo[i][0].knum;
            weight[i] = kpoint->weight_k[i] * delta_gauss(e_ref - energy[knum], epsilon);
        }
    }
}

void PHON_NS::Integration::calc_weight_smearing(const int nk,
                                                const int nk_irreducible,
                                                int *map_to_irreducible_k,
                                                double *weight,
                                                double *energy,
                                                const double e_ref,
                                                const int smearing_method)
{
    int i;

    double epsilon = this->epsilon * Hz_to_kayser / time_ry;
    double invnk = 1.0 / static_cast<double>(nk);

    for (i = 0; i < nk_irreducible; ++i) weight[i] = 0.0;

    if (smearing_method == 0) {
        for (i = 0; i < nk; ++i) {
            weight[map_to_irreducible_k[i]] += delta_lorentz(e_ref - energy[i], epsilon);
        }
    } else if (smearing_method == 1) {
        for (i = 0; i < nk; ++i) {
            weight[map_to_irreducible_k[i]] += delta_gauss(e_ref - energy[i], epsilon);
        }
    }

    for (i = 0; i < nk_irreducible; ++i) weight[i] *= invnk;
}


double Integration::volume(int *klist)
{
    int i;
    double k1[3], k2[3], k3[3];
    double vol;

    for (i = 0; i < 3; ++i) {
        k1[i] = refold(kpoint->xk[klist[1]][i] - kpoint->xk[klist[0]][i]);
        k2[i] = refold(kpoint->xk[klist[2]][i] - kpoint->xk[klist[0]][i]);
        k3[i] = refold(kpoint->xk[klist[3]][i] - kpoint->xk[klist[0]][i]);
    }

    rotvec(k1, k1, system->rlavec_p, 'T');
    rotvec(k2, k2, system->rlavec_p, 'T');
    rotvec(k3, k3, system->rlavec_p, 'T');

    vol = std::abs(k1[0] * (k2[1] * k3[2] - k2[2] * k3[1])
        + k1[1] * (k2[2] * k3[0] - k2[0] * k3[2])
        + k1[2] * (k2[0] * k3[1] - k2[1] * k3[0]));

    return vol;
}

double Integration::fij(const double ei, const double ej, const double e)
{
    return (e - ej) / (ei - ej);
}

double Integration::refold(double x)
{
    if (std::abs(x) > 0.5) {
        if (x < 0.0) {
            return x + 1.0;
        } else {
            return x - 1.0;
        }
    } else {
        return x;
    }
}

void Integration::insertion_sort(double *a, int *ind, int n)
{
    int i, j;
    double tmp;

    for (i = 0; i < n; ++i) ind[i] = i;

    for (i = 1; i < n; ++i) {
        tmp = a[i];

        j = i;
        while (j > 0 && tmp < a[j - 1]) {
            a[j] = a[j - 1];
            ind[j] = ind[j - 1];
            --j;
        }
        a[j] = tmp;
        ind[j] = i;
    }
}
