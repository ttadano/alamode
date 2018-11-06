/*
phonon_velocity.cpp

Copyright (c) 2014, 2015, 2016 Terumasa Tadano

This file is distributed under the terms of the MIT license.
Please see the file 'LICENCE.txt' in the root directory 
or http://opensource.org/licenses/mit-license.php for information.
*/

#include "mpi_common.h"
#include "phonon_velocity.h"
#include "constants.h"
#include "dynamical.h"
#include "error.h"
#include "fcs_phonon.h"
#include "kpoint.h"
#include "mathfunctions.h"
#include "memory.h"
#include "system.h"
#include <complex>

using namespace PHON_NS;

Phonon_velocity::Phonon_velocity(PHON *phon): Pointers(phon)
{
    set_default_variables();
}

Phonon_velocity::~Phonon_velocity()
{
    deallocate_variables();
}

void Phonon_velocity::set_default_variables()
{
    print_velocity = false;
    phvel = nullptr;
    phvel_xyz = nullptr;

    memory->allocate(xshift_s, 27, 3);

    for (int i = 0; i < 3; ++i) xshift_s[0][i] = 0.0;
    int icell = 0;
    for (int ix = -1; ix <= 1; ++ix) {
        for (int iy = -1; iy <= 1; ++iy) {
            for (int iz = -1; iz <= 1; ++iz) {
                if (ix == 0 && iy == 0 && iz == 0) continue;

                ++icell;

                xshift_s[icell][0] = static_cast<double>(ix);
                xshift_s[icell][1] = static_cast<double>(iy);
                xshift_s[icell][2] = static_cast<double>(iz);
            }
        }
    }
}

void Phonon_velocity::deallocate_variables()
{
    if (phvel) {
        memory->deallocate(phvel);
    }
    if (phvel_xyz) {
        memory->deallocate(phvel_xyz);
    }
    if (xshift_s) {
        memory->deallocate(xshift_s);
    }
}


void Phonon_velocity::calc_group_velocity(const int kpmode)
{
    MPI_Bcast(&print_velocity, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD);
    print_velocity_xyz = false;

    if (print_velocity) {

        unsigned int nk = kpoint->nk;
        unsigned int ns = dynamical->neval;

        memory->allocate(phvel, nk, ns);

        if (kpmode == 1) {

            calc_phonon_vel_band(phvel);

        } else if (kpmode == 2) {
            print_velocity_xyz = true;
            memory->allocate(phvel_xyz, nk, ns, 3);
            calc_phonon_vel_mesh(phvel, phvel_xyz);
        }

    }
}

void Phonon_velocity::calc_phonon_vel_band(double **phvel_out)
{
    unsigned int i;
    unsigned int idiff;
    unsigned int nk = kpoint->nk;
    unsigned int n = dynamical->neval;
    double **xk_shift;
    double *xk_tmp;
    double **omega_shift, *omega_tmp;

    double h = 1.0e-4;

    std::complex<double> **evec_tmp;

    memory->allocate(evec_tmp, 1, 1);

    if (mympi->my_rank == 0) {
        std::cout << " Calculating group velocities of phonon along given k path ... ";
    }

    unsigned int ndiff = 2;
    memory->allocate(xk_shift, ndiff, 3);
    memory->allocate(omega_shift, ndiff, n);
    memory->allocate(omega_tmp, ndiff);

    memory->allocate(xk_tmp, 3);


    for (unsigned int ik = 0; ik < nk; ++ik) {

        // Represent the given kpoint in Cartesian coordinate
        rotvec(xk_tmp, kpoint->xk[ik], system->rlavec_p, 'T');

        if (ndiff == 2) {
            // central difference
            // f'(x) =~ f(x+h)-f(x-h)/2h

            for (i = 0; i < 3; ++i) {
                xk_shift[0][i] = xk_tmp[i] - h * kpoint->kvec_na[ik][i];
                xk_shift[1][i] = xk_tmp[i] + h * kpoint->kvec_na[ik][i];
            }

        } else {
            error->exit("calc_phonon_vel_band",
                        "ndiff > 2 is not supported yet.");
        }

        for (idiff = 0; idiff < ndiff; ++idiff) {

            // Move back to fractional basis

            rotvec(xk_shift[idiff], xk_shift[idiff], system->lavec_p, 'T');
            for (i = 0; i < 3; ++i) xk_shift[idiff][i] /= 2.0 * pi;

            dynamical->eval_k(xk_shift[idiff], kpoint->kvec_na[ik],
                              fcs_phonon->fc2_ext, omega_shift[idiff],
                              evec_tmp, false);

        }

        for (i = 0; i < n; ++i) {
            for (idiff = 0; idiff < ndiff; ++idiff) {
                omega_tmp[idiff] = dynamical->freq(omega_shift[idiff][i]);
            }
            phvel_out[ik][i] = diff(omega_tmp, ndiff, h);
        }
    }
    memory->deallocate(omega_tmp);
    memory->deallocate(omega_shift);
    memory->deallocate(xk_shift);
    memory->deallocate(xk_tmp);

    memory->deallocate(evec_tmp);

    if (mympi->my_rank == 0) {
        std::cout << "done!" << std::endl;
    }
}

void Phonon_velocity::calc_phonon_vel_mesh(double **phvel_out,
                                           double ***phvel3_out)
{
    unsigned int nk = kpoint->nk;
    unsigned int ns = dynamical->neval;
    double **vel;

    if (mympi->my_rank == 0) {
        std::cout << " Calculating group velocities of phonons on uniform grid ... ";
    }

    memory->allocate(vel, ns, 3);

    for (unsigned int i = 0; i < nk; ++i) {
        phonon_vel_k(kpoint->xk[i], vel);
        //        phonon_vel_k2(kpoint->xk[i],
        //                      dynamical->eval_phonon[i],
        //                      dynamical->evec_phonon[i],
        //                      vel);

        for (unsigned int j = 0; j < ns; ++j) {
            rotvec(vel[j], vel[j], system->lavec_p);
            for (unsigned int k = 0; k < 3; ++k) {
                vel[j][k] /= 2.0 * pi;
                phvel3_out[i][j][k] = vel[j][k];
            }
            phvel_out[i][j] = std::sqrt(std::pow(vel[j][0], 2)
                + std::pow(vel[j][1], 2)
                + std::pow(vel[j][2], 2));
        }
    }

    memory->deallocate(vel);

    if (mympi->my_rank == 0) {
        std::cout << "done!" << std::endl;
    }
}

void Phonon_velocity::phonon_vel_k(const double *xk_in,
                                   double **vel_out)
{
    unsigned int j;
    unsigned int idiff;
    unsigned int n = dynamical->neval;
    double **xk_shift;
    std::complex<double> **evec_tmp;
    double **omega_shift, *omega_tmp;
    double **kvec_na_tmp;
    double h = 1.0e-4;

    unsigned int ndiff = 2;

    memory->allocate(omega_shift, ndiff, n);
    memory->allocate(xk_shift, ndiff, 3);
    memory->allocate(omega_tmp, ndiff);
    memory->allocate(evec_tmp, 1, 1);
    memory->allocate(kvec_na_tmp, 2, 3);

    for (unsigned int i = 0; i < 3; ++i) {

        for (j = 0; j < 3; ++j) {
            xk_shift[0][j] = xk_in[j];
            xk_shift[1][j] = xk_in[j];
        }

        xk_shift[0][i] -= h;
        xk_shift[1][i] += h;

        // kvec_na_tmp for nonalaytic term
        for (j = 0; j < 3; ++j) {
            kvec_na_tmp[0][j] = xk_shift[0][j];
            kvec_na_tmp[1][j] = xk_shift[1][j];
        }
        rotvec(kvec_na_tmp[0], kvec_na_tmp[0], system->rlavec_p, 'T');
        rotvec(kvec_na_tmp[1], kvec_na_tmp[1], system->rlavec_p, 'T');

        double norm = std::sqrt(kvec_na_tmp[0][0] * kvec_na_tmp[0][0]
            + kvec_na_tmp[0][1] * kvec_na_tmp[0][1]
            + kvec_na_tmp[0][2] * kvec_na_tmp[0][2]);

        if (norm > eps) {
            for (j = 0; j < 3; ++j) kvec_na_tmp[0][j] /= norm;
        }
        norm = std::sqrt(kvec_na_tmp[1][0] * kvec_na_tmp[1][0]
            + kvec_na_tmp[1][1] * kvec_na_tmp[1][1]
            + kvec_na_tmp[1][2] * kvec_na_tmp[1][2]);

        if (norm > eps) {
            for (j = 0; j < 3; ++j) kvec_na_tmp[1][j] /= norm;
        }

        for (idiff = 0; idiff < ndiff; ++idiff) {

            dynamical->eval_k(xk_shift[idiff], kvec_na_tmp[0],
                              fcs_phonon->fc2_ext, omega_shift[idiff],
                              evec_tmp, false);

        }

        for (j = 0; j < n; ++j) {
            for (idiff = 0; idiff < ndiff; ++idiff) {
                omega_tmp[idiff] = dynamical->freq(omega_shift[idiff][j]);
            }
            vel_out[j][i] = diff(omega_tmp, ndiff, h);
        }
    }

    memory->deallocate(xk_shift);
    memory->deallocate(omega_shift);
    memory->deallocate(omega_tmp);
    memory->deallocate(evec_tmp);
    memory->deallocate(kvec_na_tmp);
}

double Phonon_velocity::diff(double *f,
                             const unsigned int n,
                             double h) const
{
    double df;

    if (n == 2) {
        df = (f[1] - f[0]) / (2.0 * h);
    } else {
        error->exit("diff", "Numerical differentiation of n > 2 is not supported yet.");
    }

    return df;
}

void Phonon_velocity::phonon_vel_k2(const double *xk_in,
                                    const double *omega_in,
                                    std::complex<double> **evec_in,
                                    double **vel_out)
{
    unsigned int i, j, l, m;
    unsigned int icrd;
    int nmode = 3 * system->natmin;

    std::complex<double> ***ddyn;
    std::complex<double> ctmp;
    std::complex<double> **vel_tmp;
    std::complex<double> ***mat_tmp;
    std::complex<double> czero(0.0, 0.0);
    std::vector<int> smallgroup_k;
    double **eval_tmp;

    if (dynamical->nonanalytic) {

        error->exit("phonon_vel_k2",
                    "Sorry. Analytic calculation of \
            group velocity is not supported for NONANALYTIC>0.");

    }

    memory->allocate(ddyn, 3, nmode, nmode);
    memory->allocate(vel_tmp, 3, nmode);
    calc_derivative_dynmat_k(xk_in, fcs_phonon->fc2_ext, ddyn);

    bool do_diagonalize = false;


    if (do_diagonalize) {
        // Detect degeneracy at the given k
        double tol_omega = 1.0e-7; // Approximately equal to 0.01 cm^{-1}

        std::vector<int> degeneracy_at_k;

        degeneracy_at_k.clear();

        double omega_prev = omega_in[0];
        int ideg = 1;

        for (i = 1; i < nmode; ++i) {
            double omega_now = omega_in[i];

            if (std::abs(omega_now - omega_prev) < tol_omega) {
                ++ideg;
            } else {
                degeneracy_at_k.push_back(ideg);
                ideg = 1;
                omega_prev = omega_now;
            }
        }
        degeneracy_at_k.push_back(ideg);

        int is = 0;

        for (i = 0; i < degeneracy_at_k.size(); ++i) {
            ideg = degeneracy_at_k[i];

            if (ideg == 1) {

                // When the branch is non-degenerate, the velocity can be calculated 
                // from the diagonal element of e^{*} * DDYN * e.

                for (icrd = 0; icrd < 3; ++icrd) {
                    vel_tmp[icrd][is] = czero;

                    for (l = 0; l < nmode; ++l) {
                        ctmp = czero;
                        for (m = 0; m < nmode; ++m) {
                            ctmp += ddyn[icrd][l][m] * evec_in[is][m];
                        }
                        vel_tmp[icrd][is] += std::conj(evec_in[is][l]) * ctmp;
                    }
                    vel_tmp[icrd][is] /= 2.0 * omega_in[is];
                }

            } else if (ideg > 1) {

                // When the branch is degenerated with two or more branches, 
                // we have to construct a MxM matrix and diagonalize it to obtain 
                // group velocities.

                memory->allocate(mat_tmp, 3, ideg, ideg);
                memory->allocate(eval_tmp, 3, ideg);

                for (icrd = 0; icrd < 3; ++icrd) {

                    for (j = 0; j < ideg; ++j) {
                        for (unsigned int k = 0; k < ideg; ++k) {
                            mat_tmp[icrd][j][k] = czero;

                            for (l = 0; l < nmode; ++l) {
                                ctmp = czero;
                                for (m = 0; m < nmode; ++m) {
                                    ctmp += ddyn[icrd][l][m] * evec_in[j + is][m];
                                }
                                mat_tmp[icrd][j][k] += std::conj(evec_in[k + is][l]) * ctmp;
                            }
                        }
                    }
                    // Diagonalize the matrix here

                    diagonalize_hermite_mat(ideg, mat_tmp[icrd], eval_tmp[icrd]);

                    for (j = 0; j < ideg; ++j) {
                        vel_tmp[icrd][j + is] = eval_tmp[icrd][j] / (2.0 * omega_in[j + is]);
                    }
                }

                memory->deallocate(mat_tmp);
                memory->deallocate(eval_tmp);

            } else {
                error->exit("phonon_vel_k2", "This cannot happen.");
            }

            is += ideg;
        }
    } else {

        for (icrd = 0; icrd < 3; ++icrd) {

            for (j = 0; j < nmode; ++j) {
                vel_tmp[icrd][j] = czero;

                for (l = 0; l < nmode; ++l) {
                    ctmp = czero;
                    for (m = 0; m < nmode; ++m) {
                        ctmp += ddyn[icrd][l][m] * evec_in[j][m];
                    }
                    vel_tmp[icrd][j] += std::conj(evec_in[j][l]) * ctmp;
                }
            }
            for (j = 0; j < nmode; ++j) {
                vel_tmp[icrd][j] /= 2.0 * omega_in[j];
            }
        }
    }


    for (icrd = 0; icrd < 3; ++icrd) {
        for (i = 0; i < nmode; ++i) {
            vel_out[i][icrd] = vel_tmp[icrd][i].real();
        }
    }


    if (ddyn) {
        memory->deallocate(ddyn);
    }
    if (vel_tmp) {
        memory->deallocate(vel_tmp);
    }

    double symmetrizer_k[3][3];

    kpoint->get_small_group_k(xk_in, smallgroup_k, symmetrizer_k);

    for (i = 0; i < nmode; ++i) {
        rotvec(vel_out[i], vel_out[i], symmetrizer_k, 'T');
    }
}


void Phonon_velocity::calc_derivative_dynmat_k(const double *xk_in,
                                               const std::vector<FcsClassExtent> &fc2_in,
                                               std::complex<double> ***ddyn_out) const
{
    int i, j, k;

    int nmode = 3 * system->natmin;

    double vec[3];
    std::complex<double> im(0.0, 1.0);

    for (k = 0; k < 3; ++k) {
        for (i = 0; i < nmode; ++i) {
            for (j = 0; j < nmode; ++j) {
                ddyn_out[k][i][j] = std::complex<double>(0.0, 0.0);
            }
        }
    }

    for (const auto &it : fc2_in) {

        unsigned int atm1_p = it.atm1;
        unsigned int atm2_s = it.atm2;
        unsigned int xyz1 = it.xyz1;
        unsigned int xyz2 = it.xyz2;
        unsigned int icell = it.cell_s;

        unsigned int atm1_s = system->map_p2s[atm1_p][0];
        unsigned int atm2_p = system->map_s2p[atm2_s].atom_num;

        for (i = 0; i < 3; ++i) {
            vec[i] = system->xr_s[atm2_s][i] + xshift_s[icell][i]
                - system->xr_s[system->map_p2s[atm2_p][0]][i];
        }

        rotvec(vec, vec, system->lavec_s);
        rotvec(vec, vec, system->rlavec_p);

        double phase = vec[0] * xk_in[0] + vec[1] * xk_in[1] + vec[2] * xk_in[2];

        for (k = 0; k < 3; ++k) {
            ddyn_out[k][3 * atm1_p + xyz1][3 * atm2_p + xyz2]
                += it.fcs_val * std::exp(im * phase) * vec[k] / std::sqrt(
                    system->mass[atm1_s] * system->mass[atm2_s]);
        }

    }

    for (k = 0; k < 3; ++k) {
        for (i = 0; i < nmode; ++i) {
            for (j = 0; j < nmode; ++j) {
                ddyn_out[k][i][j] *= std::complex<double>(0.0, 1.0);
            }
        }
    }
}


void Phonon_velocity::diagonalize_hermite_mat(const int n,
                                              std::complex<double> **mat_in,
                                              double *eval_out) const
{
    std::complex<double> *mat_1D;
    int LWORK = (2 * n - 1) * 10;
    int INFO;
    std::complex<double> *WORK;
    double *RWORK;
    char JOBZ = 'N';
    char UPLO = 'U';
    int n_ = n;

    memory->allocate(mat_1D, n * n);
    memory->allocate(RWORK, 3 * n - 2);
    memory->allocate(WORK, LWORK);

    int k = 0;
    for (int j = 0; j < n; ++j) {
        for (int i = 0; i < n; ++i) {
            mat_1D[k++] = mat_in[i][j];
        }
    }

    zheev_(&JOBZ, &UPLO, &n_, mat_1D, &n_, eval_out, WORK, &LWORK, RWORK, &INFO);

    memory->deallocate(RWORK);
    memory->deallocate(WORK);
    memory->deallocate(mat_1D);
}
