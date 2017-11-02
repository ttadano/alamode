/*
 selfenergy.cpp

 Copyright (c) 2014 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#include "mpi_common.h"
#include "selfenergy.h"
#include "constants.h"
#include "relaxation.h"
#include "dynamical.h"
#include "kpoint.h"
#include "memory.h"
#include "thermodynamics.h"
#include "mathfunctions.h"
#include "integration.h"

using namespace PHON_NS;

Selfenergy::Selfenergy(PHON *phon): Pointers(phon)
{
    im = std::complex<double>(0.0, 1.0);
}

Selfenergy::~Selfenergy()
{
}

void Selfenergy::setup_selfenergy()
{
    nk = kpoint->nk;
    ns = dynamical->neval;
    epsilon = integration->epsilon;
}


void Selfenergy::mpi_reduce_complex(unsigned int N,
                                    std::complex<double> *in_mpi,
                                    std::complex<double> *out)
{
#ifdef MPI_COMPLEX16
    MPI_Reduce(&in_mpi[0], &out[0], N, MPI_COMPLEX16, MPI_SUM, 0, MPI_COMM_WORLD);
#else
    unsigned int i;
    double *ret_mpi_re, *ret_mpi_im;
    double *ret_re, *ret_im;

    memory->allocate(ret_mpi_re, N);
    memory->allocate(ret_mpi_im, N);
    memory->allocate(ret_im, N);
    memory->allocate(ret_re, N);

    for (i = 0; i < N; ++i) {
        ret_mpi_re[i] = in_mpi[i].real();
        ret_mpi_im[i] = in_mpi[i].imag();
    }
    MPI_Reduce(&ret_mpi_re[0], &ret_re[0], N, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&ret_mpi_im[0], &ret_im[0], N, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    for (i = 0; i < N; ++i) {
        out[i] = ret_re[i] + im * ret_im[i];
    }
    memory->deallocate(ret_mpi_re);
    memory->deallocate(ret_mpi_im);
    memory->deallocate(ret_re);
    memory->deallocate(ret_im);
#endif
}

void Selfenergy::selfenergy_tadpole(const unsigned int N,
                                    double *T,
                                    const double omega,
                                    const unsigned int knum,
                                    const unsigned int snum,
                                    std::complex<double> *ret)
{
    unsigned int i;
    unsigned int ik2;
    unsigned int is1, is2;
    unsigned int arr_cubic1[3], arr_cubic2[3];
    std::complex<double> v3_tmp1, v3_tmp2;
    std::complex<double> *ret_mpi, *ret_tmp;

    double T_tmp;
    double n2;
    double omega1, omega2;
    double factor;

    arr_cubic1[0] = ns * kpoint->knum_minus[knum] + snum;
    arr_cubic1[1] = ns * knum + snum;

    memory->allocate(ret_mpi, N);
    memory->allocate(ret_tmp, N);

    for (i = 0; i < N; ++i) ret[i] = std::complex<double>(0.0, 0.0);

    for (is1 = 0; is1 < ns; ++is1) {
        arr_cubic1[2] = is1;
        arr_cubic2[0] = is1;
        omega1 = dynamical->eval_phonon[0][is1];

        if (omega1 < eps8) continue;

        v3_tmp1 = relaxation->V3(arr_cubic1);

        for (i = 0; i < N; ++i) ret_mpi[i] = std::complex<double>(0.0, 0.0);

        for (ik2 = mympi->my_rank; ik2 < nk; ik2 += mympi->nprocs) {
            for (is2 = 0; is2 < ns; ++is2) {
                arr_cubic2[1] = ns * ik2 + is2;
                arr_cubic2[2] = ns * kpoint->knum_minus[ik2] + is2;

                v3_tmp2 = relaxation->V3(arr_cubic2);
                omega2 = dynamical->eval_phonon[ik2][is2];

                if (omega2 < eps8) continue;

                for (i = 0; i < N; ++i) {
                    T_tmp = T[i];
                    n2 = thermodynamics->fB(omega2, T_tmp);
                    //   ret[i] += v3_tmp1 * v3_tmp2 * (2.0*n2 + 1.0) / omega1; 
                    ret_mpi[i] += v3_tmp2 * (2.0 * n2 + 1.0);
                }
            }
        }
        mpi_reduce_complex(N, ret_mpi, ret_tmp);

        for (i = 0; i < N; ++i) {
            ret[i] += ret_tmp[i] * v3_tmp1 / omega1;
        }
    }

    factor = -1.0 / (static_cast<double>(nk) * std::pow(2.0, 3));
    for (i = 0; i < N; ++i) ret[i] *= factor;

    memory->deallocate(ret_tmp);
    memory->deallocate(ret_mpi);
}

void Selfenergy::selfenergy_a(const unsigned int N,
                              double *T,
                              const double omega,
                              const unsigned int knum,
                              const unsigned int snum,
                              std::complex<double> *ret)
{
    /*

    Diagram (a)  
    Matrix elements that appear : V3^2
    Computational cost          : O(N_k * N^2)

    */

    unsigned int i;
    unsigned int ik1, ik2;
    unsigned int is1, is2;
    unsigned int arr_cubic[3];
    double xk_tmp[3];
    double v3_tmp;
    std::complex<double> omega_shift;
    std::complex<double> omega_sum[2];
    std::complex<double> *ret_mpi;

    unsigned int nkx = kpoint->nkx;
    unsigned int nky = kpoint->nky;
    unsigned int nkz = kpoint->nkz;

    unsigned int iloc, jloc, kloc;

    double T_tmp;
    double n1, n2;
    double omega1, omega2;
    double factor;

    arr_cubic[0] = ns * kpoint->knum_minus[knum] + snum;

    omega_shift = omega + im * epsilon;

    memory->allocate(ret_mpi, N);

    for (i = 0; i < N; ++i) ret_mpi[i] = std::complex<double>(0.0, 0.0);

    for (ik1 = mympi->my_rank; ik1 < nk; ik1 += mympi->nprocs) {

        xk_tmp[0] = kpoint->xk[knum][0] - kpoint->xk[ik1][0];
        xk_tmp[1] = kpoint->xk[knum][1] - kpoint->xk[ik1][1];
        xk_tmp[2] = kpoint->xk[knum][2] - kpoint->xk[ik1][2];

        iloc = (nint(xk_tmp[0] * static_cast<double>(nkx) + static_cast<double>(2 * nkx))) % nkx;
        jloc = (nint(xk_tmp[1] * static_cast<double>(nky) + static_cast<double>(2 * nky))) % nky;
        kloc = (nint(xk_tmp[2] * static_cast<double>(nkz) + static_cast<double>(2 * nkz))) % nkz;

        ik2 = kloc + nkz * jloc + nky * nkz * iloc;

        for (is1 = 0; is1 < ns; ++is1) {

            arr_cubic[1] = ns * ik1 + is1;
            omega1 = dynamical->eval_phonon[ik1][is1];

            for (is2 = 0; is2 < ns; ++is2) {

                arr_cubic[2] = ns * ik2 + is2;
                omega2 = dynamical->eval_phonon[ik2][is2];

                v3_tmp = std::norm(relaxation->V3(arr_cubic));

                omega_sum[0] = 1.0 / (omega_shift + omega1 + omega2) - 1.0 / (omega_shift - omega1 - omega2);
                omega_sum[1] = 1.0 / (omega_shift + omega1 - omega2) - 1.0 / (omega_shift - omega1 + omega2);

                for (i = 0; i < N; ++i) {
                    T_tmp = T[i];
                    n1 = thermodynamics->fB(omega1, T_tmp);
                    n2 = thermodynamics->fB(omega2, T_tmp);

                    ret_mpi[i] += v3_tmp * ((1.0 + n1 + n2) * omega_sum[0] + (n2 - n1) * omega_sum[1]);
                }
            }
        }
    }

    factor = 1.0 / (static_cast<double>(nk) * std::pow(2.0, 4));
    for (i = 0; i < N; ++i) ret_mpi[i] *= factor;

    mpi_reduce_complex(N, ret_mpi, ret);

    memory->deallocate(ret_mpi);
}

void Selfenergy::selfenergy_b(const unsigned int N,
                              double *T,
                              const double omega,
                              const unsigned int knum,
                              const unsigned int snum,
                              std::complex<double> *ret)
{
    /*
    Diagram (b)
    Matrix elements that appear : V4
    Computational cost          : O(N_k * N)
    Note                        : This give rise to the phonon frequency-shift only.
    */

    unsigned int i;
    unsigned int ik1;
    unsigned int is1;
    unsigned int arr_quartic[4];

    double omega1;
    double n1;
    double factor;

    std::complex<double> v4_tmp;
    std::complex<double> *ret_mpi;

    memory->allocate(ret_mpi, N);

    for (i = 0; i < N; ++i) ret_mpi[i] = std::complex<double>(0.0, 0.0);

    arr_quartic[0] = ns * kpoint->knum_minus[knum] + snum;
    arr_quartic[3] = ns * knum + snum;


    for (ik1 = mympi->my_rank; ik1 < nk; ik1 += mympi->nprocs) {
        for (is1 = 0; is1 < ns; ++is1) {

            arr_quartic[1] = ns * ik1 + is1;
            arr_quartic[2] = ns * kpoint->knum_minus[ik1] + is1;

            omega1 = dynamical->eval_phonon[ik1][is1];
            if (omega1 < eps8) continue;

            v4_tmp = relaxation->V4(arr_quartic);

            for (i = 0; i < N; ++i) {
                n1 = thermodynamics->fB(omega1, T[i]);
                ret_mpi[i] += v4_tmp * (2.0 * n1 + 1.0);
            }
        }
    }

    factor = -1.0 / (static_cast<double>(nk) * std::pow(2.0, 3));
    for (i = 0; i < N; ++i) ret_mpi[i] *= factor;

    mpi_reduce_complex(N, ret_mpi, ret);

    memory->deallocate(ret_mpi);
}


void Selfenergy::selfenergy_c(const unsigned int N,
                              double *T,
                              const double omega,
                              const unsigned int knum,
                              const unsigned int snum,
                              std::complex<double> *ret)
{
    /* 

    Diagram (c)
    Matrix elements that appear : V4^2
    Computational cost          : O(N_k^2 * N^3) <-- about N_k * N times that of Diagram (a)

    */

    unsigned int i;
    unsigned int ik1, ik2, ik3;
    unsigned int is1, is2, is3;
    unsigned int arr_quartic[4];
    unsigned int iloc, jloc, kloc;

    unsigned int nkx = kpoint->nkx;
    unsigned int nky = kpoint->nky;
    unsigned int nkz = kpoint->nkz;

    double xk_tmp[3];
    double v4_tmp;
    double omega1, omega2, omega3;
    double n1, n2, n3;
    double n12, n23, n31;
    double T_tmp;
    double factor;

    std::complex<double> omega_shift;
    std::complex<double> omega_sum[4];
    std::complex<double> *ret_mpi;

    memory->allocate(ret_mpi, N);

    omega_shift = omega + im * epsilon;

    for (i = 0; i < N; ++i) ret_mpi[i] = std::complex<double>(0.0, 0.0);

    arr_quartic[0] = ns * kpoint->knum_minus[knum] + snum;

    for (ik1 = mympi->my_rank; ik1 < nk; ik1 += mympi->nprocs) {
        for (ik2 = 0; ik2 < nk; ++ik2) {

            xk_tmp[0] = kpoint->xk[knum][0] - kpoint->xk[ik1][0] - kpoint->xk[ik2][0];
            xk_tmp[1] = kpoint->xk[knum][1] - kpoint->xk[ik1][1] - kpoint->xk[ik2][1];
            xk_tmp[2] = kpoint->xk[knum][2] - kpoint->xk[ik1][2] - kpoint->xk[ik2][2];

            iloc = (nint(xk_tmp[0] * static_cast<double>(nkx) + static_cast<double>(2 * nkx))) % nkx;
            jloc = (nint(xk_tmp[1] * static_cast<double>(nky) + static_cast<double>(2 * nky))) % nky;
            kloc = (nint(xk_tmp[2] * static_cast<double>(nkz) + static_cast<double>(2 * nkz))) % nkz;

            ik3 = kloc + nkz * jloc + nky * nkz * iloc;

            for (is1 = 0; is1 < ns; ++is1) {

                arr_quartic[1] = ns * ik1 + is1;
                omega1 = dynamical->eval_phonon[ik1][is1];

                for (is2 = 0; is2 < ns; ++is2) {

                    arr_quartic[2] = ns * ik2 + is2;
                    omega2 = dynamical->eval_phonon[ik2][is2];

                    for (is3 = 0; is3 < ns; ++is3) {

                        arr_quartic[3] = ns * ik3 + is3;
                        omega3 = dynamical->eval_phonon[ik3][is3];

                        v4_tmp = std::norm(relaxation->V4(arr_quartic));

                        omega_sum[0]
                            = 1.0 / (omega_shift - omega1 - omega2 - omega3)
                            - 1.0 / (omega_shift + omega1 + omega2 + omega3);
                        omega_sum[1]
                            = 1.0 / (omega_shift - omega1 - omega2 + omega3)
                            - 1.0 / (omega_shift + omega1 + omega2 - omega3);
                        omega_sum[2]
                            = 1.0 / (omega_shift + omega1 - omega2 - omega3)
                            - 1.0 / (omega_shift - omega1 + omega2 + omega3);
                        omega_sum[3]
                            = 1.0 / (omega_shift - omega1 + omega2 - omega3)
                            - 1.0 / (omega_shift + omega1 - omega2 + omega3);

                        for (i = 0; i < N; ++i) {
                            T_tmp = T[i];

                            n1 = thermodynamics->fB(omega1, T_tmp);
                            n2 = thermodynamics->fB(omega2, T_tmp);
                            n3 = thermodynamics->fB(omega3, T_tmp);

                            n12 = n1 * n2;
                            n23 = n2 * n3;
                            n31 = n3 * n1;

                            ret_mpi[i] += v4_tmp
                                * ((n12 + n23 + n31 + n1 + n2 + n3 + 1.0) * omega_sum[0]
                                    + (n31 + n23 + n3 - n12) * omega_sum[1]
                                    + (n12 + n31 + n1 - n23) * omega_sum[2]
                                    + (n23 + n12 + n2 - n31) * omega_sum[3]);
                        }
                    }
                }
            }
        }
    }

    factor = 1.0 / (std::pow(static_cast<double>(nk), 2) * std::pow(2.0, 5) * 3.0);
    for (i = 0; i < N; ++i) ret_mpi[i] *= factor;

    mpi_reduce_complex(N, ret_mpi, ret);

    memory->deallocate(ret_mpi);
}


void Selfenergy::selfenergy_d(const unsigned int N,
                              double *T,
                              const double omega,
                              const unsigned int knum,
                              const unsigned int snum,
                              std::complex<double> *ret)
{
    /*

    Diagram (d)
    Matrix elements that appear : V3^2 V4
    Computational cost          : O(N_k^2 * N^4)
    Note                        : 2 3-point vertexes and 1 4-point vertex.

    */

    unsigned int i;
    unsigned int ik1, ik2, ik3, ik4;
    unsigned int is1, is2, is3, is4;
    unsigned int arr_cubic1[3], arr_cubic2[3];
    unsigned int arr_quartic[4];
    unsigned int iloc, jloc, kloc;
    unsigned int nkx = kpoint->nkx;
    unsigned int nky = kpoint->nky;
    unsigned int nkz = kpoint->nkz;

    double xk_tmp[3];
    double n1, n2, n3, n4;
    double omega1, omega2, omega3, omega4;
    double T_tmp;
    double factor;

    std::complex<double> v3_tmp1, v3_tmp2, v4_tmp;
    std::complex<double> v_prod;
    std::complex<double> omega_shift;
    std::complex<double> omega_sum[4];
    std::complex<double> *ret_mpi;

    memory->allocate(ret_mpi, N);

    omega_shift = omega + im * epsilon;

    for (i = 0; i < N; ++i) ret[i] = std::complex<double>(0.0, 0.0);

    arr_cubic1[0] = ns * kpoint->knum_minus[knum] + snum;
    arr_cubic2[2] = ns * knum + snum;

    for (ik1 = mympi->my_rank; ik1 < nk; ik1 += mympi->nprocs) {

        xk_tmp[0] = kpoint->xk[knum][0] - kpoint->xk[ik1][0];
        xk_tmp[1] = kpoint->xk[knum][1] - kpoint->xk[ik1][1];
        xk_tmp[2] = kpoint->xk[knum][2] - kpoint->xk[ik1][2];

        iloc = (nint(xk_tmp[0] * static_cast<double>(nkx) + static_cast<double>(2 * nkx))) % nkx;
        jloc = (nint(xk_tmp[1] * static_cast<double>(nky) + static_cast<double>(2 * nky))) % nky;
        kloc = (nint(xk_tmp[2] * static_cast<double>(nkz) + static_cast<double>(2 * nkz))) % nkz;

        ik2 = kloc + nkz * jloc + nky * nkz * iloc;

        for (ik3 = 0; ik3 < nk; ++ik3) {

            xk_tmp[0] = kpoint->xk[knum][0] - kpoint->xk[ik3][0];
            xk_tmp[1] = kpoint->xk[knum][1] - kpoint->xk[ik3][1];
            xk_tmp[2] = kpoint->xk[knum][2] - kpoint->xk[ik3][2];

            iloc = (nint(xk_tmp[0] * static_cast<double>(nkx) + static_cast<double>(2 * nkx))) % nkx;
            jloc = (nint(xk_tmp[1] * static_cast<double>(nky) + static_cast<double>(2 * nky))) % nky;
            kloc = (nint(xk_tmp[2] * static_cast<double>(nkz) + static_cast<double>(2 * nkz))) % nkz;

            ik4 = kloc + nkz * jloc + nky * nkz * iloc;

            for (is1 = 0; is1 < ns; ++is1) {

                omega1 = dynamical->eval_phonon[ik1][is1];

                arr_cubic2[0] = ns * kpoint->knum_minus[ik1] + is1;
                arr_quartic[0] = ns * ik1 + is1;

                for (is2 = 0; is2 < ns; ++is2) {

                    omega2 = dynamical->eval_phonon[ik2][is2];

                    arr_cubic2[1] = ns * kpoint->knum_minus[ik2] + is2;
                    arr_quartic[1] = ns * ik2 + is2;

                    v3_tmp2 = relaxation->V3(arr_cubic2);

                    for (is3 = 0; is3 < ns; ++is3) {

                        omega3 = dynamical->eval_phonon[ik3][is3];

                        arr_cubic1[1] = ns * ik3 + is3;
                        arr_quartic[2] = ns * kpoint->knum_minus[ik3] + is3;

                        for (is4 = 0; is4 < ns; ++is4) {

                            omega4 = dynamical->eval_phonon[ik4][is4];

                            arr_cubic1[2] = ns * ik4 + is4;
                            arr_quartic[3] = ns * kpoint->knum_minus[ik4] + is4;

                            v3_tmp1 = relaxation->V3(arr_cubic1);
                            v4_tmp = relaxation->V4(arr_quartic);

                            v_prod = v3_tmp1 * v3_tmp2 * v4_tmp;

                            omega_sum[0]
                                = 1.0 / (omega_shift + omega1 + omega2)
                                - 1.0 / (omega_shift - omega1 - omega2);
                            omega_sum[1]
                                = 1.0 / (omega_shift + omega1 - omega2)
                                - 1.0 / (omega_shift - omega1 + omega2);
                            omega_sum[2]
                                = 1.0 / (omega_shift + omega3 + omega4)
                                - 1.0 / (omega_shift - omega3 - omega4);
                            omega_sum[3]
                                = 1.0 / (omega_shift + omega3 - omega4)
                                - 1.0 / (omega_shift - omega3 + omega4);

                            for (i = 0; i < N; ++i) {
                                T_tmp = T[i];

                                n1 = thermodynamics->fB(omega1, T_tmp);
                                n2 = thermodynamics->fB(omega2, T_tmp);
                                n3 = thermodynamics->fB(omega3, T_tmp);
                                n4 = thermodynamics->fB(omega4, T_tmp);

                                ret_mpi[i] += v_prod
                                    * ((1.0 + n1 + n2) * omega_sum[0] + (n2 - n1) * omega_sum[1])
                                    * ((1.0 + n3 + n4) * omega_sum[2] + (n4 - n3) * omega_sum[3]);
                            }
                        }
                    }
                }
            }
        }
    }

    factor = -1.0 / (std::pow(static_cast<double>(nk), 2) * std::pow(2.0, 7));
    for (i = 0; i < N; ++i) ret_mpi[i] *= factor;

    mpi_reduce_complex(N, ret_mpi, ret);

    memory->deallocate(ret_mpi);
}

void Selfenergy::selfenergy_e(const unsigned int N,
                              double *T,
                              const double omega,
                              const unsigned int knum,
                              const unsigned int snum,
                              std::complex<double> *ret)
{
    /*

    Diagram (e)
    Matrix elements that appear : V3^2 V4
    Computational cost          : O(N_k^2 * N^4)
    Note                        : Double pole appears when omega1 = omega2.

    */

    unsigned int i;
    unsigned int ik1, ik2, ik3, ik4;
    unsigned int is1, is2, is3, is4;
    unsigned int iloc, jloc, kloc;
    unsigned int arr_cubic1[3], arr_cubic2[3];
    unsigned int arr_quartic[4];
    unsigned int nkx = kpoint->nkx;
    unsigned int nky = kpoint->nky;
    unsigned int nkz = kpoint->nkz;

    int ip1, ip4;

    double T_tmp;
    double factor;
    double omega1, omega2, omega3, omega4;
    double dp1, dp4;
    double dp1_inv;
    double n1, n2, n3, n4;
    double xk_tmp[3];
    double D12[2];
    double T_inv;

    std::complex<double> v3_tmp1, v3_tmp2, v4_tmp;
    std::complex<double> v_prod;
    std::complex<double> omega_shift;
    std::complex<double> omega_sum;
    std::complex<double> omega_sum14[4], omega_sum24[4];
    std::complex<double> omega_prod[6];
    std::complex<double> *prod_tmp;
    std::complex<double> *ret_mpi;

    memory->allocate(ret_mpi, N);
    memory->allocate(prod_tmp, N);

    omega_shift = omega + im * epsilon;

    for (i = 0; i < N; ++i) ret_mpi[i] = std::complex<double>(0.0, 0.0);

    arr_cubic1[0] = ns * kpoint->knum_minus[knum] + snum;
    arr_cubic2[2] = ns * knum + snum;

    for (ik1 = mympi->my_rank; ik1 < nk; ik1 += mympi->nprocs) {

        ik2 = ik1;

        xk_tmp[0] = kpoint->xk[knum][0] - kpoint->xk[ik1][0];
        xk_tmp[1] = kpoint->xk[knum][1] - kpoint->xk[ik1][1];
        xk_tmp[2] = kpoint->xk[knum][2] - kpoint->xk[ik1][2];

        iloc = (nint(xk_tmp[0] * static_cast<double>(nkx) + static_cast<double>(2 * nkx))) % nkx;
        jloc = (nint(xk_tmp[1] * static_cast<double>(nky) + static_cast<double>(2 * nky))) % nky;
        kloc = (nint(xk_tmp[2] * static_cast<double>(nkz) + static_cast<double>(2 * nkz))) % nkz;

        ik4 = kloc + nkz * jloc + nky * nkz * iloc;

        for (ik3 = 0; ik3 < nk; ++ik3) {

            for (is1 = 0; is1 < ns; ++is1) {

                omega1 = dynamical->eval_phonon[ik1][is1];

                arr_cubic1[1] = ns * ik1 + is1;
                arr_quartic[0] = ns * kpoint->knum_minus[ik1] + is1;

                for (is2 = 0; is2 < ns; ++is2) {

                    omega2 = dynamical->eval_phonon[ik2][is2];

                    arr_cubic2[0] = ns * kpoint->knum_minus[ik2] + is2;
                    arr_quartic[3] = ns * ik2 + is2;

                    if (std::abs(omega1 - omega2) < eps) {

                        for (is3 = 0; is3 < ns; ++is3) {

                            omega3 = dynamical->eval_phonon[ik3][is3];

                            arr_quartic[1] = ns * ik3 + is3;
                            arr_quartic[2] = ns * kpoint->knum_minus[ik3] + is3;

                            v4_tmp = relaxation->V4(arr_quartic);

                            for (is4 = 0; is4 < ns; ++is4) {

                                omega4 = dynamical->eval_phonon[ik4][is4];

                                arr_cubic1[2] = ns * ik4 + is4;
                                arr_cubic2[1] = ns * kpoint->knum_minus[ik4] + is4;

                                v3_tmp1 = relaxation->V3(arr_cubic1);
                                v3_tmp2 = relaxation->V3(arr_cubic2);

                                v_prod = v3_tmp1 * v3_tmp2 * v4_tmp;

                                for (i = 0; i < N; ++i) prod_tmp[i] = std::complex<double>(0.0, 0.0);

                                for (ip1 = 1; ip1 >= -1; ip1 -= 2) {
                                    dp1 = static_cast<double>(ip1) * omega1;
                                    dp1_inv = 1.0 / dp1;

                                    for (ip4 = 1; ip4 >= -1; ip4 -= 2) {
                                        dp4 = static_cast<double>(ip4) * omega4;

                                        omega_sum = 1.0 / (omega_shift + dp1 + dp4);

                                        for (i = 0; i < N; ++i) {
                                            T_tmp = T[i];

                                            n1 = thermodynamics->fB(dp1, T_tmp);
                                            n4 = thermodynamics->fB(dp4, T_tmp);


                                            if (std::abs(T_tmp) < eps) {
                                                //special treatment for T = 0
                                                // This is valid since beta always appears as a product beta*n
                                                // which is zero when T = 0.
                                                T_inv = 0.0;
                                            } else {
                                                T_inv = 1.0 / (thermodynamics->T_to_Ryd * T_tmp);
                                            }

                                            prod_tmp[i] += static_cast<double>(ip4) * omega_sum
                                                * ((1.0 + n1 + n4) * omega_sum
                                                    + (1.0 + n1 + n4) * dp1_inv + n1 * (1.0 + n1) * T_inv);
                                        }
                                    }
                                }

                                for (i = 0; i < N; ++i) {
                                    T_tmp = T[i];

                                    n3 = thermodynamics->fB(omega3, T_tmp);
                                    ret_mpi[i] += v_prod * (2.0 * n3 + 1.0) * prod_tmp[i];
                                }
                            }
                        }

                    } else {

                        D12[0] = 1.0 / (omega1 + omega2) - 1.0 / (omega1 - omega2);
                        D12[1] = 1.0 / (omega1 + omega2) + 1.0 / (omega1 - omega2);

                        for (is3 = 0; is3 < ns; ++is3) {

                            omega3 = dynamical->eval_phonon[ik3][is3];

                            arr_quartic[1] = ns * ik3 + is3;
                            arr_quartic[2] = ns * kpoint->knum_minus[ik3] + is3;

                            v4_tmp = relaxation->V4(arr_quartic);

                            for (is4 = 0; is4 < ns; ++is4) {

                                omega4 = dynamical->eval_phonon[ik4][is4];

                                arr_cubic1[2] = ns * ik4 + is4;
                                arr_cubic2[1] = ns * kpoint->knum_minus[ik4] + is4;

                                v3_tmp1 = relaxation->V3(arr_cubic1);
                                v3_tmp2 = relaxation->V3(arr_cubic2);

                                v_prod = v3_tmp1 * v3_tmp2 * v4_tmp;

                                omega_sum14[0] = 1.0 / (omega_shift + omega1 + omega4);
                                omega_sum14[1] = 1.0 / (omega_shift + omega1 - omega4);
                                omega_sum14[2] = 1.0 / (omega_shift - omega1 + omega4);
                                omega_sum14[3] = 1.0 / (omega_shift - omega1 - omega4);

                                omega_sum24[0] = 1.0 / (omega_shift + omega2 + omega4);
                                omega_sum24[1] = 1.0 / (omega_shift + omega2 - omega4);
                                omega_sum24[2] = 1.0 / (omega_shift - omega2 + omega4);
                                omega_sum24[3] = 1.0 / (omega_shift - omega2 - omega4);

                                omega_prod[0] = D12[0] * (omega_sum14[0] - omega_sum14[1]);
                                omega_prod[1] = D12[0] * (omega_sum14[2] - omega_sum14[3]);
                                omega_prod[2] = D12[1] * (omega_sum24[0] - omega_sum24[1]);
                                omega_prod[3] = D12[1] * (omega_sum24[2] - omega_sum24[3]);
                                omega_prod[4] = (omega_sum14[1] - omega_sum14[3])
                                    * (omega_sum24[1] - omega_sum24[3]);
                                omega_prod[5] = (omega_sum14[0] - omega_sum14[2])
                                    * (omega_sum24[0] - omega_sum24[2]);

                                for (i = 0; i < N; ++i) {
                                    T_tmp = T[i];

                                    n1 = thermodynamics->fB(omega1, T_tmp);
                                    n2 = thermodynamics->fB(omega2, T_tmp);
                                    n3 = thermodynamics->fB(omega3, T_tmp);
                                    n4 = thermodynamics->fB(omega4, T_tmp);

                                    ret_mpi[i] += v_prod * (2.0 * n3 + 1.0)
                                        * ((1.0 + n1) * omega_prod[0] + n1 * omega_prod[1]
                                            + (1.0 + n2) * omega_prod[2] + n2 * omega_prod[3]
                                            + (1.0 + n4) * omega_prod[4] + n4 * omega_prod[5]);

                                    /*
                                    ret[i] *= v3_tmp1 * v3_tmp2 * v4_tmp * (2.0 * n3 + 1.0) * (2.0 * omega2) / (omega1 * omega1 - omega2 * omega2)
                                    * ((1.0 + n1 + n4) * (1.0 / (omega - omega1 - omega4 + im * epsilon) - 1.0 / (omega + omega1 + omega4 + im * epsilon)) 
                                    + (n4 - n1) * (1.0 / (omega - omega1 + omega4 + im * epsilon) - 1.0 / (omega + omega1 - omega4 + im * epsilon)));
                                    */
                                }
                            }
                        }

                    }
                }
            }
        }
    }

    factor = -1.0 / (std::pow(static_cast<double>(nk), 2) * std::pow(2.0, 6));
    //	factor = -1.0 / (std::pow(static_cast<double>(nk), 2) * std::pow(2.0, 7));
    for (i = 0; i < N; ++i) ret_mpi[i] *= factor;

    mpi_reduce_complex(N, ret_mpi, ret);

    memory->deallocate(prod_tmp);
    memory->deallocate(ret_mpi);
}

void Selfenergy::selfenergy_f(const unsigned int N,
                              double *T,
                              const double omega,
                              const unsigned int knum,
                              const unsigned int snum,
                              std::complex<double> *ret)
{
    /*
    Diagram (f)
    Matrix elements that appear : V3^4
    Computational cost          : O(N_k^2 * N^5)
    Note                        : Computationally expensive & double pole when omega1 = omega5.
    */

    unsigned int i;
    unsigned int ik1, ik2, ik3, ik4, ik5;
    unsigned int is1, is2, is3, is4, is5;
    unsigned int arr_cubic1[3], arr_cubic2[3], arr_cubic3[3], arr_cubic4[3];
    unsigned int iloc, jloc, kloc;
    unsigned int nkx = kpoint->nkx;
    unsigned int nky = kpoint->nky;
    unsigned int nkz = kpoint->nkz;

    int ip1, ip2, ip3, ip4, ip5;

    double omega1, omega2, omega3, omega4, omega5;
    double n1, n2, n3, n4, n5;
    double xk_tmp[3];
    double dp1, dp2, dp3, dp4, dp5;
    double T_tmp;
    double dp1_inv;
    double factor;
    double D15, D134, D345;
    double T_inv;

    std::complex<double> omega_sum[3];
    std::complex<double> v3_tmp1, v3_tmp2, v3_tmp3, v3_tmp4;
    std::complex<double> v3_prod;
    std::complex<double> omega_shift;
    std::complex<double> *ret_mpi;

    memory->allocate(ret_mpi, N);

    omega_shift = omega + im * epsilon;

    for (i = 0; i < N; ++i) ret_mpi[i] = std::complex<double>(0.0, 0.0);

    arr_cubic1[0] = ns * kpoint->knum_minus[knum] + snum;
    arr_cubic4[2] = ns * knum + snum;

    for (ik1 = mympi->my_rank; ik1 < nk; ik1 += mympi->nprocs) {

        ik5 = ik1;

        xk_tmp[0] = kpoint->xk[knum][0] - kpoint->xk[ik1][0];
        xk_tmp[1] = kpoint->xk[knum][1] - kpoint->xk[ik1][1];
        xk_tmp[2] = kpoint->xk[knum][2] - kpoint->xk[ik1][2];

        iloc = (nint(xk_tmp[0] * static_cast<double>(nkx) + static_cast<double>(2 * nkx))) % nkx;
        jloc = (nint(xk_tmp[1] * static_cast<double>(nky) + static_cast<double>(2 * nky))) % nky;
        kloc = (nint(xk_tmp[2] * static_cast<double>(nkz) + static_cast<double>(2 * nkz))) % nkz;

        ik2 = kloc + nkz * jloc + nky * nkz * iloc;

        for (ik3 = 0; ik3 < nk; ++ik3) {

            xk_tmp[0] = kpoint->xk[ik1][0] - kpoint->xk[ik3][0];
            xk_tmp[1] = kpoint->xk[ik1][1] - kpoint->xk[ik3][1];
            xk_tmp[2] = kpoint->xk[ik1][2] - kpoint->xk[ik3][2];

            iloc = (nint(xk_tmp[0] * static_cast<double>(nkx) + static_cast<double>(2 * nkx))) % nkx;
            jloc = (nint(xk_tmp[1] * static_cast<double>(nky) + static_cast<double>(2 * nky))) % nky;
            kloc = (nint(xk_tmp[2] * static_cast<double>(nkz) + static_cast<double>(2 * nkz))) % nkz;

            ik4 = kloc + nkz * jloc + nky * nkz * iloc;

            for (is1 = 0; is1 < ns; ++is1) {

                omega1 = dynamical->eval_phonon[ik1][is1];

                arr_cubic1[1] = ns * ik1 + is1;
                arr_cubic2[0] = ns * kpoint->knum_minus[ik1] + is1;

                for (is2 = 0; is2 < ns; ++is2) {

                    omega2 = dynamical->eval_phonon[ik2][is2];

                    arr_cubic1[2] = ns * ik2 + is2;
                    arr_cubic4[1] = ns * kpoint->knum_minus[ik2] + is2;

                    v3_tmp1 = relaxation->V3(arr_cubic1);

                    for (is5 = 0; is5 < ns; ++is5) {

                        omega5 = dynamical->eval_phonon[ik5][is5];

                        arr_cubic3[2] = ns * ik5 + is5;
                        arr_cubic4[0] = ns * kpoint->knum_minus[ik5] + is5;

                        v3_tmp4 = relaxation->V3(arr_cubic4);

                        for (is3 = 0; is3 < ns; ++is3) {

                            omega3 = dynamical->eval_phonon[ik3][is3];

                            arr_cubic2[1] = ns * ik3 + is3;
                            arr_cubic3[0] = ns * kpoint->knum_minus[ik3] + is3;

                            for (is4 = 0; is4 < ns; ++is4) {

                                omega4 = dynamical->eval_phonon[ik4][is4];

                                arr_cubic2[2] = ns * ik4 + is4;
                                arr_cubic3[1] = ns * kpoint->knum_minus[ik4] + is4;

                                v3_tmp2 = relaxation->V3(arr_cubic2);
                                v3_tmp3 = relaxation->V3(arr_cubic3);

                                v3_prod = v3_tmp1 * v3_tmp2 * v3_tmp3 * v3_tmp4;

                                if (std::abs(omega1 - omega5) < eps) {

                                    for (ip1 = 1; ip1 >= -1; ip1 -= 2) {
                                        dp1 = static_cast<double>(ip1) * omega1;
                                        dp1_inv = 1.0 / dp1;

                                        for (ip2 = 1; ip2 >= -1; ip2 -= 2) {
                                            dp2 = static_cast<double>(ip2) * omega2;
                                            omega_sum[0] = 1.0 / (omega_shift + dp1 + dp2);

                                            for (ip3 = 1; ip3 >= -1; ip3 -= 2) {
                                                dp3 = static_cast<double>(ip3) * omega3;

                                                for (ip4 = 1; ip4 >= -1; ip4 -= 2) {
                                                    dp4 = static_cast<double>(ip4) * omega4;

                                                    D134 = 1.0 / (dp1 + dp3 + dp4);
                                                    omega_sum[1] = 1.0 / (omega_shift + dp2 + dp3 + dp4);

                                                    for (i = 0; i < N; ++i) {
                                                        T_tmp = T[i];

                                                        n1 = thermodynamics->fB(dp1, T_tmp);
                                                        n2 = thermodynamics->fB(dp2, T_tmp);
                                                        n3 = thermodynamics->fB(dp3, T_tmp);
                                                        n4 = thermodynamics->fB(dp4, T_tmp);

                                                        if (std::abs(T_tmp) < eps) {
                                                            T_inv = 0.0;
                                                        } else {
                                                            T_inv = 1.0 / (thermodynamics->T_to_Ryd * T_tmp);
                                                        }

                                                        ret_mpi[i]
                                                            += v3_prod * static_cast<double>(ip2 * ip3 * ip4)
                                                            * (omega_sum[1]
                                                                * (n2 * omega_sum[0]
                                                                    * ((1.0 + n3 + n4) * omega_sum[0] + (1.0 + n2 + n4) * dp1_inv)
                                                                    + (1.0 + n3) * (1.0 + n4) * D134 * (D134 + dp1_inv))
                                                                + (1.0 + n1) * (1.0 + n3 + n4) * D134
                                                                * omega_sum[0] * (omega_sum[0] + D134 + dp1_inv + n1 * T_inv));
                                                    }
                                                }
                                            }
                                        }
                                    }

                                } else {

                                    for (ip1 = 1; ip1 >= -1; ip1 -= 2) {
                                        dp1 = static_cast<double>(ip1) * omega1;

                                        for (ip5 = 1; ip5 >= -1; ip5 -= 2) {
                                            dp5 = static_cast<double>(ip5) * omega5;

                                            D15 = 1.0 / (dp1 - dp5);

                                            for (ip2 = 1; ip2 >= -1; ip2 -= 2) {
                                                dp2 = static_cast<double>(ip2) * omega2;

                                                omega_sum[0] = 1.0 / (omega_shift + dp1 + dp2);
                                                omega_sum[1] = 1.0 / (omega_shift + dp5 + dp2);

                                                for (ip3 = 1; ip3 >= -1; ip3 -= 2) {
                                                    dp3 = static_cast<double>(ip3) * omega3;

                                                    for (ip4 = 1; ip4 >= -1; ip4 -= 2) {
                                                        dp4 = static_cast<double>(ip4) * omega4;

                                                        D134 = 1.0 / (dp1 + dp3 + dp4);
                                                        D345 = 1.0 / (dp5 + dp3 + dp4);
                                                        omega_sum[2] = 1.0 / (omega_shift + dp2 + dp3 + dp4);

                                                        for (i = 0; i < N; ++i) {
                                                            T_tmp = T[i];

                                                            n1 = thermodynamics->fB(dp1, T_tmp);
                                                            n2 = thermodynamics->fB(dp2, T_tmp);
                                                            n3 = thermodynamics->fB(dp3, T_tmp);
                                                            n4 = thermodynamics->fB(dp4, T_tmp);
                                                            n5 = thermodynamics->fB(dp5, T_tmp);

                                                            ret_mpi[i]
                                                                += v3_prod * static_cast<double>(ip1 * ip2 * ip3 * ip4 * ip5)
                                                                * ((1.0 + n3 + n4)
                                                                    * (-(1.0 + n1 + n2) * D15 * D134 * omega_sum[0]
                                                                        + (1.0 + n5 + n2) * D15 * D345 * omega_sum[1])
                                                                    + (1.0 + n2 + n3 + n4 + n2 * n3 + n3 * n4 + n4 * n2)
                                                                    * D15 * (D345 - D134) * omega_sum[2]);
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    factor = 1.0 / (std::pow(static_cast<double>(nk), 2) * std::pow(2.0, 7));
    for (i = 0; i < N; ++i) ret_mpi[i] *= factor;

    mpi_reduce_complex(N, ret_mpi, ret);

    memory->deallocate(ret_mpi);
}


void Selfenergy::selfenergy_g(const unsigned int N,
                              double *T,
                              const double omega,
                              const unsigned int knum,
                              const unsigned int snum,
                              std::complex<double> *ret)
{
    /* 
    Diagram (g)
    Matrix elements that appear : V3^2 V4
    Computational cost          : O(N_k^2 * N^4)
    */

    unsigned int i;
    unsigned int ik1, ik2, ik3, ik4;
    unsigned int is1, is2, is3, is4;
    unsigned int iloc, jloc, kloc;

    unsigned int nkx = kpoint->nkx;
    unsigned int nky = kpoint->nky;
    unsigned int nkz = kpoint->nkz;

    unsigned int arr_quartic[4], arr_cubic1[3], arr_cubic2[3];

    int ip1, ip2, ip3, ip4;

    double omega1, omega2, omega3, omega4;
    double dp1, dp2, dp3, dp4;
    double n1, n2, n3, n4;
    double D124;

    double xk_tmp[3];
    double T_tmp;
    double factor;

    std::complex<double> omega_shift;
    std::complex<double> omega_sum[2];

    std::complex<double> v3_tmp1, v3_tmp2, v4_tmp;
    std::complex<double> v_prod;
    std::complex<double> *ret_mpi;

    memory->allocate(ret_mpi, N);

    omega_shift = omega + im * epsilon;

    for (i = 0; i < N; ++i) ret_mpi[i] = std::complex<double>(0.0, 0.0);

    arr_quartic[0] = ns * kpoint->knum_minus[knum] + snum;
    arr_cubic2[2] = ns * knum + snum;

    for (ik1 = mympi->my_rank; ik1 < nk; ik1 += mympi->nprocs) {

        for (ik2 = 0; ik2 < nk; ++ik2) {

            xk_tmp[0] = kpoint->xk[knum][0] - kpoint->xk[ik1][0] - kpoint->xk[ik2][0];
            xk_tmp[1] = kpoint->xk[knum][1] - kpoint->xk[ik1][1] - kpoint->xk[ik2][1];
            xk_tmp[2] = kpoint->xk[knum][2] - kpoint->xk[ik1][2] - kpoint->xk[ik2][2];

            iloc = (nint(xk_tmp[0] * static_cast<double>(nkx) + static_cast<double>(2 * nkx))) % nkx;
            jloc = (nint(xk_tmp[1] * static_cast<double>(nky) + static_cast<double>(2 * nky))) % nky;
            kloc = (nint(xk_tmp[2] * static_cast<double>(nkz) + static_cast<double>(2 * nkz))) % nkz;

            ik3 = kloc + nkz * jloc + nky * nkz * iloc;

            xk_tmp[0] = kpoint->xk[knum][0] - kpoint->xk[ik3][0];
            xk_tmp[1] = kpoint->xk[knum][1] - kpoint->xk[ik3][1];
            xk_tmp[2] = kpoint->xk[knum][2] - kpoint->xk[ik3][2];

            iloc = (nint(xk_tmp[0] * static_cast<double>(nkx) + static_cast<double>(2 * nkx))) % nkx;
            jloc = (nint(xk_tmp[1] * static_cast<double>(nky) + static_cast<double>(2 * nky))) % nky;
            kloc = (nint(xk_tmp[2] * static_cast<double>(nkz) + static_cast<double>(2 * nkz))) % nkz;

            ik4 = kloc + nkz * jloc + nky * nkz * iloc;

            for (is1 = 0; is1 < ns; ++is1) {
                omega1 = dynamical->eval_phonon[ik1][is1];

                arr_quartic[1] = ns * ik1 + is1;
                arr_cubic1[0] = ns * kpoint->knum_minus[ik1] + is1;

                for (is2 = 0; is2 < ns; ++is2) {
                    omega2 = dynamical->eval_phonon[ik2][is2];

                    arr_quartic[2] = ns * ik2 + is2;
                    arr_cubic1[1] = ns * kpoint->knum_minus[ik2] + is2;

                    for (is3 = 0; is3 < ns; ++is3) {
                        omega3 = dynamical->eval_phonon[ik3][is3];

                        arr_quartic[3] = ns * ik3 + is3;
                        arr_cubic2[0] = ns * kpoint->knum_minus[ik3] + is3;

                        v4_tmp = relaxation->V4(arr_quartic);

                        for (is4 = 0; is4 < ns; ++is4) {
                            omega4 = dynamical->eval_phonon[ik4][is4];

                            arr_cubic1[2] = ns * ik4 + is4;
                            arr_cubic2[1] = ns * kpoint->knum_minus[ik4] + is4;

                            v3_tmp1 = relaxation->V3(arr_cubic1);
                            v3_tmp2 = relaxation->V3(arr_cubic2);

                            v_prod = v4_tmp * v3_tmp1 * v3_tmp2;

                            for (ip1 = 1; ip1 >= -1; ip1 -= 2) {
                                dp1 = static_cast<double>(ip1) * omega1;
                                for (ip2 = 1; ip2 >= -1; ip2 -= 2) {
                                    dp2 = static_cast<double>(ip2) * omega2;
                                    for (ip3 = 1; ip3 >= -1; ip3 -= 2) {
                                        dp3 = static_cast<double>(ip3) * omega3;

                                        omega_sum[1] = 1.0 / (omega_shift + dp1 + dp2 + dp3);

                                        for (ip4 = 1; ip4 >= -1; ip4 -= 2) {
                                            dp4 = static_cast<double>(ip4) * omega4;

                                            omega_sum[0] = 1.0 / (omega_shift + dp3 + dp4);
                                            D124 = 1.0 / (dp1 + dp2 - dp4);

                                            for (i = 0; i < N; ++i) {
                                                T_tmp = T[i];

                                                n1 = thermodynamics->fB(dp1, T_tmp);
                                                n2 = thermodynamics->fB(dp2, T_tmp);
                                                n3 = thermodynamics->fB(dp3, T_tmp);
                                                n4 = thermodynamics->fB(dp4, T_tmp);

                                                ret_mpi[i]
                                                    += v_prod * static_cast<double>(ip1 * ip2 * ip3 * ip4) * D124
                                                    * ((1.0 + n1 + n2 + n3 + n4 + n1 * n3 + n1 * n4 + n2 * n3 + n2 * n4) * omega_sum[0]
                                                        - (1.0 + n1 + n2 + n3 + n1 * n2 + n2 * n3 + n1 * n3) * omega_sum[1]);

                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    factor = -1.0 / (std::pow(static_cast<double>(nk), 2) * std::pow(2.0, 6));
    for (i = 0; i < N; ++i) ret_mpi[i] *= factor;

    mpi_reduce_complex(N, ret_mpi, ret);

    memory->deallocate(ret_mpi);
}

void Selfenergy::selfenergy_h(const unsigned int N,
                              double *T,
                              const double omega,
                              const unsigned int knum,
                              const unsigned int snum,
                              std::complex<double> *ret)
{
    /*
    Diagram (h)
    Matrix elements that appear : V3^4
    Computational cost          : O(N_k^2 * N^5)
    Note                        : The most complicated diagram.
    */

    unsigned int i;
    unsigned int ik1, ik2, ik3, ik4, ik5;
    unsigned int is1, is2, is3, is4, is5;
    unsigned int arr_cubic1[3], arr_cubic2[3], arr_cubic3[3], arr_cubic4[3];
    unsigned int nkx = kpoint->nkx;
    unsigned int nky = kpoint->nky;
    unsigned int nkz = kpoint->nkz;

    unsigned int iloc, jloc, kloc;

    int ip1, ip2, ip3, ip4, ip5;

    double T_tmp;
    double xk_tmp[3];
    double factor;
    double omega1, omega2, omega3, omega4, omega5;
    double dp1, dp2, dp3, dp4, dp5;
    double n1, n2, n3, n4, n5;
    double D1, D2, D1_inv, D2_inv, D12_inv;
    double N12, N35, N34;
    double N_prod[4];

    std::complex<double> v3_tmp1, v3_tmp2, v3_tmp3, v3_tmp4;
    std::complex<double> v_prod;
    std::complex<double> omega_shift;
    std::complex<double> omega_sum[4];
    std::complex<double> *ret_mpi;

    memory->allocate(ret_mpi, N);

    omega_shift = omega + im * epsilon;

    for (i = 0; i < N; ++i) ret_mpi[i] = std::complex<double>(0.0, 0.0);

    arr_cubic1[0] = ns * kpoint->knum_minus[knum] + snum;
    arr_cubic4[2] = ns * knum + snum;

    for (ik1 = mympi->my_rank; ik1 < nk; ik1 += mympi->nprocs) {

        xk_tmp[0] = kpoint->xk[knum][0] - kpoint->xk[ik1][0];
        xk_tmp[1] = kpoint->xk[knum][1] - kpoint->xk[ik1][1];
        xk_tmp[2] = kpoint->xk[knum][2] - kpoint->xk[ik1][2];

        iloc = (nint(xk_tmp[0] * static_cast<double>(nkx) + static_cast<double>(2 * nkx))) % nkx;
        jloc = (nint(xk_tmp[1] * static_cast<double>(nky) + static_cast<double>(2 * nky))) % nky;
        kloc = (nint(xk_tmp[2] * static_cast<double>(nkz) + static_cast<double>(2 * nkz))) % nkz;

        ik2 = kloc + nkz * jloc + nky * nkz * iloc;

        for (ik3 = 0; ik3 < nk; ++ik3) {

            xk_tmp[0] = kpoint->xk[ik1][0] - kpoint->xk[ik3][0];
            xk_tmp[1] = kpoint->xk[ik1][1] - kpoint->xk[ik3][1];
            xk_tmp[2] = kpoint->xk[ik1][2] - kpoint->xk[ik3][2];

            iloc = (nint(xk_tmp[0] * static_cast<double>(nkx) + static_cast<double>(2 * nkx))) % nkx;
            jloc = (nint(xk_tmp[1] * static_cast<double>(nky) + static_cast<double>(2 * nky))) % nky;
            kloc = (nint(xk_tmp[2] * static_cast<double>(nkz) + static_cast<double>(2 * nkz))) % nkz;

            ik5 = kloc + nkz * jloc + nky * nkz * iloc;

            xk_tmp[0] = kpoint->xk[knum][0] - kpoint->xk[ik5][0];
            xk_tmp[1] = kpoint->xk[knum][1] - kpoint->xk[ik5][1];
            xk_tmp[2] = kpoint->xk[knum][2] - kpoint->xk[ik5][2];

            iloc = (nint(xk_tmp[0] * static_cast<double>(nkx) + static_cast<double>(2 * nkx))) % nkx;
            jloc = (nint(xk_tmp[1] * static_cast<double>(nky) + static_cast<double>(2 * nky))) % nky;
            kloc = (nint(xk_tmp[2] * static_cast<double>(nkz) + static_cast<double>(2 * nkz))) % nkz;

            ik4 = kloc + nkz * jloc + nky * nkz * iloc;


            for (is1 = 0; is1 < ns; ++is1) {
                omega1 = dynamical->eval_phonon[ik1][is1];

                arr_cubic1[1] = ns * ik1 + is1;
                arr_cubic2[0] = ns * kpoint->knum_minus[ik1] + is1;

                for (is2 = 0; is2 < ns; ++is2) {
                    omega2 = dynamical->eval_phonon[ik2][is2];

                    arr_cubic1[2] = ns * ik2 + is2;
                    arr_cubic3[0] = ns * kpoint->knum_minus[ik2] + is2;

                    v3_tmp1 = relaxation->V3(arr_cubic1);

                    for (is3 = 0; is3 < ns; ++ is3) {
                        omega3 = dynamical->eval_phonon[ik3][is3];

                        arr_cubic2[1] = ns * ik3 + is3;
                        arr_cubic3[1] = ns * kpoint->knum_minus[ik3] + is3;

                        for (is4 = 0; is4 < ns; ++is4) {
                            omega4 = dynamical->eval_phonon[ik4][is4];

                            arr_cubic3[2] = ns * ik4 + is4;
                            arr_cubic4[0] = ns * kpoint->knum_minus[ik4] + is4;

                            v3_tmp3 = relaxation->V3(arr_cubic3);

                            for (is5 = 0; is5 < ns; ++is5) {
                                omega5 = dynamical->eval_phonon[ik5][is5];

                                arr_cubic2[2] = ns * ik5 + is5;
                                arr_cubic4[1] = ns * kpoint->knum_minus[ik5] + is5;

                                v3_tmp2 = relaxation->V3(arr_cubic2);
                                v3_tmp4 = relaxation->V3(arr_cubic4);

                                v_prod = v3_tmp1 * v3_tmp2 * v3_tmp3 * v3_tmp4;

                                for (ip1 = 1; ip1 >= -1; ip1 -= 2) {
                                    dp1 = static_cast<double>(ip1) * omega1;

                                    for (ip2 = 1; ip2 >= -1; ip2 -= 2) {
                                        dp2 = static_cast<double>(ip2) * omega2;
                                        omega_sum[0] = 1.0 / (omega_shift + dp1 - dp2);

                                        for (ip3 = 1; ip3 >= -1; ip3 -= 2) {
                                            dp3 = static_cast<double>(ip3) * omega3;

                                            for (ip4 = 1; ip4 >= -1; ip4 -= 2) {
                                                dp4 = static_cast<double>(ip4) * omega4;

                                                D2 = dp4 - dp3 - dp2;
                                                D2_inv = 1.0 / D2;
                                                omega_sum[3] = 1.0 / (omega_shift + dp1 + dp3 - dp4);

                                                for (ip5 = 1; ip5 >= -1; ip5 -= 2) {
                                                    dp5 = static_cast<double>(ip5) * omega5;

                                                    D1 = dp5 - dp3 - dp1;
                                                    D1_inv = 1.0 / D1;
                                                    D12_inv = D1_inv * D2_inv;

                                                    omega_sum[1] = 1.0 / (omega_shift - dp4 + dp5);
                                                    omega_sum[2] = 1.0 / (omega_shift - dp2 - dp3 + dp5);

                                                    for (i = 0; i < N; ++i) {
                                                        T_tmp = T[i];

                                                        n1 = thermodynamics->fB(dp1, T_tmp);
                                                        n2 = thermodynamics->fB(dp2, T_tmp);
                                                        n3 = thermodynamics->fB(dp3, T_tmp);
                                                        n4 = thermodynamics->fB(dp4, T_tmp);
                                                        n5 = thermodynamics->fB(dp5, T_tmp);

                                                        N12 = n1 - n2;
                                                        N34 = n3 - n4;
                                                        N35 = n3 - n5;

                                                        N_prod[0] = N12 * (1.0 + n3);
                                                        N_prod[1] = (1.0 + n2 + n3) * (1.0 + n5) - (1.0 + n1 + n3) * (1.0 + n4);
                                                        N_prod[2] = ((1.0 + n2) * N35 - n3 * (1.0 + n5));
                                                        N_prod[3] = -((1.0 + n1) * N34 - n3 * (1.0 + n4));

                                                        ret_mpi[i]
                                                            += v_prod * static_cast<double>(ip1 * ip2 * ip3 * ip4 * ip5)
                                                            * (D12_inv
                                                                * (N_prod[0] * omega_sum[0]
                                                                    + N_prod[1] * omega_sum[1]
                                                                    + N_prod[2] * omega_sum[2]
                                                                    + N_prod[3] * omega_sum[3])
                                                                + N12 * ((1.0 + n5) * D1_inv - (1.0 + n4) * D2_inv)
                                                                * omega_sum[0] * omega_sum[1]);
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    factor = 1.0 / (std::pow(static_cast<double>(nk), 2) * std::pow(2.0, 7));
    for (i = 0; i < N; ++i) ret_mpi[i] *= factor;

    mpi_reduce_complex(N, ret_mpi, ret);

    memory->deallocate(ret_mpi);
}


void Selfenergy::selfenergy_i(const unsigned int N,
                              double *T,
                              const double omega,
                              const unsigned int knum,
                              const unsigned int snum,
                              std::complex<double> *ret)
{
    /* 

    Diagram (i)
    Matrix elements that appear : V3^2 V4
    Computational cost          : O(N_k^2 * N^4)
    Note                        : Double pole when omega2 = omega4. 
    : No frequency dependence.

    */

    unsigned int i;
    unsigned int ik1, ik2, ik3, ik4;
    unsigned int is1, is2, is3, is4;
    unsigned int arr_quartic[4];
    unsigned int arr_cubic1[3], arr_cubic2[3];
    unsigned int iloc, jloc, kloc;
    unsigned int nkx = kpoint->nkx;
    unsigned int nky = kpoint->nky;
    unsigned int nkz = kpoint->nkz;

    int ip1, ip2, ip3, ip4;

    double omega1, omega2, omega3, omega4;
    double n1, n2, n3, n4;
    double dp1, dp2, dp3, dp4;
    double D24, D123, D134;
    double dp2_inv;
    double T_tmp;
    double factor;
    double xk_tmp[3];
    double N_prod[2];
    double T_inv;

    std::complex<double> v4_tmp, v3_tmp1, v3_tmp2;
    std::complex<double> v_prod;
    std::complex<double> *ret_mpi;

    memory->allocate(ret_mpi, N);

    for (i = 0; i < N; ++i) ret_mpi[i] = std::complex<double>(0.0, 0.0);

    arr_quartic[0] = ns * kpoint->knum_minus[knum] + snum;
    arr_quartic[3] = ns * knum + snum;

    for (ik1 = mympi->my_rank; ik1 < nk; ik1 += mympi->nprocs) {
        for (ik2 = 0; ik2 < nk; ++ik2) {

            ik4 = ik2;
            xk_tmp[0] = kpoint->xk[ik2][0] - kpoint->xk[ik1][0];
            xk_tmp[1] = kpoint->xk[ik2][1] - kpoint->xk[ik1][1];
            xk_tmp[2] = kpoint->xk[ik2][2] - kpoint->xk[ik1][2];

            iloc = (nint(xk_tmp[0] * static_cast<double>(nkx) + static_cast<double>(2 * nkx))) % nkx;
            jloc = (nint(xk_tmp[1] * static_cast<double>(nky) + static_cast<double>(2 * nky))) % nky;
            kloc = (nint(xk_tmp[2] * static_cast<double>(nkz) + static_cast<double>(2 * nkz))) % nkz;

            ik3 = kloc + nkz * jloc + nky * nkz * iloc;

            for (is2 = 0; is2 < ns; ++is2) {
                omega2 = dynamical->eval_phonon[ik2][is2];

                arr_quartic[1] = ns * ik2 + is2;
                arr_cubic2[0] = ns * kpoint->knum_minus[ik2] + is2;

                for (is4 = 0; is4 < ns; ++is4) {
                    omega4 = dynamical->eval_phonon[ik4][is4];

                    arr_quartic[2] = ns * kpoint->knum_minus[ik4] + is4;
                    arr_cubic1[2] = ns * ik4 + is4;

                    v4_tmp = relaxation->V4(arr_quartic);

                    if (std::abs(omega2 - omega4) < eps) {

                        for (is3 = 0; is3 < ns; ++is3) {
                            omega3 = dynamical->eval_phonon[ik3][is3];

                            arr_cubic1[1] = ns * kpoint->knum_minus[ik3] + is3;
                            arr_cubic2[2] = ns * ik3 + is3;

                            for (is1 = 0; is1 < ns; ++is1) {
                                omega1 = dynamical->eval_phonon[ik1][is1];

                                arr_cubic1[0] = ns * kpoint->knum_minus[ik1] + is1;
                                arr_cubic2[1] = ns * ik1 + is1;

                                v3_tmp1 = relaxation->V3(arr_cubic1);
                                v3_tmp2 = relaxation->V3(arr_cubic2);

                                v_prod = v4_tmp * v3_tmp1 * v3_tmp2;

                                for (ip1 = 1; ip1 >= -1; ip1 -= 2) {
                                    dp1 = static_cast<double>(ip1) * omega1;

                                    for (ip2 = 1; ip2 >= -1; ip2 -= 2) {
                                        dp2 = static_cast<double>(ip2) * omega2;

                                        dp2_inv = 1.0 / dp2;

                                        for (ip3 = 1; ip3 >= -1; ip3 -= 2) {
                                            dp3 = static_cast<double>(ip3) * omega3;

                                            D123 = 1.0 / (dp1 + dp2 + dp3);

                                            for (i = 0; i < N; ++i) {
                                                T_tmp = T[i];

                                                n1 = thermodynamics->fB(dp1, T_tmp);
                                                n2 = thermodynamics->fB(dp2, T_tmp);
                                                n3 = thermodynamics->fB(dp3, T_tmp);

                                                N_prod[0] = (1.0 + n1) * (1.0 + n3) + n2 * (1.0 + n2 + n3);
                                                N_prod[1] = n2 * (1.0 + n2) * (1.0 + n2 + n3);

                                                if (std::abs(T_tmp) < eps) {
                                                    T_inv = 0.0;
                                                } else {
                                                    T_inv = 1.0 / (thermodynamics->T_to_Ryd * T_tmp);
                                                }

                                                ret_mpi[i]
                                                    += v_prod * static_cast<double>(ip1 * ip3)
                                                    * (D123 * (N_prod[0] * D123 + N_prod[1] * T_inv + N_prod[0] * dp2_inv));
                                            }
                                        }
                                    }
                                }
                            }
                        }

                    } else {
                        for (is3 = 0; is3 < ns; ++is3) {
                            omega3 = dynamical->eval_phonon[ik3][is3];

                            arr_cubic1[1] = ns * kpoint->knum_minus[ik3] + is3;
                            arr_cubic2[2] = ns * ik3 + is3;

                            for (is1 = 0; is1 < ns; ++is1) {
                                omega1 = dynamical->eval_phonon[ik1][is1];

                                arr_cubic1[0] = ns * kpoint->knum_minus[ik1] + is1;
                                arr_cubic2[1] = ns * ik1 + is1;

                                v3_tmp1 = relaxation->V3(arr_cubic1);
                                v3_tmp2 = relaxation->V3(arr_cubic2);

                                v_prod = v4_tmp * v3_tmp1 * v3_tmp2;

                                for (ip1 = 1; ip1 >= -1; ip1 -= 2) {
                                    dp1 = static_cast<double>(ip1) * omega1;

                                    for (ip2 = 1; ip2 >= -1; ip2 -= 2) {
                                        dp2 = static_cast<double>(ip2) * omega2;

                                        for (ip3 = 1; ip3 >= -1; ip3 -= 2) {

                                            dp3 = static_cast<double>(ip3) * omega3;
                                            D123 = 1.0 / (dp1 - dp2 + dp3);

                                            for (ip4 = 1; ip4 >= -1; ip4 -= 2) {
                                                dp4 = static_cast<double>(ip4) * omega4;

                                                D24 = 1.0 / (dp2 - dp4);
                                                D134 = 1.0 / (dp1 + dp3 - dp4);

                                                for (i = 0; i < N; ++i) {
                                                    T_tmp = T[i];

                                                    n1 = thermodynamics->fB(dp1, T_tmp);
                                                    n2 = thermodynamics->fB(dp2, T_tmp);
                                                    n3 = thermodynamics->fB(dp3, T_tmp);
                                                    n4 = thermodynamics->fB(dp4, T_tmp);

                                                    ret_mpi[i]
                                                        += v_prod * static_cast<double>(ip1 * ip2 * ip3 * ip4)
                                                        * ((1.0 + n1 + n3) * D24 * (n4 * D134 - n2 * D123)
                                                            + D123 * D134 * n1 * n3);
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    factor = -1.0 / (std::pow(static_cast<double>(nk), 2) * std::pow(2.0, 7));
    for (i = 0; i < N; ++i) ret_mpi[i] *= factor;

    mpi_reduce_complex(N, ret_mpi, ret);

    memory->deallocate(ret_mpi);
}


void Selfenergy::selfenergy_j(const unsigned int N,
                              double *T,
                              const double omega,
                              const unsigned int knum,
                              const unsigned int snum,
                              std::complex<double> *ret)
{
    /*

    Diagram (j)
    Matrix elements that appear : V4^2
    Computational cost          : O(N_k^2 * N^3)
    Note                        : Double pole when omega1 = omega3

    */

    unsigned int i;
    unsigned int ik1, ik2, ik3;
    unsigned int is1, is2, is3;
    unsigned int arr_quartic1[4], arr_quartic2[4];
    unsigned int nkx = kpoint->nkx;
    unsigned int nky = kpoint->nky;
    unsigned int nkz = kpoint->nkz;

    double T_tmp;
    double factor;
    double n1, n2, n3;
    double omega1, omega2, omega3;
    double omega1_inv;
    double D13[2];
    double T_inv;

    std::complex<double> v4_tmp1, v4_tmp2;
    std::complex<double> v_prod;
    std::complex<double> *ret_mpi;

    memory->allocate(ret_mpi, N);

    for (i = 0; i < N; ++i) ret_mpi[i] = std::complex<double>(0.0, 0.0);

    arr_quartic1[0] = ns * kpoint->knum_minus[knum] + snum;
    arr_quartic1[3] = ns * knum + snum;

    for (ik1 = mympi->my_rank; ik1 < nk; ik1 += mympi->nprocs) {

        ik3 = ik1;

        for (ik2 = 0; ik2 < nk; ++ik2) {

            for (is1 = 0; is1 < ns; ++is1) {
                omega1 = dynamical->eval_phonon[ik1][is1];

                arr_quartic1[1] = ns * ik1 + is1;
                arr_quartic2[0] = ns * kpoint->knum_minus[ik1] + is1;


                for (is3 = 0; is3 < ns; ++is3) {
                    omega3 = dynamical->eval_phonon[ik1][is3];

                    arr_quartic1[2] = ns * kpoint->knum_minus[ik3] + is3;
                    arr_quartic2[3] = ns * ik3 + is3;

                    v4_tmp1 = relaxation->V4(arr_quartic1);

                    if (std::abs(omega1 - omega3) < eps) {
                        omega1_inv = 1.0 / omega1;

                        for (is2 = 0; is2 < ns; ++is2) {
                            omega2 = dynamical->eval_phonon[ik2][is2];

                            arr_quartic2[1] = ns * ik2 + is2;
                            arr_quartic2[2] = ns * kpoint->knum_minus[ik2] + is2;

                            v4_tmp2 = relaxation->V4(arr_quartic2);

                            v_prod = v4_tmp1 * v4_tmp2;

                            for (i = 0; i < N; ++i) {
                                T_tmp = T[i];

                                n1 = thermodynamics->fB(omega1, T_tmp);
                                n2 = thermodynamics->fB(omega2, T_tmp);

                                if (std::abs(T_tmp) < eps) {
                                    T_inv = 0.0;
                                } else {
                                    T_inv = 1.0 / (thermodynamics->T_to_Ryd * T_tmp);
                                }

                                ret_mpi[i]
                                    += v_prod * (2.0 * n2 + 1.0)
                                    * (-2.0 * (1.0 + n1) * n1 * T_inv
                                        - (2.0 * n1 + 1.0) * omega1_inv);
                            }
                        }
                    } else {

                        D13[0] = 1.0 / (omega1 - omega3);
                        D13[1] = 1.0 / (omega1 + omega3);

                        for (is2 = 0; is2 < ns; ++is2) {
                            omega2 = dynamical->eval_phonon[ik2][is2];

                            arr_quartic2[1] = ns * ik2 + is2;
                            arr_quartic2[2] = ns * kpoint->knum_minus[ik2] + is2;

                            v4_tmp2 = relaxation->V4(arr_quartic2);

                            v_prod = v4_tmp1 * v4_tmp2;

                            for (i = 0; i < N; ++i) {
                                T_tmp = T[i];

                                n1 = thermodynamics->fB(omega1, T_tmp);
                                n2 = thermodynamics->fB(omega2, T_tmp);
                                n3 = thermodynamics->fB(omega3, T_tmp);

                                ret_mpi[i]
                                    += v_prod * 2.0
                                    * ((n1 - n3) * D13[0] - (1.0 + n1 + n3) * D13[1]);
                            }
                        }
                    }
                }
            }
        }
    }

    factor = -1.0 / (std::pow(static_cast<double>(nk), 2) * std::pow(2.0, 6));

    for (i = 0; i < N; ++i) ret_mpi[i] *= factor;

    mpi_reduce_complex(N, ret_mpi, ret);

    memory->deallocate(ret_mpi);
}
