/*
 phonon_velocity.cpp

 Copyright (c) 2014 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#include "mpi_common.h"
#include <complex>
#include <iomanip>
#include "fcs_phonon.h"
#include "phonon_velocity.h"
#include "kpoint.h"
#include "memory.h"
#include "dynamical.h"
#include "system.h"
#include "error.h"
#include "write_phonons.h"
#include "constants.h"
#include "mathfunctions.h"

using namespace PHON_NS;

Phonon_velocity::Phonon_velocity(PHON *phon): Pointers(phon){}

Phonon_velocity::~Phonon_velocity(){
    if (print_velocity) {
        memory->deallocate(phvel);
    }
}

void Phonon_velocity::setup_velocity()
{
    MPI_Bcast(&print_velocity, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD);
}

void Phonon_velocity::calc_group_velocity(const int kpmode)
{
    if (print_velocity) {

        unsigned int nk = kpoint->nk;
        unsigned int ns = dynamical->neval;

        memory->allocate(phvel, nk, ns);

        if (kpmode == 1) {
        
            calc_phonon_vel_band(phvel);

        } else if (kpmode == 2) {

            calc_phonon_vel_mesh(phvel);
        }

    }
}

void Phonon_velocity::calc_phonon_vel_band(double **phvel_out)
{
    unsigned int i;
    unsigned int ik, idiff;
    unsigned int nk = kpoint->nk;
    unsigned int n = dynamical->neval;
    unsigned int ndiff;
    double **xk_shift;
    double *xk_tmp;
    double **omega_shift, *omega_tmp;

    double h = 1.0e-4;

    std::complex<double> **evec_tmp;

    memory->allocate(evec_tmp, 1, 1);

    if (mympi->my_rank == 0) {
        std::cout << " Calculating group velocities of phonon along given k path ... ";
    }

    ndiff = 2;
    memory->allocate(xk_shift, ndiff, 3);
    memory->allocate(omega_shift, ndiff, n);
    memory->allocate(omega_tmp, ndiff);

    memory->allocate(xk_tmp, 3);


    for (ik = 0; ik < nk; ++ik){

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
            error->exit("calc_phonon_vel_band", "ndiff > 2 is not supported yet.");
        }

        for (idiff = 0; idiff < ndiff; ++idiff){

            // Move back to fractional basis

            rotvec(xk_shift[idiff], xk_shift[idiff], system->lavec_p, 'T');
            for (i = 0; i < 3; ++i) xk_shift[idiff][i] /= 2.0 * pi;

                dynamical->eval_k(xk_shift[idiff], kpoint->kvec_na[ik], fcs_phonon->fc2_ext, omega_shift[idiff], evec_tmp, false);

        }

        for (i = 0; i < n; ++i){
            for (idiff = 0; idiff < ndiff; ++idiff){
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

void Phonon_velocity::calc_phonon_vel_mesh(double **phvel_out)
{
    unsigned int i, j, k;
    unsigned int nk = kpoint->nk;
    unsigned int ns = dynamical->neval;
    double **vel;

    if (mympi->my_rank == 0) {
        std::cout << " Calculating group velocities of phonons at uniform grid ... ";
    }

    memory->allocate(vel, ns, 3);

    for (i = 0; i < nk; ++i) {
        phonon_vel_k(kpoint->xk[i], vel);

        for (j = 0; j < ns; ++j){
            rotvec(vel[j], vel[j], system->lavec_p, 'T');
            for (k = 0; k < 3; ++k) vel[j][k] /= 2.0 * pi;
            phvel_out[i][j] = std::sqrt(std::pow(vel[j][0], 2) + std::pow(vel[j][1], 2) + std::pow(vel[j][2], 2));
        }
    }

    memory->deallocate(vel);

    if (mympi->my_rank == 0) {
        std::cout << "done!" << std::endl;
    }
}

void Phonon_velocity::phonon_vel_k(double *xk_in, double **vel_out)
{
    unsigned int i, j;
    unsigned int idiff;
    unsigned int ndiff;
    unsigned int n = dynamical->neval;
    double **xk_shift;
    std::complex<double> **evec_tmp; 
    double **omega_shift, *omega_tmp;
    double **kvec_na_tmp;
    double norm;
    double h = 1.0e-4;

    ndiff = 2;

    memory->allocate(omega_shift, ndiff, n); 
    memory->allocate(xk_shift, ndiff, 3);
    memory->allocate(omega_tmp, ndiff);
    memory->allocate(evec_tmp, 1, 1);
    memory->allocate(kvec_na_tmp, 2, 3);

    for (i = 0; i < 3; ++i){

        for (j = 0; j < 3; ++j){
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

        norm = std::sqrt(kvec_na_tmp[0][0] * kvec_na_tmp[0][0] + kvec_na_tmp[0][1] * kvec_na_tmp[0][1] + kvec_na_tmp[0][2] * kvec_na_tmp[0][2]);
        if (norm > eps) {
            for (j = 0; j < 3; ++j) kvec_na_tmp[0][j] /= norm;
        }
        norm = std::sqrt(kvec_na_tmp[1][0] * kvec_na_tmp[1][0] + kvec_na_tmp[1][1] * kvec_na_tmp[1][1] + kvec_na_tmp[1][2] * kvec_na_tmp[1][2]);
        if (norm > eps) {
            for (j = 0; j < 3; ++j) kvec_na_tmp[1][j] /= norm;
        }

        for (idiff = 0; idiff < ndiff; ++idiff) {
           
                dynamical->eval_k(xk_shift[idiff], kvec_na_tmp[0], fcs_phonon->fc2_ext, omega_shift[idiff], evec_tmp, false);
           
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

double Phonon_velocity::diff(double *f, const unsigned int n, double h)
{
    double df;

    if (n == 2) {
        df = (f[1] - f[0]) / (2.0 * h);
    } else {
        error->exit("diff", "Numerical differentiation of n > 2 is not supported yet.");
    }

    return df;
}
