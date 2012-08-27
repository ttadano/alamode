#include "phonon_velocity.h"
#include "kpoint.h"
#include "memory.h"
#include "dynamical.h"
#include "system.h"
#include "error.h"
#include "write_phonons.h"
#include "../alm_c++/constants.h"
#include <complex>
#include <iomanip>


using namespace PHON_NS;

Phonon_velocity::Phonon_velocity(PHON *phon): Pointers(phon){}

Phonon_velocity::~Phonon_velocity(){}

void Phonon_velocity::calc_phonon_vel_band()
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

    std::cout << "Calculating group velocities of phonon along given k-path ..." << std::endl;

    memory->allocate(phvel, nk, n);

    ndiff = 2;
    memory->allocate(xk_shift, ndiff, 3);
    memory->allocate(omega_shift, ndiff, n);
    memory->allocate(omega_tmp, ndiff);

    memory->allocate(xk_tmp, 3);


    for (ik = 0; ik < nk; ++ik){

        // Represent the given kpoint in cartesian coordinate
        system->rotvec(xk_tmp, kpoint->xk[ik], system->rlavec_p, 'T');

        if (ndiff == 2) {
            // central difference
            // f'(x) =~ f(x+h)-f(x-h)/2h

            for (i = 0; i < 3; ++i) {
                xk_shift[0][i] = xk_tmp[i] - h * kpoint->kpoint_direction[ik][i];
                xk_shift[1][i] = xk_tmp[i] + h * kpoint->kpoint_direction[ik][i];
            }

        } else {
            error->exit("calc_phonon_vel_band", "ndiff > 2 is not supported yet.");
        }

        for (idiff = 0; idiff < ndiff; ++idiff){

            // Move back to fractional basis

            system->rotvec(xk_shift[idiff], xk_shift[idiff], system->lavec_p, 'T');
            for (i = 0; i < 3; ++i) xk_shift[idiff][i] /= 2.0 * pi;

            dynamical->eval_k(xk_shift[idiff], omega_shift[idiff], evec_tmp, false); 
        }

        for (i = 0; i < n; ++i){
            for (idiff = 0; idiff < ndiff; ++idiff){
                omega_tmp[idiff] = freq(omega_shift[idiff][i]);
            }
            phvel[ik][i] = diff(omega_tmp, ndiff, h);
        }
    }
    memory->deallocate(omega_tmp);
    memory->deallocate(omega_shift);
    memory->deallocate(xk_shift);
    memory->deallocate(xk_tmp);

    memory->deallocate(evec_tmp);

    std::cout << "..done!" << std::endl;
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
    double h = 1.0e-4;

    ndiff = 2;

    memory->allocate(omega_shift, ndiff, n); 
    memory->allocate(xk_shift, ndiff, 3);
    memory->allocate(omega_tmp, ndiff);

    for (i = 0; i < 3; ++i){

        for (j = 0; j < 3; ++j){
            xk_shift[0][j] = xk_in[j];
            xk_shift[1][j] = xk_in[j];
        }

        xk_shift[0][i] -= h;
        xk_shift[1][i] += h;

        for (idiff = 0; idiff < ndiff; ++idiff) {
            dynamical->eval_k(xk_shift[idiff], omega_shift[idiff], evec_tmp, false);
        }

        for (j = 0; j < n; ++j){
            for (idiff = 0; idiff < ndiff; ++idiff){
                omega_tmp[idiff] = freq(omega_shift[idiff][j]);
            }
            vel_out[j][i] = diff(omega_tmp, ndiff, h);
        }
    }

    memory->deallocate(xk_shift);
    memory->deallocate(omega_shift);
    memory->deallocate(omega_tmp);
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

double Phonon_velocity::freq(const double x) 
{
    if (x >= 0.0) {
        return std::sqrt(x);
    } else {
        return std::sqrt(-x);
    }
}
