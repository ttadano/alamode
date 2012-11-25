#include "mpi_common.h"
#include <iostream>
#include <iomanip>
#include "dynamical.h"
#include "fcs_phonon.h"
#include "error.h"
#include "kpoint.h"
#include "memory.h"
#include "phonon_thermodynamics.h"
#include "phonon_velocity.h"
#include "pointers.h"
#include "system.h"
#include "../alm_c++/constants.h"

using namespace PHON_NS;

Phonon_thermodynamics::Phonon_thermodynamics(PHON *phon): Pointers(phon) {
    T_to_Ryd = k_Boltzmann / Ryd;
}

Phonon_thermodynamics::~Phonon_thermodynamics(){};

double Phonon_thermodynamics::Cv(const double omega, const double T)
{
    double x;

    if (std::abs(T) < eps || omega == 0.0) {
        return 0.0;
    } else {
        x = omega / (T_to_Ryd * T);
        return k_Boltzmann * std::pow(x/(2.0 * sinh(0.5*x)), 2);
    }
}

double Phonon_thermodynamics::Cv_classical(const double omega, const double T)
{
    double x;

    if (std::abs(T) < eps || omega == 0.0) {
        return 0.0;
    } else {
        x = omega / (T_to_Ryd * T);
        return k_Boltzmann * std::pow(x, 2) * std::exp(-x);
    }
}

double Phonon_thermodynamics::fB(const double omega, const double T)
{
    double x;

    if (std::abs(T) < eps){
        return 0.0;
    } else {
        x = omega / (T_to_Ryd * T);
        return 1.0 / (std::exp(x) - 1.0);
    }
}

double Phonon_thermodynamics::fC(const double omega, const double T)
{
    double x;

    if (std::abs(T) < eps) {
        return 0.0;
    } else {
        x = omega / (T_to_Ryd * T);
        return std::exp(-x);
    }
}

void Phonon_thermodynamics::test_fB(const double T)
{

    unsigned int i, j;
    unsigned int nk = kpoint->nk;
    unsigned int ns = dynamical->neval;

    double omega;

    for (i = 0; i < nk; ++i){
        for (j = 0; j < ns; ++j){
            omega = dynamical->eval_phonon[i][j];
            std::cout << "omega = " << omega / (T_to_Ryd * T) << " ,fB = " << fB(omega, T) << std::endl;
        }
    }
}

double Phonon_thermodynamics::Cv_tot(const double T)
{
    unsigned int ik, is;
    unsigned int nk = kpoint->nk;
    unsigned int ns = dynamical->neval;
    double omega;

    double ret = 0.0;

    for (ik = 0; ik < nk; ++ik){
        for (is = 0; is < ns; ++is){
            omega = dynamical->eval_phonon[ik][is];
            if (omega < 0.0) continue;
            ret += Cv(omega, T);
        }
    }
    return ret / static_cast<double>(nk);
}

double Phonon_thermodynamics::Cv_Debye(const double T, const double TD)
{
    unsigned int natmin = system->natmin;
    unsigned int i;
    double d_theta, theta, theta_max;
    unsigned int ntheta;

    double x, y;
    double ret;
    d_theta = 0.001;

    if (TD < eps)  {
        error->exit("Cv_Debye", "Too small TD");
    }
    if (T < 0.0) {
        error->exit("Cv_Debye", "Negative T");
    }

    if (T < eps) {
        return 0.0;
    } else {
        x = TD / T;
        theta_max = atan(x);

        ret = 0.0;
        ntheta = static_cast<unsigned int>(theta_max/d_theta);

        for (i = 0; i < ntheta; ++i){
            theta = static_cast<double>(i) * d_theta;
            y = tan(theta);
            if (y > eps){
                ret += std::pow(y, 4)*std::exp(y) / std::pow((std::exp(y) - 1.0) * cos(theta), 2);
            }
        }
        y = tan(theta_max);
        ret += 0.5 * std::pow(y, 4)*std::exp(y) / std::pow((std::exp(y) - 1.0) * cos(theta_max), 2);

        return 9.0 * static_cast<double>(natmin) * k_Boltzmann * ret * d_theta/ std::pow(x, 3);
    }
}

void Phonon_thermodynamics::Debye_T(const double T, double &TD)
{    
    double TD_old;
    double diff_C;
    double fdegfree = 1.0 / static_cast<double>(3.0 * system->natmin);

    if (T > eps) {

        do {   
            //       std::cout << "T = " << T << " , TD = " << TD << std::endl;
            diff_C = fdegfree * (Cv_tot(T) - Cv_Debye(T, TD)) / k_Boltzmann;

            TD_old = TD;
            TD = TD_old - diff_C * 10.0;
        } while(std::abs(diff_C) > 1.0e-5);
    }
}

double Phonon_thermodynamics::Internal_Energy(const double T)
{
    unsigned int ik, is;
    unsigned int nk = kpoint->nk;
    unsigned int ns = dynamical->neval;
    double omega;

    double ret = 0.0;

    for (ik = 0; ik < nk; ++ik){
        for (is = 0; is < ns; ++is){
            omega = dynamical->eval_phonon[ik][is];
            if (omega <= 0.0) {
                // exactly zero is anomalous
                continue;
            } else {
                ret += omega * coth_T(omega, T);
            }
        }
    }
    return ret / static_cast<double>(nk);
}

double Phonon_thermodynamics::coth_T(const double omega, const double T)
{
    if (T < eps) {
        // if T = 0.0 and omega > 0, coth(hbar*omega/(2*kB*T)) = 1.0
        return 1.0;
    } else {
        double x = 0.5 * omega / (T_to_Ryd * T);
        return 1.0 + 2.0 / (std::exp(2.0 * x) - 1.0);
    }
}

void Phonon_thermodynamics::calc_gruneisen()
{
    unsigned int is, ik;
    unsigned int i, j;
    unsigned int ns = dynamical->neval;
    unsigned int nk = kpoint->nk;

    memory->allocate(gruneisen, nk, ns);

    std::complex<double> **dfc2_reciprocal;

    memory->allocate(dfc2_reciprocal, ns, ns);

    for (ik = 0; ik < nk; ++ik){
        for (is = 0; is < ns; ++is){

            calc_dfc2_reciprocal(dfc2_reciprocal, kpoint->xk[ik]);

            gruneisen[ik][is] = std::complex<double>(0.0, 0.0);

       /*     for (i = 0; i < ns; ++i){
                for (j = 0; j < ns; ++j){
                    std::cout << " " << dfc2_reciprocal[i][j];
                }
                std::cout << std::endl;
            }
            std::cout << std::endl;*/

            for (i = 0; i < ns; ++i){
                for (j = 0; j < ns; ++j){
                    gruneisen[ik][is] += std::conj(dynamical->evec_phonon[ik][is][i]) * dfc2_reciprocal[i][j] * dynamical->evec_phonon[ik][is][j];
                    // gruneisen[ik][is] += std::conj(dynamical->evec_phonon[ik][is][i])  * dynamical->evec_phonon[ik][is][j];
                }
            }

          /*  std::cout << "ik = " << ik << " , is = " << is << std::endl;
            std::cout << "gamma = " << gruneisen[ik][is] << std::endl;
*/

            //     gruneisen[ik][is] /= - 6.0 * std::pow(dynamical->evec_phonon[ik][is]);
        }
    }

}

void Phonon_thermodynamics::calc_dfc2_reciprocal(std::complex<double> **dphi2, double *xk_in)
{
    unsigned int i, j;
    unsigned int icrd, jcrd, itran;
    unsigned int ns = dynamical->neval;
    unsigned int natmin = system->natmin;
    unsigned int ntran = system->ntran;

    unsigned int atm_p1, atm_p2, atm_s2;

    double phase;
    double vec[3];

    std::complex<double> exp_phase;
    std::complex<double> im(0.0, 1.0);
    std::complex<double> ctmp[3][3];


    for (i = 0; i < ns; ++i){
        for (j = 0; j < ns; ++j){
            dphi2[i][j] = std::complex<double>(0.0, 0.0);
        }
    }

    for (i = 0; i < natmin; ++i){

        atm_p1 = system->map_p2s[i][0];

        for (j = 0; j < natmin; ++j){

            atm_p2 = system->map_p2s[j][0];

            for (icrd = 0; icrd < 3; ++icrd){
                for (jcrd = 0; jcrd < 3; ++jcrd){
                    ctmp[icrd][jcrd] = std::complex<double>(0.0, 0.0);
                }
            }

            for (itran = 0; itran < ntran; ++itran){

                atm_s2 = system->map_p2s[j][itran];

                //  std::cout << "atm_p1 = " << atm_p1 << " , atm_p2 = " << atm_p2 << " , atm_s2 = " << atm_s2 << std::endl;

                for (icrd = 0; icrd < 3; ++icrd){
                    if (system->cell_dimension[icrd] == 1) {
                        vec[icrd] = system->xr_s[atm_p1][icrd] - system->xr_s[atm_s2][icrd];
                        if (std::abs(vec[icrd]) < 0.5) {
                            vec[icrd] = 0.0;
                        } else {
                            if (system->xr_s[atm_p1][icrd] < 0.5) {
                                vec[icrd] = -1.0;
                            } else {
                                vec[icrd] = 1.0;
                            }
                        }
                    } else if (system->cell_dimension[icrd] == 2){
                        vec[icrd] = system->xr_s[atm_p2][icrd] - system->xr_s[atm_s2][icrd];
                        vec[icrd] = dynamical->fold(vec[icrd]);
                        if (std::abs(system->xr_s[atm_p1][icrd] - system->xr_s[atm_s2][icrd]) > 0.5) vec[icrd] *= -1.0;
                    } else {
                        vec[icrd] = system->xr_s[atm_p2][icrd] - system->xr_s[atm_s2][icrd];
                        vec[icrd] = dynamical->fold(vec[icrd]);
                    }
                }

                std::cout << "# i = " << i << " , j = " << j << " , itran = " << itran;
                std::cout << " : " << vec[0] << " " << vec[1] << " " << vec[2] << std::endl;

                system->rotvec(vec, vec, system->lavec_s);
                system->rotvec(vec, vec, system->rlavec_p);

                phase = vec[0] * xk_in[0] + vec[1] * xk_in[1] + vec[2] * xk_in[2];
                exp_phase = std::exp(im * phase);
                /*

                if (std::abs(exp_phase.real()) < 1.0e-14) {
                exp_phase = std::complex<double>(0.0, exp_phase.imag());
                }
                if (std::abs(exp_phase.imag()) < 1.0e-14) {
                exp_phase = std::complex<double>(exp_phase.real(), 0.0);
                }*/


                for (icrd = 0; icrd < 3; ++icrd){
                    for (jcrd = 0; jcrd < 3; ++jcrd){
                        //        ctmp[icrd][jcrd] += fcs_phonon->dfc2[i][atm_s2][icrd][jcrd] * std::exp(im * phase);
                        ctmp[icrd][jcrd] += fcs_phonon->fc2[i][atm_s2][icrd][jcrd] * exp_phase;

                    }
                }


            }

            if (i != j) {
                std::cout << "i = " << i << " , j = " << j << std::endl;
                std::cout << "xk = " << xk_in[0] << " " << xk_in[1] << " " << xk_in[2] << std::endl;
                std::cout << "ctmp = " << ctmp[0][0] << std::endl;
            }

            for (icrd = 0; icrd < 3; ++icrd){
                for (jcrd = 0; jcrd < 3; ++jcrd){
                    dphi2[3 * i + icrd][3 * j + jcrd] = ctmp[icrd][jcrd];
                    //    dphi2[3 * i + icrd][3 * j + jcrd] = ctmp[icrd][jcrd] / std::sqrt(system->mass[atm_p1] * system->mass[atm_p2]);
                }
            }

        }
    }

}