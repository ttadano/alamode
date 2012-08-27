#include <fstream>
#include <iomanip>
#include <algorithm>
#include "relaxation.h"
#include "memory.h"
#include "kpoint.h"
#include "system.h"
#include "dynamical.h"
#include "fcs_phonon.h"
#include "../alm_c++/constants.h"
#include "phonon_thermodynamics.h"
#include "phonon_velocity.h"
#include "parsephon.h"
#include "error.h"

using namespace PHON_NS;

Relaxation::Relaxation(PHON *phon): Pointers(phon) {
    im = std::complex<double>(0.0, 1.0);
    epsilon = 2.5/Hz_to_kayser*time_ry;
//    epsilon = 2.5*1.0e-1/Hz_to_kayser*time_ry;
//    std::cout << "epsilon = " << epsilon << std::endl;
}

Relaxation::~Relaxation(){};

void Relaxation::setup_relaxation()
{
    unsigned int nk = kpoint->nk;
    unsigned int nband = dynamical->neval;

    memory->allocate(V, 1);
    memory->allocate(self_E, nk*nband);
    memory->allocate(tau, nk*nband);

    unsigned int i, j, k;

    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j){
            mat_convert[i][j] = 0.0;
            for (k = 0; k < 3; ++k){
                mat_convert[i][j] += system->lavec_s[i][k] * system->rlavec_p[k][j]; 
            }
        }
    }
}

void Relaxation::finish_relaxation()
{
    V[0].clear();
    memory->deallocate(V);
    memory->deallocate(self_E);
    memory->deallocate(tau);
}

void Relaxation::calc_ReciprocalV()
{
    // Calculate V_{ks, k's', k''s''}^(3) for Self-energy calculation

    unsigned int nband = dynamical->neval;
    unsigned int nk = kpoint->nk;
    unsigned int i;
    unsigned int k1, k2, k3;
    unsigned int b1, b2, b3;
    unsigned int ks_tmp[3];

    double xk_tmp[3], xk_norm;
    std::complex<double> prod;

    std::vector<FcsClass>::iterator it;

    std::cout << std::endl;
    std::cout << "Calculating force constants in reciprocal space .." << std::endl;

    double omega[3], omega_prod;

    for (k1 = 0; k1 < nk; ++k1){
        for (k2 = 0; k2 < nk; ++k2){
            for (k3 = 0; k3 < nk; ++k3){

                xk_tmp[0] = kpoint->xk[k1][0] + kpoint->xk[k2][0] + kpoint->xk[k3][0];
                xk_tmp[1] = kpoint->xk[k1][1] + kpoint->xk[k2][1] + kpoint->xk[k3][1];
                xk_tmp[2] = kpoint->xk[k1][2] + kpoint->xk[k2][2] + kpoint->xk[k3][2];

                for (i = 0; i < 3; ++i){
                    xk_tmp[i] = std::fmod(xk_tmp[i], 1.0);
                }
                xk_norm = std::pow(xk_tmp[0], 2) + std::pow(xk_tmp[1], 2) + std::pow(xk_tmp[2], 2);

                // If the momentum-conservation is not satisfied, we skip the loop.

                if (std::sqrt(xk_norm) > eps15) continue;

                for (b1 = 0; b1 < nband; ++b1){
                    for (b2 = 0; b2 < nband; ++b2){
                        for (b3 = 0; b3 < nband; ++b3){

                            ks_tmp[0] = nband * k1 + b1;
                            ks_tmp[1] = nband * k2 + b2;
                            ks_tmp[2] = nband * k3 + b3;

                            if (ks_tmp[0] > ks_tmp[1] || ks_tmp[1] > ks_tmp[2]) continue;

                            prod = V3(ks_tmp[0], ks_tmp[1], ks_tmp[2]);

                            // Add to list
                            if (std::abs(prod) > eps15) {
                                V[0].push_back(ReciprocalVs(prod, ks_tmp, 3));
                            }
                        }
                    }
                }
            }
        }
    }
    std::cout << "Done !" << std::endl;
    std::cout << "Number of nonzero V's: " << V[0].size() << std::endl;
}

std::complex<double> Relaxation::V3(const unsigned int ks1, const unsigned int ks2, const unsigned int ks3)
{
    unsigned int i;
    unsigned int ind[3], atmn[3], crdn[3];
    unsigned int k1, k2, k3;
    unsigned int b1, b2, b3;
    unsigned int atm_p1, atm_s1, atm_s2;
    unsigned int nband = dynamical->neval;

    double mass_prod, phase, omega_prod;
    double vec1[3], vec2[3];
    double omega[3];

    std::vector<FcsClass>::iterator it;
    std::complex<double> ret, tmp;

    ret = (0.0, 0.0);

    k1 = ks1 / nband;
    k2 = ks2 / nband;
    k3 = ks3 / nband;

    b1 = ks1 % nband;
    b2 = ks2 % nband;
    b3 = ks3 % nband;
    
    omega[0] = freq2(dynamical->eval_phonon[k1][b1]);
    omega[1] = freq2(dynamical->eval_phonon[k2][b2]);
    omega[2] = freq2(dynamical->eval_phonon[k3][b3]);
    omega_prod = omega[0] * omega[1] * omega[2];

    for (it = fcs_phonon->force_constant[1].begin(); it != fcs_phonon->force_constant[1].end(); ++it){

        for (i = 0; i < 3; ++i) {
            ind[i] = (*it).elems[i];
            atmn[i] = ind[i] / 3;
            crdn[i] = ind[i] % 3;
        }
        
        atm_p1 = system->map_p2s[0][system->map_s2p[atmn[0]].tran_num];
        atm_s1 = system->map_p2s[0][system->map_s2p[atmn[1]].tran_num];
        atm_s2 = system->map_p2s[0][system->map_s2p[atmn[2]].tran_num];

        for (i = 0; i < 3; ++i){
            vec1[i] = system->xr_s[atm_s1][i] - system->xr_s[atm_p1][i];
            vec2[i] = system->xr_s[atm_s2][i] - system->xr_s[atm_p1][i];
            vec1[i] = dynamical->fold(vec1[i]);
            vec2[i] = dynamical->fold(vec2[i]);
        }

        system->rotvec(vec1, vec1, system->lavec_s);
        system->rotvec(vec1, vec1, system->rlavec_p);
        system->rotvec(vec2, vec2, system->lavec_s);
        system->rotvec(vec2, vec2, system->rlavec_p);
/*
        system->rotvec(vec1, vec1, mat_convert);
        system->rotvec(vec2, vec2, mat_convert);
*/
        phase = vec1[0] * kpoint->xk[k2][0] 
              + vec1[1] * kpoint->xk[k2][1] 
              + vec1[2] * kpoint->xk[k2][2]
              + vec2[0] * kpoint->xk[k3][0] 
              + vec2[1] * kpoint->xk[k3][1]  
              + vec2[2] * kpoint->xk[k3][2];

        mass_prod = system->mass[atmn[0]] * system->mass[atmn[1]] * system->mass[atmn[2]];

        ret += (*it).fcs_val * std::exp(im * phase) / std::sqrt(mass_prod)
            * dynamical->evec_phonon[k1][b1][3 * system->map_s2p[atmn[0]].atom_num + crdn[0]]
            * dynamical->evec_phonon[k2][b2][3 * system->map_s2p[atmn[1]].atom_num + crdn[1]]
            * dynamical->evec_phonon[k3][b3][3 * system->map_s2p[atmn[2]].atom_num + crdn[2]];

    }

    return ret/std::sqrt(omega_prod);
}

void Relaxation::calc_selfenergy()
{
    double Tmin, Tmax, dT;

    Tmin = 1.0;
    Tmax = 2.0;
    dT= 1.0;

    unsigned int NT= static_cast<unsigned int>((Tmax - Tmin) / dT);

    unsigned int i;

    unsigned int nband = dynamical->neval;
    unsigned int nk = kpoint->nk;

    unsigned int nks = nband*nk;

    unsigned int iks;
    double T;

    std::string file_selfenergy;
    std::ofstream ofs_selfenergy;

    file_selfenergy = input->job_title + ".selfE";
    ofs_selfenergy.open(file_selfenergy.c_str(), std::ios::out);
    if(!ofs_selfenergy) error->exit("write_selfenergy", "cannot open file_selfenergy");

    ofs_selfenergy << "#Temperature, k-point, branch, Re[Sigma] [cm^-1], Im[Sigma]^[cm^-1]" << std::endl;
    ofs_selfenergy.setf(std::ios::scientific);

    for (i = 0; i < NT; ++i){
        T = Tmin + dT * static_cast<double>(i);
        calc_selfenergy_at_T(T);
        for (iks = 0; iks < nks; ++iks){
            ofs_selfenergy << std::setw(5) << T;
            ofs_selfenergy << std::setw(5) << iks / nband;
            ofs_selfenergy << std::setw(5) << iks % nband;
            ofs_selfenergy << std::setw(15) << relaxation->self_E[iks].real()/time_ry*Hz_to_kayser; 
            ofs_selfenergy << std::setw(15) << relaxation->self_E[iks].imag()/time_ry*Hz_to_kayser;
            ofs_selfenergy << std::endl;
        }
        ofs_selfenergy << std::endl;
    }
    ofs_selfenergy.close();

}

void Relaxation::calc_selfenergy_at_T(const double T)
{
    unsigned int i;
    unsigned int nk = kpoint->nk;
    unsigned int nband = dynamical->neval;
    unsigned int nks;
    double v_norm;

    unsigned int *ind, *knum, *snum;
    double *omega, omega_prod;
    double n1, n2;

    double fcell = 1.0 / static_cast<double>(system->ntran);

    unsigned int iks;

    memory->allocate(ind, 3);
    memory->allocate(knum, 3);
    memory->allocate(snum, 3);
    memory->allocate(omega, 3);

    nks = nk * nband;

    for (i = 0; i < nks; ++i) self_E[i] = 0.0;

    for (std::vector<ReciprocalVs>::iterator it = V[0].begin(); it != V[0].end(); ++it){
        ReciprocalVs obj = *it;
        v_norm = std::norm(obj.v);

        do {
            for (i = 0; i < 3; ++i){
                ind[i] = obj.ks[i];
                knum[i] = ind[i] / nband;
                snum[i] = ind[i] % nband;
                omega[i] = freq2(dynamical->eval_phonon[knum[i]][snum[i]]);
            }

            n1 = phonon_thermodynamics->fB(omega[1], T) + phonon_thermodynamics->fB(omega[2], T) + 1.0;
            n2 = phonon_thermodynamics->fB(omega[1], T) - phonon_thermodynamics->fB(omega[2], T);

            iks = nband * kpoint->knum_minus[knum[0]] + snum[0];
            if (iks >= 3 && iks <= 5) {
                std::cout << "knum = " << iks;
                std::cout << " ,-k1s1, k2s2, k3s3 = ";
                std::cout << std::setw(5) << knum[0] << std::setw(5) << snum[0];
                std::cout << std::setw(5) << knum[1] << std::setw(5) << snum[1];
                std::cout << std::setw(5) << knum[2] << std::setw(5) << snum[2];
                std::cout << std::setw(15) << obj.v.real();
                std::cout << std::setw(15) << obj.v.imag() << std::endl;
            }
            self_E[nband * kpoint->knum_minus[knum[0]] + snum[0]] 
                += std::pow(0.5, 4) * fcell * v_norm
                * ( n1 / (omega[0] + omega[1] + omega[2] + im * epsilon)
                    - n1 / (omega[0] - omega[1] - omega[2] + im * epsilon) 
                    + n2 / (omega[0] - omega[1] + omega[2] + im * epsilon)
                    - n2 / (omega[0] + omega[1] - omega[2] + im * epsilon));

        } while (std::next_permutation(obj.ks.begin(), obj.ks.end()));

    }
}

std::complex<double> Relaxation::delta_lorentz(const double omega)
{
    return 1.0 / (omega - im*epsilon);
}

void Relaxation::test_delta(const double T)
{
    unsigned int nk = kpoint->nk;
    unsigned int ns = dynamical->neval;

    int i, j;
    double omega;

    for (i = -1000; i <= 1000; ++i){
        omega = static_cast<double>(i) * 0.1;
        std::cout << "omega = " << omega << ", delta(omega)= " << delta_lorentz(omega).real() << " " << delta_lorentz(omega).imag() << std::endl;
    }
}

double Relaxation::freq2(const double x) 
{
    if (std::abs(x) < eps15) return 0.0;
    if (x >= 0.0) {
        return std::sqrt(x);
    } else {
        return std::sqrt(-x);
    }
}
