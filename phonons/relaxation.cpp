#include "relaxation.h"
#include "memory.h"
#include "kpoint.h"
#include "system.h"
#include "dynamical.h"
#include "fcs_phonon.h"
#include "../alm_c++/constants.h"
#include <iomanip>
#include <algorithm>
#include "phonon_thermodynamics.h"
#include "phonon_velocity.h"

using namespace PHON_NS;

Relaxation::Relaxation(PHON *phon): Pointers(phon) {
    im = std::complex<double>(0.0, 1.0);
    epsilon = 1.0e-8;
}
Relaxation::~Relaxation(){};


void Relaxation::setup_relaxation()
{
    unsigned int nk = kpoint->nk;
    unsigned int nband = dynamical->neval;

    memory->allocate(V, 1);
    memory->allocate(tau, nk*nband);
}

void Relaxation::calc_ReciprocalV()
{
    // Calculate V_{ks, k's', k''s''}^(3) for Self-energy calculation

    unsigned int nband = dynamical->neval;
    unsigned int natmin = system->natmin;
    unsigned int nat = system->nat;
    unsigned int nk = kpoint->nk;

    unsigned int i;
    unsigned int *atmn, *ind, *ind_mod;
    unsigned int *crdn;

    unsigned int k1, k2, k3;
    unsigned int b1, b2, b3;

    double mass_prod, omega_prod;

    std::complex<double> prod, prod_tmp;


    std::cout << std::endl;
    std::cout << "Calculating force constants in reciprocal space .." << std::endl;

    memory->allocate(atmn, 3);
    memory->allocate(crdn, 3);
    memory->allocate(ind, 3);
    memory->allocate(ind_mod, 3);

    std::set<FcsClass>::iterator location;
    std::vector<FcsClass>::iterator it;

    double xk_tmp[3], xk_norm;

    std::complex<double> exp_sum;

    std::vector<unsigned int> ks_tmp;


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

                            ks_tmp.clear();

                            ks_tmp.push_back(nband * k1 + b1);
                            ks_tmp.push_back(nband * k2 + b2);
                            ks_tmp.push_back(nband * k3 + b3);

                            if (ks_tmp[0] > ks_tmp[1] || ks_tmp[1] > ks_tmp[2]) continue;

                            prod = (0.0, 0.0);

                            for (it = fcs_phonon->force_constant[1].begin(); it != fcs_phonon->force_constant[1].end(); ++it){

                                for (i = 0; i < 3; ++i) {
                                    ind[i] = (*it).elems[i];
                                    atmn[i] = ind[i] / 3;
                                    crdn[i] = ind[i] % 3;
                                }

                                mass_prod = 1.0 / std::sqrt(system->mass[atmn[0]] * system->mass[atmn[1]] * system->mass[atmn[2]]);

                                prod += mass_prod * (*it).fcs_val 
                                    * dynamical->evec_phonon[k1][b1][3 * system->map_s2p[atmn[0]].atom_num + crdn[0]]
                                * dynamical->evec_phonon[k2][b2][3 * system->map_s2p[atmn[1]].atom_num + crdn[1]]
                                * dynamical->evec_phonon[k3][b3][3 * system->map_s2p[atmn[2]].atom_num + crdn[2]];
                            }

                            if (std::norm(prod) > eps12){
                               // omega_prod = dynamical->eval_phonon[k1][b1] * dynamical->eval_phonon[k2][b2] * dynamical->eval_phonon[k3][b3];
                               // prod *= std::pow(0.5, 1.5) * std::sqrt(fcell / omega_prod);

                                V[0].push_back(ReciprocalVs(prod, ks_tmp));
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


void Relaxation::calc_selfenergy(const double T)
{
    unsigned int i;
    unsigned int nk = kpoint->nk;
    unsigned int nband = dynamical->neval;
    unsigned int nks;
    double v_norm;

    unsigned int *ind, *knum, *snum;
    double *omega, omega_prod;

    
    double fcell = 1.0 / static_cast<double>(system->ntran);

    memory->allocate(ind, 3);
    memory->allocate(knum, 3);
    memory->allocate(snum, 3);
    memory->allocate(omega, 3);

    nks = nk * nband;

    for (i = 0; i < nks; ++i) tau[i] = 0.0;

    for (std::vector<ReciprocalVs>::iterator it = V[0].begin(); it != V[0].end(); ++it){
        ReciprocalVs obj = *it;
        v_norm = std::norm(obj.v);

        do {

            for (i = 0; i < 3; ++i){
                ind[i] = obj.ks[i];
                knum[i] = ind[i] / nband;
                snum[i] = ind[i] % nband;
                omega[i] = phonon_velocity->freq(dynamical->eval_phonon[knum[i]][snum[i]]);
            }

            omega_prod = omega[0] * omega[1] * omega[2];
            
            tau[ind[0]] += std::pow(0.5, 3) * fcell / omega_prod
                * (v_norm
                * ( -(phonon_thermodynamics->fB(omega[1], T) + phonon_thermodynamics->fB(omega[2], T) + 1.0) * delta_lorentz(omega[0] + omega[1] + omega[2])
                    +(phonon_thermodynamics->fB(omega[1], T) + phonon_thermodynamics->fB(omega[2], T) + 1.0) * delta_lorentz(omega[0] - omega[1] - omega[2])
                    -(phonon_thermodynamics->fB(omega[1], T) - phonon_thermodynamics->fB(omega[2], T)) * delta_lorentz(omega[0] - omega[1] + omega[2])
                    +(phonon_thermodynamics->fB(omega[1], T) - phonon_thermodynamics->fB(omega[2], T)) * delta_lorentz(omega[0] + omega[1] - omega[2]))).real();

        } while (std::next_permutation(obj.ks.begin(), obj.ks.end()));

    }

    for (i = 0; i < nks; ++i){
    std::cout << "i = " << i << " " << tau[i] << std::endl;
    }
}

std::complex<double> Relaxation::delta_lorentz(const double omega)
{
    return 1.0 / (omega - im*epsilon);
}