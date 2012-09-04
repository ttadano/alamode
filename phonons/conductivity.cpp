#include "conductivity.h"
#include "memory.h"
#include "kpoint.h"
#include "dynamical.h"
#include "relaxation.h"
#include "phonon_velocity.h"
#include "phonon_thermodynamics.h"
#include "integration.h"
#include <fstream>
#include <iomanip>
#include "parsephon.h"
#include "error.h"
#include "write_phonons.h"
#include "../alm_c++/constants.h"
#include <iostream>
#include "system.h"

using namespace PHON_NS;

Conductivity::Conductivity(PHON *phon): Pointers(phon) {}
Conductivity::~Conductivity(){};

void Conductivity::setup_kl()
{
    nk = kpoint->nk;
    ns = dynamical->neval;

    memory->allocate(vel, nk, ns, 3);
    memory->allocate(tau, nk, ns);

    unsigned int i, j, k;

    for (i = 0; i < nk; ++i){
        phonon_velocity->phonon_vel_k(kpoint->xk[i], vel[i]);

        // Generate phonon velocity in cartesian coordinate
        for (j = 0; j < ns; ++j){
            system->rotvec(vel[i][j], vel[i][j], system->lavec_p, 'T');
            for (k = 0; k < 3; ++k) vel[i][j][k] /= 2.0 * pi;
        }
    }
/*
    double vtmp[3][3];
    double **func;
    unsigned int ik, is;

    memory->allocate(func, nk, ns);
    for (ik = 0; ik < nk; ++ik){
        for (is = 0; is < ns; ++is){
            func[ik][is] = 10.0 * static_cast<double>(ik + kpoint->knum_minus[ik]) + 1.0 * static_cast<double>(is);
        }
    }
    
    for(i = 0; i < 3; ++i){
        for (j = 0; j < 3; ++j){
            vtmp[i][j] = 0.0;
        }
    }

    for (i = 0; i < nk; ++i){
        for (j = 0; j < ns; ++j){
            std::cout << "#" << std::setw(5) << i << std::setw(5) << j << ":";
            for(unsigned int k = 0; k < 3; ++k){
                std::cout << std::setw(15) << vel[i][j][k];
            }
            std::cout << std::endl;

            for (unsigned int mu = 0; mu < 3; ++mu){
                for (unsigned int nu = 0; nu < 3; ++nu){
                    vtmp[mu][nu] += vel[i][j][mu] * vel[i][j][nu] * func[i][j];
                }
            }
        }

    }

    for(i = 0; i < 3; ++i){
        for (j = 0; j < 3; ++j){
            std::cout << std::setw(15) << vtmp[i][j]; 
        }
        std::cout << std::endl;
    }

    error->exit("conductivity", "KSK");
*/

    std::cout.setf(std::ios::fixed);

    std::cout << " Tmin = " << std::setw(10) << system->Tmin; 
    std::cout << " Tmax = " << std::setw(10) << system->Tmax; 
    std::cout << " dT   = " << std::setw(10) << system->dT; 
    std::cout << std::endl;

    std::cout.unsetf(std::ios::fixed);
}
void Conductivity::finish_kl()
{
    memory->deallocate(vel);
    memory->deallocate(tau);
}

void Conductivity::calc_kl()
{
    unsigned int iT, i, j;
    double T;

    double Tmax = system->Tmax;
    double Tmin = system->Tmin;
    double dT = system->dT;

    std::string file_kl;
    std::ofstream ofs_kl;

    unsigned int NT= static_cast<unsigned int>((Tmax - Tmin) / dT);

    file_kl = input->job_title + ".kl";
    ofs_kl.open(file_kl.c_str(), std::ios::out);
    if(!ofs_kl) error->exit("calc_kl", "cannot open file_kl");

    ofs_kl << "# Temperature [K], Thermal Conductivity (xx, xy, xz, yx, yy, yz, zx, zy, zz) [W/mK]" << std::endl;

    for (iT = 0; iT <= NT; ++iT){
        T = Tmin + dT * static_cast<double>(iT);
        calc_kl_at_T(T);

        ofs_kl << std::setw(5) << T;
        for (i = 0; i < 3; ++i){
            for (j = 0; j < 3; ++j){
                ofs_kl << std::setw(15) << kl[i][j];
            }
        }
        ofs_kl << std::endl;
    }

    ofs_kl.close();
}

void Conductivity::calc_kl_at_T(const double T)
{
    unsigned int i, j;
    unsigned int is, ik;
    unsigned int jk;
    double omega;
    unsigned int knum;
    double tau_tmp;

/*
    for (ik = 0; ik < kpoint->kpIBZ.size(); ++ik){

        for (is = 0; is < ns; ++is) {
            omega = dynamical->eval_phonon[kpoint->kpIBZ[ik].knum][is];
            tau[ik][is] = 1.0 / (2.0 * relaxation->selfenergy(T, omega, kpoint->kpIBZ[ik].knum, is).imag());
        }
    }


    jk = 0;

    for (ik = 0; ik < kpoint->nk_equiv.size(); ++ik){
        for (is = 0; is < ns; ++is){
            omega = dynamical->eval_phonon[kpoint->kpIBZ[jk].knum][is];
            tau[jk][is] = 1.0 / (2.0 * relaxation->selfenergy(T, omega, kpoint->kpIBZ[jk].knum, is).imag());
        }
        jk += kpoint->nk_equiv[ik];
    }
*/

    for (i = 0; i < 3; ++i){
        for (j = 0; j < 3; ++j){

            kl[i][j] = 0.0;
        }
    }

    jk = 0;

    for (ik = 0; ik < kpoint->nk_equiv.size(); ++ik){

        knum = kpoint->kpIBZ[jk].knum;

        for (is = 0; is < ns; ++is){

            omega = dynamical->eval_phonon[knum][is];
            tau[knum][is] = 1.0 / (2.0 * relaxation->selfenergy(T, omega, knum, is).imag());

            for (i = 0; i < 3; ++i){
                for (j = 0; j < 3; ++j){
                    kl[i][j] += kpoint->weight_k[ik] * phonon_thermodynamics->Cv(omega, T) * vel[knum][is][i] * vel[knum][is][j] * tau[knum][is];
                }
            }
        }

        jk += kpoint->nk_equiv[ik];
    }

    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            kl[i][j] /= Bohr_in_Angstrom * 1.0e-10 * time_ry * system->volume_p;
        }
    }
}
