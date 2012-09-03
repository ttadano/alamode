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
Conductivity::~Conductivity(){
    memory->deallocate(vel);
    memory->deallocate(tau);
};

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
    double omega;

    for (ik = 0; ik < nk; ++ik){
        for (is = 0; is < ns; ++is){
            omega = phonon_velocity->freq(dynamical->eval_phonon[ik][is]);
            tau[ik][is] = 1.0 / (2.0 * relaxation->selfenergy(T, omega, ik, is).imag());
        }
    }

    for (i = 0; i < 3; ++i){
        for (j = 0; j < 3; ++j){

            kl[i][j] = 0.0;

            for (ik = 0; ik < nk; ++ik){
                for (is = 0; is < ns; ++is){

                    omega = phonon_velocity->freq(dynamical->eval_phonon[ik][is]);

                    kl[i][j] += phonon_thermodynamics->Cv(omega, T) * vel[ik][is][i] * vel[ik][is][j] * tau[ik][is];
                }
            }
            kl[i][j] /= Bohr_in_Angstrom * 1.0e-10 * time_ry * static_cast<double>(nk) * system->volume_p;
        }
    }

}
