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


using namespace PHON_NS;

Conductivity::Conductivity(PHON *phon): Pointers(phon) {
}

Conductivity::~Conductivity(){};

void Conductivity::setup_kl()
{
    nk = kpoint->nk;
    ns = dynamical->neval;

    memory->allocate(tau, nk, ns);
    memory->allocate(vel, nk, ns, 3);
    memory->allocate(func, nk);

    unsigned int i;

    for (i = 0; i < nk; ++i){
        phonon_velocity->phonon_vel_k(kpoint->xk[i], vel[i]);
    }
}

void Conductivity::gen_tau(const double T)
{
    unsigned int i, j;

    relaxation->calc_selfenergy_at_T(T);

    unsigned int k = 0;

    for (i = 0; i < nk; ++i){
        for (j = 0; j < ns; ++j){
            tau[i][j] = 1.0 / (2.0 * relaxation->self_E[k++].imag());
            tau[i][j] = 1.0;
        }
    }
}

void Conductivity::calc_kl()
{
    unsigned int iT, i, j;
    double T;

    std::string file_kl;
    std::ofstream ofs_kl;

    unsigned int NT= static_cast<unsigned int>((Tmax - Tmin) / dT);

    file_kl = input->job_title + ".kl";
    ofs_kl.open(file_kl.c_str(), std::ios::out);
    if(!ofs_kl) error->exit("calc_kl", "cannot open file_kl");

    ofs_kl << "# Temperature [K], Thermal Conductivity (xx, xy, xz, yx, yy, yz, zx, zy, zz)" << std::endl;

    for (iT = 0; iT < NT; ++iT){
        T = Tmin + dT * static_cast<double>(iT);
        gen_tau(T);
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

    double **eval;

    unsigned int iomega;

    double omega_min = 0.0;
    double omega_max = 600.0;
    double domega = 1.0;

    unsigned int nomega = static_cast<unsigned int>((omega_max - omega_min) / domega);

    memory->allocate(eval, ns, nk);

    for (ik = 0; ik < nk; ++ik){
        for (is = 0; is < ns; ++is){
            eval[is][ik] = writes->in_kayser(dynamical->eval_phonon[ik][is]);
        }
    }
    std::cout << "Cv, v_i, v_j, tau" << std::endl;
    for (i = 0; i < 3; ++i){

        for (j = 0; j < 3; ++j){

            kl[i][j] = 0.0;

            for (is = 0; is < ns; ++is){
                std::cout << "is = " << is << std::endl;
                for (ik = 0; ik < nk; ++ik){
                    func[ik] =  phonon_thermodynamics->Cv(phonon_velocity->freq(dynamical->eval_phonon[ik][is]), T) 
                        * vel[ik][is][i] * vel[ik][is][j] * tau[ik][is];
                    std::cout << "ik = " << ik << " " << phonon_thermodynamics->Cv(phonon_velocity->freq(dynamical->eval_phonon[ik][is]), T) << " ";
                    std::cout << vel[ik][is][i] << " " << vel[ik][is][j] << " " << tau[ik][is] << std::endl;
                }
                for (iomega = 0; iomega < nomega; ++iomega){
                    kl[i][j] += domega*integration->do_tetrahedron(eval[is], func, omega_min + domega*static_cast<double>(iomega));
                }
            }
            kl[i][j] /= Bohr_in_Angstrom * 1.0e-10 * time_ry * static_cast<double>(nk);
        }

    }
}