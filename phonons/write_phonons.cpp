#include "write_phonons.h"
#include "system.h"
#include "dynamical.h"
#include "kpoint.h"
#include "parsephon.h"
#include "error.h"
#include "phonon_dos.h"
#include "phonon_velocity.h"
#include "../alm_c++/constants.h"
#include "memory.h"
#include <iomanip>
#include <fstream>

using namespace PHON_NS;

Writes::Writes(PHON *phon): Pointers(phon){
    Ry_to_kayser = std::pow(Hz_to_kayser, 2) / (amu_ry * std::pow(time_ry, 2));
};

Writes::~Writes(){};

void Writes::write_phonon_info()
{
    if (nbands < 0 || nbands > 3 * system->natmin) {
        std::cout << "WARNING: nbands < 0 or nbands > 3 * natmin" << std::endl;
        std::cout << "All modes will be printed." << std::endl;    
        nbands =  3 * system->natmin;
    }

    if(kpoint->kpoint_mode == 1){
        write_phonon_bands();
        write_phonon_vel();
    }

    if(dos->flag_dos) {
        write_phonon_dos();
    }

    if(writeanime) {
        write_mode_anime();
    }

    if(dynamical->eigenvectors) {
        write_eigenvectors();
    }
}

void Writes::write_phonon_bands()
{
    std::ofstream ofs_bands;

    file_bands = input->job_title + ".bands";
    ofs_bands.open(file_bands.c_str(), std::ios::out);
    if(!ofs_bands) error->exit("write_phonon_bands", "cannot open file_bands");

    unsigned int i, j;

    unsigned int nk = kpoint->nk;

    double *kaxis = kpoint->kaxis;
    double **eval = dynamical->eval_phonon;

    ofs_bands << "# k-axis, Eigenvalues [cm^-1]" << std::endl;
    ofs_bands.setf(std::ios::fixed);

    for (i = 0; i < nk; ++i){
        ofs_bands << std::setw(8) << kaxis[i];
        for (j = 0; j < nbands; ++j){
            ofs_bands << std::setw(12) << in_kayser(eval[i][j]);
        }
        ofs_bands << std::endl;
    }

    ofs_bands.close();
}

void Writes::write_phonon_vel()
{
    std::ofstream ofs_vel;

    file_vel = input->job_title + ".phvel";
    ofs_vel.open(file_vel.c_str(), std::ios::out);
    if(!ofs_vel) error->exit("write_phonon_vel", "cannot open file_vel");

        unsigned int i, j;

    unsigned int nk = kpoint->nk;

    double *kaxis = kpoint->kaxis;
    double **eval = dynamical->eval_phonon;

    ofs_vel << "# k-axis, Velocity [Ry Bohr]" << std::endl;
    ofs_vel.setf(std::ios::fixed);

    for (i = 0; i < nk; ++i){
        ofs_vel << std::setw(8) << kaxis[i];
        for (j = 0; j < nbands; ++j){
            ofs_vel << std::setw(12) << phonon_velocity->phvel[i][j];
        }
        ofs_vel << std::endl;
    }

    ofs_vel.close();
}


void Writes::write_phonon_dos()
{
    int i;
    std::ofstream ofs_dos;

    file_bands = input->job_title + ".dos";
    ofs_dos.open(file_bands.c_str(), std::ios::out);
    if(!ofs_dos) error->exit("write_phonon_dos", "cannot open file_dos");

    ofs_dos << "# Energy [cm^-1], population" << std::endl;
    ofs_dos.setf(std::ios::scientific);

    for (i = 0; i < dos->n_energy; ++i){
        ofs_dos << std::setw(15) << dos->energy_dos[i] << std::setw(15) << dos->dos_phonon[i] << std::endl;
    } 

    ofs_dos.close();
}

void Writes::write_mode_anime()
{
    std::ofstream ofs_anime;

    file_anime = input->job_title + ".axsf";
    ofs_anime.open(file_anime.c_str(), std::ios::out);
    if(!ofs_anime) error->exit("write_mode_anime", "cannot open file_anime");

    ofs_anime.setf(std::ios::scientific);

    unsigned int i, j, k;
    unsigned int natmin = system->natmin;
    unsigned int nk = kpoint->nk;

    double force_factor = 100.0;

    double **xmod;
    std::string *kd_tmp;

    memory->allocate(xmod, natmin, 3);
    memory->allocate(kd_tmp, natmin);

    ofs_anime << "ANIMSTEPS " << nbands * nk << std::endl;
    ofs_anime << "CRYSTAL" << std::endl;
    ofs_anime << "PRIMVEC" << std::endl;

    for (i = 0; i < 3; ++i){
        for (j = 0; j < 3; ++j){
            ofs_anime << std::setw(15) << system->lavec_p[j][i]*Bohr_in_Angstrom;
        }
        ofs_anime << std::endl;
    }

    for (i = 0; i < natmin; ++i){
        k = system->map_p2s[i][0];
        for (j = 0; j < 3; ++j){
            xmod[i][j] = system->xc[k][j];
        }
        // system->rotvec(system->lavec_p, xmod[i], xmod[i]);

        for (j = 0; j < 3; ++j){
            xmod[i][j] *= Bohr_in_Angstrom;
        }
        kd_tmp[i] = system->symbol_kd[system->kd[k]];
    }

    unsigned int ik, imode;
    double norm;
    std::complex<double> evec_tmp;
    unsigned int m;
    i = 0;

    for (ik = 0; ik < nk; ++ik){
        for (imode = 0; imode < nbands; ++imode){
            ofs_anime << "PRIMCOORD " << std::setw(10) << i + 1 << std::endl;
            ofs_anime << std::setw(10) << natmin << std::setw(10) << 1 << std::endl;
            norm = 0.0;

            for (j = 0; j < 3 * natmin; ++j){
                evec_tmp = dynamical->evec_phonon[ik][imode][j];
                norm += std::pow(evec_tmp.real(), 2) + std::pow(evec_tmp.imag(), 2);
            }

            norm *= force_factor / static_cast<double>(natmin);

            for (j = 0; j < natmin; ++j){

                m = system->map_p2s[j][0];

                ofs_anime << std::setw(10) << kd_tmp[j];

                for (k = 0; k < 3; ++k){
                    ofs_anime << std::setw(15) << xmod[j][k];
                }
                for (k = 0; k < 3; ++k){
                    ofs_anime << std::setw(15) << dynamical->evec_phonon[ik][imode][3 * j + k].real() / (std::sqrt(system->mass[m]) * norm);
                }
                ofs_anime << std::endl;
            }

            ++i;
        }
    }

    memory->deallocate(xmod);
    memory->deallocate(kd_tmp);

    ofs_anime.close();
}

void Writes::write_eigenvectors()
{
    std::ofstream ofs_evec;
    file_evec = input->job_title + ".evec";
    ofs_evec.open(file_evec.c_str(), std::ios::out);
    if(!ofs_evec) error->exit("write_eigenvectors", "cannot open file_evec");

    ofs_evec.setf(std::ios::scientific);
    unsigned int i, j, k;

    ofs_evec << "Lattice vectors of the primitive lattice" << std::endl;

    for (i = 0; i < 3; ++i){
        for (j = 0; j < 3; ++j){
            ofs_evec << std::setw(15) << system->lavec_p[j][i];
        }
        ofs_evec << std::endl;
    }

    ofs_evec << std::endl;

    for (i = 0; i < 3; ++i){
        for (j = 0; j < 3; ++j){
            ofs_evec << std::setw(15) << system->rlavec_p[i][j];
        }
        ofs_evec << std::endl;
    }

    unsigned int nk = kpoint->nk;
    unsigned int neval = dynamical->neval;
    ofs_evec << "Modes and k-points information below" << std::endl;
    ofs_evec << std::setw(10) << nbands;
    ofs_evec << std::setw(10) << nk << std::endl;

    for (i = 0; i < nk; ++i){
        ofs_evec << "#" << std::setw(10) << i + 1;
        for (j = 0; j < 3; ++j){
            ofs_evec << std::setw(15) << kpoint->xk[i][j];
        }
        ofs_evec << std::endl;
        for (j = 0; j < nbands; ++j){
            ofs_evec << std::setw(15) << dynamical->eval_phonon[i][j] << std::endl;

            for (k = 0; k < neval; ++k){
                ofs_evec << std::setw(15) << real(dynamical->evec_phonon[i][j][k]);
                ofs_evec << std::setw(15) << imag(dynamical->evec_phonon[i][j][k]) << std::endl;
            }
            ofs_evec << std::endl;
        }
        ofs_evec << std::endl;
    }
    ofs_evec.close();
}

double Writes::in_kayser(const double x)
{
    double val = x;
    val *= Ry_to_kayser;

    if(val < 0.0){
        return -std::sqrt(-val);
    } else {
        return std::sqrt(val);
    }
}
