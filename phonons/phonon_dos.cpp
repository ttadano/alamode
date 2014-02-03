#include "mpi_common.h"
#include "phonon_dos.h"
#include "kpoint.h"
#include "constants.h"
#include "error.h"
#include "system.h"
#include "memory.h"
#include "dynamical.h"
#include "write_phonons.h"
#include "integration.h"
#include <algorithm>
#include <vector>
#include <fstream>
#include <iomanip>
#include "parsephon.h"

using namespace PHON_NS;

Dos::Dos(PHON *phon): Pointers(phon){}

Dos::~Dos(){
    if(flag_dos) {
        memory->deallocate(energy_dos);
        memory->deallocate(dos_phonon);
        if(dynamical->eigenvectors) {
            memory->deallocate(pdos_phonon);
        }
    }
}

void Dos::setup()
{
    // This function must not called before dynamica->setup_dynamical()

    int i;

    MPI_Bcast(&emin, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&emax, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&delta_e, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if(kpoint->kpoint_mode == 2) {
        flag_dos = true;
    } else {
        flag_dos = false;
    }

    if(flag_dos && delta_e < eps12) error->exit("dos_setup()", "Too small delta_e");

    if(flag_dos) {
        n_energy = static_cast<int>((emax - emin) / delta_e);
        memory->allocate(energy_dos, n_energy);
        memory->allocate(dos_phonon, n_energy);

        for (i = 0; i < n_energy; ++i){
            energy_dos[i] = emin + delta_e * static_cast<double>(i);
        }

        if (dynamical->eigenvectors) {
            memory->allocate(pdos_phonon, system->natmin, n_energy);
        }
    }
}

void Dos::calc_dos()
{
    int i;
    unsigned int j, k;
    unsigned int nk = kpoint->nk;
    unsigned int neval = dynamical->neval;
    double **eval;

    double *dos_local;

    memory->allocate(eval, neval, nk);
    memory->allocate(dos_local, n_energy);

    for (j = 0; j < nk; ++j){
        for (k = 0; k < neval; ++k){
            eval[k][j] = writes->in_kayser(dynamical->eval_phonon[j][k]);
        }
    }

    for (i = 0; i < n_energy; ++i){
        dos_local[i] = 0.0;

        for (j = mympi->my_rank; j < neval; j += mympi->nprocs) {
            dos_local[i] += integration->dos_integration(eval[j], energy_dos[i]);
        }
    }

    MPI_Reduce(&dos_local[0], &dos_phonon[0], n_energy, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    memory->deallocate(dos_local);

    if (dynamical->eigenvectors) {

        // Calculate atom projected phonon-DOS

        unsigned int ik, imode, iat, icrd;
        unsigned int natmin = system->natmin;

        double **proj;
        double **pdos_local;
        memory->allocate(proj, neval, nk);
        memory->allocate(pdos_local, natmin, n_energy);

        for (iat = 0; iat < natmin; ++iat){

            for (imode = 0; imode < neval; ++imode){
                for (ik = 0; ik < nk; ++ik){
                    proj[imode][ik] = 0.0;

                    for (icrd = 0; icrd < 3; ++icrd){
                        proj[imode][ik] += std::norm(dynamical->evec_phonon[ik][imode][3 * iat + icrd]);
                    }
                }
            }

            for (i = 0; i < n_energy; ++i){
                pdos_local[iat][i] = 0.0;

                for (imode = mympi->my_rank; imode < neval; imode += mympi->nprocs) {
                    pdos_local[iat][i] += integration->do_tetrahedron(eval[imode], proj[imode], energy_dos[i]);
                }
            }            
        }

        MPI_Reduce(&pdos_local[0][0], &pdos_phonon[0][0], natmin*n_energy, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

        memory->deallocate(proj);
        memory->deallocate(pdos_local);
    }

    memory->deallocate(eval);
}

void Dos::calc_tdos()
{
    unsigned int nk = kpoint->nk;
    unsigned int ns = dynamical->neval;
    unsigned int i;
    unsigned int ik, jk;
    unsigned int is, js;
    unsigned int kcount;

    double k_tmp[3], xk_tmp[3];
    double xk_norm;

    double **e_tmp;    
    double *energy_dos, **tdos;

    n_energy = static_cast<int>((emax - emin) / delta_e);

    memory->allocate(e_tmp, 4, nk);
    memory->allocate(energy_dos, n_energy);
    memory->allocate(tdos, 4, n_energy);

    k_tmp[0] = 0.0; k_tmp[1] = 0.0; k_tmp[2] = 0.0;

    for (i = 0; i < n_energy; ++i){
        energy_dos[i] = emin + delta_e * static_cast<double>(i);
        for (unsigned int j = 0; j < 4; ++j) tdos[j][i] = 0.0;
    }

    for (is = 0; is < ns; ++is){
        for (js = 0; js < ns; ++js){

            kcount = 0;

            for (ik = 0; ik < nk; ++ik) {
                for (jk = 0; jk < nk; ++jk){

                    xk_tmp[0] = k_tmp[0] + kpoint->xk[ik][0] + kpoint->xk[jk][0];
                    xk_tmp[1] = k_tmp[1] + kpoint->xk[ik][1] + kpoint->xk[jk][1];
                    xk_tmp[2] = k_tmp[2] + kpoint->xk[ik][2] + kpoint->xk[jk][2];

                    for (i = 0; i < 3; ++i)  xk_tmp[i] = std::fmod(xk_tmp[i], 1.0);
                    xk_norm = std::pow(xk_tmp[0], 2) + std::pow(xk_tmp[1], 2) + std::pow(xk_tmp[2], 2);
                    if (std::sqrt(xk_norm) > eps15) continue; 

                    e_tmp[0][kcount] = - writes->in_kayser(dynamical->eval_phonon[ik][is] + dynamical->eval_phonon[jk][js]);
                    e_tmp[1][kcount] = writes->in_kayser(dynamical->eval_phonon[ik][is] + dynamical->eval_phonon[jk][js]);
                    e_tmp[2][kcount] = writes->in_kayser(dynamical->eval_phonon[ik][is] - dynamical->eval_phonon[jk][js]);
                    e_tmp[3][kcount] = - writes->in_kayser(dynamical->eval_phonon[ik][is] - dynamical->eval_phonon[jk][js]);

                    ++kcount;
                }
            }

            for (i = 0; i < n_energy; ++i){
                for (unsigned int j = 0; j < 4; ++j) tdos[j][i] += integration->dos_integration(e_tmp[j], energy_dos[i]);
            }

            std::cout << "kcount = " << kcount << std::endl;
        }
    }

    std::string file_tdos;
    std::ofstream ofs_tdos;

    file_tdos = input->job_title + ".tdos";

    ofs_tdos.open(file_tdos.c_str(), std::ios::out);

    ofs_tdos << "# Two-phonon DOS at" << std::setw(15) << k_tmp[0] << std::setw(15) << k_tmp[1] << std::setw(15) << k_tmp[2] << std::endl;
    for (i = 0; i < n_energy; ++i) {
        ofs_tdos << std::setw(15) << energy_dos[i];
        ofs_tdos << std::setw(15) << tdos[0][i];
        ofs_tdos << std::setw(15) << tdos[1][i];
        ofs_tdos << std::setw(15) << tdos[2][i];
        ofs_tdos << std::setw(15) << tdos[3][i] << std::endl;
    }

    ofs_tdos.close();
    memory->deallocate(e_tmp);
    memory->deallocate(tdos);
    memory->deallocate(energy_dos);
}
