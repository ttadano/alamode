/*
phonon_dos.cpp

Copyright (c) 2014 Terumasa Tadano

This file is distributed under the terms of the MIT license.
Please see the file 'LICENCE.txt' in the root directory 
or http://opensource.org/licenses/mit-license.php for information.
*/

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
#include "symmetry_core.h"

using namespace PHON_NS;

Dos::Dos(PHON *phon): Pointers(phon){}

Dos::~Dos(){
    if (flag_dos) {
        memory->deallocate(energy_dos);
        memory->deallocate(dos_phonon);
        if (projected_dos) {
            memory->deallocate(pdos_phonon);
        }
        if (two_phonon_dos) {
            memory->deallocate(dos2_phonon);
        }
        memory->deallocate(kmap_irreducible);
    }
}

void Dos::setup()
{
    // This function must not be called before dynamica->setup_dynamical()

    int i;

    MPI_Bcast(&emin, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&emax, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&delta_e, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&projected_dos, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD);
    MPI_Bcast(&two_phonon_dos, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD);

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

        if (projected_dos) {
            memory->allocate(pdos_phonon, system->natmin, n_energy);
        }

        if (two_phonon_dos) {
            int n_energy2 = static_cast<int>((emax * 2.0 - emin) / delta_e);
            memory->allocate(dos2_phonon, n_energy2, 4);
        }


        int ***symmetry_tmp;

        memory->allocate(kmap_irreducible, kpoint->nk);
        memory->allocate(symmetry_tmp, symmetry->nsym, 3, 3);

        for (i = 0; i < symmetry->nsym; ++i) {
            for (int j = 0; j < 3; ++j) {
                for (int k = 0; k < 3; ++k) {
                    symmetry_tmp[i][j][k] = symmetry->SymmList[i].symop[3 * j + k];
                }
            }
        }

        kpoint->generate_irreducible_kmap(kmap_irreducible, nk_irreducible, k_irreducible,
            kpoint->nkx, kpoint->nky, kpoint->nkz, 
            kpoint->xk, symmetry->nsym, symmetry_tmp);

        memory->deallocate(symmetry_tmp);
    }
}

void Dos::calc_dos_all()
{
    int i;
    unsigned int j, k;
    unsigned int nk = kpoint->nk;
    unsigned int neval = dynamical->neval;
    double **eval;
    double *weight;

    memory->allocate(eval, neval, nk);

    for (j = 0; j < nk; ++j){
        for (k = 0; k < neval; ++k){
            eval[k][j] = writes->in_kayser(dynamical->eval_phonon[j][k]);
        }
    }

    calc_dos(nk_irreducible, kmap_irreducible, eval, n_energy, energy_dos,
        dos_phonon, neval, integration->ismear, kpoint->kpoint_irred_all);


    if (projected_dos) {
        calc_atom_projected_dos(nk, eval, n_energy, energy_dos,
            pdos_phonon, neval, system->natmin, integration->ismear, dynamical->evec_phonon);        
    }
    memory->deallocate(eval);

    if (two_phonon_dos) {
        int n_energy2 = static_cast<int>((emax * 2.0 - emin) / delta_e);
        double *energy2;

        memory->allocate(energy2, n_energy2);

        for (i = 0; i < n_energy2; ++i) {
            energy2[i] = emin + delta_e * static_cast<double>(i);
        }

        calc_two_phonon_dos(n_energy2, energy2, dos2_phonon, integration->ismear, kpoint->kpoint_irred_all);       
        memory->deallocate(energy2);
    }
}

void PHON_NS::Dos::calc_dos(const unsigned int nk_irreducible, int *map_k, double **eval, 
                            const unsigned int n, double *energy, double *ret, const unsigned int neval, const int smearing_method,
                            std::vector<std::vector<KpointList> > &kpinfo)
{
    int i, j, k;
    double *weight;

    if (mympi->my_rank == 0) std::cout << " Calculating phonon DOS ...";

#pragma omp parallel private (weight, k)
    {
        memory->allocate(weight, nk_irreducible);

#pragma omp for
        for (i = 0; i < n; ++i) {

            dos_phonon[i] = 0.0;

            for (k = 0; k < neval; ++k) {
                if (smearing_method == -1) {
                    integration->calc_weight_tetrahedron(nk_irreducible, map_k, weight, eval[k], energy[i]);
                } else {
                    integration->calc_weight_smearing(kpinfo, weight, eval[k], energy[i], smearing_method);
                }

                for (j = 0; j < nk_irreducible; ++j) {
                    ret[i] += weight[j];
                }
            }
        }
        memory->deallocate(weight);
    }

    if (mympi->my_rank == 0) std::cout << " done." << std::endl;
}

void PHON_NS::Dos::calc_atom_projected_dos(const unsigned int nk, double **eval, const unsigned int n, double *energy,
                                           double **ret, const unsigned int neval, const unsigned int natmin, const int smearing_method,
                                           std::complex<double> ***evec)
{
    // Calculate atom projected phonon-DOS

    int i;
    unsigned int j, k;
    unsigned int ik, imode, iat, icrd;
    int *kmap_identity;
    double *weight;
    double **proj;

    if (mympi->my_rank == 0) std::cout << " PDOS = 1 : Calculating atom-projected phonon DOS ...";

    memory->allocate(kmap_identity, nk);
    memory->allocate(proj, neval, nk);


    for (i = 0; i < nk; ++i) kmap_identity[i] = i;

    for (iat = 0; iat < natmin; ++iat){

        for (imode = 0; imode < neval; ++imode){
            for (i = 0; i < nk; ++i){

                proj[imode][i] = 0.0;

                for (icrd = 0; icrd < 3; ++icrd){
                    proj[imode][i] += std::norm(evec[i][imode][3 * iat + icrd]);
                }
            }
        }

#pragma omp parallel private (weight, k, j)
        {
            memory->allocate(weight, nk);

#pragma omp for
            for (i = 0; i < n; ++i){
                ret[iat][i] = 0.0;

                for (k = 0; k < neval; ++k) {
                    if (smearing_method == -1) {
                        integration->calc_weight_tetrahedron(nk, kmap_identity, weight, eval[k], energy[i]);
                    } else {
                        integration->calc_weight_smearing(nk, nk, kmap_identity, weight, eval[k], energy[i], smearing_method);
                    }

                    for (j = 0; j < nk; ++j) {
                        ret[iat][i] += proj[k][j] * weight[j];
                    }
                }
            }            

            memory->deallocate(weight);
        }
    }
    memory->deallocate(proj);
    memory->deallocate(kmap_identity);

    if (mympi->my_rank == 0) std::cout << " done." << std::endl;
}


void Dos::calc_two_phonon_dos(const unsigned int n, double *energy, double **ret, const int smearing_method,
                              std::vector<std::vector<KpointList> > kpinfo)
{
    int i, j;
    int is, js, ik, jk;
    int k;

    int knum;

    unsigned int nk = kpoint->nk;
    unsigned int ns = dynamical->neval;

    int *kmap_identity;

    double multi;
    double **e_tmp;    
    double **tdos_q;

    double **weight;
    double emax2 = 2.0 * emax;

    double xk_tmp[3];

    int loc;

    if (mympi->my_rank == 0) {
        std::cout << " TDOS = 1 : Calculating two-phonon DOS with doubled EMAX" << std::endl;
        std::cout << "            This can take a while ... ";
    }

    memory->allocate(tdos_q, n, 4);
    memory->allocate(kmap_identity, nk);
    memory->allocate(e_tmp, 4, nk);
    memory->allocate(weight, n, nk);

    for (i = 0; i < nk; ++i) kmap_identity[i] = i;

    for (i = 0; i < n; ++i) {
        for (j = 0; j < 4; ++j) {
            ret[i][j] = 0.0;
        }
    }

    for (is = 0; is < ns; ++is) {
        for (js = 0; js < ns; ++js) {

            for (ik = 0; ik < kpinfo.size(); ++ik) {

                multi = static_cast<double>(kpinfo[ik].size()) / static_cast<double>(nk);
                knum = kpinfo[ik][0].knum;

                for (i = 0; i < n; ++i) {
                    for (j = 0; j < 4; ++j) {
                        tdos_q[i][j] = 0.0;
                    }
                }

#pragma omp parallel for private (j, loc, xk_tmp) 
                for (jk = 0; jk < nk; ++jk) {

                    for (j = 0; j < 3; ++j) xk_tmp[j] = kpoint->xk[knum][j] + kpoint->xk[jk][j];
                    loc = kpoint->get_knum(xk_tmp[0], xk_tmp[1], xk_tmp[2]);

                    e_tmp[0][jk] = - writes->in_kayser(dynamical->eval_phonon[jk][is] + dynamical->eval_phonon[loc][js]);
                    e_tmp[1][jk] = - e_tmp[0][jk];
                    e_tmp[2][jk] =   writes->in_kayser(dynamical->eval_phonon[jk][is] - dynamical->eval_phonon[loc][js]);
                    e_tmp[3][jk] = - e_tmp[2][jk];
                }

                if (smearing_method == -1) {

                    for (j = 1; j < 4; ++j) {
#pragma omp parallel for private (k)
                        for (i = 0; i < n; ++i) {
                            integration->calc_weight_tetrahedron(nk, kmap_identity, weight[i], e_tmp[j], energy[i]);
                            for (k = 0; k < nk; ++k) {
                                tdos_q[i][j] += weight[i][k];
                            }
                        }
                    }

                } else {

                    for (j = 0; j < 4; ++j) {
#pragma omp parallel for private (k)
                        for (i = 0; i < n; ++i) {
                            integration->calc_weight_smearing(nk, nk, kmap_identity, weight[i], e_tmp[j], energy[i], smearing_method);
                            for (k = 0; k < nk; ++k) {
                                tdos_q[i][j] += weight[i][k];
                            }
                        }
                    }
                }

                for (i = 0; i < n; ++i) {
                    for (j = 0; j < 4; ++j) {
                        ret[i][j] += multi * tdos_q[i][j];
                    }
                }
            }
        }
    }
    memory->deallocate(e_tmp);
    memory->deallocate(weight);
    memory->deallocate(kmap_identity);
    memory->deallocate(tdos_q);

    if (mympi->my_rank == 0) {
        std::cout << "done." << std::endl;
    }  
}
