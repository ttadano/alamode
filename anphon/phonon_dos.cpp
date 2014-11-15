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
        if (scattering_phase_space) {
            memory->deallocate(sps3_mode);
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
    MPI_Bcast(&scattering_phase_space, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD);

    if (kpoint->kpoint_mode == 2) {
        flag_dos = true;
    } else {
        flag_dos = false;
    }

    if (flag_dos && delta_e < eps12) error->exit("dos_setup()", "Too small delta_e");

    if (flag_dos) {
        n_energy = static_cast<int>((emax - emin) / delta_e);
        memory->allocate(energy_dos, n_energy);
        memory->allocate(dos_phonon, n_energy);

        for (i = 0; i < n_energy; ++i) {
            energy_dos[i] = emin + delta_e * static_cast<double>(i);
        }

        if (projected_dos) {
            memory->allocate(pdos_phonon, system->natmin, n_energy);
        }

        if (two_phonon_dos) {
            memory->allocate(dos2_phonon, kpoint->nk_reduced, n_energy, 4);
        }

        if (scattering_phase_space) {
            memory->allocate(sps3_mode, kpoint->nk_reduced, dynamical->neval);
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
        calc_two_phonon_dos(n_energy, energy_dos, dos2_phonon, integration->ismear, kpoint->kpoint_irred_all);       
    }

    if (scattering_phase_space) {
        calc_total_scattering_phase_space(dynamical->eval_phonon, integration->ismear, 
            kpoint->kpoint_irred_all, sps3_mode, total_sps3);
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

            ret[i] = 0.0;

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


void Dos::calc_two_phonon_dos(const unsigned int n, double *energy, double ***ret, const int smearing_method,
                              std::vector<std::vector<KpointList> > kpinfo)
{
    int i, j;
    int is, js, ik, jk;
    int k;
    int ib;
    int knum;

    unsigned int nk = kpoint->nk;
    unsigned int ns = dynamical->neval;
    unsigned int nk_reduced = kpoint->nk_reduced;

    int ns2 = ns * ns;

    int *kmap_identity;

    double multi;
    double **e_tmp;    
    double **weight;
    double emax2 = 2.0 * emax;

    double xk_tmp[3];

    int loc;
    int *k_pair;

    if (mympi->my_rank == 0) {
        std::cout << " TDOS = 1 : Calculating two-phonon DOS for all irreducible k points." << std::endl;
        std::cout << "            This may take a while ... ";
    }

    memory->allocate(kmap_identity, nk);
    memory->allocate(e_tmp, 2, nk);
    memory->allocate(weight, n, nk);
    memory->allocate(k_pair, nk);


    for (i = 0; i < nk; ++i) kmap_identity[i] = i;

    for (ik = 0; ik < nk_reduced; ++ik) {

        knum = kpinfo[ik][0].knum;
   //     multi = static_cast<double>(kpinfo[ik].size()) / static_cast<double>(nk);

        for (jk = 0; jk < nk; ++jk) {
            for (i = 0; i < 3; ++i) xk_tmp[i] = kpoint->xk[knum][i] + kpoint->xk[jk][i];
            k_pair[jk] = kpoint->get_knum(xk_tmp[0], xk_tmp[1], xk_tmp[2]);
        }

        for (i = 0; i < n; ++i) {
            for (j = 0; j < 2; ++j) {
                ret[ik][i][j] = 0.0;
            }
        }

        for (ib = 0; ib < ns2; ++ib) {

            is = ib / ns;
            js = ib % ns;

#pragma omp parallel for private(loc)
            for (jk = 0; jk < nk; ++jk) {
                loc = k_pair[jk];
                e_tmp[0][jk] = writes->in_kayser(dynamical->eval_phonon[jk][is] + dynamical->eval_phonon[loc][js]);
                e_tmp[1][jk] = writes->in_kayser(dynamical->eval_phonon[jk][is] - dynamical->eval_phonon[loc][js]);
            }


            if (smearing_method == -1) {

                for (j = 0; j < 2; ++j) {
#pragma omp parallel for private(k)
                    for (i = 0; i < n; ++i) {
                        integration->calc_weight_tetrahedron(nk, kmap_identity, weight[i], e_tmp[j], energy[i]);
                        for (k = 0; k < nk; ++k) {
                            ret[ik][i][j] += weight[i][k];
                        }
                    }
                }

            } else {

                for (j = 0; j < 2; ++j) {
#pragma omp parallel for private(k)
                    for (i = 0; i < n; ++i) {
                        integration->calc_weight_smearing(nk, nk, kmap_identity, weight[i], e_tmp[j], energy[i], smearing_method);
                        for (k = 0; k < nk; ++k) {
                            ret[ik][i][j] += weight[i][k];
                        }
                    }
                }
            }
        }
    }

    memory->deallocate(e_tmp);
    memory->deallocate(weight);
    memory->deallocate(kmap_identity);
    memory->deallocate(k_pair);

    if (mympi->my_rank == 0) {
        std::cout << "done." << std::endl;
    }  
}

void Dos::calc_total_scattering_phase_space(double **omega, const int smearing_method, 
                                            std::vector<std::vector<KpointList> > kpinfo, double **ret_mode, double &ret)
{
    int i, j, k;
    int is, ik;
    int knum;

    unsigned int nk = kpoint->nk;
    unsigned int ns = dynamical->neval;
    int ns2 = ns * ns;
    int ib;

    int *kmap_identity;

    double multi;
    double omega0;
    double sps_tmp1, sps_tmp2;
    double sps_sum1, sps_sum2;

    if (mympi->my_rank == 0) {
        std::cout << " SPS = 1 : Calculating three-phonon scattering phase space ..." << std::endl;
    }

    memory->allocate(kmap_identity, nk);

    for (i = 0; i < nk; ++i) kmap_identity[i] = i;

    ret = 0.0;
    sps_sum1 = 0.0;
    sps_sum2 = 0.0;

    for (ik = 0; ik < kpinfo.size(); ++ik) {

        knum = kpinfo[ik][0].knum;
        multi = static_cast<double>(kpinfo[ik].size()) / static_cast<double>(nk);

        for (is = 0; is < ns; ++is) {

            omega0 = writes->in_kayser(omega[knum][is]);
  //          omega0 = omega[knum][is];

            sps_tmp1 = 0.0;
            sps_tmp2 = 0.0;

#pragma omp parallel 
            {
                double **e_tmp;    
                double *weight;
                int js, ks;
                int jk, loc;
                double xk_tmp[3];

                memory->allocate(weight, nk);
                memory->allocate(e_tmp, 2, nk);

#pragma omp for private(i, j), reduction(+: sps_tmp1, sps_tmp2)
                for (ib = 0; ib < ns2; ++ib) {

                    js = ib / ns;
                    ks = ib % ns; 

                    for (jk = 0; jk < nk; ++jk) {

                        for (i = 0; i < 3; ++i) xk_tmp[i] = kpoint->xk[knum][i] + kpoint->xk[jk][i];
                        loc = kpoint->get_knum(xk_tmp[0], xk_tmp[1], xk_tmp[2]);

                        e_tmp[0][jk] =  writes->in_kayser(omega[jk][js] + omega[loc][ks]);
                        e_tmp[1][jk] =  writes->in_kayser(omega[jk][js] - omega[loc][ks]);

                       //   e_tmp[0][jk] =  omega[jk][js] + omega[loc][ks];
                      //    e_tmp[1][jk] =  omega[jk][js] - omega[loc][ks];
                    }

                    if (smearing_method == -1) {

                        integration->calc_weight_tetrahedron(nk, kmap_identity, weight, e_tmp[0], omega0);
                        for (j = 0; j < nk; ++j) sps_tmp1 += weight[j];
                        integration->calc_weight_tetrahedron(nk, kmap_identity, weight, e_tmp[1], omega0);
                        for (j = 0; j < nk; ++j) sps_tmp2 += weight[j];


                    } else {

                        integration->calc_weight_smearing(nk, nk, kmap_identity, weight, e_tmp[0], omega0,
                            smearing_method);
                        for (j = 0; j < nk; ++j) sps_tmp1 += weight[j];
                        integration->calc_weight_smearing(nk, nk, kmap_identity, weight, e_tmp[1], omega0,
                            smearing_method);
                        for (j = 0; j < nk; ++j) sps_tmp2 += weight[j];

                    }

                }
                memory->deallocate(e_tmp);
                memory->deallocate(weight);
            }
            sps_sum1 += multi * sps_tmp1;
            sps_sum2 += multi * sps_tmp2;

            ret_mode[ik][is] = (sps_tmp1 + 2.0 * sps_tmp2) / (3.0 * static_cast<double>(std::pow(ns,3)));
        }   
    }

    memory->deallocate(kmap_identity);

    ret = (sps_sum1 + 2.0 * sps_sum2) / (3.0 * static_cast<double>(std::pow(ns, 3)));

    if (mympi->my_rank == 0) {
        std::cout << "done." << std::endl;
    }  
}
