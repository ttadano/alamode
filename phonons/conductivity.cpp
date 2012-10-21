#include "mpi_common.h"
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

    if (mympi->my_rank == 0) {
        std::cout.setf(std::ios::fixed);
        std::cout << std::endl;
        std::cout << " Tmin = " << std::setw(10) << system->Tmin; 
        std::cout << " Tmax = " << std::setw(10) << system->Tmax; 
        std::cout << " dT   = " << std::setw(10) << system->dT; 
        std::cout << std::endl;

        std::cout.unsetf(std::ios::fixed);
    }
}

void Conductivity::finish_kl()
{
    memory->deallocate(vel);
    memory->deallocate(tau);
}

void Conductivity::calc_kl()
{
    unsigned int iT, i, j, k;

    double Tmax = system->Tmax;
    double Tmin = system->Tmin;
    double dT = system->dT;

    std::string file_kl;
    std::ofstream ofs_kl;

    unsigned int NT = static_cast<unsigned int>((Tmax - Tmin) / dT);

    double ***kl_T_loc;
    double ***kl_T;
    double *T_arr;

    memory->allocate(T_arr, NT);
    memory->allocate(kl_T_loc, NT, 3, 3);
    memory->allocate(kl_T, NT, 3, 3);

    for (iT = 0; iT < NT; ++iT){
        T_arr[iT] = Tmin  + static_cast<double>(iT)*dT;
    }

    if (mympi->my_rank == 0) {
        file_kl = input->job_title + ".kl";
        ofs_kl.open(file_kl.c_str(), std::ios::out);
        if(!ofs_kl) error->exit("calc_kl", "cannot open file_kl");

        ofs_kl << "# Temperature [K], Thermal Conductivity (xx, xy, xz, yx, yy, yz, zx, zy, zz) [W/mK]" << std::endl;
    }

   
    unsigned int nk_g, nk_l;
    unsigned int *k_local;
    unsigned int *nk_equiv_l;
    unsigned int **k_reduced_l;
    double *wk_l;
    unsigned int ik_loc;


    nk_g = kpoint->nk_reduced;
    nk_l = nk_g / mympi->nprocs;

    int nrem = nk_g - nk_l * mympi->nprocs;
    if (mympi->my_rank + 1 <= nrem) ++nk_l;

    memory->allocate(k_local, nk_l);
    memory->allocate(wk_l, nk_l);
    memory->allocate(nk_equiv_l, nk_l);
    memory->allocate(k_reduced_l, nk_l, kpoint->nequiv_max);

    ik_loc = 0;

    for (i = 0; i < nk_g; ++i) {
    if (i % mympi->nprocs == mympi->my_rank){
    k_local[ik_loc] = kpoint->k_reduced[i][0];
    nk_equiv_l[ik_loc] = kpoint->nk_equiv_arr[i];
    for (k = 0; k < nk_equiv_l[ik_loc]; ++k){
    k_reduced_l[ik_loc][k] = kpoint->k_reduced[i][k];
    }
    ++ik_loc;
    }
    }

    for (i = 0; i < nk_l; ++i){
    wk_l[i] = static_cast<double>(nk_equiv_l[i]) / static_cast<double>(kpoint->nk);
    }

    calc_kl_mpi(nk_l, k_local, wk_l, nk_equiv_l, k_reduced_l, NT, T_arr, kl_T_loc);
    MPI_Reduce(&kl_T_loc[0][0][0], &kl_T[0][0][0], 9*NT, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    memory->deallocate(k_local);
    memory->deallocate(wk_l);
    memory->deallocate(nk_equiv_l);
    memory->deallocate(k_reduced_l);


/*
    unsigned int **ks_local;
    unsigned int nks_g, nks_l;
    double *wks_local;
    unsigned int iks_local;
    unsigned int *nk_equiv_local;
    unsigned int **k_reduced_local;

    nks_g = kpoint->nk_reduced * ns;
    nks_l = nks_g / mympi->nprocs;

    int nrem = nks_g - nks_l * mympi->nprocs;
    if (mympi->my_rank + 1 <= nrem) ++nks_l;

    memory->allocate(ks_local, nks_l, 2);
    memory->allocate(wks_local, nks_l);
    memory->allocate(nk_equiv_local, nks_l);
    memory->allocate(k_reduced_local, nks_l, kpoint->nequiv_max);

    iks_local = 0;
    for (i = 0; i < nks_g; ++i) {
        if (i % mympi->nprocs == mympi->my_rank) {
            ks_local[iks_local][0] = kpoint->k_reduced[i / ns][0];
            ks_local[iks_local][1] = i % ns;
            nk_equiv_local[iks_local] = kpoint->nk_equiv_arr[i / ns];
            for (k = 0; k < nk_equiv_local[iks_local]; ++k){
                k_reduced_local[iks_local][k] = kpoint->k_reduced[i / ns][k];
            }
            ++iks_local;
        }
    }

    for (i = 0; i < nks_l; ++i){
        wks_local[i] = static_cast<double>(nk_equiv_local[i]) / static_cast<double>(kpoint->nk);
    }

    calc_kl_mpi2(nks_l, ks_local, wks_local, nk_equiv_local, k_reduced_local, NT, T_arr, kl_T_loc);
    MPI_Reduce(&kl_T_loc[0][0][0], &kl_T[0][0][0], 9*NT, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    */

    for (iT = 0; iT < NT; ++iT){
        if (mympi->my_rank == 0) {
            ofs_kl << std::setw(5) << T_arr[iT];
            for (i = 0; i < 3; ++i){
                for (j = 0; j < 3; ++j){
                    ofs_kl << std::setw(15) << kl_T[iT][i][j] / (Bohr_in_Angstrom * 1.0e-10 * time_ry * system->volume_p);
                }
            }
            ofs_kl << std::endl;
        }
    }


    if (mympi->my_rank == 0) ofs_kl.close();

    memory->deallocate(kl_T_loc);
    memory->deallocate(kl_T);
    memory->deallocate(T_arr);
}

void Conductivity::calc_kl_mpi(const unsigned int nk_local, unsigned int *k_local, double *weight, 
    unsigned int *nk_equivalent, unsigned int **k_equivalent, const unsigned int NT, double *temperature, double ***kl)
{
    unsigned int i, j, k;
    unsigned int iT;
    unsigned int ik, is;
    unsigned int knum, ktmp;
    double omega;
    double vv_tmp;
    double *damping;

    memory->allocate(damping, NT);

    for (iT = 0; iT < NT; ++iT){
        for (i = 0; i < 3; ++i){
            for (j = 0; j < 3; ++j){
                kl[iT][i][j] = 0.0;
            }
        }
    }

    for (ik = 0; ik < nk_local; ++ik){

        knum = k_local[ik];

        for (is = 0; is < ns; ++is){

            omega = dynamical->eval_phonon[knum][is];

            if (relaxation->ksum_mode == 0) {
                relaxation->calc_damping(NT, temperature, omega, knum, is, damping);
            } else if (relaxation->ksum_mode == -1) {
                relaxation->calc_damping_tetra(NT, temperature, omega, knum, is, damping);
            }

            for (i = 0; i < 3; ++i){
                for (j = 0; j < 3; ++j){

                    vv_tmp = 0.0;

                    for (k = 0; k < nk_equivalent[ik]; ++k){
                        ktmp = k_equivalent[ik][k];
                        vv_tmp += vel[ktmp][is][i] * vel[ktmp][is][j];
                    }
                    vv_tmp /= static_cast<double>(nk_equivalent[ik]);

                    for (iT = 0; iT < NT; ++iT){
                        kl[iT][i][j] += weight[ik] * phonon_thermodynamics->Cv(omega, temperature[iT]) * vv_tmp / (2.0 * damping[iT]);
                    }
                }
            }

        }
    }

    memory->deallocate(damping);
}

void Conductivity::calc_kl_mpi2(const unsigned int nks_local, unsigned int **ks_local, double *weight, 
    unsigned int *nk_equivalent, unsigned int **k_equivalent, const unsigned int NT, double *temperature, double ***kl)
{
    unsigned int i, j, k;
    unsigned int iT;
    unsigned int ik, is;
    unsigned int knum, snum, ktmp;
    double omega;
    double vv_tmp;
    double *damping;

    memory->allocate(damping, NT);

    for (iT = 0; iT < NT; ++iT){
        for (i = 0; i < 3; ++i){
            for (j = 0; j < 3; ++j){
                kl[iT][i][j] = 0.0;
            }
        }
    }

    for (ik = 0; ik < nks_local; ++ik){

        knum = ks_local[ik][0];
        snum = ks_local[ik][1];

        omega = dynamical->eval_phonon[knum][snum];

        if (relaxation->ksum_mode == 0) {
            relaxation->calc_damping(NT, temperature, omega, knum, snum, damping);
        } else if (relaxation->ksum_mode == -1) {
            relaxation->calc_damping_tetra(NT, temperature, omega, knum, snum, damping);
        }

        for (i = 0; i < 3; ++i){
            for (j = 0; j < 3; ++j){

                vv_tmp = 0.0;

                for (k = 0; k < nk_equivalent[ik]; ++k){
                    ktmp = k_equivalent[ik][k];
                    vv_tmp += vel[ktmp][snum][i] * vel[ktmp][snum][j];
                }
                vv_tmp /= static_cast<double>(nk_equivalent[ik]);

                for (iT = 0; iT < NT; ++iT){
                    kl[iT][i][j] += weight[ik] * phonon_thermodynamics->Cv(omega, temperature[iT]) * vv_tmp / (2.0 * damping[iT]);
                }
            }

        }
    }

    memory->deallocate(damping);
}



void Conductivity::calc_kl_at_T(const double T, double kl[3][3])
{
    unsigned int i, j;
    unsigned int is, ik;
    unsigned int jk, kk;
    double omega;
    unsigned int knum;

    std::complex<double> tmp1;
    unsigned int ktmp;

    unsigned int ikIBZ = 0, nsame = 0;

    kk = 0;

    for (i = 0; i < 3; ++i){
        for (j = 0; j < 3; ++j){
            kl[i][j] = 0.0;
        }
    }

    double vv_tmp;

    jk = 0;

    for (ik = 0; ik < kpoint->nk_equiv.size(); ++ik){

        knum = kpoint->kpIBZ[jk].knum;

        for (is = 0; is < ns; ++is){

            omega = dynamical->eval_phonon[knum][is];
            //tau[knum][is] = 1.0 / (2.0 * relaxation->self_E[ns * knum + is].imag());
            tau[knum][is] = 1.0 / (2.0 * relaxation->selfenergy(T, omega, knum, is).imag());
            //          tau[ik][is] = 1.0 / (2.0 * relaxation->self_tetra(T, omega, ik, is));

            for (i = 0; i < 3; ++i){
                for (j = 0; j < 3; ++j){

                    vv_tmp = 0.0;

                    for (kk = 0; kk < kpoint->nk_equiv[ik]; ++kk){
                        ktmp = kpoint->kpIBZ[jk + kk].knum;
                        vv_tmp += vel[ktmp][is][i] * vel[ktmp][is][j];
                    }
                    vv_tmp /= static_cast<double>(kpoint->nk_equiv[ik]);

                    kl[i][j] += kpoint->weight_k[ik] * phonon_thermodynamics->Cv(omega, T) * vv_tmp * tau[knum][is];
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
