#include "mpi_common.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include "conductivity.h"
#include "dynamical.h"
#include "error.h"
#include "integration.h"
#include "kpoint.h"
#include "memory.h"
#include "parsephon.h"
#include "phonon_thermodynamics.h"
#include "phonon_velocity.h"
#include "relaxation.h"
#include "system.h"
#include "write_phonons.h"
#include "../alm_c++/constants.h"

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

        // Generate phonon velocity in Cartesian coordinate
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

        if (use_classical_Cv == 1) {
            std::cout << "Heat capacity will be replaced by kB (classical limit)" << std::endl;
        }
    }

    MPI_Bcast(&use_classical_Cv, 1, MPI_INT, 0, MPI_COMM_WORLD);
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

    double ***kl_T_loc, ***kl_T;
    double *T_arr;
    double *wks_local;

    unsigned int nks_g, nks_l;
    unsigned int iks_local;
    unsigned int **ks_local;
    unsigned int *nk_equiv_local;
    unsigned int **k_reduced_local;

    memory->allocate(T_arr, NT);
    memory->allocate(kl_T_loc, NT, 3, 3);
    memory->allocate(kl_T, NT, 3, 3);

    // Distribute (k,s) to individual MPI threads
    nks_g = kpoint->nk_reduced * ns;
    nks_l = nks_g / mympi->nprocs;

    int nrem = nks_g - nks_l * mympi->nprocs;
    if (mympi->my_rank + 1 <= nrem) ++nks_l;

    memory->allocate(ks_local, nks_l, 2);
    memory->allocate(wks_local, nks_l);
    memory->allocate(nk_equiv_local, nks_l);
    memory->allocate(k_reduced_local, nks_l, kpoint->nequiv_max);

    for (iT = 0; iT < NT; ++iT){
        T_arr[iT] = Tmin  + static_cast<double>(iT)*dT;
    }

    if (mympi->my_rank == 0) {
        file_kl = input->job_title + ".kl";
        ofs_kl.open(file_kl.c_str(), std::ios::out);
        if(!ofs_kl) error->exit("calc_kl", "cannot open file_kl");

        ofs_kl << "# Temperature [K], Thermal Conductivity (xx, xy, xz, yx, yy, yz, zx, zy, zz) [W/mK]" << std::endl;

        std::cout << " The number of total (k,s) = " << nks_g << std::endl;
    }

    std::cout << " #RANK = " << std::setw(5) << mympi->my_rank << " , nks = " << std::setw(5) << nks_l << std::endl;

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

    calc_kl_mpi(nks_l, ks_local, wks_local, nk_equiv_local, k_reduced_local, NT, T_arr, kl_T_loc);
    MPI_Reduce(&kl_T_loc[0][0][0], &kl_T[0][0][0], 9*NT, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

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

    memory->deallocate(ks_local);
    memory->deallocate(wks_local);
    memory->deallocate(nk_equiv_local);
    memory->deallocate(k_reduced_local);

    memory->deallocate(kl_T_loc);
    memory->deallocate(kl_T);
    memory->deallocate(T_arr);
}

void Conductivity::calc_kl_mpi(const unsigned int nks_local, unsigned int **ks_local, double *weight, 
    unsigned int *nk_equivalent, unsigned int **k_equivalent, const unsigned int NT, double *temperature, double ***kl)
{
    unsigned int i, j, k;
    unsigned int iT;
    unsigned int ik;
    unsigned int knum, snum, ktmp;
    double omega;
    double vv_tmp;
    double *damping;

    std::string file_tau;
    std::stringstream str_rank;

    str_rank << mympi->my_rank;
    file_tau = input->job_title + ".tau_" + str_rank.str();

    std::ofstream ofs_tau;

    ofs_tau.open(file_tau.c_str(), std::ios::out);
    

    //MPI_File thefile;
    //MPI_Offset bufsize = NT;
    //MPI_Offset disp = mympi->my_rank*bufsize*sizeof(double);

    //MPI_File_open(MPI_COMM_WORLD, const_cast<char*>(file_tau.c_str()), MPI_MODE_WRONLY, MPI_INFO_NULL, &thefile);
    //MPI_File_set_view(thefile, disp, MPI_DOUBLE, MPI_DOUBLE, "native", MPI_INFO_NULL);

    std::cout << file_tau << std::endl;

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

        if (relaxation->ksum_mode == 0 || relaxation->ksum_mode == 1) {
            relaxation->calc_damping(NT, temperature, omega, knum, snum, damping);
        } else if (relaxation->ksum_mode == -1) {
            relaxation->calc_damping_tetra(NT, temperature, omega, knum, snum, damping);
        }

       // MPI_File_write(thefile, damping, bufsize, MPI_DOUBLE, MPI_STATUS_IGNORE);

        ofs_tau << "# Relaxation time [ps] of a phonon at xk = ";
        for (i = 0; i < 3; ++i) {
            ofs_tau << std::setw(15) << kpoint->xk[knum][i];
        }
        ofs_tau << std::endl;
        ofs_tau << "# Branch = " << snum << std::endl;

        for (i = 0; i < NT; ++i){
            ofs_tau << std::setw(8) << knum << std::setw(8) << snum;
            ofs_tau << std::setw(15) << writes->in_kayser(omega);
            ofs_tau << std::setw(15) << temperature[i] << std::setw(15) << time_ry / ( 2.0 * damping[i]) * 1.0e+12 << std::endl;
        }
        ofs_tau << std::endl;

        if (mympi->my_rank == 0) {
            std::cout <<  "ELEMENT " << std::setw(5) << ik + 1 << " done." << std::endl;
        }

        for (i = 0; i < 3; ++i){
            for (j = 0; j < 3; ++j){

                vv_tmp = 0.0;

                for (k = 0; k < nk_equivalent[ik]; ++k){
                    ktmp = k_equivalent[ik][k];
                    vv_tmp += vel[ktmp][snum][i] * vel[ktmp][snum][j];
                }
                vv_tmp /= static_cast<double>(nk_equivalent[ik]);

                if (use_classical_Cv == 1) {
                    for (iT = 0; iT < NT; ++iT){
                        kl[iT][i][j] += weight[ik] * phonon_thermodynamics->Cv_classical(omega, temperature[iT]) * vv_tmp / (2.0 * damping[iT]);
                    }
                } else {
                    for (iT = 0; iT < NT; ++iT){
                        kl[iT][i][j] += weight[ik] * phonon_thermodynamics->Cv(omega, temperature[iT]) * vv_tmp / (2.0 * damping[iT]);
                    }
                }
            }

        }
    }

   // MPI_File_close(&thefile);
    ofs_tau.close();

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
