#include "mpi_common.h"
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <vector>
#include "relaxation.h"
#include "memory.h"
#include "kpoint.h"
#include "system.h"
#include "dynamical.h"
#include "fcs_phonon.h"
#include "symmetry_core.h"
#include "../alm_c++/constants.h"
#include "phonon_thermodynamics.h"
#include "phonon_velocity.h"
#include "parsephon.h"
#include "error.h"
#include "conductivity.h"
#include "write_phonons.h"
#include "integration.h"
#include "parsephon.h"
#include "phonon_dos.h"
#include <omp.h>
#include "timer.h"

using namespace PHON_NS;

Relaxation::Relaxation(PHON *phon): Pointers(phon) {
    im = std::complex<double>(0.0, 1.0);
}

Relaxation::~Relaxation(){};

void Relaxation::setup_relaxation()
{
    nk = kpoint->nk;
    ns = dynamical->neval;
    nks = ns*nk;

    memory->allocate(V, 1);
    memory->allocate(self_E, nks);

    unsigned int i, j, k;

    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j){
            mat_convert[i][j] = 0.0;
            for (k = 0; k < 3; ++k){
                mat_convert[i][j] += system->rlavec_p[i][k] * system->lavec_s[k][j]; 
            }
        }
    }

    memory->allocate(relvec, system->nat, system->nat, 3);
    memory->allocate(invsqrt_mass_p, system->natmin);

    if (mympi->my_rank == 0) {
        double vec[3];
        memory->allocate(vec_s, system->ntran, 3);

        for (i = 0; i < system->ntran; ++i){
            for (j = 0; j < 3; ++j){
                vec_s[i][j] = system->xr_s[system->map_p2s[0][i]][j];
            }
        }

        for (i = 0; i < system->nat; ++i){
            for (j = 0; j < system->nat; ++j){

                for (k = 0; k < 3; ++k){

                    if (system->cell_dimension[k] == 1) {

                        vec[k] = system->xr_s[i][k] - system->xr_s[j][k];
                        /*
                        if (vec[k] <= -0.5) {
                        vec[k] = -1.0;
                        } else if (vec[k] > 0.5){
                        vec[k] = 1.0;
                        } else {
                        vec[k] = 0.0;
                        }
                        */
                        if (std::abs(vec[k]) < 0.5) {
                            vec[k] = 0.0;
                        } else {
                            if (system->xr_s[i][k] < 0.5) {
                                vec[k] = -1.0;
                            } else {
                                vec[k] = 1.0;
                            }
                        }

                    } else if (system->cell_dimension[k] == 2) {

                        vec[k] = system->xr_s[system->map_p2s[0][system->map_s2p[i].tran_num]][k] 
                        - system->xr_s[system->map_p2s[0][system->map_s2p[j].tran_num]][k];  
                        vec[k] = dynamical->fold(vec[k]);
                        if (std::abs(system->xr_s[i][k] - system->xr_s[j][k]) > 0.5) vec[k] *= -1.0;

                    } else {

                        vec[k] = system->xr_s[system->map_p2s[0][system->map_s2p[i].tran_num]][k] 
                        - system->xr_s[system->map_p2s[0][system->map_s2p[j].tran_num]][k];  
                        vec[k] = dynamical->fold(vec[k]);

                    }

                    relvec[i][j][k] = vec[k];
                }
                system->rotvec(relvec[i][j], relvec[i][j], mat_convert);
            }
        }
        memory->deallocate(vec_s);
    }
    MPI_Bcast(&relvec[0][0][0], 3*system->nat*system->nat, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    for (i = 0; i < system->natmin; ++i){
        invsqrt_mass_p[i] = std::sqrt(1.0 / system->mass[system->map_p2s[i][0]]);
    }

    memory->allocate(vec_for_v3, fcs_phonon->force_constant[1].size(), 3, 2);
    memory->allocate(invmass_for_v3, fcs_phonon->force_constant[1].size());

    j = 0;
    unsigned int atom_num[3];

    for (std::vector<FcsClass>::const_iterator it = fcs_phonon->force_constant[1].begin(); it != fcs_phonon->force_constant[1].end(); ++it) {

        for (i = 0; i < 3; ++i) atom_num[i] = system->map_p2s[(*it).elems[i].atom][(*it).elems[i].cell];

        for (i = 0; i < 3; ++i) {
            vec_for_v3[j][i][0] = relvec[atom_num[1]][atom_num[0]][i];
            vec_for_v3[j][i][1] = relvec[atom_num[2]][atom_num[0]][i];
        }

        invmass_for_v3[j] =  invsqrt_mass_p[(*it).elems[0].atom] * invsqrt_mass_p[(*it).elems[1].atom] * invsqrt_mass_p[(*it).elems[2].atom];
        ++j;     
    }

    // Tentative modification for tuning
    memory->allocate(cexp_phase, nk, fcs_phonon->force_constant[1].size(), 2);

    for (i = 0; i < nk; ++i){
        for (j = 0; j < fcs_phonon->force_constant[1].size(); ++j){
         cexp_phase[i][j][0] =  std::exp(im * (vec_for_v3[j][0][0] * kpoint->xk[i][0] + vec_for_v3[j][1][0] * kpoint->xk[i][1] + vec_for_v3[j][2][0] * kpoint->xk[i][2]));
         cexp_phase[i][j][1] =  std::exp(im * (vec_for_v3[j][0][1] * kpoint->xk[i][0] + vec_for_v3[j][1][1] * kpoint->xk[i][1] + vec_for_v3[j][2][1] * kpoint->xk[i][2]));
        }
    }

    // For tetrahedron method
    memory->allocate(e_tmp, 4, nk);
    memory->allocate(f_tmp, 4, nk);

    if (mympi->my_rank == 0) {
        double domega_min;

        domega_min = 10000.0;

        domega_min = std::min<double>(domega_min, writes->in_kayser(dynamical->eval_phonon[1][0]));
        domega_min = std::min<double>(domega_min, writes->in_kayser(dynamical->eval_phonon[kpoint->nkz][0]));
        domega_min = std::min<double>(domega_min, writes->in_kayser(dynamical->eval_phonon[kpoint->nkz*kpoint->nky][0]));

        std::cout << std::endl;
        std::cout << "Estimated minimum energy difference (cm^-1) = " << domega_min << std::endl;
        std::cout << "Given epsilon (cm^-1) = " << epsilon << std::endl << std::endl;

        if (ksum_mode == 0) {
            std::cout << "Lorentzian broadning will be used." << std::endl;
        } else if (ksum_mode == 1) {
            std::cout << "Gaussian broadning will be used." << std::endl;
        } else if (ksum_mode == -1) {
            std::cout << "Tetrahedron method will be used." << std::endl;
        } else {
            error->exit("setup_relaxation", "Invalid ksum_mode");
        }
    }

    modify_eigenvectors();

    epsilon *= time_ry / Hz_to_kayser;

    MPI_Bcast(&ksum_mode, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&epsilon, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    /*
    for (unsigned int ik = 0; ik < nk; ++ik){
    for (unsigned int is = 0; is < ns; ++is){ 
    std::cout << "ik = " << ik << "is = " << is << " " << dynamical->eval_phonon[ik][is] << " " << dynamical->eval_phonon[kpoint->knum_minus[ik]][is] << std::endl;
    for (unsigned int js = 0; js < ns; ++js) {
    if (std::abs(dynamical->evec_phonon[ik][is][js] - std::conj(dynamical->evec_phonon[kpoint->knum_minus[ik]][is][js])) > eps10) {
    std::cout << "js = " << js;
    std::cout << " " << dynamical->evec_phonon[ik][is][js];
    std::cout << " " << dynamical->evec_phonon[kpoint->knum_minus[ik]][is][js] << std::endl;
    std::cout << std::endl;
    }
    }
    }
    }
    */
    /*
    unsigned int tmp[3];
    std::complex<double> v3tmp1, v3tmp2;

    for (unsigned int ik = 0; ik < nk; ++ik){
    for (unsigned int jk = 0; jk < nk; ++jk){
    for (unsigned int kk = 0; kk < nk; ++kk){
    std::cout << "k_orig=" << ik << "," << jk << "," << kk << std::endl;
    std::cout << "k_minus=" << kpoint->knum_minus[ik] << "," << kpoint->knum_minus[jk] << "," << kpoint->knum_minus[kk] << std::endl;
    for (unsigned int is = 0; is < ns; ++is){
    for (unsigned int js = 0; js < ns; ++js){
    for (unsigned int ks = 0; ks < ns; ++ks){

    tmp[0] = ns * ik + is;
    tmp[1] = ns * jk + js;
    tmp[2] = ns * kk + ks;

    v3tmp1 = V3new(tmp);

    tmp[0] = ns * kpoint->knum_minus[ik] + is;
    tmp[1] = ns * kpoint->knum_minus[jk] + js; tmp[2] = ns * kpoint->knum_minus[kk] + ks;

    v3tmp2 = V3new(tmp);

    if (std::abs(v3tmp1 - conj(v3tmp2)) > eps10) {
    std::cout << "(" << ik << "," << is;
    std::cout << ";" << jk << "," << js;
    std::cout << ";" << kk << "," << ks;
    std::cout << ") = " << v3tmp1 << v3tmp2 << std::endl;
    }
    }
    }
    }
    }
    }
    }
    error->exit("hoge", "fuga");
    */
}

void Relaxation::finish_relaxation()
{
    V[0].clear();
    memory->deallocate(V);
    memory->deallocate(self_E);

    memory->deallocate(relvec);
    memory->deallocate(invsqrt_mass_p);
    memory->deallocate(vec_for_v3);
    memory->deallocate(invmass_for_v3);
}

void Relaxation::calc_ReciprocalV()
{
    // Calculate V_{ks, k's', k''s''}^(3) for Self-energy calculation

    int ik;
    int k1, k2, k3;
    unsigned int i;
    unsigned int b1, b2, b3;
    unsigned int nkp;

    double xk_tmp[3], xk_norm;

    bool is_necessary;

    StructKS ks_tmp;
    std::vector<StructKS> kslist;
    std::vector<FcsClass>::iterator it;

    std::cout << std::endl;
    std::cout << "Calculating force constants in reciprocal space .." << std::endl;

    nkp = 0;

    kslist.clear();

    for (k1 = 0; k1 < nk; ++k1){
        for (k2 = k1; k2 < nk; ++k2){
            for (k3 = k2; k3 < nk; ++k3){

                is_necessary 
                    = (kpoint->kpset_uniq.find(kpoint->knum_minus[k1]) != kpoint->kpset_uniq.end()) 
                    || (kpoint->kpset_uniq.find(kpoint->knum_minus[k2]) != kpoint->kpset_uniq.end())
                    || (kpoint->kpset_uniq.find(kpoint->knum_minus[k3]) != kpoint->kpset_uniq.end());

                if(!is_necessary) continue;

                xk_tmp[0] = kpoint->xk[k1][0] + kpoint->xk[k2][0] + kpoint->xk[k3][0];
                xk_tmp[1] = kpoint->xk[k1][1] + kpoint->xk[k2][1] + kpoint->xk[k3][1];
                xk_tmp[2] = kpoint->xk[k1][2] + kpoint->xk[k2][2] + kpoint->xk[k3][2];

                for (i = 0; i < 3; ++i){
                    xk_tmp[i] = std::fmod(xk_tmp[i], 1.0);
                }
                xk_norm = std::pow(xk_tmp[0], 2) + std::pow(xk_tmp[1], 2) + std::pow(xk_tmp[2], 2);

                // If the momentum-conservation is not satisfied, we skip the loop.

                if (std::sqrt(xk_norm) > eps15) continue;

                for (b1 = 0; b1 < ns; ++b1){
                    for (b2 = 0; b2 < ns; ++b2){
                        for (b3 = 0; b3 < ns; ++b3){

                            ks_tmp.ks1 = ns * k1 + b1;
                            ks_tmp.ks2 = ns * k2 + b2;
                            ks_tmp.ks3 = ns * k3 + b3;

                            if (ks_tmp.ks1 > ks_tmp.ks2 || ks_tmp.ks2 > ks_tmp.ks3) continue;

                            kslist.push_back(ks_tmp);

                        }
                    }
                }
            }
        }
    }

    nkp = kslist.size();

#pragma omp parallel
    {
        unsigned int ks_arr[3];
        std::complex<double> prod;

#pragma omp for schedule(static)
        for (ik = 0; ik < nkp; ++ik){

            ks_arr[0] = kslist[ik].ks1;
            ks_arr[1] = kslist[ik].ks2;
            ks_arr[2] = kslist[ik].ks3;
            prod = V3(ks_arr[0], ks_arr[1], ks_arr[2]);

            // Add to list
            if (std::abs(prod) > eps) {
#pragma omp critical
                V[0].push_back(ReciprocalVs(prod, ks_arr, 3));
            }
        }
    }

    kslist.clear();
    std::cout << "Done !" << std::endl;
    std::cout << "Number of nonzero V's: " << V[0].size() << std::endl;
}

std::complex<double> Relaxation::V3(const unsigned int ks1, const unsigned int ks2, const unsigned int ks3)
{
    unsigned int i;
    unsigned int k1, k2, k3;
    unsigned int b1, b2, b3;
    int ielem;

    double omega_prod;
    double omega[3];

    //  std::vector<FcsClass>::iterator it;

    double ret_re = 0.0;
    double ret_im = 0.0;

    k1 = ks1 / ns;
    k2 = ks2 / ns;
    k3 = ks3 / ns;

    b1 = ks1 % ns;
    b2 = ks2 % ns;
    b3 = ks3 % ns;

    omega[0] = dynamical->eval_phonon[k1][b1];
    omega[1] = dynamical->eval_phonon[k2][b2];
    omega[2] = dynamical->eval_phonon[k3][b3];
    omega_prod = omega[0] * omega[1] * omega[2];

#pragma omp parallel
    {
        double vec1[3], vec2[3];
        double invsqrt_mass_prod, phase; 
        std::complex<double> tmp;
        tmp = std::complex<double>(0.0, 0.0);

#pragma omp parallel for reduction(+: ret_re, ret_im)
        for (ielem = 0; ielem < fcs_phonon->force_constant[1].size(); ++ielem){

            for (i = 0; i < 3; ++i){
                vec1[i] = relvec[fcs_phonon->force_constant[1][ielem].elems[1].cell][fcs_phonon->force_constant[1][ielem].elems[0].cell][i];
                vec2[i] = relvec[fcs_phonon->force_constant[1][ielem].elems[2].cell][fcs_phonon->force_constant[1][ielem].elems[0].cell][i];
            }

            phase = vec1[0] * kpoint->xk[k2][0] + vec1[1] * kpoint->xk[k2][1] + vec1[2] * kpoint->xk[k2][2]
            + vec2[0] * kpoint->xk[k3][0] + vec2[1] * kpoint->xk[k3][1] + vec2[2] * kpoint->xk[k3][2];

            invsqrt_mass_prod = invsqrt_mass_p[fcs_phonon->force_constant[1][ielem].elems[0].atom]
            * invsqrt_mass_p[fcs_phonon->force_constant[1][ielem].elems[1].atom] 
            * invsqrt_mass_p[fcs_phonon->force_constant[1][ielem].elems[2].atom];


            tmp = fcs_phonon->force_constant[1][ielem].fcs_val * std::exp(im * phase) * invsqrt_mass_prod
                * dynamical->evec_phonon[k1][b1][3 * fcs_phonon->force_constant[1][ielem].elems[0].atom + fcs_phonon->force_constant[1][ielem].elems[0].xyz]
            * dynamical->evec_phonon[k2][b2][3 * fcs_phonon->force_constant[1][ielem].elems[1].atom + fcs_phonon->force_constant[1][ielem].elems[1].xyz]
            * dynamical->evec_phonon[k3][b3][3 * fcs_phonon->force_constant[1][ielem].elems[2].atom + fcs_phonon->force_constant[1][ielem].elems[2].xyz];
            ret_re += tmp.real();
            ret_im += tmp.imag();
        }
    }

    return (ret_re + im * ret_im) / std::sqrt(omega_prod);
}

std::complex<double> Relaxation::V3new(const unsigned int ks[3])
{
    unsigned int i;
    unsigned int kn[3];
    unsigned int sn[3];

    int ielem;

    double omega[3];

    double ret_re = 0.0;
    double ret_im = 0.0;

    for (i = 0; i < 3; ++i){
        kn[i] = ks[i] / ns;
        sn[i] = ks[i] % ns;
        omega[i] = dynamical->eval_phonon[kn[i]][sn[i]];
    }

#pragma omp parallel 
    {
        std::complex<double> ctmp;

#pragma omp for reduction(+: ret_re, ret_im)
        for (ielem = 0; ielem < fcs_phonon->force_constant[1].size(); ++ielem) {

            ctmp = dynamical->evec_phonon[kn[0]][sn[0]][3 * fcs_phonon->force_constant[1][ielem].elems[0].atom + fcs_phonon->force_constant[1][ielem].elems[0].xyz] 
            * dynamical->evec_phonon[kn[1]][sn[1]][3 * fcs_phonon->force_constant[1][ielem].elems[1].atom + fcs_phonon->force_constant[1][ielem].elems[1].xyz]
            * dynamical->evec_phonon[kn[2]][sn[2]][3 * fcs_phonon->force_constant[1][ielem].elems[2].atom + fcs_phonon->force_constant[1][ielem].elems[2].xyz]
            * fcs_phonon->force_constant[1][ielem].fcs_val * invmass_for_v3[ielem] * cexp_phase[kn[1]][ielem][0] * cexp_phase[kn[2]][ielem][1];

            ret_re += ctmp.real();
            ret_im += ctmp.imag();
        }
    }

    return (ret_re + im * ret_im) / std::sqrt(omega[0] * omega[1] * omega[2]);
}

std::complex<double> Relaxation::V3new2(const unsigned int ks[3])
{
    unsigned int i;
    unsigned int kn[3];
    unsigned int sn[3];

    int ielem;

    double omega[3];

    double ret_re = 0.0;
    double ret_im = 0.0;

    for (i = 0; i < 3; ++i){
        kn[i] = ks[i] / ns;
        sn[i] = ks[i] % ns;
        omega[i] = dynamical->eval_phonon[kn[i]][sn[i]];
    }

#pragma omp parallel 
    {
        double invsqrt_mass_prod;
        double phase;
        std::complex<double> ctmp;
        double vec1[3], vec2[3];
        unsigned int atom_num[3];

#pragma omp for reduction(+: ret_re, ret_im)
        for (ielem = 0; ielem < fcs_phonon->force_constant[1].size(); ++ielem) {

            for (i = 0; i < 3; ++i) atom_num[i] = system->map_p2s[fcs_phonon->force_constant[1][ielem].elems[i].atom][fcs_phonon->force_constant[1][ielem].elems[i].cell];
            for (i = 0; i < 3; ++i) {
                vec1[i] = relvec[atom_num[1]][atom_num[0]][i];
                vec2[i] = relvec[atom_num[2]][atom_num[0]][i];
            }

            phase = 0.0;
            ctmp = std::complex<double>(1.0, 0.0);
            invsqrt_mass_prod = 1.0;

            for (i = 0; i < 3; ++i){
                phase += vec1[i] * kpoint->xk[kn[1]][i] + vec2[i] * kpoint->xk[kn[2]][i];
                invsqrt_mass_prod *= invsqrt_mass_p[fcs_phonon->force_constant[1][ielem].elems[i].atom];
                ctmp *= dynamical->evec_phonon[kn[i]][sn[i]][3 * fcs_phonon->force_constant[1][ielem].elems[i].atom + fcs_phonon->force_constant[1][ielem].elems[i].xyz];
            }

            ctmp *=  fcs_phonon->force_constant[1][ielem].fcs_val * std::exp(im * phase) * invsqrt_mass_prod;

            ret_re += ctmp.real();
            ret_im += ctmp.imag();
        }
    }

    return (ret_re + im * ret_im) / std::sqrt(omega[0] * omega[1] * omega[2]);
}

std::complex<double> Relaxation::selfenergy(const double T, const double omega, const unsigned int knum, const unsigned int snum)
{
    std::complex<double> ret(0.0, 0.0);

    unsigned int ik, jk;
    unsigned int is, js;

    unsigned int i;
    unsigned int knum_inv;
    double ret_re, ret_im;

    knum_inv = kpoint->knum_minus[knum];

#pragma omp parallel private(ik, jk, is, js)
    {
        int iks2;
        unsigned int iks, jks;
        unsigned int arr[3];
        double xk_tmp[3], xk_norm;
        double omega_inner[2];
        double n1, n2;

        std::complex<double> ctmp;

        ret_re = 0.0;
        ret_im = 0.0;
        arr[0] = ns * kpoint->knum_minus[knum] + snum;

#pragma omp for reduction (+:ret_re, ret_im)
        for (iks2 = 0; iks2 < nks*nks; ++iks2){
            iks = iks2 / nks;
            jks = iks2 % nks;
            ik = iks / ns;
            is = iks % ns;
            jk = jks / ns;
            js = jks % ns;

            arr[1] = iks;
            arr[2] = jks;

            xk_tmp[0] = kpoint->xk[knum_inv][0] + kpoint->xk[ik][0] + kpoint->xk[jk][0];
            xk_tmp[1] = kpoint->xk[knum_inv][1] + kpoint->xk[ik][1] + kpoint->xk[jk][1];
            xk_tmp[2] = kpoint->xk[knum_inv][2] + kpoint->xk[ik][2] + kpoint->xk[jk][2];

            for (i = 0; i < 3; ++i)  xk_tmp[i] = std::fmod(xk_tmp[i], 1.0);
            xk_norm = std::pow(xk_tmp[0], 2) + std::pow(xk_tmp[1], 2) + std::pow(xk_tmp[2], 2);
            if (std::sqrt(xk_norm) > eps15) continue; 

            omega_inner[0] = dynamical->eval_phonon[ik][is];
            omega_inner[1] = dynamical->eval_phonon[jk][js];
            n1 = phonon_thermodynamics->fB(omega_inner[0], T) + phonon_thermodynamics->fB(omega_inner[1], T) + 1.0;
            n2 = phonon_thermodynamics->fB(omega_inner[0], T) - phonon_thermodynamics->fB(omega_inner[1], T);

            //            ctmp =  std::norm(V3new2(arr))
            ctmp =  std::norm(V3new(arr))
                * ( n1 / (omega + omega_inner[0] + omega_inner[1] + im * epsilon)
                - n1 / (omega - omega_inner[0] - omega_inner[1] + im * epsilon) 
                + n2 / (omega - omega_inner[0] + omega_inner[1] + im * epsilon)
                - n2 / (omega + omega_inner[0] - omega_inner[1] + im * epsilon));


            ret_re += ctmp.real();
            ret_im += ctmp.imag();
        }
    }

    return (ret_re + im * ret_im) * std::pow(0.5, 4) / static_cast<double>(nk);
}

std::complex<double> Relaxation::selfenergy2(const double T, const double omega, const unsigned int knum, const unsigned int snum)
{
    std::complex<double> ret(0.0, 0.0);

    unsigned int ik, jk;
    unsigned int is, js;

    unsigned int i;
    unsigned int knum_inv;

    knum_inv = kpoint->knum_minus[knum];

    unsigned int arr[3];
    double xk_tmp[3], xk_norm;
    double omega_inner[2];
    double n1, n2;

    std::complex<double> ctmp;

    arr[0] = ns * kpoint->knum_minus[knum] + snum;

    for (ik = 0; ik < nk; ++ik){
        for (jk = 0; jk < nk; ++jk){

            xk_tmp[0] = kpoint->xk[knum_inv][0] + kpoint->xk[ik][0] + kpoint->xk[jk][0];
            xk_tmp[1] = kpoint->xk[knum_inv][1] + kpoint->xk[ik][1] + kpoint->xk[jk][1];
            xk_tmp[2] = kpoint->xk[knum_inv][2] + kpoint->xk[ik][2] + kpoint->xk[jk][2];

            for (i = 0; i < 3; ++i){
                xk_tmp[i] = std::fmod(xk_tmp[i], 1.0);
            }
            xk_norm = std::pow(xk_tmp[0], 2) + std::pow(xk_tmp[1], 2) + std::pow(xk_tmp[2], 2);


            if (std::sqrt(xk_norm) > eps15) continue; 

            for (is = 0; is < ns; ++is){
                for (js = 0; js < ns; ++js){

                    arr[1] = ns * ik + is;
                    arr[2] = ns * jk + js;
                    omega_inner[0] = dynamical->eval_phonon[ik][is];
                    omega_inner[1] = dynamical->eval_phonon[jk][js];
                    n1 = phonon_thermodynamics->fB(omega_inner[0], T) + phonon_thermodynamics->fB(omega_inner[1], T) + 1.0;
                    n2 = phonon_thermodynamics->fB(omega_inner[0], T) - phonon_thermodynamics->fB(omega_inner[1], T);

                    ctmp =  std::norm(V3new(arr))
                        * ( n1 / (omega + omega_inner[0] + omega_inner[1] + im * epsilon)
                        - n1 / (omega - omega_inner[0] - omega_inner[1] + im * epsilon) 
                        + n2 / (omega - omega_inner[0] + omega_inner[1] + im * epsilon)
                        - n2 / (omega + omega_inner[0] - omega_inner[1] + im * epsilon));

                    ret += ctmp;
                }
            }
        }
    }

    return ret * std::pow(0.5, 4) / static_cast<double>(nk); 
}

void Relaxation::calc_damping(const unsigned int N, double *T, const double omega, const unsigned int knum, const unsigned int snum, double *ret)
{
    unsigned int i;
    unsigned int ik, jk;
    unsigned int is, js; 
    unsigned int iks, jks, iks2;
    unsigned int arr[3];

    double T_tmp;
    double n1, n2;
    double v3_tmp;
    double xk_tmp[3], xk_norm;
    double omega_inner[2];

    for (i = 0; i < N; ++i) ret[i] = 0.0;

    arr[0] = ns * kpoint->knum_minus[knum] + snum;

    unsigned int nkx = kpoint->nkx;
    unsigned int nky = kpoint->nky;
    unsigned int nkz = kpoint->nkz;

    int iloc, jloc, kloc;

    double time_tmp;

    for (ik = 0; ik < nk; ++ik) {

        xk_tmp[0] = kpoint->xk[knum][0] - kpoint->xk[ik][0];
        xk_tmp[1] = kpoint->xk[knum][1] - kpoint->xk[ik][1];
        xk_tmp[2] = kpoint->xk[knum][2] - kpoint->xk[ik][2];

        iloc = (kpoint->nint(xk_tmp[0]*static_cast<double>(nkx) + static_cast<double>(2*nkx))) % nkx;
        jloc = (kpoint->nint(xk_tmp[1]*static_cast<double>(nky) + static_cast<double>(2*nky))) % nky;
        kloc = (kpoint->nint(xk_tmp[2]*static_cast<double>(nkz) + static_cast<double>(2*nkz))) % nkz;

        jk = kloc + nkz * jloc + nky * nkz * iloc;

        for (is = 0; is < ns; ++is){
            for (js = 0; js < ns; ++js){

                arr[1] = ns * ik + is;
                arr[2] = ns * jk + js;

                omega_inner[0] = dynamical->eval_phonon[ik][is];
                omega_inner[1] = dynamical->eval_phonon[jk][js];

                v3_tmp = std::norm(V3new(arr));

                for (i = 0; i < N; ++i) {
                    T_tmp = T[i];
                    n1 = phonon_thermodynamics->fB(omega_inner[0], T_tmp) + phonon_thermodynamics->fB(omega_inner[1], T_tmp) + 1.0;
                    n2 = phonon_thermodynamics->fB(omega_inner[0], T_tmp) - phonon_thermodynamics->fB(omega_inner[1], T_tmp);

                    ret[i] += v3_tmp 
                        * (- n1 * delta_lorentz(omega + omega_inner[0] + omega_inner[1])
                        + n1 * delta_lorentz(omega - omega_inner[0] - omega_inner[1])
                        - n2 * delta_lorentz(omega - omega_inner[0] + omega_inner[1])
                        + n2 * delta_lorentz(omega + omega_inner[0] - omega_inner[1]));
                }
            }
        }
    }

    for (i = 0; i < N; ++i) ret[i] *=  pi * std::pow(0.5, 4) / static_cast<double>(nk);
}


void Relaxation::calc_damping_tetra(const unsigned int N, double *T, const double omega, const unsigned int knum, const unsigned int snum, double *ret)
{
    unsigned int i, j;
    unsigned int is, js, ik, jk;
    unsigned int ks_tmp[3];
    unsigned int kcount;

    double xk_tmp[3], xk_norm;
    double n1, n2;
    double *v3_tmp;
    double **omega_inner;

    unsigned int nkx = kpoint->nkx;
    unsigned int nky = kpoint->nky;
    unsigned int nkz = kpoint->nkz;

    int iloc, jloc, kloc;

    for (i = 0; i < N; ++i) ret[i] = 0.0;

    ks_tmp[0] = ns * kpoint->knum_minus[knum] + snum;

    memory->allocate(v3_tmp, nk);
    memory->allocate(omega_inner, nk, 2);

    for (is = 0; is < ns; ++is){
        for (js = 0; js < ns; ++js){

            for (ik = 0; ik < nk; ++ik) {

                xk_tmp[0] = kpoint->xk[knum][0] - kpoint->xk[ik][0];
                xk_tmp[1] = kpoint->xk[knum][1] - kpoint->xk[ik][1];
                xk_tmp[2] = kpoint->xk[knum][2] - kpoint->xk[ik][2];

                iloc = (kpoint->nint(xk_tmp[0]*static_cast<double>(nkx) + static_cast<double>(2*nkx))) % nkx;
                jloc = (kpoint->nint(xk_tmp[1]*static_cast<double>(nky) + static_cast<double>(2*nky))) % nky;
                kloc = (kpoint->nint(xk_tmp[2]*static_cast<double>(nkz) + static_cast<double>(2*nkz))) % nkz;

                jk = kloc + nkz * jloc + nky * nkz * iloc;

                ks_tmp[1] = ik * ns + is;
                ks_tmp[2] = jk * ns + js;

                //    if (kcount == nk) error->exitall("calc_damping_tetra", "The number of element is larger than nk");

                omega_inner[ik][0] = dynamical->eval_phonon[ik][is];
                omega_inner[ik][1] = dynamical->eval_phonon[jk][js];

                v3_tmp[ik] = std::norm(V3new(ks_tmp));

                // e_tmp[0][kcount] = -omega_inner[kcount][0] - omega_inner[kcount][1];
                e_tmp[1][ik] = omega_inner[ik][0] + omega_inner[ik][1];
                e_tmp[2][ik] = omega_inner[ik][0] - omega_inner[ik][1];
                e_tmp[3][ik] = -omega_inner[ik][0] + omega_inner[ik][1];
            }

            for (j = 0; j < N; ++j){
                for (i = 0; i < nk; ++i){
                    n1 = phonon_thermodynamics->fB(omega_inner[i][0], T[j]) + phonon_thermodynamics->fB(omega_inner[i][1], T[j]) + 1.0;
                    n2 = phonon_thermodynamics->fB(omega_inner[i][0], T[j]) - phonon_thermodynamics->fB(omega_inner[i][1], T[j]);
                    // f_tmp[0][i] = -v3_tmp[i] * n1;
                    f_tmp[1][i] = v3_tmp[i] * n1;
                    f_tmp[2][i] = -v3_tmp[i] * n2;
                    f_tmp[3][i] = v3_tmp[i] * n2;
                }

                for (i = 1; i < 4; ++i) {
                    ret[j] += integration->do_tetrahedron(e_tmp[i], f_tmp[i], omega);
                }

            }
        }
    }

    for (i = 0; i < N; ++i) {
        ret[i] *=  pi * std::pow(0.5, 4);
    }

    memory->deallocate(v3_tmp);
    memory->deallocate(omega_inner);
}

double Relaxation::self_tetra(const double T, const double omega, const unsigned int knum, const unsigned int snum)
{
    double ret = 0.0;

    unsigned int i;
    unsigned int is, js, ik, jk;
    unsigned int knum_inv;
    unsigned int ks_tmp[3], kcount;

    double xk_tmp[3], xk_norm;
    double v3_tmp, omega_inner[2];
    double n1, n2;    

    knum_inv= kpoint->knum_minus[knum];

    ks_tmp[0] = knum_inv * ns + snum;

    for (is = 0; is < ns; ++is){
        for (js = 0; js < ns; ++js){

            kcount = 0;

            for (ik = 0; ik < nk; ++ik) {
                for (jk = 0; jk < nk; ++jk){

                    xk_tmp[0] = kpoint->xk[knum_inv][0] + kpoint->xk[ik][0] + kpoint->xk[jk][0];
                    xk_tmp[1] = kpoint->xk[knum_inv][1] + kpoint->xk[ik][1] + kpoint->xk[jk][1];
                    xk_tmp[2] = kpoint->xk[knum_inv][2] + kpoint->xk[ik][2] + kpoint->xk[jk][2];

                    for (i = 0; i < 3; ++i)  xk_tmp[i] = std::fmod(xk_tmp[i], 1.0);
                    xk_norm = std::pow(xk_tmp[0], 2) + std::pow(xk_tmp[1], 2) + std::pow(xk_tmp[2], 2);
                    if (std::sqrt(xk_norm) > eps15) continue; 

                    ks_tmp[1] = ik * ns + is;
                    ks_tmp[2] = jk * ns + js;

                    v3_tmp = std::norm(V3new(ks_tmp));

                    omega_inner[0] = dynamical->eval_phonon[ik][is];
                    omega_inner[1] = dynamical->eval_phonon[jk][js];

                    n1 = phonon_thermodynamics->fB(omega_inner[0], T) + phonon_thermodynamics->fB(omega_inner[1], T) + 1.0;
                    n2 = phonon_thermodynamics->fB(omega_inner[0], T) - phonon_thermodynamics->fB(omega_inner[1], T);

                    // e_tmp[0][kcount] = -omega_inner[0] - omega_inner[1];
                    e_tmp[1][kcount] = omega_inner[0] + omega_inner[1];
                    e_tmp[2][kcount] = omega_inner[0] - omega_inner[1];
                    e_tmp[3][kcount] = -omega_inner[0] + omega_inner[1];

                    // f_tmp[0][kcount] = -v3_tmp * n1;
                    f_tmp[1][kcount] = v3_tmp * n1;
                    f_tmp[2][kcount] = -v3_tmp * n2;
                    f_tmp[3][kcount] = v3_tmp * n2;

                    ++kcount;
                }
            }

            for (unsigned int j = 1; j < 4; ++j) {
                ret += integration->do_tetrahedron(e_tmp[j], f_tmp[j], omega);
            }
        }
    }
    return ret * pi * std::pow(0.5, 4);
}

void Relaxation::calc_selfenergy_at_T(const double T)
{
    unsigned int i;
    int ielem;

    for (i = 0; i < nks; ++i) {
        self_E[i] = std::complex<double>(0.0, 0.0);
    }

#pragma omp parallel
    {
        unsigned int ind[3], knum[3], snum[3];
        double omega[3];
        double n1, n2; 

        std::complex<double> ctmp;

#pragma omp for
        for (ielem = 0; ielem < V[0].size(); ++ielem){

            do {
                for (i = 0; i < 3; ++i){
                    ind[i] = V[0][ielem].ks[i];
                    knum[i] = ind[i] / ns;
                    snum[i] = ind[i] % ns;
                    omega[i] = dynamical->eval_phonon[knum[i]][snum[i]];
                }

                n1 = phonon_thermodynamics->fB(omega[1], T) + phonon_thermodynamics->fB(omega[2], T) + 1.0;
                n2 = phonon_thermodynamics->fB(omega[1], T) - phonon_thermodynamics->fB(omega[2], T);

                ctmp = std::norm(V[0][ielem].v)
                    * ( n1 / (omega[0] + omega[1] + omega[2] + im * epsilon)
                    - n1 / (omega[0] - omega[1] - omega[2] + im * epsilon) 
                    + n2 / (omega[0] - omega[1] + omega[2] + im * epsilon)
                    - n2 / (omega[0] + omega[1] - omega[2] + im * epsilon));

                self_E[ns * kpoint->knum_minus[knum[0]] + snum[0]] += ctmp;

            } while (std::next_permutation(V[0][ielem].ks.begin(), V[0][ielem].ks.end()));
        }
    }

    for (i = 0; i < nks; ++i){
        self_E[i] *= std::pow(0.5, 4) / static_cast<double>(nk);
    }
}

void Relaxation::modify_eigenvectors()
{
    bool *flag_done;
    unsigned int ik;
    unsigned int is, js;
    unsigned int nk_inv;
    std::complex<double> *evec_tmp;

    if (mympi->my_rank == 0) {
        std::cout << "**********      NOTICE      **********" << std::endl;
        std::cout << "For the brevity of the calculation, " << std::endl;
        std::cout << "phonon eigenvectors will be modified" << std::endl;
        std::cout << "so that e_{-ks}^{mu} = (e_{ks}^{mu})^{*}. " << std::endl;
    }

    memory->allocate(flag_done, nk);
    memory->allocate(evec_tmp, ns);



    for (ik = 0; ik < nk; ++ik) flag_done[ik] = false;

    for (ik = 0; ik < nk; ++ik){

        if (!flag_done[ik]) {

            nk_inv = kpoint->knum_minus[ik];   

            for (is = 0; is < ns; ++is){
                for (js = 0; js < ns; ++js){
                    evec_tmp[js] = dynamical->evec_phonon[ik][is][js];
                }

                for (js = 0; js < ns; ++js){
                    dynamical->evec_phonon[nk_inv][is][js] = std::conj(evec_tmp[js]);
                }
            }

            flag_done[ik] = true;
            flag_done[nk_inv] = true;
        }
    }

    memory->deallocate(flag_done);
    memory->deallocate(evec_tmp);

    MPI_Barrier(MPI_COMM_WORLD);
    if (mympi->my_rank == 0) {
        std::cout << "done !" << std::endl;
        std::cout << "**************************************" << std::endl;
    }
}

double Relaxation::delta_lorentz(const double omega)
{
    return epsilon / (omega*omega + epsilon*epsilon) / pi;
}

void Relaxation::calc_selfenergy()
{
    double Tmin = system->Tmin;
    double Tmax = system->Tmax;
    double dT = system->dT;

    unsigned int NT = static_cast<unsigned int>((Tmax - Tmin) / dT);
    unsigned int i, j;

    double *T_arr, T;

    std::string file_selfenergy;
    std::ofstream ofs_selfenergy;

    std::string file_test;
    std::ofstream ofs_test;

    unsigned int knum, snum;
    double omega;
    double *damping;

    memory->allocate(T_arr, NT);
    memory->allocate(damping, NT);

    for (i = 0; i < NT; ++i) T_arr[i] = Tmin + dT * static_cast<double>(i);

    /*
    file_test = input->job_title + ".test";
    ofs_test.open(file_test.c_str(), std::ios::out);
    if(!ofs_test) error->exit("write_selfenergy", "cannot open file_test");

    for (i = 0; i < dos->n_energy; ++i){
    omega_tmp = dos->emin + dos->delta_e * static_cast<double>(i);
    ofs_test << std::setw(15) << omega_tmp;
    omega_tmp *= time_ry / Hz_to_kayser;
    for (j = 3; j < 4; ++j){
    ctmp = selfenergy(T, omega_tmp, 0, j);
    ofs_test << std::setw(15) << ctmp.real()/time_ry*Hz_to_kayser;
    ofs_test << std::setw(15) << ctmp.imag()/time_ry*Hz_to_kayser;
    }
    ofs_test << std::endl;
    }

    ofs_test.close();
    */    
    if (mympi->my_rank == 0) {
        file_test = input->job_title + ".damp_T";
        ofs_test.open(file_test.c_str(), std::ios::out);
        if(!ofs_test) error->exit("write_selfenergy", "cannot open file_test");

        knum = 1;
        snum = 0;

        omega = dynamical->eval_phonon[knum][snum];

        ofs_test << "# Damping function [cm] of a phonon at xk = ";
        for (i = 0; i < 3; ++i) {
            ofs_test << std::setw(15) << kpoint->xk[knum][i];
        }
        ofs_test << std::endl;
        ofs_test << "# Branch = " << snum << std::endl;

        if (relaxation->ksum_mode == 0) {
            relaxation->calc_damping(NT, T_arr, omega, knum, snum, damping);
        } else if (relaxation->ksum_mode == -1) {
            relaxation->calc_damping_tetra(NT, T_arr, omega, knum, snum, damping);
        }

        for (i = 0; i < NT; ++i) {
            T = Tmin + dT * static_cast<double>(i);
            ofs_test  << std::setw(5) << T;
            ofs_test << std::setw(15) << 2.0 * damping[i]/time_ry*Hz_to_kayser;
            ofs_test << std::setw(15) << time_ry / (2.0 * damping[i]) * 1.0e+12;
            ofs_test << std::endl;
        }
        ofs_test.close();
    }
  //  error->exitall("hoge", "tomare!");
}
