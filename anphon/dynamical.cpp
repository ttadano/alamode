/*
dynamical.cpp

Copyright (c) 2014 Terumasa Tadano

This file is distributed under the terms of the MIT license.
Please see the file 'LICENCE.txt' in the root directory 
or http://opensource.org/licenses/mit-license.php for information.
*/

#include "mpi_common.h"
#include "dynamical.h"
#include "system.h"
#include "memory.h"
#include "kpoint.h"
#include "fcs_phonon.h"
#include <complex>
#include <vector>
#include "constants.h"
#include "fcs_phonon.h"
#include <iomanip>
#include <fstream>
#include "timer.h"
#include "error.h"
#include "symmetry_core.h"
#include "mathfunctions.h"
#include "write_phonons.h"
#include "phonon_dos.h"
#include "gruneisen.h"

using namespace PHON_NS;

Dynamical::Dynamical(PHON *phon): Pointers(phon){}
Dynamical::~Dynamical(){}

void Dynamical::setup_dynamical(std::string mode)
{
    int i;
    int ix, iy, iz;
    int icell = 0;

    neval = 3 * system->natmin;
    UPLO = 'U';

    if (mympi->my_rank == 0) {
        std::cout << std::endl << std::endl;
        std::cout << " ------------------------------------------------------------" << std::endl << std::endl;
        if (nonanalytic == 1) {
            std::cout << std::endl;
            std::cout << "  NONANALYTIC = 1 : Non-analytic part of the dynamical matrix will be included " << std::endl;
            std::cout << "                    by the Parlinski's method." << std::endl;
            std::cout << "                    The damping factor for the non-analytic term : " << na_sigma << std::endl;
            std::cout << std::endl;
        } else if (nonanalytic == 2) {
            std::cout << std::endl;
            std::cout << "  NONANALYTIC = 2 : Non-analytic part of the dynamical matrix will be included " << std::endl;
            std::cout << "                    by the mixed-space approach." << std::endl;
            std::cout << std::endl;
        }
    }

    memory->allocate(xshift_s, 27, 3);

    for (i = 0; i < 3; ++i) xshift_s[0][i] = 0.0;

    for (ix = -1; ix <= 1; ++ix) {
        for (iy = -1; iy <= 1; ++iy) {
            for (iz = -1; iz <= 1; ++iz) {
                if (ix == 0 && iy == 0 && iz == 0) continue;

                ++icell;

                xshift_s[icell][0] = static_cast<double>(ix);
                xshift_s[icell][1] = static_cast<double>(iy);
                xshift_s[icell][2] = static_cast<double>(iz);
            }
        }
    }

    if (mympi->my_rank == 0) {
        eigenvectors = false;

        if (phon->mode == "RTA") {
            eigenvectors = true;
        } else {
            if (print_eigenvectors || writes->print_msd || writes->writeanime 
                || dos->projected_dos || gruneisen->print_gruneisen) {
                eigenvectors = true;
            }
        }
    }


    MPI_Bcast(&eigenvectors, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD);
    MPI_Bcast(&nonanalytic, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD);

    if (nonanalytic) {

        memory->allocate(borncharge, system->natmin, 3, 3);

        if (mympi->my_rank == 0) load_born();

        MPI_Bcast(&dielec[0][0], 9, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&borncharge[0][0][0], 9*system->natmin, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&na_sigma, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        memory->allocate(mindist_list, system->natmin, system->nat);
        prepare_mindist_list(mindist_list);

//         for (i = 0; i < system->natmin; ++i) {
//             for (int j = 0; j < system->nat; ++j) {
//                 std::cout << std::setw(5) << i + 1;
//                 std::cout << std::setw(5) << j + 1 << std::endl;
// 
//                 for (int k = 0; k < mindist_list[i][j].size(); ++k) {
//                     std::cout << std::setw(5) << mindist_list[i][j][k] << std::endl;
//                 }
//             }
//         }
    }
}

void Dynamical::prepare_mindist_list(std::vector<int> **mindist_out)
{
    unsigned int i, j, k;
    unsigned int icell;
    unsigned int nneib = 27;

    double dist_tmp;
    double ***xcrd;
    double vec[3];

    unsigned int iat;
    unsigned int nat = system->nat;
    unsigned int natmin = system->natmin;
    int isize, jsize, ksize;

    std::vector<DistWithCell> **distall;

    memory->allocate(distall, natmin, nat);
    memory->allocate(xcrd, nneib, nat, 3);

    for (i = 0; i < nat; ++i){
        for (j = 0; j < 3; ++j){
            xcrd[0][i][j] = system->xr_s[i][j];
        }
    }
    icell = 0;
    for (isize = -1; isize <= 1; ++isize){
        for (jsize = -1; jsize <= 1; ++jsize){
            for (ksize = -1; ksize <= 1; ++ksize){

                if (isize == 0 && jsize == 0 && ksize == 0) continue;

                ++icell;
                for (i = 0; i < nat; ++i){
                    xcrd[icell][i][0] = system->xr_s[i][0] + static_cast<double>(isize);
                    xcrd[icell][i][1] = system->xr_s[i][1] + static_cast<double>(jsize);
                    xcrd[icell][i][2] = system->xr_s[i][2] + static_cast<double>(ksize);
                }
            }
        }
    }

    for (icell = 0; icell < nneib; ++icell) {
        for (i = 0; i < nat; ++i) {
            rotvec(xcrd[icell][i], xcrd[icell][i], system->lavec_s);
        }
    }

    for (i = 0; i < natmin; ++i) {
        iat = system->map_p2s[i][0];
        for (j = 0; j < nat; ++j) {
            distall[i][j].clear();
            for (icell = 0; icell < nneib; ++icell) {

                dist_tmp = distance(xcrd[0][iat], xcrd[icell][j]);
                distall[i][j].push_back(DistWithCell(icell, dist_tmp));
            }
            std::sort(distall[i][j].begin(), distall[i][j].end());
        }
    }



    // Construct pairs of minimum distance.

    double dist_min;
    for (i = 0; i < natmin; ++i) {
        for (j = 0; j < nat; ++j) {
            mindist_out[i][j].clear();

            dist_min = distall[i][j][0].dist;
            for (std::vector<DistWithCell>::const_iterator it = distall[i][j].begin(); it != distall[i][j].end(); ++it) {
                if (std::abs((*it).dist - dist_min) < eps8) {
                    mindist_out[i][j].push_back((*it).cell);
                }
            }
        }
    }

    memory->deallocate(distall);
    memory->deallocate(xcrd);
}

double Dynamical::distance(double *x1, double *x2)
{
    double dist;    
    dist = std::pow(x1[0] - x2[0], 2) + std::pow(x1[1] - x2[1], 2) + std::pow(x1[2] - x2[2], 2);
    dist = std::sqrt(dist);

    return dist;
}

void Dynamical::eval_k(double *xk_in, double *kvec_in, std::vector<FcsClassExtent> fc2_ext,
                       double *eval_out, std::complex<double> **evec_out, bool require_evec) {

    // Calculate phonon energy for the specific k-point given in fractional basis

    unsigned int i, j;
    std::complex<double> **dymat_k;

    memory->allocate(dymat_k, neval, neval);

    calc_analytic_k(xk_in, fc2_ext, dymat_k);

    if (nonanalytic) {

        std::complex<double> **dymat_na_k;

        memory->allocate(dymat_na_k, neval, neval);
        
        if (nonanalytic == 1) {
            calc_nonanalytic_k(xk_in, kvec_in, dymat_na_k);
        } else if (nonanalytic == 2) {
            calc_nonanalytic_k2(xk_in, kvec_in, fc2_ext, dymat_na_k);
        }        

        // Add non-analytic correction

        for (i = 0; i < neval; ++i) {
            for (j = 0; j < neval; ++j) {
                dymat_k[i][j] += dymat_na_k[i][j];
            }
        }
        memory->deallocate(dymat_na_k);
    }

    char JOBZ;
    int INFO, LWORK;
    double *RWORK;
    std::complex<double> *WORK;

    LWORK = (2 * neval - 1) * 10;
    memory->allocate(RWORK, 3*neval - 2);
    memory->allocate(WORK, LWORK);

    std::complex<double> *amat;
    memory->allocate(amat, neval * neval);

    unsigned int k = 0;
    int n = dynamical->neval;
    for(i = 0; i < neval; ++i){
        for (j = 0; j < neval; ++j){
            amat[k++] = dymat_k[i][j];
        }
    }

    memory->deallocate(dymat_k);

    if (require_evec) {
        JOBZ = 'V';
    } else {
        JOBZ = 'N';
    }

    // Perform diagonalization
    zheev_(&JOBZ, &UPLO, &n, amat, &n, eval_out, WORK, &LWORK, RWORK, &INFO);

    if (eigenvectors && require_evec){
        k = 0;
        for(i = 0; i < neval; ++i){
            for (j = 0; j < neval; ++j){
                evec_out[i][j] = amat[k++];
            }
        }
    }

    memory->deallocate(RWORK);
    memory->deallocate(WORK);
    memory->deallocate(amat);
}

void Dynamical::calc_analytic_k(double *xk_in, std::vector<FcsClassExtent> fc2_in, std::complex<double> **dymat_out)
{
    int i, j;
    unsigned int atm1_s, atm2_s;
    unsigned int atm1_p, atm2_p;
    unsigned int xyz1, xyz2;
    unsigned int icell;

    double vec[3];
    double phase;
    std::complex<double> im(0.0, 1.0);
    std::complex<double> **ctmp;

    memory->allocate(ctmp, 3*system->natmin, 3*system->natmin);

    for (i = 0; i < 3*system->natmin; ++i) {
        for (j = 0; j < 3*system->natmin; ++j) {
            dymat_out[i][j] = std::complex<double>(0.0, 0.0);
        }
    }

    for (std::vector<FcsClassExtent>::const_iterator it = fc2_in.begin(); it != fc2_in.end(); ++it) {

        atm1_p = (*it).atm1;
        atm2_s = (*it).atm2;
        xyz1 = (*it).xyz1;
        xyz2 = (*it).xyz2;
        icell = (*it).cell_s;

        atm1_s = system->map_p2s[atm1_p][0];
        atm2_p = system->map_s2p[atm2_s].atom_num;

        for (i = 0; i < 3; ++i) {
            vec[i] = system->xr_s[atm2_s][i] + xshift_s[icell][i]
                   - system->xr_s[system->map_p2s[atm2_p][0]][i] ;
        }

        rotvec(vec, vec, system->lavec_s);
        rotvec(vec, vec, system->rlavec_p);

        phase = vec[0] * xk_in[0] + vec[1] * xk_in[1] + vec[2] * xk_in[2];

        dymat_out[3 * atm1_p + xyz1][3 * atm2_p + xyz2] 
        += (*it).fcs_val * std::exp(im * phase) / std::sqrt(system->mass[atm1_s] * system->mass[atm2_s]);
    }
}

void Dynamical::calc_nonanalytic_k(double *xk_in, double *kvec_na_in, std::complex<double> **dymat_na_out)
{
    // Calculate the non-analytic part of dynamical matrices 
    // by Parlinski's method.

    unsigned int i, j;
    unsigned int iat, jat;
    unsigned int atm_p1, atm_p2;
    unsigned int natmin = system->natmin;
    double kepsilon[3];
    double kz1[3], kz2[3];
    double denom, norm2;
    double born_tmp[3][3];
    double xk_tmp[3], xdiff[3];
    double factor, phase;
    std::complex<double> im(0.0, 1.0);


    for (i = 0; i < neval; ++i) {
        for (j = 0; j < neval; ++j) {
            dymat_na_out[i][j] = std::complex<double>(0.0, 0.0);
        }
    }

    rotvec(kepsilon, kvec_na_in, dielec);
    denom = kvec_na_in[0] * kepsilon[0] + kvec_na_in[1] * kepsilon[1] + kvec_na_in[2] * kepsilon[2];

    if (denom > eps) {

        for (iat = 0; iat < natmin; ++iat) {
            atm_p1 = system->map_p2s[iat][0];

            for (i = 0; i <3; ++i) {
                for (j = 0; j < 3; ++j) {
                    born_tmp[i][j] = borncharge[iat][i][j];
                }
            }

            rotvec(kz1, kvec_na_in, born_tmp, 'T');

            for (jat = 0; jat < natmin; ++jat) {
                atm_p2 = system->map_p2s[jat][0];


                for (i = 0; i <3; ++i) {
                    for (j = 0; j < 3; ++j) {
                        born_tmp[i][j] = borncharge[jat][i][j];
                    }
                }

                rotvec(kz2, kvec_na_in, born_tmp, 'T');

                for (i = 0; i < 3; ++i) {
                    for (j = 0; j < 3; ++j) {

                        dymat_na_out[3 * iat + i][3 * jat + j] 
                        = kz1[i] * kz2[j] / (denom * std::sqrt(system->mass[atm_p1] * system->mass[atm_p2]));

                    }
                }
            }
        }
    }

    rotvec(xk_tmp, xk_in, system->rlavec_p, 'T');
    norm2 = xk_tmp[0] * xk_tmp[0] + xk_tmp[1] * xk_tmp[1] + xk_tmp[2] * xk_tmp[2];

    factor = 8.0 * pi / system->volume_p * std::exp(-norm2 / std::pow(na_sigma, 2));

    for (i = 0; i < neval; ++i) {
        for (j = 0; j < neval; ++j) {
            dymat_na_out[i][j] *= factor;
        }
    }

    // Multiply an additional phase factor for the non-analytic term.

    for (iat = 0; iat < natmin; ++iat) {
        for (jat = 0; jat < natmin; ++jat) {

            for (i = 0; i < 3; ++i) {
                xdiff[i] = system->xr_s[system->map_p2s[iat][0]][i]
                         - system->xr_s[system->map_p2s[jat][0]][i];
            }

            rotvec(xdiff, xdiff, system->lavec_s);
            rotvec(xdiff, xdiff, system->rlavec_p);

            phase = xk_in[0] * xdiff[0] + xk_in[1] * xdiff[1] + xk_in[2] * xdiff[2];

            for (i = 0; i < 3; ++i) {
                for (j = 0; j < 3; ++j) {
                    dymat_na_out[3 * iat + i][3 * jat + j] *= exp(im * phase);
                }
            }
        }
    }
}

void Dynamical::calc_nonanalytic_k2(double *xk_in, double *kvec_na_in, 
                                    std::vector<FcsClassExtent> fc2_in, std::complex<double> **dymat_na_out)
{
    // Calculate the non-analytic part of dynamical matrices 
    // by the mixed-space approach.

    unsigned int i, j, k;
    unsigned int iat, jat;
    unsigned int atm_p1, atm_p2, atm_s2;
    unsigned int natmin = system->natmin;
    unsigned int itran, icell, cell;
    double kepsilon[3];
    double kz1[3], kz2[3];
    double denom, norm2;
    double born_tmp[3][3];
    double xk_tmp[3], vec[3];
    double factor, phase;
    std::complex<double> im(0.0, 1.0);
    std::complex<double> exp_phase, exp_phase_tmp;


    for (i = 0; i < neval; ++i) {
        for (j = 0; j < neval; ++j) {
            dymat_na_out[i][j] = std::complex<double>(0.0, 0.0);
        }
    }

    rotvec(kepsilon, kvec_na_in, dielec);
    denom = kvec_na_in[0] * kepsilon[0] + kvec_na_in[1] * kepsilon[1] + kvec_na_in[2] * kepsilon[2];

    if (denom > eps) {

        for (iat = 0; iat < natmin; ++iat) {
            atm_p1 = system->map_p2s[iat][0];

            for (i = 0; i <3; ++i) {
                for (j = 0; j < 3; ++j) {
                    born_tmp[i][j] = borncharge[iat][i][j];
                }
            }

            rotvec(kz1, kvec_na_in, born_tmp, 'T');

            for (jat = 0; jat < natmin; ++jat) {
                atm_p2 = system->map_p2s[jat][0];


                for (i = 0; i <3; ++i) {
                    for (j = 0; j < 3; ++j) {
                        born_tmp[i][j] = borncharge[jat][i][j];
                    }
                }

                rotvec(kz2, kvec_na_in, born_tmp, 'T');

                exp_phase = std::complex<double>(0.0, 0.0);

                for (i = 0; i < system->ntran; ++i) {

                    exp_phase_tmp = std::complex<double>(0.0, 0.0);
                    atm_s2 = system->map_p2s[jat][i];

                    // Average over mirror atoms

                    for (j = 0; j < mindist_list[iat][atm_s2].size(); ++j) {
                        cell = mindist_list[iat][atm_s2][j];

                        for (k = 0; k < 3; ++k) {
                            vec[k] = system->xr_s[system->map_p2s[jat][i]][k] + xshift_s[cell][k]
                                   - system->xr_s[atm_p2][k];
                        }

                        rotvec(vec, vec, system->lavec_s);
                        rotvec(vec, vec, system->rlavec_p);

                        phase = vec[0] * xk_in[0] + vec[1] * xk_in[1] + vec[2] * xk_in[2];

                        exp_phase_tmp += std::exp(im * phase);
                    }
                    exp_phase += exp_phase_tmp / static_cast<double>(mindist_list[iat][atm_s2].size());
                }
                exp_phase /= static_cast<double>(system->ntran);

              //  std::cout << std::setw(15) << exp_phase.real() << std::setw(15) << exp_phase.imag() << std::endl; 


                for (i = 0; i < 3; ++i) {
                    for (j = 0; j < 3; ++j) {
                        dymat_na_out[3 * iat + i][3 * jat + j] 
                        = kz1[i] * kz2[j] / (denom * std::sqrt(system->mass[atm_p1] * system->mass[atm_p2])) * exp_phase;
                    }
                }
            }
        }
    }
    
    factor = 8.0 * pi / system->volume_p;
    
    for (i = 0; i < neval; ++i) {
        for (j = 0; j < neval; ++j) {
            dymat_na_out[i][j] *= factor;
         }
     }
}


void Dynamical::diagonalize_dynamical_all()
{
    int ik;
    unsigned int i;
    unsigned int is;
    unsigned int nk = kpoint->nk;
    bool require_evec; 

    if (mympi->my_rank == 0) {
        std::cout << std::endl << " Diagonalizing dynamical matrices for all k points ... ";
    }

    memory->allocate(eval_phonon, nk, neval);
    if (eigenvectors) {
        require_evec = true;
        memory->allocate(evec_phonon, nk, neval, neval);
    } else {
        require_evec = false;
        memory->allocate(evec_phonon, nk, 1, 1);
    }

    // Calculate phonon eigenvalues and eigenvectors for all k-points

#pragma omp parallel for private (is)
    for (ik = 0; ik < nk; ++ik){

        eval_k(kpoint->xk[ik], kpoint->kvec_na[ik], fcs_phonon->fc2_ext, eval_phonon[ik], evec_phonon[ik], require_evec);

        // Phonon energy is the square-root of the eigenvalue 
        for (is = 0; is < neval; ++is){
            eval_phonon[ik][is] = freq(eval_phonon[ik][is]);
        }
    }

    if (mympi->my_rank == 0) {
        std::cout << "done !" << std::endl;
    }

    if (kpoint->kpoint_mode == 2 && eigenvectors) {
        //		modify_eigenvectors();
        //		modify_eigenvectors_sym();
    }

    // #ifdef _DEBUG
    // 
    // 
    //     // Check if D(k) - D(-k)^{*} is satisfied when KPMODE = 2.
    // 
    //     if (kpoint->kpoint_mode == 2) {
    //         int i, j;
    //         int nk_minus;
    //         std::complex<double> **dmat1, **dmat2;
    // 
    //         double res;
    //         // 	double *eval1, *eval2;
    //         // 		std::complex<double> **evec1, **evec2;
    //         double *xk1, *xk2;
    // 
    //         memory->allocate(dmat1, neval, neval);
    //         memory->allocate(dmat2, neval, neval);
    // 
    //         memory->allocate(xk1, 3);
    //         memory->allocate(xk2, 3);
    // 
    // 
    //         for (ik = 0; ik < nk; ++ik) {
    // 
    //             for (i = 0; i < 3; ++i) {
    //                 xk1[i] = kpoint->xk[ik][i];
    //                 xk2[i] = -kpoint->xk[ik][i];
    //             }
    // 
    //             nk_minus = kpoint->get_knum(xk2[0], xk2[1], xk2[2]);
    // 
    //             for (i = 0; i < 3; ++i) {
    //                 xk2[i] = kpoint->xk[nk_minus][i];
    //             }
    // 
    //             // 		std::cout << "ik = " << ik << std::endl;
    //             // 		std::cout << "xk = " << xk[0] << " " << xk[1] << " " << xk[2] << std::endl;
    // 
    //             calc_analytic_k(xk1, fcs_phonon->fc2, dmat1);	  
    //             calc_analytic_k(xk2, fcs_phonon->fc2, dmat2);
    // 
    //             res = 0.0;
    // 
    //             for (i = 0; i < neval; ++i) {
    //                 for (j = 0; j < neval; ++j) {
    // 
    //                     res += std::norm(dmat1[i][j] - std::conj(dmat2[i][j]));
    //                 }
    //             }
    // 
    //             if (std::sqrt(res)/static_cast<double>(neval) > eps12) {
    //                 std::cout << "ik = " << ik << std::endl;
    //                 error->warn("diagonalize_dynamical_all", "D(k) = D(-k)^{*} is not satisfied. This might imply a bug.");
    //             }
    // 
    // 
    //             // 	  std::cout << "D(k) - D(-k)^{*}" << std::endl;
    //             // 	  for (int icrd = 0; icrd < neval; ++icrd){
    //             // 		  for (int jcrd = 0; jcrd < neval; ++jcrd){
    //             // 			  std::cout << "(" << std::setw(15) << std::scientific << real(dmat1[icrd][jcrd] - std::conj(dmat2[icrd][jcrd]));
    //             // 			  std::cout << std::setw(15) << std::scientific << imag(dmat1[icrd][jcrd] - std::conj(dmat2[icrd][jcrd])) << ") ";
    //             // 		  }
    //             // 		  std::cout << std::endl;
    //             // 	  }
    //             // 	  std::cout << "Compare eigenvectors" << std::endl;
    //             // 	  for (int icrd = 0; icrd < neval; ++icrd) {
    //             // 		 for (int jcrd = 0; jcrd < neval; ++jcrd) {
    //             // 			 std::cout << evec1[icrd][jcrd] << " " << evec2[icrd][jcrd] << " " << evec1[icrd][jcrd] - std::conj(evec2[icrd][jcrd]) << std::endl;
    //             // 		 }
    //             // 		 std::cout << std::endl;
    //             // 	}
    //             // 	}
    // 
    //         }
    //         memory->deallocate(dmat1);
    //         memory->deallocate(dmat2);
    //         memory->deallocate(xk1);
    //         memory->deallocate(xk2);
    // 
    //     }
    // 
    //     // 	 	 	std::complex<double> prod;
    //     // 	 	 	std::cout << "orthogonality check 1" << std::endl;
    //     // 	 	 	for (ik = 0; ik < nk; ++ik) {
    //     // 	 	 		
    //     // 	 	 
    //     // 	 	 		for (int is1 = 0; is1 < neval; ++is1) {
    //     // 	 	 			for (int is2 = 0; is2 < neval; ++is2) {
    //     // 	 	 
    //     // 	 	 				prod = std::complex<double>(0.0, 0.0);
    //     // 	 	 
    //     // 	 	 				for (int i = 0; i < neval; ++i) {
    //     // 	 	 					prod += std::conj(evec_phonon[ik][is1][i]) * evec_phonon[ik][is2][i];
    //     // 	 	 				}
    //     // 	 	 				if (std::norm(prod) > eps) {
    //     // 	 	 					std::cout << "ik = " << ik << " is1, is2 =" << is1 << " " << is2 << " prod = " << prod << std::endl;
    //     // 	 	 				}
    //     // 	 	 			}
    //     // 	 	 		}
    //     // 	 	 	}
    //     // 	 	 	std::cout << "orthogonality check 2" << std::endl;
    //     // 	 	 
    //     // 	 	 	for (ik = 0; ik < nk; ++ik) {
    //     // 	 	 				for (int i = 0; i < neval; ++i) {
    //     // 	 	 					for (int j = 0; j < neval; ++j) {
    //     // 	 	 
    //     // 	 	 						prod = std::complex<double>(0.0, 0.0);
    //     // 	 	 
    //     // 	 	 						for (int is1 = 0; is1 < neval; ++is1) {
    //     // 	 	 							prod += std::conj(evec_phonon[ik][is1][i]) * evec_phonon[ik][is1][j];
    //     // 	 	 						}
    //     // 	 	 						if (std::norm(prod) > eps) {
    //     // 	 	 							std::cout << "ik = " << ik << " i, j =" << i << " " << j << " prod = " << prod << std::endl;
    //     // 	 	 						}
    //     // 	 	 
    //     // 	 	 					}
    //     // 	 	 		}
    //     // 	 		}
    // 
    //     //  	double xk_tmp[3];
    //     //  	std::complex<double> **dymat;
    //     //  
    //     //  	memory->allocate(dymat, neval, neval);
    //     //  
    //     //  	for (ik = 0; ik < kpoint->nk_reduced; ++ik) {
    //     //  		for (unsigned int i = 0; i < kpoint->nk_equiv[ik]; ++i) {
    //     //  			std::cout << "ik = " << ik + 1 << " i = " << i + 1 << " kp = " << kpoint->k_reduced[ik][i];
    //     //  	
    //     //  			for (int j = 0; j < 3; ++j) {
    //     //  				xk_tmp[j] = kpoint->xk[kpoint->k_reduced[ik][i]][j];
    //     //  				std::cout << std::setw(15) << xk_tmp[j];
    //     //  			}
    //     //  			std::cout << std::endl;
    //     //  
    //     //  		  if (fcs_phonon->is_fc2_ext) {
    //     //  			  calc_analytic_k(xk_tmp, fcs_phonon->fc2_ext, dymat);
    //     //  		  } else {
    //     //  			  calc_analytic_k(xk_tmp, fcs_phonon->fc2, dymat);
    //     //  		  }
    //     //  		
    //     //  		  for (int j = 0; j < neval; ++j) {
    //     //  			  for (int k = 0; k < neval; ++k) {
    //     //  				  if (std::abs(dymat[j][k].real()) < eps12) {
    //     //  					  dymat[j][k] = std::complex<double>(0.0, dymat[j][k].imag());
    //     //  				  }
    //     //  				  if (std::abs(dymat[j][k].imag()) < eps12) {
    //     //  					  dymat[j][k] = std::complex<double>(dymat[j][k].real(), 0.0);
    //     //  				  }
    //     //  
    //     //  				  std::cout << std::setw(15) << dymat[j][k].real() << std::setw(15) << dymat[j][k].imag();
    //     //  		  }
    //     //  			  std::cout << std::endl;
    //     //  		}
    //     //  		  std::cout << std::endl;
    //     //  	}
    //     //  	}
    // #endif
}

void Dynamical::modify_eigenvectors()
{
    bool *flag_done;
    unsigned int ik;
    unsigned int is, js;
    unsigned int nk_inv;
    std::complex<double> *evec_tmp;

    unsigned int nk = kpoint->nk;
    unsigned int ns = neval;

    if (mympi->my_rank == 0) {
        std::cout << " **********      NOTICE      ********** " << std::endl;
        std::cout << " For the brevity of the calculation, " << std::endl;
        std::cout << " phonon eigenvectors will be modified" << std::endl;
        std::cout << " so that e_{-ks}^{mu} = (e_{ks}^{mu})^{*}. " << std::endl;
    }

    memory->allocate(flag_done, nk);
    memory->allocate(evec_tmp, ns);

    for (ik = 0; ik < nk; ++ik) flag_done[ik] = false;

    for (ik = 0; ik < nk; ++ik){

        if (!flag_done[ik]) {

            nk_inv = kpoint->knum_minus[ik];   

            for (is = 0; is < ns; ++is){
                for (js = 0; js < ns; ++js){
                    evec_tmp[js] = evec_phonon[ik][is][js];
                    //	evec_tmp[js] = 0.5 * (std::conj(evec_phonon[ik][is][js]) + evec_phonon[nk_inv][is][js]);
                }

                for (js = 0; js < ns; ++js){
                    evec_phonon[nk_inv][is][js] = std::conj(evec_tmp[js]);
                    //		evec_phonon[ik][is][js] = evec_tmp[js];
                    //		evec_phonon[nk_inv][is][js] = std::conj(evec_tmp[js]);
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
        std::cout << " done !" << std::endl;
        std::cout << " **************************************" << std::endl;
    }
}


void Dynamical::load_born()
{
    // Read the dielectric tensor and born effective charges from file_born

    unsigned int i, j, k;
    double sum_born[3][3];
    double res;
    std::ifstream ifs_born;

    ifs_born.open(file_born.c_str(), std::ios::in);
    if (!ifs_born) error->exit("load_born", "cannot open file_born");

    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            ifs_born >> dielec[i][j];
        }
    }

    for (i = 0; i < system->natmin; ++i) {
        for (j = 0; j < 3; ++j) {
            for (k = 0; k < 3; ++k) {
                ifs_born >> borncharge[i][j][k];
            }
        }
    }
    std::cout << "  Dielectric constants and Born effective charges are read from " 
        << file_born << "." << std::endl << std::endl;
    std::cout << "  Dielectric constant tensor in Cartesian coordinate : " << std::endl;
    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            std::cout << std::setw(15) << dielec[i][j];
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
    std::cout << "  Born effective charge tensor in Cartesian coordinate" << std::endl;
    for (i = 0; i < system->natmin; ++i) {
        std::cout << "  Atom" << std::setw(5) << i + 1 << "(" 
            << std::setw(3) << system->symbol_kd[system->kd[system->map_p2s[i][0]]] << ") :" << std::endl;

        for (j = 0; j < 3; ++j) {
            for (k = 0; k < 3; ++k) {
                std::cout << std::setw(15) << borncharge[i][j][k];
            }
            std::cout << std::endl;
        }
    }

    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            sum_born[i][j] = 0.0;
        }
    }

    for (i = 0; i < system->natmin; ++i) {
        for (j = 0; j < 3; ++j) {
            for (k = 0; k < 3; ++k) {
                sum_born[j][k] += borncharge[i][j][k];
            }
        }
    }

    res = 0.0;
    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            res += std::pow(sum_born[i][j], 2);
        }
    }

    if (res > eps10) {
        std::cout << std::endl;
        std::cout << "  WARNING: Born effective charges do not satisfy the acoustic sum rule." << std::endl;
        std::cout << "           The born effective charges will be modified as follows." << std::endl;

        for (i = 0; i < system->natmin; ++i) {
            for (j = 0; j < 3; ++j) {
                for (k = 0; k < 3; ++k) {
                    borncharge[i][j][k] -= sum_born[j][k] / static_cast<double>(system->natmin);
                }
            }
        }
        std::cout << std::endl;
        std::cout << "  New Born effective charge tensor in Cartesian coordinate." << std::endl;
        for (i = 0; i < system->natmin; ++i) {
            std::cout << "  Atom" << std::setw(5) << i + 1 << "(" 
                << std::setw(3) << system->symbol_kd[system->kd[system->map_p2s[i][0]]] << ") :" << std::endl;

            for (j = 0; j < 3; ++j) {
                for (k = 0; k < 3; ++k) {
                    std::cout << std::setw(15) << borncharge[i][j][k];
                }
                std::cout << std::endl;
            }
        }
    }
}

double Dynamical::fold(double x)
{
    return x - static_cast<double>(nint(x));
}

double Dynamical::freq(const double x) 
{
    if (std::abs(x) < eps) {
        // Special treatment to avoid the divergence of computation.
        return eps15;
    } else {
        if (x > 0.0) {
            return std::sqrt(x);
        } else {
            return -std::sqrt(-x);
        }
    }
}


void PHON_NS::Dynamical::finish_dynamical()
{
    memory->deallocate(xshift_s);

    if (kpoint->kpoint_mode < 3) {
        memory->deallocate(evec_phonon);
        memory->deallocate(eval_phonon);
    }

    if (nonanalytic) {
        memory->deallocate(borncharge);
    }
}


// void Dynamical::eval_k(double *xk_in, double *kvec_in, double ****fc2_in, 
//                        double *eval_out, std::complex<double> **evec_out, bool require_evec) {
// 
//                            // Calculate phonon energy for the specific k-point given in fractional basis
// 
//                            unsigned int i, j;
// 
//                            std::complex<double> **dymat_k;
//                            double **dymat_na_k;
//                            std::complex<double> **dymat_na_mod;
// 
//                            memory->allocate(dymat_k, neval, neval);
// 
//                            calc_analytic_k(xk_in, fc2_in, dymat_k);
// 
//                            if (nonanalytic) {
// 
//                                memory->allocate(dymat_na_k, neval, neval);
//                                memory->allocate(dymat_na_mod, neval, neval);
// 
//                                calc_nonanalytic_k(xk_in, kvec_in, dymat_na_k);
// 
//                                double xdiff[3];
//                                double phase;
//                                std::complex<double> im(0.0, 1.0);
//                                unsigned int icrd, jcrd;
// 
//                                // Multiply a phase factor for the non-analytic term.
//                                for (i = 0; i < system->natmin; ++i) {
//                                    for (j = 0; j < system->natmin; ++j) {
// 
//                                        for (icrd = 0; icrd < 3; ++icrd) {
//                                            xdiff[icrd] = system->xr_s[system->map_p2s[i][0]][icrd] 
//                                            - system->xr_s[system->map_p2s[j][0]][icrd];
//                                        }
// 
//                                        rotvec(xdiff, xdiff, system->lavec_s);
//                                        rotvec(xdiff, xdiff, system->rlavec_p);
// 
//                                        phase = xk_in[0] * xdiff[0] + xk_in[1] * xdiff[1] + xk_in[2] * xdiff[2];
// 
//                                        for (icrd = 0; icrd < 3; ++icrd) {
//                                            for (jcrd = 0; jcrd < 3; ++jcrd) {
//                                                dymat_na_mod[3 * i + icrd][3 * j + jcrd] = dymat_na_k[3 * i + icrd][3 * j + jcrd] 
//                                                * exp(im * phase);
//                                            }
//                                        }
//                                    }
//                                }
// 
//                                for (i = 0; i < neval; ++i) {
//                                    for (j = 0; j < neval; ++j) {
//                                        dymat_k[i][j] += dymat_na_mod[i][j];
//                                    }
//                                }
//                                memory->deallocate(dymat_na_k);
//                                memory->deallocate(dymat_na_mod);
//                            }
// 
// 
//                            // Hermitize the dynamical matrix
// 
//                            std::complex<double> **dymat_tmp, **dymat_transpose;
// 
//                            memory->allocate(dymat_tmp, neval, neval);
//                            memory->allocate(dymat_transpose, neval, neval);
// 
//                            for (i = 0; i < neval; ++i){
//                                for (j = 0; j < neval; ++j){
//                                    dymat_tmp[i][j] = dymat_k[i][j];
//                                    dymat_transpose[i][j] = dymat_k[j][i];
//                                }
//                            }
//                            for (i = 0; i < neval; ++i){
//                                for (j = 0; j < neval; ++j){
//                                    dymat_k[i][j] = 0.5 * (dymat_tmp[i][j] + std::conj(dymat_transpose[i][j]));
//                                }
//                            }
// 
//                            memory->deallocate(dymat_tmp);
//                            memory->deallocate(dymat_transpose);
// 
// 
//                            char JOBZ;
//                            int INFO, LWORK;
//                            double *RWORK;
//                            std::complex<double> *WORK;
// 
//                            LWORK = (2 * neval - 1) * 10;
//                            memory->allocate(RWORK, 3*neval - 2);
//                            memory->allocate(WORK, LWORK);
// 
//                            std::complex<double> *amat;
//                            memory->allocate(amat, neval * neval);
// 
//                            unsigned int k = 0;
//                            int n = dynamical->neval;
//                            for(i = 0; i < neval; ++i){
//                                for (j = 0; j < neval; ++j){
//                                    amat[k++] = dymat_k[i][j];
//                                }
//                            }
// 
//                            memory->deallocate(dymat_k);
// 
//                            if (require_evec) {
//                                JOBZ = 'V';
//                            } else {
//                                JOBZ = 'N';
//                            }
// 
//                            // Perform diagonalization
//                            zheev_(&JOBZ, &UPLO, &n, amat, &n, eval_out, WORK, &LWORK, RWORK, &INFO);
// 
//                            if (eigenvectors && require_evec){
//                                k = 0;
//                                for(i = 0; i < neval; ++i){
//                                    for (j = 0; j < neval; ++j){
//                                        evec_out[i][j] = amat[k++];
//                                    }
//                                }
//                            }
// 
//                            memory->deallocate(RWORK);
//                            memory->deallocate(WORK);
//                            memory->deallocate(amat);
// }

// void Dynamical::calc_analytic_k(double *xk_in, double ****fc2_in, std::complex<double> **dymat_out) {
// 
//     unsigned int i, j;
//     unsigned int icrd, jcrd;
//     unsigned int itran;
//     unsigned int ntran = system->ntran;
//     unsigned int natmin = system->natmin;
//     unsigned int atm_p1, atm_p2, atm_s2;
// 
//     double phase;
//     double vec[3];
//     std::complex<double> ctmp[3][3];
//     std::complex<double> im(0.0, 1.0);
//     std::complex<double> exp_phase;
// 
// 
//     for (i = 0; i < natmin; ++i){
// 
//         atm_p1 = system->map_p2s[i][0];
// 
//         for (j = 0; j < natmin; ++j){
// 
//             atm_p2 = system->map_p2s[j][0];
// 
//             for(icrd = 0; icrd < 3; ++icrd){
//                 for(jcrd = 0; jcrd < 3; ++jcrd){
//                     ctmp[icrd][jcrd] = std::complex<double>(0.0, 0.0);
//                 }
//             }
// 
//             for(itran = 0; itran < ntran; ++itran){
//                 atm_s2 =system->map_p2s[j][itran];
// 
//                 for(icrd = 0; icrd < 3; ++icrd){
//                     if (system->cell_dimension[icrd] == 1) {
//                         vec[icrd] = system->xr_s[atm_p1][icrd] - system->xr_s[atm_s2][icrd];
//                         if (std::abs(vec[icrd]) < 0.5) {
//                             vec[icrd] = 0.0;
//                         } else {
//                             if (system->xr_s[atm_p1][icrd] < 0.5) {
//                                 vec[icrd] = 1.0;
//                             } else {
//                                 vec[icrd] = -1.0;
//                             }
//                         }
//                     } else if (system->cell_dimension[icrd] == 2) {
//                         vec[icrd] = system->xr_s[atm_p2][icrd] - system->xr_s[atm_s2][icrd];
//                         vec[icrd] = fold(vec[icrd]);
//                         if (std::abs(system->xr_s[atm_p1][icrd] - system->xr_s[atm_s2][icrd]) > 0.5) vec[icrd] *= -1.0;
//                     } else {
//                         vec[icrd] = system->xr_s[atm_p1][icrd] - system->xr_s[atm_s2][icrd];
//                         vec[icrd] = fold(vec[icrd]);
// 
//                         vec[icrd] += system->xr_s[atm_p2][icrd] - system->xr_s[atm_p1][icrd];
//                     }
//                 }
// 
//                 rotvec(vec, vec, system->lavec_s);
//                 rotvec(vec, vec, system->rlavec_p);
// 
//                 phase = vec[0] * xk_in[0] + vec[1] * xk_in[1] + vec[2] * xk_in[2];
//                 exp_phase = std::exp(-im * phase);
// 
// 
// 
//                 for (icrd = 0; icrd < 3; ++icrd){
//                     for (jcrd = 0; jcrd < 3; ++jcrd){
//                         ctmp[icrd][jcrd] += fc2_in[i][atm_s2][icrd][jcrd] * exp_phase;
//                     }
//                 }
//             }
// 
//             for (icrd = 0; icrd < 3; ++icrd){
//                 for (jcrd = 0; jcrd < 3; ++jcrd){
//                     dymat_out[3 * i + icrd][3 * j + jcrd] = ctmp[icrd][jcrd] / std::sqrt(system->mass[atm_p1] * system->mass[atm_p2]);
//                 }
//             }
//         }
//     }
// }

// void Dynamical::modify_eigenvectors_sym()
// {
//     // This is under test.
// 
//     double Sk[3], S[3][3], k[3];
//     double T[3][3], S_recip[3][3];
//     double x[3], x_new[3], x_tmp[3];
//     double AA[3][3], BB[3][3];
//     double shift[3];
//     double mat_tmp[3][3];
// 
//     std::complex<double> **gamma, *gamma_1d;
//     std::complex<double> *evec_in, *evec_out;
// 
//     unsigned int natmin = system->natmin;
//     unsigned int ik;
//     unsigned int nk = kpoint->nk;
//     int ns = neval;
//     unsigned int i, j;
//     bool *flag_done;
//     int knum_sym;
//     unsigned int iat, jat, icrd, jcrd;
//     unsigned int *atom_mapped;
//     double phase;
//     double det;
// 
//     std::complex<double> exp_phase;
//     std::complex<double> im = std::complex<double>(0.0, 1.0);
//     std::complex<double> S_evec;
//     std::complex<double> alpha = std::complex<double>(1.0, 0.0);
//     std::complex<double> beta = std::complex<double>(0.0, 0.0);
// 
//     std::complex<double> **dmat, *dmat_1d;
//     std::complex<double> **dmat_sym, *dmat_sym_1d;
//     std::complex<double> *gamma_dagger;
// 
//     char TRANSA, TRANSB;
//     TRANSA = 'N';
//     TRANSB = 'N';
// 
//     memory->allocate(flag_done, nk);
//     for (ik = 0; ik < nk; ++ik) flag_done[ik] = false;
// 
//     memory->allocate(atom_mapped, system->natmin);
// 
//     for (i = 0; i < 3; ++i) {
//         for (j = 0; j < 3; ++j) {
//             AA[i][j] = system->lavec_p[i][j];
//             BB[i][j] = system->rlavec_p[i][j];
//         }
//     }
// 
//     memory->allocate(gamma, ns, ns);
//     memory->allocate(gamma_1d, ns * ns);
//     memory->allocate(evec_in, ns * ns);
//     memory->allocate(evec_out, ns * ns);
//     memory->allocate(gamma_dagger, ns * ns);
// 
//     memory->allocate(dmat, ns, ns);
//     memory->allocate(dmat_sym, ns, ns);
//     memory->allocate(dmat_1d, ns*ns);
//     memory->allocate(dmat_sym_1d, ns*ns);
// 
//     int isym = 0;
// 
//     for (ik = 0; ik < nk; ++ik) {
// 
//         for (i = 0; i < 3; ++i) k[i] = kpoint->xk[ik][i];
// 
//         isym = 0;
// 
//         for (std::vector<SymmetryOperationWithMapping>::const_iterator it = symmetry->SymmListWithMap.begin(); it != symmetry->SymmListWithMap.end(); ++it) {
// 
//             ++isym;
// 
//             // Initialize the Gamma matrix
//             for (i = 0; i < ns; ++i) {
//                 for (j = 0; j < ns; ++j) {
//                     gamma[i][j] = 0.0;
//                 }
//             }
// 
//             for (i = 0; i < 3; ++i) {
//                 for (j = 0; j < 3; ++j) {
//                     S[i][j] = (*it).rot[3 * i + j];
//                     T[i][j] = (*it).rot_real[3 * i + j];
//                     S_recip[i][j] = (*it).rot_reciprocal[3 * i + j];
//                 }
//                 shift[i] = (*it).shift[i];
//             }
//             det = S[0][0] * (S[1][1] * S[2][2] - S[2][1] * S[1][2])
//                 - S[1][0] * (S[0][1] * S[2][2] - S[2][1] * S[0][2])
//                 + S[2][0] * (S[0][1] * S[1][2] - S[1][1] * S[0][2]);
// 
//             if (std::abs(det + 1.0) < eps12) continue;
// 
//             if (ik == 0) {
//                 std::cout << "#Sym : " << std::setw(5) << isym << std::endl;
//                 std::cout << "S" << std::endl;
//                 for (i = 0; i < 3; ++i) {
//                     for (j = 0; j < 3; ++j) {
//                         std::cout << std::setw(5) << S[i][j];
//                     }
//                     std::cout << std::endl;
//                 }
//                 std::cout << "T" << std::endl;
//                 for (i = 0; i < 3; ++i) {
//                     for (j = 0; j < 3; ++j) {
//                         std::cout << std::setw(5) << T[i][j];
//                     }
//                     std::cout << std::endl;
//                 }
//                 std::cout << "S reciprocal" << std::endl;
//                 for (i = 0; i < 3; ++i) {
//                     for (j = 0; j < 3; ++j) {
//                         std::cout << std::setw(5) << S_recip[i][j];
//                     }
//                     std::cout << std::endl;
//                 }
//                 std::cout << std::endl;
// 
//                 det = S[0][0] * (S[1][1] * S[2][2] - S[2][1] * S[1][2])
//                     - S[1][0] * (S[0][1] * S[2][2] - S[2][1] * S[0][2])
//                     + S[2][0] * (S[0][1] * S[1][2] - S[1][1] * S[0][2]);
// 
//                 std::cout << "Det = " << det << std::endl;
//             }
// 
//             rotvec(Sk, k, S_recip);
//             for (i = 0; i < 3; ++i) Sk[i] = Sk[i] - nint(Sk[i]);
//             knum_sym = kpoint->get_knum(Sk[0], Sk[1], Sk[2]);
//             if (knum_sym == -1) error->exit("modify_eigenvectors_sym", "kpoint not found");
// 
//             calc_analytic_k(k, fcs_phonon->fc2_ext, dmat);
//             calc_analytic_k(Sk, fcs_phonon->fc2_ext, dmat_sym);
// 
//             if (!flag_done[knum_sym]) {
// 
//                 if (knum_sym == ik) {
//                     flag_done[knum_sym] = true;
//                     continue;
//                 }
// 
//                 if (knum_sym == kpoint->knum_minus[ik]) {
// 
//                     for (iat = 0; iat < natmin; ++iat) {
//                         jat = (*it).mapping[iat];
// 
//                         matmul3(mat_tmp, S, AA);
//                         matmul3(mat_tmp, BB, mat_tmp);
// 
//                         for (i = 0; i < 3; ++i) {
//                             x[i] = system->xr_p[system->map_p2s[iat][0]][i];
//                             x_new[i] = system->xr_p[system->map_p2s[jat][0]][i];
//                         }
// 
//                         for (i = 0; i < 3; ++i) x_tmp[i] = x[i];
//                         rotvec(x_tmp, x_tmp, mat_tmp);
// 
//                         phase = 2.0 * pi * (k[0] * x_new[0] + k[1] * x_new[1] + k[2] * x_new[2]) - (k[0] * x_tmp[0] + k[1] * x_tmp[1] + k[2] * x_tmp[2]);
//                         exp_phase = exp(im * phase);
// 
//                         for (jcrd = 0; jcrd < 3; ++jcrd) {
//                             for (icrd = 0; icrd < 3; ++icrd) {
//                                 gamma[3 * jat + jcrd][3 * iat + icrd] = -S[jcrd][icrd] * exp_phase;
//                             }
//                         }
//                     }
// 
//                 } else {
// 
//                     for (iat = 0; iat < natmin; ++iat) {
//                         jat = (*it).mapping[iat];
// 
//                         transpose3(mat_tmp, S);
//                         matmul3(mat_tmp, mat_tmp, AA);
//                         matmul3(mat_tmp, BB, mat_tmp);
// 
//                         for (i = 0; i < 3; ++i) {
//                             x[i] = system->xr_p[system->map_p2s[iat][0]][i];
//                             x_new[i] = system->xr_p[system->map_p2s[jat][0]][i];
//                         }
// 
//                         for (i = 0; i < 3; ++i) x_tmp[i] = x_new[i] - shift[i];
//                         rotvec(x_tmp, x_tmp, mat_tmp);
// 
//                         phase = k[0] * x_tmp[0] + k[1] * x_tmp[1] + k[2] * x_tmp[2] - 2.0 * pi * (k[0] * x[0] + k[1] * x[1] + k[2] * x[2]);
// 
//                         //	std::cout << "phase/2pi = " << phase / (2.0 * pi) << std::endl; 
//                         exp_phase = exp(im * phase);
// 
//                         for (jcrd = 0; jcrd < 3; ++jcrd) {
//                             for (icrd = 0; icrd < 3; ++icrd) {
//                                 gamma[3 * jat + jcrd][3 * iat + icrd] = S[jcrd][icrd] * exp_phase;
//                             }
//                         }
//                     }
//                 } // end if
// 
//                 for (i = 0; i < ns; ++i) {
//                     for (j = 0; j < ns; ++j) {
//                         evec_in[ns * i + j] = evec_phonon[ik][i][j];
//                         gamma_1d[ns * j + i] = gamma[i][j];
//                     }
//                 }
// 
//                 zgemm_(&TRANSA, &TRANSB, &ns, &ns, &ns, &alpha, gamma_1d, &ns, evec_in, &ns, &beta, evec_out, &ns);
// 
//                 for (i = 0; i < ns; ++i) {
//                     for (j = 0; j < ns; ++j) {
//                         evec_phonon[knum_sym][i][j] = evec_out[ns * i + j];
//                     }
//                 }
// 
//                 flag_done[knum_sym] = true;
//             }
//         }
//     }
// 
//     memory->deallocate(gamma);
//     memory->deallocate(gamma_1d);
//     memory->deallocate(evec_in);
//     memory->deallocate(evec_out);
//     memory->deallocate(gamma_dagger);
// 
//     memory->deallocate(dmat);
//     memory->deallocate(dmat_1d);
//     memory->deallocate(dmat_sym);
//     memory->deallocate(dmat_sym_1d);
// 
// }
