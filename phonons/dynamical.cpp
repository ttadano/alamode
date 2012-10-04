#include "mpi_common.h"
#include "dynamical.h"
#include "system.h"
#include "memory.h"
#include "kpoint.h"
#include <complex>
#include "../alm_c++/constants.h"
#include "fcs_phonon.h"
#include <iomanip>
#include "timer.h"

using namespace PHON_NS;

Dynamical::Dynamical(PHON *phon): Pointers(phon){
    eigenvectors = false;
}

Dynamical::~Dynamical(){}

void Dynamical::setup_dynamical(std::string mode)
{
    neval = 3 * system->natmin;
    UPLO = 'U';

    if (mode == "boltzmann") {
        eigenvectors = true;
    }
}

void Dynamical::eval_k(double *xk_in, double *eval_out, std::complex<double> **evec_out, bool require_evec) {

    // Calculate phonon energy for the specific k-point given in fractional basis

    unsigned int i, j;

    std::complex<double> **dymat_k;
    std::complex<double> **dymat_tmp, **dymat_transpose;

    memory->allocate(dymat_k, neval, neval);

    calc_analytic_k(dymat_k, xk_in);

    memory->allocate(dymat_tmp, neval, neval);
    memory->allocate(dymat_transpose, neval, neval);

    for (i = 0; i < neval; ++i){
        for (j = 0; j < neval; ++j){
            dymat_tmp[i][j] = dymat_k[i][j];
            dymat_transpose[i][j] = dymat_k[j][i];
        }
    }
    for (i = 0; i < neval; ++i){
        for (j = 0; j < neval; ++j){
            dymat_k[i][j] = 0.5 * (dymat_tmp[i][j] + std::conj(dymat_transpose[i][j]));
        }
    }

    memory->deallocate(dymat_tmp);
    memory->deallocate(dymat_transpose);

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

void Dynamical::calc_analytic_k(std::complex<double> **dymat_out, double *xk_in) {

    unsigned int i, j;
    unsigned int icrd, jcrd;
    unsigned int itran;
    unsigned int ntran = system->ntran;
    unsigned int natmin = system->natmin;
    unsigned int atm_p1, atm_p2, atm_s2;

    double phase;
    double vec[3];
    std::complex<double> ctmp[3][3];
    std::complex<double> im(0.0, 1.0);
    std::complex<double> exp_phase;

    for (i = 0; i < natmin; ++i){

        atm_p1 = system->map_p2s[i][0];

        for (j = 0; j < natmin; ++j){

            atm_p2 = system->map_p2s[j][0];

            for(icrd = 0; icrd < 3; ++icrd){
                for(jcrd = 0; jcrd < 3; ++jcrd){
                    ctmp[icrd][jcrd] = std::complex<double>(0.0, 0.0);
                }
            }

            for(itran = 0; itran < ntran; ++itran){
                atm_s2 =system->map_p2s[j][itran];

                for(icrd = 0; icrd < 3; ++icrd){                    
                    if (system->cell_dimension[icrd] = 1) {
                        vec[icrd] = system->xr_s[atm_p1][icrd] - system->xr_s[atm_s2][icrd];
                        if (vec[icrd] < -0.5) {
                            vec[icrd] = -1.0;
                        } else if (vec[icrd] >= 0.5){
                            vec[icrd] = 1.0;
                        } else {
                            vec[icrd] = 0.0;
                        }
                    } else if (system->cell_dimension[i] == 2) {
                         vec[icrd] = system->xr_s[atm_p1][icrd] - system->xr_s[atm_s2][icrd];
                         if (vec[icrd] < -0.5) {
                            vec[icrd] = -0.5;
                        } else if (vec[icrd] >= 0.5){
                            vec[icrd] = 0.5;
                        } else {
                            vec[icrd] = 0.0;
                        }
                    } else {
                        vec[icrd] = system->xr_s[atm_p2][icrd] - system->xr_s[atm_s2][icrd];
                        vec[icrd] = fold(vec[icrd]);
                    }
                }

                system->rotvec(vec, vec, system->lavec_s);
                system->rotvec(vec, vec, system->rlavec_p);

             //   std::cout << "r[" << atm_p1 << "] - r[" << atm_s2 << "] = " << vec[0] << " " << vec[1] << " " << vec[2] << std::endl;

                phase = vec[0] * xk_in[0] + vec[1] * xk_in[1] + vec[2] * xk_in[2];
                exp_phase = std::exp(im * phase);

                if (std::abs(exp_phase.real()) < 1.0e-14) {
                    exp_phase = std::complex<double>(0.0, exp_phase.imag());
                }
                if (std::abs(exp_phase.imag()) < 1.0e-14) {
                    exp_phase = std::complex<double>(exp_phase.real(), 0.0);
                }


                for (icrd = 0; icrd < 3; ++icrd){
                    for (jcrd = 0; jcrd < 3; ++jcrd){
                        ctmp[icrd][jcrd] += fcs_phonon->fc2[i][atm_s2][icrd][jcrd] * exp_phase;
                    }
                }
            }

            for (icrd = 0; icrd < 3; ++icrd){
                for (jcrd = 0; jcrd < 3; ++jcrd){
                    dymat_out[3 * i + icrd][3 * j + jcrd] = ctmp[icrd][jcrd] / std::sqrt(system->mass[atm_p1] * system->mass[atm_p2]);
                }
            }
        }
    }
}

void Dynamical::diagonalize_dynamical_all()
{
    unsigned int ik, is;
    unsigned int nk = kpoint->nk;
    bool require_evec; 

    if (mympi->my_rank == 0) {
        std::cout << std::endl << "Diagonalizing dynamical matrices for all k-points ..." << std::endl;
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

    for (ik = 0; ik < nk; ++ik){
        eval_k(kpoint->xk[ik], eval_phonon[ik], evec_phonon[ik], require_evec);

        // Phonon energy is the square-room of the eigenvalue 
        for (is = 0; is < neval; ++is){
            eval_phonon[ik][is] = freq(eval_phonon[ik][is]);
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    if (mympi->my_rank == 0) {
        timer->print_elapsed();
        std::cout << "done !" << std::endl;
    }
}

void Dynamical::calc_dynamical_matrix()
{
    // Calculate dynamical matrix of all k-points

    neval = 3 * system->natmin;
    memory->allocate(dymat, kpoint->nk, neval, neval);
    calc_analytic();

    std::complex<double> **dymat_tmp, **dymat_transpose;
    memory->allocate(dymat_tmp, neval, neval);
    memory->allocate(dymat_transpose, neval, neval);

    unsigned int nk = kpoint->nk;
    unsigned int i, j, ik;

    // Hermitize dynamical matrix

    for (ik = 0; ik < nk; ++ik){
        for (i = 0; i < neval; ++i){
            for (j = 0; j < neval; ++j){
                dymat_tmp[i][j] = dymat[ik][i][j];
                dymat_transpose[i][j] = dymat[ik][j][i];
            }
        }
        for (i = 0; i < neval; ++i){
            for (j = 0; j < neval; ++j){
                dymat[ik][i][j] = 0.5 * (dymat_tmp[i][j] + std::conj(dymat_transpose[i][j]));
            }
        }
    }

    memory->deallocate(dymat_tmp);
    memory->deallocate(dymat_transpose);
}

void Dynamical::calc_analytic()
{
    unsigned int i, j, itran, ik;
    unsigned int icrd, jcrd;

    unsigned int natmin = system->natmin;
    unsigned int ntran = system->ntran;

    unsigned int atm_p1, atm_p2, atm_s2;

    unsigned int nk = kpoint->nk;

    std::complex<double> ***ctmp;
    std::complex<double> im(0.0, 1.0);
    std::complex<double> exp_phase;

    double vec[3];

    double *phase;

    memory->allocate(ctmp, nk, 3, 3);
    memory->allocate(phase, nk);

    for (i = 0; i < natmin; ++i){

        atm_p1 = system->map_p2s[i][0];

        for (j = 0; j < natmin; ++j){

            atm_p2 = system->map_p2s[j][0];

            for(ik = 0; ik < nk; ++ik){
                for(icrd = 0; icrd < 3; ++icrd){
                    for(jcrd = 0; jcrd < 3; ++jcrd){
                        ctmp[ik][icrd][jcrd] = (0.0, 0.0);
                    }
                }
            }

            for(itran = 0; itran < ntran; ++itran){
                atm_s2 =system->map_p2s[j][itran];

                for(icrd = 0; icrd < 3; ++icrd){
                    vec[icrd] = system->xr_s[atm_p1][icrd] - system->xr_s[atm_s2][icrd];
                    vec[icrd] = fold(vec[icrd]);
                }

                system->rotvec(vec, vec, system->lavec_s);
                system->rotvec(vec, vec, system->rlavec_p);

                for (ik = 0; ik < nk; ++ik){
                    phase[ik] = vec[0] * kpoint->xk[ik][0] + vec[1] * kpoint->xk[ik][1] + vec[2] * kpoint->xk[ik][2];
                    exp_phase = std::exp(im * phase[ik]);

                    for (icrd = 0; icrd < 3; ++icrd){
                        for (jcrd = 0; jcrd < 3; ++jcrd){
                            ctmp[ik][icrd][jcrd] += fcs_phonon->fc2[i][atm_s2][icrd][jcrd] * exp_phase;
                        }
                    }
                }
            }

            for (ik = 0; ik < nk; ++ik){
                for (icrd = 0; icrd < 3; ++icrd){
                    for (jcrd = 0; jcrd < 3; ++jcrd){
                        dymat[ik][3 * i + icrd][3 * j + jcrd] = ctmp[ik][icrd][jcrd] / std::sqrt(system->mass[atm_p1] * system->mass[atm_p2]);
                    }
                }
            }

        }
    }

#ifdef _DEBUG
    for(ik = 0; ik < nk; ++ik){
        std::cout << "D" << std::endl;
        for (icrd = 0; icrd < 6; ++icrd){
            for (jcrd = 0; jcrd < 6; ++jcrd){
                std::cout << "(" << std::setw(15) << std::scientific << real(dymat[ik][icrd][jcrd]);
                std::cout << std::setw(15) << std::scientific << imag(dymat[ik][icrd][jcrd]) << ") ";
            }
            std::cout << std::endl;
        }
        std::cout << "D - D^{\\dagger}" << std::endl;
        for (icrd = 0; icrd < 6; ++icrd){
            for (jcrd = 0; jcrd < 6; ++jcrd){
                std::cout << "(" << std::setw(15) << std::scientific << real(dymat[ik][icrd][jcrd]) -real(dymat[ik][jcrd][icrd]) ;
                std::cout << std::setw(15) << std::scientific << imag(dymat[ik][icrd][jcrd]) + imag(dymat[ik][jcrd][icrd]) << ") ";
            }
            std::cout << std::endl;
        }
    } 
#endif

    memory->deallocate(ctmp);
    memory->deallocate(phase);
}

void Dynamical::diagonalize_dynamical()
{
    char JOBZ, UPLO = 'U';
    int INFO, LWORK;
    double *RWORK;
    std::complex<double> *WORK;

    if(eigenvectors){
        JOBZ = 'V';
    } else {
        JOBZ = 'N';
    }

    LWORK = (2 * neval - 1) * 10;
    memory->allocate(RWORK, 3*neval - 2);
    memory->allocate(WORK, LWORK);

    double *eval;
    std::complex<double> *amat;

    memory->allocate(eval, neval);
    memory->allocate(amat, neval * neval);

    int n = neval;
    unsigned int i, j;
    unsigned int ik, k;

    unsigned int nk = kpoint->nk;

    memory->allocate(eval_phonon, nk, neval);

    for (ik = 0; ik < nk; ++ik){
        k = 0;
        for(i = 0; i < neval; ++i){
            for (j = 0; j < neval; ++j){
                amat[k++] = dymat[ik][i][j];
            }
        }

        zheev_(&JOBZ, &UPLO, &n, amat, &n, eval, WORK, &LWORK, RWORK, &INFO);

        // save eigenvalues

        for (i = 0; i < neval; ++i){
            eval_phonon[ik][i] = eval[i];
        }

        // save eigenvectors if necessary

        if(eigenvectors){
            k = 0;
            for(i = 0; i < neval; ++i){
                for (j = 0; j < neval; ++j){
                    dymat[ik][i][j] = amat[k++];
                }
            }
        }
    }
    memory->deallocate(eval);
    memory->deallocate(amat);
}

double Dynamical::fold(double x)
{
    if (x >= -0.5 && x < 0.5) {
        return x;
    } else if (x < 0.0) {
        return x + 1.0;
    } else {
        return x - 1.0;
    }
}

double Dynamical::freq(const double x) 
{
    if (x >= 0.0) {
        return std::sqrt(x);
    } else {
        if (std::abs(x) < eps) {
            return std::sqrt(-x);
        } else { 
            return -std::sqrt(-x);
        }
    }
}
