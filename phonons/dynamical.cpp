#include "dynamical.h"
#include "system.h"
#include "memory.h"
#include "kpoint.h"
#include <complex>
#include "../alm_c++/constants.h"
#include "fcs_phonon.h"
#include <iomanip>

using namespace PHON_NS;

Dynamical::Dynamical(PHON *phon): Pointers(phon){}

Dynamical::~Dynamical(){}


void Dynamical::setup_dynamical()
{
    neval = 3 * system->natmin;
}

void Dynamical::eval_k(double *eval_out, double *xk_in) {

    // Calculate phonon energy for the given k-point

    int i, j;

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

    char JOBZ = 'N', UPLO = 'U';
    int INFO, LWORK;
    double *RWORK;
    std::complex<double> *WORK;

    LWORK = (2 * neval - 1) * 10;
    memory->allocate(RWORK, 3*neval - 2);
    memory->allocate(WORK, LWORK);

   // double *eval;
    std::complex<double> *amat;

   // memory->allocate(eval, neval);
    memory->allocate(amat, neval * neval);

    unsigned int k = 0;
    int n = dynamical->neval;
    for(i = 0; i < neval; ++i){
        for (j = 0; j < neval; ++j){
            amat[k++] = dymat_k[i][j];
        }
    }

    zheev_(&JOBZ, &UPLO, &n, amat, &n, eval_out, WORK, &LWORK, RWORK, &INFO);

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
                    ctmp[icrd][jcrd] = (0.0, 0.0);
                }
            }

            for(itran = 0; itran < ntran; ++itran){
                atm_s2 =system->map_p2s[j][itran];

                for(icrd = 0; icrd < 3; ++icrd){
                    vec[icrd] = system->xr_s[atm_p1][icrd] - system->xr_s[atm_s2][icrd];
                    vec[icrd] = fold(vec[icrd]);
                }

                system->rotvec(system->lavec_s, vec, vec);
                system->rotvec(system->rlavec_p, vec, vec);


                phase = vec[0] * xk_in[0] + vec[1] * xk_in[1] + vec[2] * xk_in[2];
                exp_phase = std::exp(im * phase);

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

void Dynamical::calc_dynamical_matrix()
{
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

                system->rotvec(system->lavec_s, vec, vec);
                system->rotvec(system->rlavec_p, vec, vec);

                /*  for(icrd = 0; icrd < 3; ++icrd){
                vec[icrd] /= 2.0 * pi;
                }*/

                for (ik = 0; ik < nk; ++ik){
                    // phase[ik] = 2.0 * pi * (vec[0] * kpoint->xk[ik][0] + vec[1] * kpoint->xk[ik][1] + vec[2] * kpoint->xk[ik][2]);
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
    double xabs = std::abs(x);
    if(xabs <= 0.5){
        return x;
    } else {
        if(x < 0.0) {
            return x + 1.0;
        } else {
            return x - 1.0;
        }
    }
}
