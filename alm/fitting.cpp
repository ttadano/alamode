/*
 fitting.cpp

 Copyright (c) 2014, 2015, 2016 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include <boost/lexical_cast.hpp>
#include "fitting.h"
#include "files.h"
#include "error.h"
#include "memory.h"
#include "symmetry.h"
#include "system.h"
#include "fcs.h"
#include "interaction.h"
#include "timer.h"
#include "combination.h"
#include "constants.h"
#include "constraint.h"
#include "mathfunctions.h"

#ifdef _USE_EIGEN
#include <Eigen/Dense>
#endif

#include <time.h>

#ifdef _VSL
#include "mkl_vsl.h"

#else
#include <cstdlib>
#endif

using namespace ALM_NS;


Fitting::Fitting(ALM *alm): Pointers(alm)
{
    seed = static_cast<unsigned int>(time(NULL));
#ifdef _VSL
    brng = VSL_BRNG_MT19937;
    vslNewStream(&stream, brng, seed);
#else
    std::srand(seed);
#endif
}

Fitting::~Fitting()
{
    if (alm->mode == "fitting") {
        memory->deallocate(params);
    }
}

void Fitting::fitmain()
{
    int i;
    int nat = system->nat;
    int natmin = symmetry->nat_prim;
    int ntran = symmetry->ntran;

    int ndata = system->ndata;
    int nstart = system->nstart;
    int nend = system->nend;
    int nskip = system->nskip;

    int N, M, N_new;
    int maxorder = interaction->maxorder;
    int P = constraint->P;

    int nmulti;
    int ndata_used = nend - nstart + 1;

    double **u, **f;
    double **amat, *amat_1D, *fsum;
    double *fsum_orig;
    double *param_tmp;

    amat = NULL;
    amat_1D = NULL;
    fsum = NULL;
    fsum_orig = NULL;
    param_tmp = NULL;

    std::cout << " FITTING" << std::endl;
    std::cout << " =======" << std::endl << std::endl;

    std::cout << "  Reference files" << std::endl;
    std::cout << "   Displacement: " << files->file_disp << std::endl;
    std::cout << "   Force       : " << files->file_force << std::endl;
    std::cout << std::endl;

    std::cout << "  NSTART = " << nstart << "; NEND = " << nend << std::endl;
    std::cout << "  " << ndata_used << " entries will be used for fitting."
        << std::endl << std::endl;

    // Read displacement-force training data set from files

    data_multiplier(nat, ndata, nstart, nend, ndata_used, nmulti,
                    symmetry->multiply_data, u, f,
                    files->file_disp, files->file_force);

    N = 0;
    for (i = 0; i < maxorder; ++i) {
        N += fcs->nequiv[i].size();
    }
    std::cout << "  Total Number of Parameters : "
        << N << std::endl << std::endl;

    // Calculate matrix elements for fitting

    M = 3 * natmin * ndata_used * nmulti;

    if (constraint->constraint_algebraic) {


        N_new = 0;
        for (i = 0; i < maxorder; ++i) {
            N_new += constraint->index_bimap[i].size();
        }
        std::cout << "  Total Number of Free Parameters : "
            << N_new << std::endl << std::endl;

        //   memory->allocate(amat, M, N_new);
        memory->allocate(amat_1D, N_new * M);
        memory->allocate(fsum, M);
        memory->allocate(fsum_orig, M);

        calc_matrix_elements_algebraic_constraint(M, N, N_new, nat, natmin, ndata_used,
                                                  nmulti, maxorder, u, f, amat_1D, fsum,
                                                  fsum_orig);

    } else {

        memory->allocate(amat, M, N);
        memory->allocate(fsum, M);

        calc_matrix_elements(M, N, nat, natmin, ndata_used,
                             nmulti, maxorder, u, f, amat, fsum);
    }

    memory->deallocate(u);
    memory->deallocate(f);

    // Execute fitting

    memory->allocate(param_tmp, N);

    if (nskip == 0) {

        // Fitting with singular value decomposition or QR-Decomposition

        if (constraint->constraint_algebraic) {
            fit_algebraic_constraints(N_new, M, amat_1D, fsum, param_tmp,
                                      fsum_orig, maxorder);

        } else if (constraint->exist_constraint) {
            fit_with_constraints(N, M, P, amat, fsum, param_tmp,
                                 constraint->const_mat,
                                 constraint->const_rhs);
        } else {
            fit_without_constraints(N, M, amat, fsum, param_tmp);
        }

    } else if (nskip > 0) {

        // Execute fittings consecutively with different input data.

        if (constraint->exist_constraint) {
            fit_consecutively(N, P, natmin, ndata_used, nmulti, nskip, amat, fsum,
                              constraint->const_mat,
                              constraint->const_rhs);
        } else {
            error->exit("fitmain", "nskip has to be 0 when constraint_mode = 0");
        }
    } else {

        // Execute bootstrap simulation for estimating deviations of parameters.

        if (constraint->exist_constraint) {
            fit_bootstrap(N, P, natmin, ndata_used, nmulti, amat, fsum,
                          constraint->const_mat,
                          constraint->const_rhs);

            fit_with_constraints(N, M, P, amat, fsum, param_tmp,
                                 constraint->const_mat,
                                 constraint->const_rhs);
        } else {
            error->exit("fitmain",
                        "bootstrap analysis for LSE without constraint is not supported yet");
        }
    }

    // Copy force constants to public variable "params"

    memory->allocate(params, N);

    if (constraint->constraint_algebraic) {

        for (i = 0; i < N; ++i) {
            params[i] = param_tmp[i];
        }
        memory->deallocate(fsum_orig);

    } else {

        for (i = 0; i < N; ++i) params[i] = param_tmp[i];

    }

    if (amat) {
        memory->deallocate(amat);
    }
    if (amat_1D) {
        memory->deallocate(amat_1D);
    }
    if (fsum) {
        memory->deallocate(fsum);
    }
    if (param_tmp) {
        memory->deallocate(param_tmp);
    }

    std::cout << std::endl;
    timer->print_elapsed();
    std::cout << " -------------------------------------------------------------------" << std::endl;
    std::cout << std::endl;
}

void Fitting::data_multiplier(const int nat,
                              const int ndata,
                              const int nstart,
                              const int nend,
                              const int ndata_used,
                              int &nmulti,
                              const int multiply_data,
                              double **&u,
                              double **&f,
                              const std::string file_disp,
                              const std::string file_force)
{
    int i, j, k;
    int idata, itran, isym;
    int n_mapped;
    double u_in, f_in;
    double *u_tmp, *f_tmp;
    std::vector<int> vec_data;
    unsigned int nline_f, nline_u;
    unsigned int nreq;

    std::ifstream ifs_disp, ifs_force;

    ifs_disp.open(file_disp.c_str(), std::ios::in);
    if (!ifs_disp) error->exit("openfiles", "cannot open disp file");
    ifs_force.open(file_force.c_str(), std::ios::in);
    if (!ifs_force) error->exit("openfiles", "cannot open force file");

    nreq = 3 * nat * ndata;

    memory->allocate(u_tmp, nreq);
    memory->allocate(f_tmp, nreq);

    // Read displacements from DFILE

    nline_u = 0;
    while (ifs_disp >> u_in) {
        u_tmp[nline_u++] = u_in;
        if (nline_u == nreq) break;
    }
    if (nline_u < nreq)
        error->exit("data_multiplier",
                    "The number of lines in DFILE is too small for the given NDATA = ",
                    ndata);

    // Read forces from FFILE

    nline_f = 0;
    while (ifs_force >> f_in) {
        f_tmp[nline_f++] = f_in;
        if (nline_f == nreq) break;
    }
    if (nline_f < nreq)
        error->exit("data_multiplier",
                    "The number of lines in FFILE is too small for the given NDATA = ",
                    ndata);

    // Multiply data

    if (multiply_data == 0) {

        std::cout << " MULTDAT = 0: Given displacement-force data sets will be used as is."
            << std::endl << std::endl;

        nmulti = 1;

        memory->allocate(u, ndata_used * nmulti, 3 * nat);
        memory->allocate(f, ndata_used * nmulti, 3 * nat);

        idata = 0;

        for (i = 0; i < ndata; ++i) {
            if (i < nstart - 1) continue;
            if (i > nend - 1) break;

            for (j = 0; j < nat; ++j) {
                for (k = 0; k < 3; ++k) {
                    u[idata][3 * j + k] = u_tmp[3 * nat * i + 3 * j + k];
                    f[idata][3 * j + k] = f_tmp[3 * nat * i + 3 * j + k];
                }
            }
            ++idata;
        }

    } else if (multiply_data == 1) {

        std::cout << "  MULTDAT = 1: Generate symmetrically equivalent displacement-force " << std::endl;
        std::cout << "               data sets by using pure translational operations only." << std::endl << std::endl;

        nmulti = symmetry->ntran;

        memory->allocate(u, ndata_used * nmulti, 3 * nat);
        memory->allocate(f, ndata_used * nmulti, 3 * nat);

        idata = 0;

        for (i = 0; i < ndata; ++i) {
            if (i < nstart - 1) continue;
            if (i > nend - 1) break;

            for (itran = 0; itran < symmetry->ntran; ++itran) {
                for (j = 0; j < nat; ++j) {
                    n_mapped = symmetry->map_sym[j][symmetry->symnum_tran[itran]];

                    for (k = 0; k < 3; ++k) {
                        u[idata][3 * n_mapped + k] = u_tmp[3 * nat * i + 3 * j + k];
                        f[idata][3 * n_mapped + k] = f_tmp[3 * nat * i + 3 * j + k];
                    }
                }
                ++idata;
            }
        }

    } else if (multiply_data == 2) {

        double u_rot[3], f_rot[3];

        std::cout << "  MULTDAT = 2: Generate symmetrically equivalent displacement-force" << std::endl;
        std::cout << "                data sets. (including rotational part) " << std::endl << std::endl;

        nmulti = symmetry->nsym;

        memory->allocate(u, ndata_used * nmulti, 3 * nat);
        memory->allocate(f, ndata_used * nmulti, 3 * nat);

        idata = 0;

        for (i = 0; i < ndata; ++i) {
            if (i < nstart - 1) continue;
            if (i > nend - 1) break;

#pragma omp parallel for private(j, n_mapped, k, u_rot, f_rot)
            for (isym = 0; isym < symmetry->nsym; ++isym) {
                for (j = 0; j < nat; ++j) {
                    n_mapped = symmetry->map_sym[j][isym];

                    for (k = 0; k < 3; ++k) {
                        u_rot[k] = u_tmp[3 * nat * i + 3 * j + k];
                        f_rot[k] = f_tmp[3 * nat * i + 3 * j + k];
                    }

                    rotvec(u_rot, u_rot, symmetry->SymmData[isym].rotation_cart);
                    rotvec(f_rot, f_rot, symmetry->SymmData[isym].rotation_cart);

                    for (k = 0; k < 3; ++k) {
                        u[nmulti * idata + isym][3 * n_mapped + k] = u_rot[k];
                        f[nmulti * idata + isym][3 * n_mapped + k] = f_rot[k];
                    }
                }
            }
            ++idata;
        }

    } else {
        error->exit("data_multiplier", "Unsupported MULTDAT");
    }

    memory->deallocate(u_tmp);
    memory->deallocate(f_tmp);

    ifs_disp.close();
    ifs_force.close();
}

void Fitting::fit_without_constraints(int N,
                                      int M,
                                      double **amat,
                                      double *bvec,
                                      double *param_out)
{
    int i, j;
    unsigned long k;
    int nrhs = 1, nrank, INFO, LWORK;
    int LMIN, LMAX;
    double rcond = -1.0;
    double f_square = 0.0;
    double *WORK, *S, *amat_mod, *fsum2;

    std::cout << "  Entering fitting routine: SVD without constraints" << std::endl;

    LMIN = std::min<int>(M, N);
    LMAX = std::max<int>(M, N);

    LWORK = 3 * LMIN + std::max<int>(2 * LMIN, LMAX);
    LWORK = 2 * LWORK;

    memory->allocate(WORK, LWORK);
    memory->allocate(S, LMIN);

    // transpose matrix A
    memory->allocate(amat_mod, M * N);
    memory->allocate(fsum2, LMAX);

    k = 0;
    for (j = 0; j < N; ++j) {
        for (i = 0; i < M; ++i) {
            amat_mod[k++] = amat[i][j];
        }
    }
    for (i = 0; i < M; ++i) {
        fsum2[i] = bvec[i];
        f_square += std::pow(bvec[i], 2);
    }
    for (i = M; i < LMAX; ++i) fsum2[i] = 0.0;

    std::cout << "  SVD has started ... ";

    // Fitting with singular value decomposition
    dgelss_(&M, &N, &nrhs, amat_mod, &M, fsum2, &LMAX,
            S, &rcond, &nrank, WORK, &LWORK, &INFO);

    std::cout << "finished !" << std::endl << std::endl;

    std::cout << "  RANK of the matrix = " << nrank << std::endl;
    if (nrank < N)
        error->warn("fit_without_constraints",
                    "Matrix is rank-deficient. Force constants could not be determined uniquely :(");

    if (nrank == N) {
        double f_residual = 0.0;
        for (i = N; i < M; ++i) {
            f_residual += std::pow(fsum2[i], 2);
        }
        std::cout << std::endl << "  Residual sum of squares for the solution: "
            << sqrt(f_residual) << std::endl;
        std::cout << "  Fitting error (%) : "
            << sqrt(f_residual / f_square) * 100.0 << std::endl;
    }

    for (i = 0; i < N; ++i) {
        param_out[i] = fsum2[i];
    }

    memory->deallocate(WORK);
    memory->deallocate(S);
    memory->deallocate(fsum2);
    memory->deallocate(amat_mod);
}

void Fitting::fit_with_constraints(int N,
                                   int M,
                                   int P,
                                   double **amat,
                                   double *bvec,
                                   double *param_out,
                                   double **cmat,
                                   double *dvec)
{
    int i, j;
    unsigned long k;
    int nrank;
    double f_square, f_residual;
    double *fsum2;

    std::cout << "  Entering fitting routine: QRD with constraints" << std::endl;

    memory->allocate(fsum2, M);

#ifdef _USE_EIGEN_DISABLED

    double **mat_tmp2;
    memory->allocate(mat_tmp2, M + P, N);
    for (i = 0; i < M; ++i) {
        for (j = 0; j < N; ++j) {
            mat_tmp2[i][j] = amat[i][j];
        }
    }
    for (i = 0; i < P; ++i) {
        for (j = 0; j < N; ++j) {
            mat_tmp2[M + i][j] = cmat[i][j];
        }
    }

    nrank = getRankEigen(M+P, N, mat_tmp2);
    memory->deallocate(mat_tmp2);

#else

    double *mat_tmp;

    memory->allocate(mat_tmp, (M + P) * N);

    k = 0;

    for (j = 0; j < N; ++j) {
        for (i = 0; i < M; ++i) {
            mat_tmp[k++] = amat[i][j];
        }
        for (i = 0; i < P; ++i) {
            mat_tmp[k++] = cmat[i][j];
        }
    }

    nrank = rankQRD((M + P), N, mat_tmp, eps12);
    memory->deallocate(mat_tmp);

#endif

    if (nrank != N) {
        std::cout << std::endl;
        std::cout << " **************************************************************************" << std::endl;
        std::cout << "  WARNING : rank deficient.                                                " << std::endl;
        std::cout << "  rank ( (A) ) ! = N            A: Fitting matrix     B: Constraint matrix " << std::endl;
        std::cout << "       ( (B) )                  N: The number of parameters                " << std::endl;
        std::cout << "  rank = " << nrank << " N = " << N << std::endl << std::endl;
        std::cout << "  This can cause a difficulty in solving the fitting problem properly      " << std::endl;
        std::cout << "  with DGGLSE, especially when the difference is large. Please check if    " << std::endl;
        std::cout << "  you obtain reliable force constants in the .fcs file.                    " << std::endl << std::endl;
        std::cout << "  This issue may be resolved by setting MULTDAT = 2 in the &fitting field. " << std::endl;
        std::cout << "  If not, you may need to reduce the cutoff radii and/or increase NDATA    " << std::endl;
        std::cout << "  by giving linearly-independent displacement patterns.                    " << std::endl;
        std::cout << " **************************************************************************" << std::endl;
        std::cout << std::endl;
    }

    f_square = 0.0;
    for (i = 0; i < M; ++i) {
        fsum2[i] = bvec[i];
        f_square += std::pow(bvec[i], 2);
    }
    std::cout << "  QR-Decomposition has started ...";

    double *amat_mod, *cmat_mod;
    memory->allocate(amat_mod, M * N);
    memory->allocate(cmat_mod, P * N);

    // transpose matrix A and C
    k = 0;
    for (j = 0; j < N; ++j) {
        for (i = 0; i < M; ++i) {
            amat_mod[k++] = amat[i][j];
        }
    }
    k = 0;
    for (j = 0; j < N; ++j) {
        for (i = 0; i < P; ++i) {
            cmat_mod[k++] = cmat[i][j];
        }
    }

    // Fitting

    int LWORK = P + std::min<int>(M, N) + 10 * std::max<int>(M, N);
    int INFO;
    double *WORK, *x;
    memory->allocate(WORK, LWORK);
    memory->allocate(x, N);

    dgglse_(&M, &N, &P, amat_mod, &M, cmat_mod, &P,
            fsum2, dvec, x, WORK, &LWORK, &INFO);

    std::cout << " finished. " << std::endl;

    f_residual = 0.0;
    for (i = N - P; i < M; ++i) {
        f_residual += std::pow(fsum2[i], 2);
    }
    std::cout << std::endl << "  Residual sum of squares for the solution: "
        << sqrt(f_residual) << std::endl;
    std::cout << "  Fitting error (%) : "
        << std::sqrt(f_residual / f_square) * 100.0 << std::endl;

    // copy fcs to bvec

    for (i = 0; i < N; ++i) {
        param_out[i] = x[i];
    }

    memory->deallocate(amat_mod);
    memory->deallocate(cmat_mod);
    memory->deallocate(WORK);
    memory->deallocate(x);
    memory->deallocate(fsum2);
}

void Fitting::fit_algebraic_constraints(int N,
                                        int M,
                                        double *amat,
                                        double *bvec,
                                        double *param_out,
                                        double *bvec_orig,
                                        const int maxorder)
{
    int i, j;
    unsigned long k;
    int nrhs = 1, nrank, INFO, LWORK;
    int LMIN, LMAX;
    double rcond = -1.0;
    double f_square = 0.0;
    double *WORK, *S, *amat_mod, *fsum2;

    std::cout << "  Entering fitting routine: SVD with constraints considered algebraically." << std::endl;

    LMIN = std::min<int>(M, N);
    LMAX = std::max<int>(M, N);

    LWORK = 3 * LMIN + std::max<int>(2 * LMIN, LMAX);
    LWORK = 2 * LWORK;

    memory->allocate(WORK, LWORK);
    memory->allocate(S, LMIN);

    // transpose matrix A
    //  memory->allocate(amat_mod, M * N);
    memory->allocate(fsum2, LMAX);

    //k = 0;
    //for (j = 0; j < N; ++j) {
    //    for (i = 0; i < M; ++i) {
    //        amat_mod[k++] = amat[i][j];
    //    }
    //}
    for (i = 0; i < M; ++i) {
        fsum2[i] = bvec[i];
        f_square += std::pow(bvec_orig[i], 2);
    }
    for (i = M; i < LMAX; ++i) fsum2[i] = 0.0;

    std::cout << "  SVD has started ... ";

    // Fitting with singular value decomposition
    dgelss_(&M, &N, &nrhs, amat, &M, fsum2, &LMAX,
            S, &rcond, &nrank, WORK, &LWORK, &INFO);

    std::cout << "finished !" << std::endl << std::endl;

    std::cout << "  RANK of the matrix = " << nrank << std::endl;
    if (nrank < N)
        error->warn("fit_without_constraints",
                    "Matrix is rank-deficient. Force constants could not be determined uniquely :(");

    if (nrank == N) {
        double f_residual = 0.0;
        for (i = N; i < M; ++i) {
            f_residual += std::pow(fsum2[i], 2);
        }
        std::cout << std::endl;
        std::cout << "  Residual sum of squares for the solution: "
            << sqrt(f_residual) << std::endl;
        std::cout << "  Fitting error (%) : "
            << sqrt(f_residual / f_square) * 100.0 << std::endl;
    }

    int ishift = 0;
    int iparam = 0;
    double tmp;
    int inew, iold;

    for (i = 0; i < maxorder; ++i) {
        for (j = 0; j < constraint->const_fix[i].size(); ++j) {
            param_out[constraint->const_fix[i][j].p_index_target + ishift]
                = constraint->const_fix[i][j].val_to_fix;
        }

        for (boost::bimap<int, int>::const_iterator it = constraint->index_bimap[i].begin();
             it != constraint->index_bimap[i].end(); ++it) {
            inew = (*it).left + iparam;
            iold = (*it).right + ishift;

            param_out[iold] = fsum2[inew];
        }

        for (j = 0; j < constraint->const_relate[i].size(); ++j) {
            tmp = 0.0;

            for (k = 0; k < constraint->const_relate[i][j].alpha.size(); ++k) {
                tmp += constraint->const_relate[i][j].alpha[k]
                    * param_out[constraint->const_relate[i][j].p_index_orig[k] + ishift];
            }
            param_out[constraint->const_relate[i][j].p_index_target + ishift] = -tmp;
        }

        ishift += fcs->nequiv[i].size();
        iparam += constraint->index_bimap[i].size();
    }

    memory->deallocate(WORK);
    memory->deallocate(S);
    memory->deallocate(fsum2);
    //   memory->deallocate(amat_mod);
}


void Fitting::fit_bootstrap(int N,
                            int P,
                            int natmin,
                            int ndata_used,
                            int nmulti,
                            double **amat,
                            double *bvec,
                            double **cmat,
                            double *dvec)
{
    int i, j;
    unsigned long k, l;
    int M_Start, M_End;
    int mset;
    unsigned int iboot;
    int M;

    mset = 3 * natmin * nmulti;

    M_Start = 0;
    M_End = mset * ndata_used;

    M = M_End - M_Start;


    std::string file_fcs_bootstrap;
    file_fcs_bootstrap = files->job_title + ".fcs_bootstrap";

    std::ofstream ofs_fcs_boot;
    ofs_fcs_boot.open(file_fcs_bootstrap.c_str(), std::ios::out);
    if (!ofs_fcs_boot)
        error->exit("fit_bootstrap",
                    "cannot open file_fcs_bootstrap");

    ofs_fcs_boot.setf(std::ios::scientific);

    double f_residual, f_square;
    double *fsum2;
    int INFO;
    double *WORK, *x;
    double *amat_mod, *cmat_mod;
    double *const_tmp;

    int *rnd_index;
    int iloc;

    memory->allocate(x, N);
    memory->allocate(cmat_mod, P * N);
    memory->allocate(const_tmp, P);

    std::cout << "  NSKIP < 0: Bootstrap analysis for error estimation." << std::endl;
    std::cout << "             The number of trials is NBOOT (=" << nboot << ")" << std::endl;
    std::cout << std::endl;
    std::cout << "  Relative errors and FCs are stored in file: "
        << file_fcs_bootstrap << std::endl;

    ofs_fcs_boot << "# Relative Error(%), FCs ...";

    for (i = 0; i < interaction->maxorder; ++i) {
        ofs_fcs_boot << std::setw(10) << fcs->nequiv[i].size();
    }
    ofs_fcs_boot << std::endl;

    memory->allocate(fsum2, M);
    memory->allocate(amat_mod, N * M);
    memory->allocate(rnd_index, ndata_used);

    int LWORK = P + std::min<int>(M, N) + 100 * std::max<int>(M, N);
    memory->allocate(WORK, LWORK);

    for (iboot = 0; iboot < nboot; ++iboot) {
#ifdef _VSL
        // Use Intel MKL VSL if available
        viRngUniform(VSL_METHOD_IUNIFORM_STD, stream, ndata_used, rnd_index, 0, ndata_used);
#else
        for (i = 0; i < ndata_used; ++i) {
            // random number uniformly distributed in [0, ndata_used)
            rnd_index[i] = std::rand() % ndata_used;
        }
#endif

        f_square = 0.0;
        k = 0;
        for (i = 0; i < ndata_used; ++i) {
            iloc = rnd_index[i];
            for (j = iloc * mset; j < (iloc + 1) * mset; ++j) {
                fsum2[k++] = bvec[j];
                f_square += std::pow(bvec[j], 2);
            }
        }
        l = 0;
        for (j = 0; j < N; ++j) {
            for (i = 0; i < ndata_used; ++i) {
                iloc = rnd_index[i];
                for (k = iloc * mset; k < (iloc + 1) * mset; ++k) {
                    amat_mod[l++] = amat[k][j];
                }
            }
        }

        k = 0;
        for (j = 0; j < N; ++j) {
            for (i = 0; i < P; ++i) {
                cmat_mod[k++] = cmat[i][j];
            }
        }

        for (i = 0; i < P; ++i) {
            const_tmp[i] = dvec[i];
        }

        dgglse_(&M, &N, &P, amat_mod, &M, cmat_mod, &P,
                fsum2, const_tmp, x, WORK, &LWORK, &INFO);

        f_residual = 0.0;
        for (i = N - P; i < M; ++i) {
            f_residual += std::pow(fsum2[i], 2);
        }

        ofs_fcs_boot << 100.0 * std::sqrt(f_residual / f_square);

        for (i = 0; i < N; ++i) {
            ofs_fcs_boot << std::setw(15) << x[i];
        }
        ofs_fcs_boot << std::endl;
    }
    ofs_fcs_boot.close();

    memory->deallocate(x);
    memory->deallocate(fsum2);
    memory->deallocate(cmat_mod);
    memory->deallocate(const_tmp);
    memory->deallocate(amat_mod);
    memory->deallocate(rnd_index);
    memory->deallocate(WORK);

    std::cout << "  Bootstrap analysis finished." << std::endl;
    std::cout << "  Normal fitting will be performed" << std::endl;
}


void Fitting::fit_consecutively(int N,
                                int P,
                                const int natmin,
                                const int ndata_used,
                                const int nmulti,
                                const int nskip,
                                double **amat,
                                double *bvec,
                                double **cmat,
                                double *dvec)
{
    int i, j;
    unsigned long k;
    int iend;
    int M_Start, M_End;
    int mset;
    int M;

    mset = 3 * natmin * nmulti;

    M_Start = 0;


    std::string file_fcs_sequence;
    file_fcs_sequence = files->job_title + ".fcs_sequence";

    std::ofstream ofs_fcs_seq;
    ofs_fcs_seq.open(file_fcs_sequence.c_str(), std::ios::out);
    if (!ofs_fcs_seq)
        error->exit("fit_consecutively",
                    "cannot open file_fcs_sequence");

    ofs_fcs_seq.setf(std::ios::scientific);

    double f_residual, f_square;
    double *fsum2;
    int INFO;
    double *WORK, *x;
    double *amat_mod, *cmat_mod;
    double *const_tmp;

    memory->allocate(x, N);
    memory->allocate(cmat_mod, P * N);
    memory->allocate(const_tmp, P);

    std::cout << "  NSKIP > 0: Fitting will be performed consecutively" << std::endl;
    std::cout << "             with variously changing NEND as NEND = NSTART + i*NSKIP" << std::endl;
    std::cout << std::endl;
    std::cout << "  Relative errors and FCs will be stored in the file "
        << file_fcs_sequence << std::endl;

    ofs_fcs_seq << "# Relative Error(%), FCS...";

    for (i = 0; i < interaction->maxorder; ++i) {
        ofs_fcs_seq << std::setw(10) << fcs->nequiv[i].size();
    }
    ofs_fcs_seq << std::endl;

    for (iend = 1; iend <= ndata_used; iend += nskip) {

        M_End = mset * iend;
        M = M_End - M_Start;

        memory->allocate(fsum2, M);

        f_square = 0.0;
        j = 0;
        for (i = M_Start; i < M_End; ++i) {
            fsum2[j++] = bvec[i];
            f_square += std::pow(bvec[i], 2);
        }

        memory->allocate(amat_mod, M * N);

        // Transpose matrix A
        k = 0;
        for (j = 0; j < N; ++j) {
            for (i = M_Start; i < M_End; ++i) {
                amat_mod[k++] = amat[i][j];
            }
        }

        k = 0;
        for (j = 0; j < N; ++j) {
            for (i = 0; i < P; ++i) {
                cmat_mod[k++] = cmat[i][j];
            }
        }

        for (i = 0; i < P; ++i) {
            const_tmp[i] = dvec[i];
        }

        // Fitting

        int LWORK = P + std::min<int>(M, N) + 100 * std::max<int>(M, N);
        memory->allocate(WORK, LWORK);

        dgglse_(&M, &N, &P, amat_mod, &M, cmat_mod, &P,
                fsum2, const_tmp, x, WORK, &LWORK, &INFO);

        memory->deallocate(amat_mod);
        memory->deallocate(WORK);

        f_residual = 0.0;
        for (i = N - P; i < M; ++i) {
            f_residual += std::pow(fsum2[i], 2);
        }

        ofs_fcs_seq << 100.0 * std::sqrt(f_residual / f_square);

        for (i = 0; i < N; ++i) {
            ofs_fcs_seq << std::setw(15) << x[i];
        }
        ofs_fcs_seq << std::endl;
        memory->deallocate(fsum2);
    }

    for (i = 0; i < N; ++i) {
        bvec[i] = x[i];
    }

    memory->deallocate(cmat_mod);
    memory->deallocate(const_tmp);
    memory->deallocate(x);

    ofs_fcs_seq.close();

    std::cout << "  Consecutive fitting finished." << std::endl;
}

void Fitting::calc_matrix_elements(const int M,
                                   const int N,
                                   const int nat,
                                   const int natmin,
                                   const int ndata_fit,
                                   const int nmulti,
                                   const int maxorder,
                                   double **u,
                                   double **f,
                                   double **amat,
                                   double *bvec)
{
    int i, j;
    int irow;
    int ncycle;

    std::cout << "  Calculation of matrix elements for direct fitting started ... ";
    for (i = 0; i < M; ++i) {
        for (j = 0; j < N; ++j) {
            amat[i][j] = 0.0;
        }
        bvec[i] = 0.0;
    }

    ncycle = ndata_fit * nmulti;

#ifdef _OPENMP
#pragma omp parallel private(irow, i, j)
#endif
    {
        int *ind;
        int mm, order, iat, k;
        int im, idata, iparam;
        double amat_tmp;

        memory->allocate(ind, maxorder + 1);

#ifdef _OPENMP
#pragma omp for schedule(guided)
#endif
        for (irow = 0; irow < ncycle; ++irow) {

            // generate r.h.s vector B
            for (i = 0; i < natmin; ++i) {
                iat = symmetry->map_p2s[i][0];
                for (j = 0; j < 3; ++j) {
                    im = 3 * i + j + 3 * natmin * irow;
                    bvec[im] = f[irow][3 * iat + j];
                }
            }

            // generate l.h.s. matrix A

            idata = 3 * natmin * irow;
            iparam = 0;

            for (order = 0; order < maxorder; ++order) {

                mm = 0;

                for (std::vector<int>::iterator iter = fcs->nequiv[order].begin();
                     iter != fcs->nequiv[order].end(); ++iter) {
                    for (i = 0; i < *iter; ++i) {
                        ind[0] = fcs->fc_table[order][mm].elems[0];
                        k = idata + inprim_index(fcs->fc_table[order][mm].elems[0]);
                        amat_tmp = 1.0;
                        for (j = 1; j < order + 2; ++j) {
                            ind[j] = fcs->fc_table[order][mm].elems[j];
                            amat_tmp *= u[irow][fcs->fc_table[order][mm].elems[j]];
                        }
                        amat[k][iparam] -= gamma(order + 2, ind) * fcs->fc_table[order][mm].sign * amat_tmp;
                        ++mm;
                    }
                    ++iparam;
                }
            }
        }

        memory->deallocate(ind);

    }

    std::cout << "done!" << std::endl << std::endl;
}


void Fitting::calc_matrix_elements_algebraic_constraint(const int M,
                                                        const int N,
                                                        const int N_new,
                                                        const int nat,
                                                        const int natmin,
                                                        const int ndata_fit,
                                                        const int nmulti,
                                                        const int maxorder,
                                                        double **u,
                                                        double **f,
                                                        double *amat,
                                                        double *bvec,
                                                        double *bvec_orig)
{
    int i, j;
    int irow;
    int ncycle;

    std::cout << "  Calculation of matrix elements for direct fitting started ... ";

    ncycle = ndata_fit * nmulti;
    int natmin3 = 3 * natmin;


#ifdef _OPENMP
#pragma omp parallel for private(j)
#endif
    for (i = 0; i < M; ++i) {
        //        for (j = 0; j < N_new; ++j) {
        //           amat[i][j] = 0.0;
        //       }
        bvec[i] = 0.0;
        bvec_orig[i] = 0.0;
    }

#ifdef _OPENMP
#pragma omp parallel private(irow, i, j)
#endif
    {
        int *ind;
        int mm, order, iat, k;
        int im, idata, iparam;
        int ishift;
        int iold, inew;
        double amat_tmp;
        double **amat_orig;
        double **amat_mod;

        memory->allocate(ind, maxorder + 1);
        memory->allocate(amat_orig, natmin3, N);
        memory->allocate(amat_mod, natmin3, N_new);

#ifdef _OPENMP
#pragma omp for schedule(guided)
#endif
        for (irow = 0; irow < ncycle; ++irow) {

            // generate r.h.s vector B
            for (i = 0; i < natmin; ++i) {
                iat = symmetry->map_p2s[i][0];
                for (j = 0; j < 3; ++j) {
                    im = 3 * i + j + natmin3 * irow;
                    bvec[im] = f[irow][3 * iat + j];
                    bvec_orig[im] = f[irow][3 * iat + j];
                }
            }

            for (i = 0; i < natmin3; ++i) {
                for (j = 0; j < N; ++j) {
                    amat_orig[i][j] = 0.0;
                }
                for (j = 0; j < N_new; ++j) {
                    amat_mod[i][j] = 0.0;
                }
            }

            // generate l.h.s. matrix A

            idata = natmin3 * irow;
            iparam = 0;

            for (order = 0; order < maxorder; ++order) {

                mm = 0;

                for (std::vector<int>::iterator iter = fcs->nequiv[order].begin();
                     iter != fcs->nequiv[order].end(); ++iter) {
                    for (i = 0; i < *iter; ++i) {
                        ind[0] = fcs->fc_table[order][mm].elems[0];
                        k = inprim_index(ind[0]);

                        amat_tmp = 1.0;
                        for (j = 1; j < order + 2; ++j) {
                            ind[j] = fcs->fc_table[order][mm].elems[j];
                            amat_tmp *= u[irow][fcs->fc_table[order][mm].elems[j]];
                        }
                        amat_orig[k][iparam] -= gamma(order + 2, ind) * fcs->fc_table[order][mm].sign * amat_tmp;
                        ++mm;
                    }
                    ++iparam;
                }
            }

            ishift = 0;
            iparam = 0;

            for (order = 0; order < maxorder; ++order) {

                for (i = 0; i < constraint->const_fix[order].size(); ++i) {

                    for (j = 0; j < natmin3; ++j) {
                        bvec[j + idata] -= constraint->const_fix[order][i].val_to_fix
                            * amat_orig[j][ishift + constraint->const_fix[order][i].p_index_target];
                    }
                }

                for (boost::bimap<int, int>::const_iterator it = constraint->index_bimap[order].begin();
                     it != constraint->index_bimap[order].end(); ++it) {
                    inew = (*it).left + iparam;
                    iold = (*it).right + ishift;

                    for (j = 0; j < natmin3; ++j) {
                        amat_mod[j][inew] = amat_orig[j][iold];
                    }
                }

                for (i = 0; i < constraint->const_relate[order].size(); ++i) {

                    iold = constraint->const_relate[order][i].p_index_target + ishift;

                    for (j = 0; j < constraint->const_relate[order][i].alpha.size(); ++j) {

                        inew = constraint->index_bimap[order].right.at(
                                                                 constraint->const_relate[order][i].p_index_orig[j])
                            + iparam;
                        for (k = 0; k < natmin3; ++k) {
                            amat_mod[k][inew] -= amat_orig[k][iold] * constraint->const_relate[order][i].alpha[j];
                        }
                    }
                }

                ishift += fcs->nequiv[order].size();
                iparam += constraint->index_bimap[order].size();
            }

            for (i = 0; i < natmin3; ++i) {
                for (j = 0; j < N_new; ++j) {
                    //           amat[i + idata][j] = amat_mod[i][j];
                    // Transpose the matrix for later use of lapack
                    amat[natmin3 * ncycle * j + i + idata] = amat_mod[i][j];
                }
            }

        }

        memory->deallocate(ind);
        memory->deallocate(amat_orig);
        memory->deallocate(amat_mod);
    }

    std::cout << "done!" << std::endl << std::endl;
}


int Fitting::inprim_index(const int n)
{
    int in;
    int atmn = n / 3;
    int crdn = n % 3;

    for (int i = 0; i < symmetry->nat_prim; ++i) {
        if (symmetry->map_p2s[i][0] == atmn) {
            in = 3 * i + crdn;
            break;
        }
    }
    return in;
}

double Fitting::gamma(const int n, const int *arr)
{
    int *arr_tmp, *nsame;
    int i;
    int ind_front, nsame_to_front;

    memory->allocate(arr_tmp, n);
    memory->allocate(nsame, n);

    for (i = 0; i < n; ++i) {
        arr_tmp[i] = arr[i];
        nsame[i] = 0;
    }

    ind_front = arr[0];
    nsame_to_front = 1;

    interaction->insort(n, arr_tmp);

    int nuniq = 1;
    int iuniq = 0;

    nsame[0] = 1;

    for (i = 1; i < n; ++i) {
        if (arr_tmp[i] == arr_tmp[i - 1]) {
            ++nsame[iuniq];
        } else {
            ++nsame[++iuniq];
            ++nuniq;
        }

        if (arr[i] == ind_front) ++nsame_to_front;
    }

    int denom = 1;

    for (i = 0; i < nuniq; ++i) {
        denom *= factorial(nsame[i]);
    }

    memory->deallocate(arr_tmp);
    memory->deallocate(nsame);

    return static_cast<double>(nsame_to_front) / static_cast<double>(denom);
}

int Fitting::factorial(const int n)
{
    if (n == 1 || n == 0) {
        return 1;
    } else {
        return n * factorial(n - 1);
    }
}

#ifdef _USE_EIGEN_DISABLED
int Fitting::getRankEigen(const int m, const int n, double **mat)
{
    using namespace Eigen;

    MatrixXd mat_tmp(m, n);

    int i, j;

    for (i = 0; i < m; ++i) {
        for (j = 0; j < n; ++j) {
            mat_tmp(i,j) = mat[i][j];
        }
    }
    ColPivHouseholderQR<MatrixXd> qr(mat_tmp);
    return qr.rank();
}
#endif

int Fitting::rankQRD(const int m,
                     const int n,
                     double *mat,
                     const double tolerance)
{
    // Return the rank of matrix mat revealed by the column pivoting QR decomposition
    // The matrix mat is destroyed.

    int m_ = m;
    int n_ = n;

    int LDA = m_;

    int LWORK = 10 * n_;
    int INFO;
    int *JPVT;
    double *WORK, *TAU;

    int nmin = std::min<int>(m_, n_);

    memory->allocate(JPVT, n_);
    memory->allocate(WORK, LWORK);
    memory->allocate(TAU, nmin);

    for (int i = 0; i < n_; ++i) JPVT[i] = 0;

    dgeqp3_(&m_, &n_, mat, &LDA, JPVT, TAU, WORK, &LWORK, &INFO);

    memory->deallocate(JPVT);
    memory->deallocate(WORK);
    memory->deallocate(TAU);

    if (std::abs(mat[0]) < eps) return 0;

    double **mat_tmp;
    memory->allocate(mat_tmp, m_, n_);

    unsigned long k = 0;

    for (int j = 0; j < n_; ++j) {
        for (int i = 0; i < m_; ++i) {
            mat_tmp[i][j] = mat[k++];
        }
    }

    int nrank = 0;
    for (int i = 0; i < nmin; ++i) {
        if (std::abs(mat_tmp[i][i]) > tolerance * std::abs(mat[0])) ++nrank;
    }

    memory->deallocate(mat_tmp);

    return nrank;
}

int Fitting::rankSVD(const int m,
                     const int n,
                     double *mat,
                     const double tolerance)
{
    int i;
    int m_ = m;
    int n_ = n;

    int LWORK = 10 * m;
    int INFO;
    int *IWORK;
    int ldu = 1, ldvt = 1;
    double *s, *WORK;
    double u[1], vt[1];

    int nmin = std::min<int>(m, n);

    memory->allocate(IWORK, 8 * nmin);
    memory->allocate(WORK, LWORK);
    memory->allocate(s, nmin);

    char mode[] = "N";

    dgesdd_(mode, &m_, &n_, mat, &m_, s, u, &ldu, vt, &ldvt,
            WORK, &LWORK, IWORK, &INFO);

    int rank = 0;
    for (i = 0; i < nmin; ++i) {
        if (s[i] > s[0] * tolerance) ++rank;
    }

    memory->deallocate(WORK);
    memory->deallocate(IWORK);
    memory->deallocate(s);

    return rank;
}

int Fitting::rankSVD2(const int m_in,
                      const int n_in,
                      double **mat,
                      const double tolerance)
{
    // Reveal the rank of matrix mat without destroying the matrix elements

    int i, j, k;
    double *arr;

    int m = m_in;
    int n = n_in;

    memory->allocate(arr, m * n);

    k = 0;

    for (j = 0; j < n; ++j) {
        for (i = 0; i < m; ++i) {
            arr[k++] = mat[i][j];
        }
    }

    int LWORK = 10 * m;
    int INFO;
    int *IWORK;
    int ldu = 1, ldvt = 1;
    double *s, *WORK;
    double u[1], vt[1];

    int nmin = std::min<int>(m, n);

    memory->allocate(IWORK, 8 * nmin);
    memory->allocate(WORK, LWORK);
    memory->allocate(s, nmin);

    char mode[] = "N";

    dgesdd_(mode, &m, &n, arr, &m, s, u, &ldu, vt, &ldvt,
            WORK, &LWORK, IWORK, &INFO);

    int rank = 0;
    for (i = 0; i < nmin; ++i) {
        if (s[i] > s[0] * tolerance) ++rank;
    }

    memory->deallocate(IWORK);
    memory->deallocate(WORK);
    memory->deallocate(s);
    memory->deallocate(arr);

    return rank;
}
