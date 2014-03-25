#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include <set>
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


Fitting::Fitting(ALM *alm): Pointers(alm){
    seed = (unsigned int) time(NULL);
#ifdef _VSL
    brng = VSL_BRNG_MT19937;
    vslNewStream(&stream, brng, seed);
#else
    std::srand(seed);
#endif
}

Fitting::~Fitting() {
    if (alm->mode == "fitting") {
        memory->deallocate(params);
    }
}

void Fitting::fitmain()
{
    int nat = system->nat;
    int natmin = symmetry->natmin;
    int ntran = symmetry->ntran;

    int i;
    int ndata = system->ndata;
    int nstart = system->nstart;
    int nend = system->nend;
    int nskip = system->nskip;

    int N, M;
    int maxorder = interaction->maxorder;
    int P = constraint->P;

    int M_Start, M_End;

    int nmulti;
    int ndata_used = nend - nstart + 1;

    double **amat, *fsum;


    std::cout << " FITTING" << std::endl;
    std::cout << " =======" << std::endl << std::endl;

    std::cout << "  Reference files" << std::endl;
    std::cout << "   Displacement: " << files->file_disp << std::endl;
    std::cout << "   Force       : " << files->file_force << std::endl;
    std::cout << std::endl;

    std::cout << "  NSTART = " << nstart << "; NEND = " << nend << std::endl;
    std::cout << "  " << nend - nstart + 1 << " entries will be used for fitting." << std::endl << std::endl;

    data_multiplier(nat, ndata, nstart, nend, ndata_used, nmulti,symmetry->multiply_data);

    N = 0;
    for(i = 0; i < maxorder; ++i){
        N += fcs->ndup[i].size();
    }
    std::cout << "  Total Number of Parameters : " << N << std::endl << std::endl;

    M = 3 * natmin * ndata_used * nmulti;

    memory->allocate(amat, M, N);
    memory->allocate(fsum, M);

    // Calculate matrix elements for fitting

    calc_matrix_elements(M, N, nat, natmin, ndata_used, nmulti, maxorder, amat, fsum);

    /*
    // Calculate Hessian Matrix of Amat

    std::cout << "Calculating Covariance Matrix for error estimation" << std::endl;
    memory->allocate(varcovar, N, N);
    calc_covariance(M, N);
    timer->print_elapsed(); */

    M_Start = 0;
    M_End = 3 * natmin * ndata_used * nmulti;

    if (nskip == 0){

        // Fitting with singular value decomposition or QR-Decomposition

        if (constraint->exist_constraint) {
            fit_with_constraints(N, M_Start, M_End, P, amat, fsum, constraint->const_mat, constraint->const_rhs);
        } else {
            fit_without_constraints(N, M_Start, M_End, amat, fsum);
        }

    } else if (nskip > 0) {

        // Execute fittings consecutively with different input data.

        if (constraint->exist_constraint) {
            fit_consecutively(N, P, natmin, ndata_used, nmulti, nskip, amat, fsum, constraint->const_mat, constraint->const_rhs);
        } else {
            error->exit("fitmain", "nskip has to be 0 when constraint_mode = 0");
        }
    } else {

        // Execute bootstrap simulation for estimating deviations of parameters.

        if (constraint->exist_constraint) {
            fit_bootstrap(N, P, natmin, ndata_used, nmulti, amat, fsum, constraint->const_mat, constraint->const_rhs);
            fit_with_constraints(N, M_Start, M_End, P, amat, fsum, constraint->const_mat, constraint->const_rhs);
        } else {
            error->exit("fitmain", "bootstrap analysis for LSE without constraint is not supported yet");
        }
    }

    // Copy force constants to public variable "params"

    memory->allocate(params, N);

    for(i = 0; i < N; ++i) params[i] = fsum[i];

    memory->deallocate(amat);
    memory->deallocate(fsum);

    std::cout << std::endl;
    timer->print_elapsed();
    std::cout << " --------------------------------------------------------------" << std::endl;
    std::cout << std::endl;

}

void Fitting::data_multiplier(const int nat, const int ndata, const int nstart, const int nend, 
                              const int ndata_used, int &nmulti, const int multiply_data) 
{
    int i, j, k;
    int idata, itran, isym;
    int n_mapped;
    double u_rot[3], f_rot[3];
    double ***u_tmp, ***f_tmp;

    memory->allocate(u_tmp, ndata, nat, 3);
    memory->allocate(f_tmp, ndata, nat, 3);

    files->ifs_disp.seekg(0, std::ios::beg);
    files->ifs_force.seekg(0, std::ios::beg);

    for (i = 0; i < ndata; ++i) {
        for (j = 0; j < nat; ++j) {
            files->ifs_disp >> u_tmp[i][j][0] >> u_tmp[i][j][1] >> u_tmp[i][j][2];
            files->ifs_force >> f_tmp[i][j][0] >> f_tmp[i][j][1] >> f_tmp[i][j][2];
        }
    }

    if (multiply_data == 0) {

        std::cout << " MULTDAT = 0: Given displacement-force data sets will be used as is." << std::endl << std::endl;

        nmulti = 1;

        memory->allocate(u, ndata_used * nmulti, 3 * nat);
        memory->allocate(f, ndata_used * nmulti, 3 * nat);

        idata = 0;

        for (i = 0; i < ndata; ++i) {
            if (i < nstart - 1) continue;
            if (i > nend - 1) break;

            for (j = 0; j < nat; ++j) {
                for (k = 0; k < 3; ++k) {
                    u[idata][3 * j + k] = u_tmp[i][j][k];
                    f[idata][3 * j + k] = f_tmp[i][j][k];
                }
            }
            ++idata;
        }

    } else if (multiply_data == 1) {

        std::cout << "  MULTDAT = 1: Generate symmetrically equivalent displacement-force data sets " << std::endl;
        std::cout << "               by using pure translational operations only." << std::endl << std::endl;

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
                        u[idata][3 * n_mapped + k] = u_tmp[i][j][k];
                        f[idata][3 * n_mapped + k] = f_tmp[i][j][k];
                    }
                }
                ++idata;
            }
        }

    } else if (multiply_data == 2) {

        std::cout << "  MULTDAT = 2: Generate symmetrically equivalent displacement-force data sets." << std::endl;
        std::cout << "               (including rotational part) " << std::endl << std::endl;

        nmulti = symmetry->nsym;

        memory->allocate(u, ndata_used * nmulti, 3 * nat);
        memory->allocate(f, ndata_used * nmulti, 3 * nat);

        idata = 0;

        for (i = 0; i < ndata; ++i) {
            if (i < nstart - 1) continue;
            if (i > nend - 1) break;

            for (isym = 0; isym < symmetry->nsym; ++isym) {
                for (j = 0; j < nat; ++j) {
                    n_mapped = symmetry->map_sym[j][isym];

                    rotvec(u_rot, u_tmp[i][j], symmetry->symrel[isym]);
                    rotvec(f_rot, f_tmp[i][j], symmetry->symrel[isym]);

                    for (k = 0; k < 3; ++k) {
                        u[idata][3 * n_mapped + k] = u_rot[k];
                        f[idata][3 * n_mapped + k] = f_rot[k];
                    }
                }
                ++idata;
            }
        }

    } else {
        error->exit("data_multiplier", "Unsupported MULTDAT");
    }

    memory->deallocate(u_tmp);
    memory->deallocate(f_tmp);
}

void Fitting::fit_without_constraints(int N, int M_Start, int M_End, double **amat, double *bvec)
{
    int i, j, k;
    int nrhs = 1, nrank, INFO, LWORK;
    double *WORK, *S, *amat_mod;
    double rcond = -1.0;
    double f_square = 0.0;
    double *fsum2;

    int M = M_End - M_Start;

    std::cout << "  Entering fitting routine: SVD without constraints" << std::endl;

    LWORK = 3 * std::min<int>(M, N) + std::max<int>(2*std::min<int>(M, N), std::max<int>(M, N));
    LWORK = 2 * LWORK;

    memory->allocate(WORK, LWORK);
    memory->allocate(S, N);

    // transpose matrix A
    memory->allocate(amat_mod, M * N);
    memory->allocate(fsum2, M);

    k = 0;
    for (j = 0; j < N; ++j){
        for(i = M_Start; i < M_End; ++i){
            amat_mod[k++] = amat[i][j];
        }
    }
    j = 0;
    for (i = M_Start; i < M_End; ++i){
        fsum2[j++] = bvec[i];
        f_square += std::pow(bvec[i], 2);
    }

    std::cout << "  SVD has started ... ";
    // fitting with singular value decomposition
    dgelss_(&M, &N, &nrhs, amat_mod, &M, fsum2, &M, S, &rcond, &nrank, WORK, &LWORK, &INFO);

    std::cout << "finished !" << std::endl << std::endl;

    std::cout << "  RANK of the matrix = " << nrank << std::endl;
    if(nrank < N) error->warn("fit_without_constraints", 
        "Matrix is rank-deficient. Force constants could not be determined uniquely :(");

    if(nrank == N) {
        double f_residual = 0.0;
        for (i = N; i < M; ++i){
            f_residual += std::pow(fsum2[i], 2);
        }
        std::cout << std::endl << "  Residual sum of squares for the solution: " << sqrt(f_residual) << std::endl;
        std::cout << "  Fitting error (%) : "<< sqrt(f_residual/f_square) * 100.0 << std::endl;
    }

    for (i = 0; i < N; ++i){
        bvec[i] = fsum2[i];
    }
    memory->deallocate(fsum2);
    memory->deallocate(amat_mod);
    memory->deallocate(WORK);
    memory->deallocate(S);
}

void Fitting::fit_with_constraints(int N, int M_Start, int M_End, int P, double **amat, double *bvec, double **cmat, double *dvec)
{
    int i, j, k;
    int nrank;
    double *mat_tmp;
    double f_square, f_residual;
    double *fsum2;

    int M = M_End - M_Start;

    std::cout << "  Entering fitting routine: QRD with constraints" << std::endl;

    memory->allocate(mat_tmp, (M + P) * N);
    memory->allocate(fsum2, M);

    k = 0;

    for(j = 0; j < N; ++j){
        for(i = M_Start; i < M_End; ++i){
            mat_tmp[k++] = amat[i][j];
        }
    }

    for(j = 0; j < N; ++j){
        for(i = 0; i < P; ++i){
            mat_tmp[k++] = cmat[i][j];
        }
    }

    nrank = rank((M+P), N, mat_tmp);
    memory->deallocate(mat_tmp);

    if(nrank != N){
        std::cout << std::endl;
        std::cout << "  WARNING: rank ( (A) ) ! = N" << std::endl;
        std::cout << "                ( (B) )      " << std::endl;
        std::cout << "  rank = " << nrank << " N = " << N << std::endl;
        std::cout << "  This must be a problem when solving equality constrained LSE problem with DGGLSE." << std::endl;
        std::cout << "  Please change the cutoff radius or {u,f} data." << std::endl << std::endl;
    }

    f_square = 0.0;
    j = 0;
    for (i = M_Start; i < M_End; ++i){
        fsum2[j++] = bvec[i];
        f_square += std::pow(bvec[i], 2);
    }

    std::cout << "  QR-Decomposition has started ...";

    double *amat_mod, *cmat_mod;
    memory->allocate(amat_mod, M * N);
    memory->allocate(cmat_mod, P * N);

    // transpose matrix A and C
    k = 0;
    for(j = 0; j < N; ++j){
        for(i = M_Start; i < M_End; ++i){
            amat_mod[k++] = amat[i][j];
        }
    }
    k = 0;
    for (j = 0; j < N; ++j){
        for(i = 0; i < P; ++i){
            cmat_mod[k++] = cmat[i][j];
        }
    }

    // Fitting

    int LWORK = P + std::min<int>(M, N) + 100 * std::max<int>(M, N);
    int INFO;
    double *WORK, *x;
    memory->allocate(WORK, LWORK);
    memory->allocate(x, N);

    dgglse_(&M, &N, &P, amat_mod, &M, cmat_mod, &P, fsum2, dvec, x, WORK, &LWORK, &INFO);

    memory->deallocate(amat_mod);
    memory->deallocate(cmat_mod);
    memory->deallocate(WORK);

    std::cout << " finished. " << std::endl;

    f_residual = 0.0;
    for (i = N - P; i < M; ++i){
        f_residual += std::pow(fsum2[i], 2);
    }
    std::cout << std::endl << "  Residual sum of squares for the solution: " << sqrt(f_residual) << std::endl;
    std::cout << "  Fitting error (%) : "<< std::sqrt(f_residual/f_square) * 100.0 << std::endl;

    // copy fcs to bvec

    for(i = 0; i < N; ++i){
        bvec[i] = x[i];
    }

    memory->deallocate(x);
    memory->deallocate(fsum2);
}

void Fitting::fit_bootstrap(int N, int P, int natmin, int ndata_used, int nmulti, double **amat, double *bvec, double **cmat, double *dvec)
{
    int i, j, k, l;
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
    if(!ofs_fcs_boot) error->exit("fit_bootstrap", "cannot open file_fcs_bootstrap");

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
    std::cout << "  Relative errors and FCs are stored in file: " << file_fcs_bootstrap << std::endl;

    ofs_fcs_boot << "# Relative Error(%), FCs ..." ;

    for (i = 0; i < interaction->maxorder; ++i){
        ofs_fcs_boot << std::setw(10) << fcs->ndup[i].size();
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
        for (i = 0; i < ndata_used; ++i){
            rnd_index[i] = std::rand() % ndata_used; // random number uniformly distributed in [0, ndata_used)
        }
#endif

        f_square = 0.0;
        k = 0;
        for (i = 0; i < ndata_used; ++i){
            iloc = rnd_index[i];
            for (j = iloc * mset; j < (iloc + 1) * mset; ++j){
                fsum2[k++] = bvec[j];
                f_square += std::pow(bvec[j], 2);
            }
        }
        l = 0;
        for(j = 0; j < N; ++j){
            for(i = 0; i < ndata_used; ++i){
                iloc = rnd_index[i];
                for (k = iloc * mset; k < (iloc + 1) * mset; ++k){
                    amat_mod[l++] = amat[k][j];
                }
            }
        }

        k = 0;
        for (j = 0; j < N; ++j){
            for(i = 0; i < P; ++i){
                cmat_mod[k++] = cmat[i][j];
            }
        }

        for (i = 0; i < P; ++i){
            const_tmp[i] = dvec[i];
        }     

        dgglse_(&M, &N, &P, amat_mod, &M, cmat_mod, &P, fsum2, const_tmp, x, WORK, &LWORK, &INFO);

        f_residual = 0.0;
        for (i = N - P; i < M; ++i){
            f_residual += std::pow(fsum2[i], 2);
        }

        ofs_fcs_boot << 100.0 * std::sqrt(f_residual/f_square);

        for (i = 0; i < N; ++i){
            ofs_fcs_boot << std::setw(15) << x[i];
        }
        ofs_fcs_boot << std::endl;
    }
    ofs_fcs_boot.close();

    std::cout << "  Bootstrap analysis finished." << std::endl;
    std::cout << "  Normal fitting will be performed" << std::endl;
}

void Fitting::fit_consecutively(int N, int P, const int natmin, const int ndata_used, const int nmulti, const int nskip, double **amat, double *bvec, double **cmat, double *dvec)
{
    int i, j, k;
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
    if(!ofs_fcs_seq) error->exit("fit_consecutively", "cannot open file_fcs_sequence");

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
    std::cout << "             with variously changing NEMD as NEND = NSTART + i*NSKIP" << std::endl;
    std::cout << std::endl;
    std::cout << "  Relative errors and FCs will be stored in the file " << file_fcs_sequence << std::endl;

    ofs_fcs_seq << "# Relative Error(%), FCS..." ;

    for (i = 0; i < interaction->maxorder; ++i){
        ofs_fcs_seq << std::setw(10) << fcs->ndup[i].size();
    }
    ofs_fcs_seq << std::endl;

    for (iend = 1; iend <= ndata_used; iend += nskip) {

        M_End = mset * iend;
        M = M_End - M_Start;

        memory->allocate(fsum2, M);

        f_square = 0.0;
        j = 0;
        for (i = M_Start; i < M_End; ++i){
            fsum2[j++] = bvec[i];
            f_square += std::pow(bvec[i], 2);
        }

        memory->allocate(amat_mod, M * N);

        // Transpose matrix A
        k = 0;
        for(j = 0; j < N; ++j){
            for(i = M_Start; i < M_End; ++i){
                amat_mod[k++] = amat[i][j];
            }
        }

        k = 0;
        for (j = 0; j < N; ++j){
            for(i = 0; i < P; ++i){
                cmat_mod[k++] = cmat[i][j];
            }
        }

        for (i = 0; i < P; ++i){
            const_tmp[i] = dvec[i];
        }

        // Fitting

        int LWORK = P + std::min<int>(M, N) + 100 * std::max<int>(M, N);
        memory->allocate(WORK, LWORK);

        dgglse_(&M, &N, &P, amat_mod, &M, cmat_mod, &P, fsum2, const_tmp, x, WORK, &LWORK, &INFO);

        memory->deallocate(amat_mod);
        memory->deallocate(WORK);

        f_residual = 0.0;
        for (i = N - P; i < M; ++i){
            f_residual += std::pow(fsum2[i], 2);
        }

        ofs_fcs_seq << 100.0 * std::sqrt(f_residual/f_square);

        for (i = 0; i < N; ++i){
            ofs_fcs_seq << std::setw(15) << x[i];
        }
        ofs_fcs_seq << std::endl;
        memory->deallocate(fsum2);
    }

    for (i = 0; i < N; ++i){
        bvec[i] = x[i];
    }

    memory->deallocate(cmat_mod);
    memory->deallocate(const_tmp);
    memory->deallocate(x);

    ofs_fcs_seq.close();

    std::cout << "  Consecutive fitting finished." << std::endl;
}

void Fitting::calc_matrix_elements(const int M, const int N, const int nat, const int natmin, const int ndata_fit, const int nmulti, const int maxorder, double **amat, double *bvec)
{
    int i, j;
    int irow;
    int ncycle;

    std::cout << "  Calculation of matrix elements for direct fitting started ... ";
    for (i = 0; i < M; ++i){
        for (j = 0; j < N; ++j){
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
        for(irow = 0; irow < ncycle; ++irow){

            // generate r.h.s vector B
            for (i = 0; i < natmin; ++i){
                iat = symmetry->map_p2s[i][0];
                for (j = 0; j < 3; ++j){
                    im = 3 * i + j + 3 * natmin * irow;
                    bvec[im] = f[irow][3 * iat + j];
                }
            }

            // generate l.h.s. matrix A

            idata = 3 * natmin * irow;
            iparam = 0;

            for(order = 0; order < maxorder; ++order){

                mm = 0;

                for(std::vector<int>::iterator iter = fcs->ndup[order].begin(); iter != fcs->ndup[order].end(); ++iter){
                    for (i = 0; i < *iter; ++i){
                        ind[0] = fcs->fc_set[order][mm].elems[0];
                        k = idata + inprim_index(fcs->fc_set[order][mm].elems[0]);
                        amat_tmp = 1.0;
                        for(j = 1; j < order + 2; ++j){
                            ind[j] = fcs->fc_set[order][mm].elems[j];
                            amat_tmp *= u[irow][fcs->fc_set[order][mm].elems[j]];
                        }
                        amat[k][iparam] -= gamma(order + 2, ind) * fcs->fc_set[order][mm].coef * amat_tmp;
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


int Fitting::inprim_index(const int n)
{
    int in;
    int atmn = n / 3;
    int crdn = n % 3;

    for (int i = 0; i < symmetry->natmin; ++i){
        if(symmetry->map_p2s[i][0] == atmn){
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

    for(i = 1; i < n; ++i){
        if(arr_tmp[i] == arr_tmp[i-1]) {
            ++nsame[iuniq];
        } else {
            ++nsame[++iuniq];
            ++nuniq;
        }

        if(arr[i] == ind_front) ++nsame_to_front;
    }

    int denom = 1;

    for(i = 0; i < nuniq; ++i){
        denom *= factorial(nsame[i]);
    }

    memory->deallocate(arr_tmp);
    memory->deallocate(nsame);

    return static_cast<double>(nsame_to_front) / static_cast<double>(denom);
}

int Fitting::factorial(const int n)
{
    if(n == 1 || n == 0){
        return 1;
    }else{
        return n * factorial(n - 1);
    }
}

//int Fitting::rank(const int m, const int n, double **mat)
//{
//  using namespace Eigen;
// 
//    MatrixXd mat_tmp(m, n);
//
//    int i, j;
//
//    for(i = 0; i < m; ++i){
//        for(j = 0; j < n; ++j){
//            mat_tmp(i,j) = mat[i][j];
//        }
//    }
//    ColPivHouseholderQR<MatrixXd> qr(mat_tmp);
//    return qr.rank();
//}


/* unsuccessful
int Fitting::getRankEigen(int m, int p, int n)
{
using namespace Eigen;
MatrixXd mat(m + p, n);


int i, j;
int k = 0;
for(i = 0; i < m; ++i){
for(j = 0; j < n; ++j){
mat(k, j) = amat[i][j];
}
++k;
}
std::cout <<"OK";

for(i = 0; i < p; ++i){
for(j = 0; j < n; ++j){
mat(k, j) = const_mat[i][j];
}
++k;
}
std::cout <<"OK";
ColPivHouseholderQR<MatrixXd> qr(mat);
return qr.rank();
}
*/

//int Fitting::rank(int m, int n, double *mat)
//{
//    int LWORK = 10 * n;
//    int INFO;
//    double *WORK, *TAU;
//    int lda = m;
//    
//    int nmin = std::min<int>(m, n);
//
//    memory->allocate(WORK, LWORK);
//    memory->allocate(TAU, n);
//
//    dgeqrf_(&m, &n, mat, TAU, WORK, &LWORK, &INFO);
//
//
//}

int Fitting::rank(int m, int n, double *mat)
{
    int i;

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

    char mode[]  = "N";

    dgesdd_(mode, &m, &n, mat, &m, s, u, &ldu, vt, &ldvt, WORK, &LWORK, IWORK, &INFO); 

    int rank = 0;
    for(i = 0; i < nmin; ++i){
        if(s[i] > eps12) ++rank;
    }

    memory->deallocate(IWORK);
    memory->deallocate(s);

    return rank;
}

int Fitting::rank2(const int m_in, const int n_in, double **mat) 
{
    // Reveal the rank of matrix mat without destroying the matrix elements

    int i, j, k;
    double *arr;

    int m = m_in;
    int n = n_in;

    memory->allocate(arr, m*n);

    k = 0;

    for (j = 0; j < n; ++j ) {
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

    char mode[]  = "N";

    dgesdd_(mode, &m, &n, arr, &m, s, u, &ldu, vt, &ldvt, WORK, &LWORK, IWORK, &INFO); 

    int rank = 0;
    for(i = 0; i < nmin; ++i){
        if(s[i] > eps12) ++rank;
    }

    memory->deallocate(IWORK);
    memory->deallocate(s);

    memory->deallocate(arr);

    return rank;
}

/*
void Fitting::calc_covariance(int m, int n)
{
Eigen::MatrixXd Atmp(m, n), Hess(n, n);
int i, j;

for (i = 0; i < m; ++i){
for (j = 0; j < n; ++j){
Atmp(i,j) = amat[i][j];
}
}

Hess = (Atmp.transpose()*Atmp).inverse();

for (i = 0; i < n; ++i){
for (j = 0; j < n; ++j){
varcovar[i][j] = Hess(i, j);
}
}

}
*/
