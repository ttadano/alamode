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
#include <Eigen/Dense>

using namespace ALM_NS;

Fitting::Fitting(ALM *alm): Pointers(alm){}

Fitting::~Fitting() {
    memory->deallocate(params);
}

void Fitting::fitmain()
{
    files->ifs_disp_sym.open(files->file_disp_sym.c_str(), std::ios::in | std::ios::binary);
    if(!files->ifs_disp_sym) error->exit("fitmain", "cannot open file disp_sym");

    files->ifs_force_sym.open(files->file_force_sym.c_str(), std::ios::in | std::ios::binary);
    if(!files->ifs_force_sym) error->exit("fitmain", "cannot open file force_sym");

    int ntran = symmetry->ntran;
    int nat = system->nat;
    int natmin = symmetry->natmin;

    int i;
    int ndata = system->ndata;
    int nstart = system->nstart;
    int nend = system->nend;
    int nskip = system->nskip;

    int N, M;
    int maxorder = interaction->maxorder;
    int P = constraint->P;

    N = 0;
    for(i = 0; i < maxorder; ++i){
        N += fcs->ndup[i].size();
    }
    std::cout << "Total Number of Parameters : " << N << std::endl;
    M = 3 * natmin * ndata * ntran;

    memory->allocate(amat, M, N);
    memory->allocate(fsum, M);

    // Calculate matrix elements for fitting

    calc_matrix_elements(M, N, nat, natmin, ntran, ndata, maxorder);
    timer->print_elapsed();

    if (nskip == 0){

        // Fitting with singular value decomposition or QR-Decomposition

        int M_Start = 3 * natmin * ntran * (nstart - 1);
        int M_End   = 3 * natmin * ntran * nend;

        if(constraint->constraint_mode == 0) {
            fit_without_constraints(N, M_Start, M_End);
        } else {
            fit_with_constraints(N, M_Start, M_End, P);
        }
    } else {
        if (constraint->constraint_mode == 0) {
            error->exit("fitmain", "nskip has to be 0 when constraint = 0");
        } else {
            fit_consecutively(N, P, natmin, ntran, ndata, nstart, nend, nskip);
        }
    }

    // Copy force constants to public valiable "params"

    memory->allocate(params, N);

    for(i = 0; i < N; ++i) params[i] = fsum[i];

    memory->deallocate(amat);
    memory->deallocate(fsum);
    timer->print_elapsed();
}

void Fitting::fit_without_constraints(int N, int M_Start, int M_End)
{
    int i, j, k;
    int nrhs = 1, nrank, INFO, LWORK;
    double *WORK, *S, *amat_mod;
    double rcond = -1.0;
    double f_square = 0.0;
    double *fsum2;

    int M = M_End - M_Start;

    std::cout << "Entering Fitting Routine: SVD without constraints" << std::endl << std::endl;

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
        fsum2[j++] = fsum[i];
        f_square += std::pow(fsum[i], 2);
    }

    // fitting with singular value decomposition
    dgelss_(&M, &N, &nrhs, amat_mod, &M, fsum2, &M, S, &rcond, &nrank, WORK, &LWORK, &INFO);

    std::cout << "Finished !" << std::endl << std::endl;

    std::cout << "RANK of the MATRIX: " << nrank << std::endl;
    if(nrank < N) error->warn("fit_without_constraints", 
        "Matrix is rank-deficient. Force constants could not be determined uniquely :(");

    if(nrank == N) {
        double f_residual = 0.0;
        for (i = N; i < M; ++i){
            f_residual += std::pow(fsum2[i], 2);
        }
        std::cout << std::endl << "Residual sum of squares for the solution: " << sqrt(f_residual) << std::endl;
        std::cout << "Fitting Error (%) : "<< sqrt(f_residual/f_square) * 100.0 << std::endl;
    }

    for (i = 0; i < N; ++i){
        fsum[i] = fsum2[i];
    }
    memory->deallocate(fsum2);
    memory->deallocate(amat_mod);
    memory->deallocate(WORK);
    memory->deallocate(S);
}

void Fitting::fit_with_constraints(int N, int M_Start, int M_End, int P)
{
    int i, j, k;
    int nrank;
    double *mat_tmp;

    //    int maxorder = interaction->maxorder;

    double f_square, f_residual;
    double *fsum2;

    int M = M_End - M_Start;

    std::cout << "Entering Fitting Routine: QRD with constraints" << std::endl << std::endl;

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
            mat_tmp[k++] = constraint->const_mat[i][j];
        }
    }

    nrank = rank((M+P), N, mat_tmp);
    memory->deallocate(mat_tmp);

    if(nrank != N){
        std::cout << std::endl;
        std::cout << "!!WARNING: rank ( (A) ) ! = N" << std::endl;
        std::cout << "                ( (B) )      " << std::endl;
        std::cout << "rank = " << nrank << " N = " << N << std::endl;
        std::cout << "This must be a problem when solving equality constrained LSE problem with DGGLSE." << std::endl;
        std::cout << "Please change the cutoff radius or {u,f} data." << std::endl << std::endl;
    }

    f_square = 0.0;
    j = 0;
    for (i = M_Start; i < M_End; ++i){
        fsum2[j++] = fsum[i];
        f_square += std::pow(fsum[i], 2);
    }

    std::cout << "QR-Decomposition Started ...";

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
            cmat_mod[k++] = constraint->const_mat[i][j];
        }
    }

    // Fitting

    int LWORK = P + std::min<int>(M, N) + 100 * std::max<int>(M, N);
    int INFO;
    double *WORK, *x;
    memory->allocate(WORK, LWORK);
    memory->allocate(x, N);

    dgglse_(&M, &N, &P, amat_mod, &M, cmat_mod, &P, fsum2, constraint->const_rhs, x, WORK, &LWORK, &INFO);

    memory->deallocate(amat_mod);
    memory->deallocate(cmat_mod);
    memory->deallocate(WORK);


    std::cout << " finished. " << std::endl;

    f_residual = 0.0;
    for (i = N - P; i < M; ++i){
        f_residual += std::pow(fsum2[i], 2);
    }
    std::cout << std::endl << "Residual sum of squares for the solution: " << sqrt(f_residual) << std::endl;
    std::cout << "Fitting Error (%) : "<< std::sqrt(f_residual/f_square) * 100.0 << std::endl;

    // copy fcs to fsum

    for(i = 0; i < N; ++i){
        fsum[i] = x[i];
    }

    memory->deallocate(x);
    memory->deallocate(fsum2);
}

void Fitting::fit_consecutively(int N, int P, const int natmin, const int ntran, const int ndata,
    const int nstart, const int nend, const int nskip)
{
    int i, j, k;
    int iend;
    int M_Start, M_End;
    int mset = 3 * natmin * ntran;

    M_Start = mset * (nstart - 1);

    int M;

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

    std::cout << "nskip != 0: Consecutive Fitting Started!!" << std::endl;
    std::cout << "Relative errors and FCS are stored in file: " << file_fcs_sequence << std::endl;

    ofs_fcs_seq << "# Relative Error(%), FCS..." ;

    for (i = 0; i < interaction->maxorder; ++i){
        ofs_fcs_seq << std::setw(10) << fcs->ndup[i].size();
    }
    ofs_fcs_seq << std::endl;

    for (iend = nstart; iend <= nend; iend += nskip){

        M_End = mset * iend;
        M = M_End - M_Start;

        memory->allocate(fsum2, M);

        f_square = 0.0;
        j = 0;
        for (i = M_Start; i < M_End; ++i){
            fsum2[j++] = fsum[i];
            f_square += std::pow(fsum[i], 2);
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
                cmat_mod[k++] = constraint->const_mat[i][j];
            }
        }

        for (i = 0; i < P; ++i){
            const_tmp[i] = constraint->const_rhs[i];
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
        fsum[i] = x[i];
    }

    memory->deallocate(cmat_mod);
    memory->deallocate(const_tmp);
    memory->deallocate(x);

    ofs_fcs_seq.close();
}


void Fitting::calc_matrix_elements(const int M, const int N, const int nat, const int natmin,
    const int ntran, const int ndata, const int maxorder)
{
    int i, j;
    int itran;
    double **u;
    double **f;
    int ncycle;
    int irow;

    std::cout << "Calculation of Matrix Elements for Direct Fitting Started ..." << std::endl;
    for (i = 0; i < M; ++i){
        for (j = 0; j < N; ++j){
            amat[i][j] = 0.0;
        }
        fsum[i] = 0.0;
    }

    ncycle = ntran * ndata;
    memory->allocate(u, ncycle, 3 * nat);
    memory->allocate(f, ncycle, 3 * nat);

    // read all displacement-force data set

    for(int data = 0; data < ndata; ++data){
        for(itran = 0; itran < ntran; ++itran){

            irow = data * ntran + itran;

            for(i = 0; i < nat; ++i){
                for(j = 0; j < 3; ++j){
                    files->ifs_disp_sym.read((char *) &u[irow][3 * i + j], sizeof(double));
                    files->ifs_force_sym.read((char *) &f[irow][3 * i + j], sizeof(double));
                }
            }
        }
    }

#pragma omp parallel private(irow, i, j)
    { 
        int *ind;
        int mm, order, iat, k;
        int im, idata, iparam;
        double amat_tmp;

        memory->allocate(ind, maxorder + 1);

#pragma omp for schedule(guided)
        for(irow = 0; irow < ncycle; ++irow){

            // generate r.h.s vector B
            for (i = 0; i < natmin; ++i){
                iat = symmetry->map_p2s[i][0];
                for (j = 0; j < 3; ++j){
                    im = 3 * i + j + 3 * natmin * irow;
                    fsum[im] = f[irow][3 * iat + j];
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

    memory->deallocate(u);
    memory->deallocate(f);

    std::cout << " Finished !" << std::endl;
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

