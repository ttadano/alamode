/*
 dense_hermitian_eigen.cpp

 Copyright (c) 2026 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory
 or http://opensource.org/licenses/mit-license.php for information.
*/

#include "dense_hermitian_eigen.h"
#include <mpi.h>
#include "error.h"
#include "ndarray.h"

extern "C"
{
    void zheev_(const char *jobz, const char *uplo, int *n, std::complex<double> *a, int *lda, double *w,
                std::complex<double> *work, int *lwork, double *rwork, int *info);
    void zheevd_(const char *jobz, const char *uplo, int *n, std::complex<double> *a, int *lda, double *w,
                 std::complex<double> *work, int *lwork, double *rwork, int *lrwork, int *iwork, int *liwork,
                 int *info);
}

using namespace PHON_NS;

void PHON_NS::solve_dense_hermitian(int n, const std::complex<double> *const *mat_in, double *eval_out,
                                    std::complex<double> **evec_out, bool compute_evec, char uplo)
{
    if (solve_dense_hermitian_info(n, mat_in, eval_out, evec_out, compute_evec, uplo) != 0) {
        exit("solve_dense_hermitian", "zheev failed to diagonalize the Hermitian matrix (INFO != 0).");
    }
}

int PHON_NS::solve_dense_hermitian_info(int n, const std::complex<double> *const *mat_in, double *eval_out,
                                        std::complex<double> **evec_out, bool compute_evec, char uplo)
{
    int INFO;
    int LWORK = (2 * n - 1) * 10;
    NDArray<std::complex<double>, 1> amat(n * n);
    NDArray<std::complex<double>, 1> WORK(LWORK);
    NDArray<double, 1> RWORK(3 * n - 2);

    unsigned int k = 0;
    for (int j = 0; j < n; ++j) {
        for (int i = 0; i < n; ++i) {
            amat[k++] = mat_in[i][j];
        }
    }

    char JOBZ = compute_evec ? 'V' : 'N';

    zheev_(&JOBZ, &uplo, &n, amat, &n, eval_out, WORK, &LWORK, RWORK, &INFO);
    if (INFO != 0) return INFO;

    if (evec_out) {
        for (int j = 0; j < n; ++j) {
            for (int i = 0; i < n; ++i) {
                evec_out[j][i] = amat[j * n + i];
            }
        }
    }
    return 0;
}

void PHON_NS::solve_dense_hermitian_dc(const Eigen::MatrixXcd &mat_in, Eigen::VectorXd &eval_out,
                                       Eigen::MatrixXcd *evec_out)
{
    int n = static_cast<int>(mat_in.rows());
    if (n == 0) {
        eval_out.resize(0);
        return;
    }
    // zheevd overwrites the matrix with the eigenvectors (column-major, as Eigen)
    Eigen::MatrixXcd amat = mat_in;
    eval_out.resize(n);
    char JOBZ = evec_out ? 'V' : 'N';
    char UPLO = 'L';
    int INFO = 0;

    // workspace query
    std::complex<double> work_query;
    double rwork_query;
    int iwork_query;
    int LWORK = -1, LRWORK = -1, LIWORK = -1;
    zheevd_(&JOBZ,
            &UPLO,
            &n,
            amat.data(),
            &n,
            eval_out.data(),
            &work_query,
            &LWORK,
            &rwork_query,
            &LRWORK,
            &iwork_query,
            &LIWORK,
            &INFO);
    if (INFO != 0) {
        exit("solve_dense_hermitian_dc", "zheevd workspace query failed (INFO != 0).");
    }
    LWORK = static_cast<int>(work_query.real());
    LRWORK = static_cast<int>(rwork_query);
    LIWORK = iwork_query;
    NDArray<std::complex<double>, 1> WORK(LWORK);
    NDArray<double, 1> RWORK(LRWORK);
    NDArray<int, 1> IWORK(LIWORK);

    zheevd_(&JOBZ, &UPLO, &n, amat.data(), &n, eval_out.data(), WORK, &LWORK, RWORK, &LRWORK, IWORK, &LIWORK, &INFO);
    if (INFO != 0) {
        exit("solve_dense_hermitian_dc", "zheevd failed to diagonalize the Hermitian matrix (INFO != 0).");
    }
    if (evec_out) {
        *evec_out = amat;
    }
}
