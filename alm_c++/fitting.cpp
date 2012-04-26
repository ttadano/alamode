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

    int N, M, P;
    int maxorder = interaction->maxorder;

    N = 0;
    for(i = 0; i < maxorder; ++i){
        N += fcs->ndup[i].size();
    }
    std::cout << "Total Number of Parameters : " << N << std::endl;
    M = 3 * natmin * ndata * ntran;

    memory->allocate(amat, M, N);
    memory->allocate(fsum, M);

    if (constraint != 0){

        // Generate constraint matrix

        int Pmax, order;

        memory->allocate(const_translation, maxorder);
        translational_invariance();

        memory->allocate(const_rotation_self, maxorder);
        memory->allocate(const_rotation_cross, maxorder);
        rotational_invariance();

        Pmax = 0;
        for (order = 0; order < maxorder; ++order){
            //   Pmax += const_translation[order].size() + const_rotation[order].size();
            Pmax += const_translation[order].size();
        }
        if(constraint == 2){
            //    Pmax -= const_translation[0].size() + const_rotation[0].size();
            Pmax -= const_translation[0].size();
            Pmax += fcs->ndup[0].size();
        }
        memory->allocate(const_mat, Pmax, N);
        memory->allocate(const_rhs, Pmax);

        calc_constraint_matrix(N, P);
        std::cout << "Total number of constraints: " << P << std::endl;
    }

    // Calculate matrix elements for fitting

    calc_matrix_elements(M, N, nat, natmin, ntran, ndata, maxorder);
    timer->print_elapsed();

    if (nskip == 0){

        // Fitting with singular value decomposition or QR-Decomposition

        int M_Start = 3 * natmin * ntran * (nstart - 1);
        int M_End   = 3 * natmin * ntran * nend;

        if(constraint == 0) {
            fit_without_constraints(N, M_Start, M_End);
        } else {
            fit_with_constraints(N, M_Start, M_End, P);
        }
    } else {
        if (constraint == 0) {
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

    if (constraint != 0){
        memory->deallocate(const_mat);
        memory->deallocate(const_rhs);
        memory->deallocate(const_translation);
        memory->deallocate(const_rotation_self);
        memory->deallocate(const_rotation_cross);
    }
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
            mat_tmp[k++] = const_mat[i][j];
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
            cmat_mod[k++] = const_mat[i][j];
        }
    }

    // Fitting

    int LWORK = P + std::min<int>(M, N) + 100 * std::max<int>(M, N);
    int INFO;
    double *WORK, *x;
    memory->allocate(WORK, LWORK);
    memory->allocate(x, N);

    dgglse_(&M, &N, &P, amat_mod, &M, cmat_mod, &P, fsum2, const_rhs, x, WORK, &LWORK, &INFO);

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
                cmat_mod[k++] = const_mat[i][j];
            }
        }

        for (i = 0; i < P; ++i){
            const_tmp[i] = const_rhs[i];
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

void Fitting::calc_constraint_matrix(const int N, int &P){

    int i, j;
    int maxorder = interaction->maxorder;
    int order;
    int icol, irow;
    int nrow, ncol;
    int *nrank, *nparam;
    double *arr_tmp;
    std::vector<Constraint> *const_vec;

    memory->allocate(nrank, maxorder);
    memory->allocate(nparam, maxorder);
    memory->allocate(const_vec, maxorder);

    std::cout << "Removing redundant constraints ...";

    using namespace Eigen;

    for(order = 0; order < maxorder; ++order){

        const_vec[order].clear();

        nrow = fcs->ndup[order].size();
        ncol = const_translation[order].size();

        nparam[order] = nrow;

        memory->allocate(arr_tmp, nrow);

        MatrixXd mat_tmp(nrow, ncol);
        icol = 0;

        for (std::set<Constraint>::iterator p = const_translation[order].begin(); p != const_translation[order].end(); ++p){
            Constraint const_now = *p;
            for (i = 0; i < nrow; ++i){
                mat_tmp(i, icol) = const_now.w_const[i];
            }
            ++icol;
        }

        FullPivLU<MatrixXd> lu_decomp(mat_tmp);
        nrank[order] = lu_decomp.rank();
        MatrixXd c_reduced = lu_decomp.image(mat_tmp);

        for(icol = 0; icol < nrank[order]; ++icol){
            for(irow = 0; irow < nrow; ++irow){
                arr_tmp[irow] = c_reduced(irow, icol);
            }

            const_vec[order].push_back(Constraint(nrow, arr_tmp));
        }
        memory->deallocate(arr_tmp);
    }

    std::cout << " done." << std::endl << std::endl;

    std::cout << "Rank of the constraint matrices for each order ..." << std::endl;
    P = 0;
    for (order = 0; order < maxorder; ++order){
        P += nrank[order];
        std::cout << std::setw(9) << interaction->str_order[order] << ": " << const_vec[order].size() << std::endl;
    }
    std::cout << std::endl;

    if(constraint == 2) {
        std::cout << "Harmonic Force Constants will be fixed to the values in the given reference file: " << fc2_file << std::endl;
        std::cout << "Constraint Matrix for Harmonic fcs will be updated." << std::endl << std::endl;
        P =  P - nrank[0] + nparam[0];
    }

    for(i = 0; i < P; ++i){
        for(j = 0; j < N; ++j){
            const_mat[i][j] = 0.0;
        }
        const_rhs[i] = 0.0;
    }

    int minorder= 0;

    irow = 0;
    icol = 0;

    if(constraint == 2){
        std::ifstream ifs_fc2;
        ifs_fc2.open(fc2_file.c_str(), std::ios::in);
        if(!ifs_fc2) error->exit("calc_constraint_matrix", "cannot open file fc2_file");

        bool is_found = false;

        int nparam_harmonic;
        std::string str_tmp;
        while(!ifs_fc2.eof())
        {
            std::getline(ifs_fc2, str_tmp);

            if(str_tmp == "##HARMONIC FORCE CONSTANTS")
            {
                ifs_fc2 >> nparam_harmonic;
                if(nparam_harmonic != nparam[0]) error->exit("calc_constraint_matrix", "Number of fc2 not the same");

                is_found = true;

                for (i = 0; i < nparam[0]; ++i){
                    const_mat[i][i] = 1.0;
                    ifs_fc2 >> const_rhs[i];
                }
                break;
            }
        }
        minorder = 1;
        irow += nparam[0];
        icol += nparam[0];
        ifs_fc2.close();
        if(!is_found) error->exit("calc_constraint_matrix", "HARMONIC FORCE CONSTANTS flag not found in the fc2_file");
    }

    for(order = minorder; order < maxorder; ++order){

        for(std::vector<Constraint>::iterator it = const_vec[order].begin(); it != const_vec[order].end(); ++it){
            Constraint const_now = *it;

            for(i = 0; i < nparam[order]; ++i){
                const_mat[irow][i + icol] = const_now.w_const[i];
            }
            ++irow;
        }
        icol += nparam[order];
    }

    memory->deallocate(nrank);
    memory->deallocate(const_vec);
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

void Fitting::translational_invariance()
{
    // create constraint matrix arising from translational invariance.
    // might be a little tricky. (hard to follow)

    int i, j;
    int iat, jat, icrd, jcrd;
    int order;
    int maxorder = interaction->maxorder;

    int *ind;
    int *intarr, *intarr_copy;
    int **xyzcomponent;

    int nxyz;
    int natmin = symmetry->natmin;
    int nat = system->nat;

    double *arr_constraint;

    std::vector<int> intlist;
    std::set<FcProperty> list_found;
    std::set<FcProperty>::iterator iter_found;

    std::cout << "Start generating constraints for translational invariance ..." << std::endl;

    memory->allocate(ind, maxorder + 1);

    for (order = 0; order < maxorder; ++order){

        std::cout << std::setw(8) << interaction->str_order[order] << " ...";

        const_translation[order].clear();
        int nparams = fcs->ndup[order].size();

        if(nparams == 0) {
            std::cout << " skipped."<< std::endl;
            continue;
        }

        // make interaction list

        list_found.clear();
        for(std::vector<FcProperty>::iterator p = fcs->fc_set[order].begin(); p != fcs->fc_set[order].end(); ++p){
            FcProperty list_tmp = *p; //using copy constructor
            for (i = 0; i < order + 2; ++i){
                ind[i] = list_tmp.elems[i];
            }
            if(list_found.find(FcProperty(order + 2, list_tmp.coef, ind, list_tmp.mother)) != list_found.end()) {
                error->exit("translational invariance", "Duplicate interaction list found");
            }
            list_found.insert(FcProperty(order + 2, list_tmp.coef, ind, list_tmp.mother));
        }

        // generate xyz component for each order

        nxyz = static_cast<int>(pow(static_cast<long double>(3), order + 1));
        memory->allocate(xyzcomponent, nxyz, order + 1);
        fcs->get_xyzcomponent(order + 1, xyzcomponent);

        memory->allocate(arr_constraint, nparams);
        memory->allocate(intarr, order + 2);
        memory->allocate(intarr_copy, order + 2);

        for(i = 0; i < natmin; ++i){

            iat = symmetry->map_p2s[i][0];

            // generate atom pairs for each order

            if(order == 0){
                for(icrd = 0; icrd < 3; ++icrd){

                    intarr[0] = 3 * iat + icrd;

                    for (jcrd = 0; jcrd < 3; ++jcrd){

                        // Reset the temporary array for next constraint
                        for (j = 0; j < nparams; ++j) arr_constraint[j] = 0.0;

                        for (jat = 0; jat < 3 * nat; jat += 3){
                            intarr[1] = jat + jcrd;

                            iter_found = list_found.find(FcProperty(order + 2, 1.0, intarr, 1));

                            // corresponding fcs found.
                            if(iter_found != list_found.end()){
                                FcProperty arrtmp = *iter_found;
                                arr_constraint[arrtmp.mother] += arrtmp.coef;
                            }
                        }
                        if(!is_allzero(nparams,arr_constraint)){
                            // add to constraint list
                            const_translation[order].insert(Constraint(nparams, arr_constraint));
                        }
                    }
                }
            } else {
                for(j = 0; j < interaction->ninter[i][order]; ++j){
                    intlist.push_back(interaction->intpairs[i][order][j]);
                }
                std::sort(intlist.begin(), intlist.end());

                CombinationWithRepetition<int> g(intlist.begin(), intlist.end(), order);
                do {
                    std::vector<int> data = g.now();

                    intarr[0] = iat;
                    intarr[1] = data[0];
                    for (unsigned int isize = 1; isize < data.size(); ++isize){
                        intarr[isize + 1] = data[isize];
                    }

                    if(!interaction->is_incutoff(order + 1, intarr)) continue;

                    for(int ixyz = 0; ixyz < nxyz; ++ixyz){
                        for(int jcrd = 0; jcrd < 3; ++jcrd){

                            // Reset the temporary array for next constraint
                            for (j = 0; j < nparams; ++j) arr_constraint[j] = 0.0;

                            for(int jat = 0; jat < 3 * nat; jat += 3){
                                intarr[order + 1] = jat / 3;

                                if(!interaction->is_incutoff(order + 2, intarr)) continue;

                                for(j = 0; j < order + 1; ++j)  intarr_copy[j] = 3 * intarr[j] + xyzcomponent[ixyz][j];
                                intarr_copy[order + 1] = jat + jcrd;

                                fcs->sort_tail(order + 2, intarr_copy);

                                iter_found = list_found.find(FcProperty(order + 2, 1.0, intarr_copy, 1));
                                if(iter_found != list_found.end()){
                                    FcProperty arrtmp = *iter_found;
                                    arr_constraint[arrtmp.mother] += arrtmp.coef;                                
                                } 

                            }
                            if(!is_allzero(nparams,arr_constraint)){
                                const_translation[order].insert(Constraint(nparams, arr_constraint));
                            }
                        }
                    }

                } while(g.next());
                intlist.clear();
            }
        }

        memory->deallocate(xyzcomponent);
        memory->deallocate(arr_constraint);
        memory->deallocate(intarr);
        memory->deallocate(intarr_copy);

        remove_redundant_rows(nparams, const_translation[order]);
        std::cout << " done." << std::endl;
    }
    memory->deallocate(ind);

    std::cout << "Finished !" << std::endl << std::endl;
    for(order = 0;  order < maxorder; ++order){
        std::cout << "Number of Constraints for" << std::setw(9) << interaction->str_order[order] << " : " << const_translation[order].size() << std::endl;
    }
    std::cout << std::endl;
}

void Fitting::rotational_invariance()
{

    // A first implemention of complicated constraints
    // 2012.4.26

    std::cout << "Start generating constraint matrix for rotational invariance..." << std::endl;

    std::ofstream ofs_constraint;
    ofs_constraint.open("CONSTRAINT", std::ios::out);

    int i, j;
    int iat, jat;
    int icrd, jcrd;
    int order;
    int maxorder = interaction->maxorder;
    int natmin = symmetry->natmin;
    int nat = system->nat;
    int mu, nu;
    int ixyz, nxyz, nxyz2;
    int mu_lambda, lambda;
    int levi_factor;

    int *ind;
    int **xyzcomponent, **xyzcomponent2;
    int *nparams, nparam_sub;
    int *interaction_index, *interaction_atom;
    int *interaction_tmp;

    double *arr_constraint;
    double *arr_constraint_self;

    std::vector<int> interaction_list, interaction_list_old, interaction_list_now;

    std::set<FcProperty> list_found;
    std::set<FcProperty> list_found_last;
    std::set<FcProperty>::iterator iter_found;

    CombinationWithRepetition<int> g;

    memory->allocate(ind, maxorder + 1);
    memory->allocate(nparams, maxorder);

    for (order = 0; order < maxorder; ++order) {
        const_rotation_self[order].clear();
        const_rotation_cross[order].clear();
    }

    for (order = 0; order < maxorder; ++order){

        nparams[order] = fcs->ndup[order].size();

        if (order == 0) {
            std::cout << "Constraint beteen " << std::setw(8) << "1st-order IFCs (which are zero) and " 
                << std::setw(8) << interaction->str_order[order] << " ..." << std::endl;
            nparam_sub = nparams[order];
        } else {
            std::cout << "Constraint between " << std::setw(8) << interaction->str_order[order - 1] << " and "
                << std::setw(8) << interaction->str_order[order] << " ..." << std::endl;
            nparam_sub = nparams[order] + nparams[order - 1];
        }

        memory->allocate(arr_constraint, nparam_sub);
        memory->allocate(arr_constraint_self, nparams[order]);
        memory->allocate(interaction_atom, order + 2);
        memory->allocate(interaction_index, order + 2);
        memory->allocate(interaction_tmp, order + 2);

        if (order > 0) {
            list_found_last = list_found;
            nxyz = static_cast<int>(pow(static_cast<double>(3), order));
            memory->allocate(xyzcomponent, nxyz, order);
            fcs->get_xyzcomponent(order, xyzcomponent);
        }

        list_found.clear();

        for (std::vector<FcProperty>::iterator p = fcs->fc_set[order].begin(); p != fcs->fc_set[order].end(); ++p){
            FcProperty list_tmp = *p; // using copy constructor
            for (i = 0; i < order + 2; ++i){
                ind[i] = list_tmp.elems[i];
            }
            list_found.insert(FcProperty(order + 2, list_tmp.coef, ind, list_tmp.mother));
        }

        for (i = 0; i < natmin; ++i){

            iat = symmetry->map_p2s[i][0];

            if (order == 0){

                 interaction_list_now.clear();
                for (j = 0; j < interaction->ninter[i][order]; ++j){
                    interaction_list_now.push_back(interaction->intpairs[i][order][j]);
                }
                std::sort(interaction_list_now.begin(), interaction_list_now.end());

                // special treatment for harmonic force constants

                for (icrd = 0; icrd < 3; ++icrd){

                    interaction_index[0] = 3 * iat + icrd;

                    for (mu = 0; mu < 3; ++mu){
                        for (nu = 0; nu < 3; ++nu){

                            if (mu == nu) continue;

                            // clear history

                            for (j = 0; j < nparam_sub; ++j) arr_constraint[j] = 0.0;

                            for (jat = 0; jat < nat; ++jat){

                                interaction_index[1] = 3 * jat + mu;
                                iter_found = list_found.find(FcProperty(order + 2, 1.0, interaction_index, 1));
                                if(iter_found != list_found.end()){
                                    FcProperty arrtmp = *iter_found;
                                    //                                    arr_constraint[arrtmp.mother] += arrtmp.coef * system->x_cartesian[jat][nu];                                
                                    arr_constraint[arrtmp.mother] += arrtmp.coef * interaction->minvec[i][jat][nu];
                                }

                                interaction_index[1] = 3 * jat + nu;
                                iter_found = list_found.find(FcProperty(order + 2, 1.0, interaction_index, 1));
                                if(iter_found != list_found.end()){
                                    FcProperty arrtmp = *iter_found;
                                    //                                    arr_constraint[arrtmp.mother] -= arrtmp.coef * system->x_cartesian[jat][mu];                             
                                    arr_constraint[arrtmp.mother] -= arrtmp.coef * interaction->minvec[i][jat][mu];                             
                                }
                            }

                            if(!is_allzero(nparam_sub,arr_constraint)){
                                // add to constraint list
                                const_rotation_self[order].insert(Constraint(nparam_sub, arr_constraint));
                            }

                        } // nu
                    }  // mu
                }
            } else {

                // constraint between force constants of different orders

                interaction_list_old.clear();
                interaction_list_now.clear();

                for (j = 0; j < interaction->ninter[i][order]; ++j){
                    interaction_list_now.push_back(interaction->intpairs[i][order][j]);
                }
                for (j = 0; j < interaction->ninter[i][order - 1]; ++j){
                    interaction_list_old.push_back(interaction->intpairs[i][order - 1][j]);
                }

                std::sort(interaction_list_now.begin(), interaction_list_now.end());
                std::sort(interaction_list_old.begin(), interaction_list_old.end());

                interaction_atom[0] = iat;

                for (icrd = 0; icrd < 3; ++icrd){

                    interaction_index[0] = 3 * iat + icrd;

                    CombinationWithRepetition<int> g_now(interaction_list_now.begin(), interaction_list_now.end(), order);
                    CombinationWithRepetition<int> g_old(interaction_list_old.begin(), interaction_list_old.end(), order);

                    // from m-th order to (m-1)-th order and
                    // from (m-1)-th order to m-th order
                    // 2-different directions to find all constraints

                    for (unsigned int direction = 0; direction < 2; ++direction){

                        if(direction == 0) {
                            g = g_now;
                        } else {
                            g = g_old;
                        }

                        // loop for the interacting pairs

                        do {
                            std::vector<int> data = g.now();

                            for (unsigned int idata = 0; idata < data.size(); ++idata)  interaction_atom[idata + 1] = data[idata];

                            for (ixyz = 0; ixyz < nxyz; ++ixyz){

                                for (j = 0; j < order + 1; ++j) interaction_index[j + 1] = 3 * interaction_atom[j + 1] + xyzcomponent[ixyz][j];

                                for (mu = 0; mu < 3; ++mu){
                                    for (nu = 0; nu < 3; ++nu){

                                        if (mu == nu) continue;

                                        // Search for a new constraint below

                                        for (j = 0; j < nparam_sub; ++j) arr_constraint[j] = 0.0;

                                        ofs_constraint << "-------------------------------------------------------------" << std::endl;
                                        ofs_constraint << "New m_1 a_1 ... m_N a_N, mu_1 ... mu_N, mu, nu" << std::endl;
                                        for (j = 0; j < order + 1; ++j){
                                            ofs_constraint << std::setw(10) << fcs->easyvizint(interaction_index[j]);
                                        }
                                        ofs_constraint << ", mu = " << mu << ", nu = " << nu << std::endl;
                                        ofs_constraint << "Entering loop for m_N+1, a_N+1" << std::endl;
                                        ofs_constraint << "********************************" << std::endl;

                                        // loop for m_{N+1}, a_{N+1}
                                        for (jat = 0; jat < nat; ++jat){ 

                                            interaction_atom[order + 1] = jat;
                                            if(!interaction->is_incutoff(order + 2, interaction_atom)) continue;

                                            interaction_index[order + 1] = 3 * jat + mu;
                                            for (j = 0; j < order + 2; ++j) interaction_tmp[j] = interaction_index[j];

                                            ofs_constraint << "m_N+1 a_N+1 = " << std::setw(5) << fcs->easyvizint(interaction_index[order + 1]);
                                            fcs->sort_tail(order + 2, interaction_tmp);

                                            iter_found = list_found.find(FcProperty(order + 2, 1.0, interaction_tmp, 1));
                                            if(iter_found != list_found.end()){
                                                FcProperty arrtmp = *iter_found;
                                                ofs_constraint << " --> Index = " << std::setw(10) << nparams[order - 1] + arrtmp.mother << ", Coefficient = " 
                                                    /* << std::setw(15) << arrtmp.coef * (system->x_cartesian[jat][nu] - system->x_cartesian[iat][nu]) << std::endl;
                                                    arr_constraint[nparams[order - 1] + arrtmp.mother] += arrtmp.coef * system->x_cartesian[jat][nu]; */
                                                    << std::setw(15) << arrtmp.coef * interaction->minvec[i][jat][nu] << std::endl;
                                                arr_constraint[nparams[order - 1] + arrtmp.mother] += arrtmp.coef * interaction->minvec[i][jat][nu];
                                            } else {
                                                ofs_constraint <<" --> Not Found!" << std::endl;
                                            }

                                            interaction_index[order + 1] = 3 * jat + nu;
                                            for (j = 0; j < order + 2; ++j) interaction_tmp[j] = interaction_index[j];

                                            ofs_constraint << "m_N+1 a_N+1 = " << std::setw(5) << fcs->easyvizint(interaction_index[order + 1]);
                                            fcs->sort_tail(order + 2, interaction_tmp);

                                            iter_found = list_found.find(FcProperty(order + 2, 1.0, interaction_tmp, 1));
                                            if(iter_found != list_found.end()){
                                                FcProperty arrtmp = *iter_found;
                                                ofs_constraint << " --> Index = " << std::setw(10) << nparams[order - 1] + arrtmp.mother << ", Coefficient = " 
                                                    /*  << std::setw(15) << -arrtmp.coef * (system->x_cartesian[jat][mu] - system->x_cartesian[iat][mu]) << std::endl;
                                                    arr_constraint[nparams[order - 1] + arrtmp.mother] -= arrtmp.coef * system->x_cartesian[jat][mu];*/
                                                    << std::setw(15) << arrtmp.coef * interaction->minvec[i][jat][mu] << std::endl;
                                                arr_constraint[nparams[order - 1] + arrtmp.mother] -= arrtmp.coef * interaction->minvec[i][jat][mu];
                                            } else {
                                                ofs_constraint << " --> Not Found!" << std::endl;
                                            }
                                        }
                                        ofs_constraint << "********************************" << std::endl;
                                        ofs_constraint << "Loop for lambda start" << std::endl;

                                        for (lambda = 0; lambda < order + 1; ++lambda){

                                            ofs_constraint << "lambda = " << lambda << std::endl;

                                            mu_lambda = interaction_index[lambda] % 3;

                                            for (jcrd = 0; jcrd < 3; ++jcrd){

                                                for (j = 0; j < order + 1; ++j) interaction_tmp[j] = interaction_index[j];

                                                interaction_tmp[lambda] = 3 * (interaction_index[lambda] / 3) + jcrd;

                                                levi_factor = 0;

                                                for (j = 0; j < 3; ++j){
                                                    levi_factor += levi_civita(j, mu, nu)*levi_civita(j, mu_lambda, jcrd); 
                                                }

                                                /*ofs_constraint << "#!" << "(" << std::setw(5) << mu << "," << std::setw(5) << mu_lambda << ")"
                                                << "(" << std::setw(5) << nu << "," << std::setw(5) << jcrd << ") - "
                                                << "(" << std::setw(5) << nu << "," << std::setw(5) << mu_lambda << ")"
                                                << "(" << std::setw(5) << mu << "," << std::setw(5) << jcrd << ")" << " = " << levi_factor << std::endl; */

                                                if(levi_factor == 0) continue;

                                                ofs_constraint << "mu_lambda = " << jcrd << ":";
                                                for (j = 0; j < order + 1; ++j){
                                                    ofs_constraint << std::setw(10) << fcs->easyvizint(interaction_tmp[j]);
                                                }

                                                fcs->sort_tail(order + 1, interaction_tmp);

                                                iter_found = list_found_last.find(FcProperty(order + 1, 1.0, interaction_tmp, 1));
                                                if(iter_found != list_found_last.end()){
                                                    FcProperty arrtmp = *iter_found;
                                                    ofs_constraint << " --> Index = " << std::setw(10) << arrtmp.mother << ", Coefficient = " 
                                                        << std::setw(15) << arrtmp.coef * static_cast<double>(levi_factor) << std::endl;
                                                    arr_constraint[arrtmp.mother] += arrtmp.coef * static_cast<double>(levi_factor);                                
                                                } else {
                                                    ofs_constraint << " --> Not Found!" << std::endl;

                                                }
                                            }
                                        }

                                        if(!is_allzero(nparam_sub,arr_constraint)){

                                            // add to constraint list

                                            if(is_allzero(nparam_sub, arr_constraint, nparams[order - 1])){
                                                const_rotation_self[order - 1].insert(Constraint(nparams[order - 1], arr_constraint));                                      
                                            } else if (is_allzero(nparams[order - 1], arr_constraint)) {
                                                const_rotation_self[order].insert(Constraint(nparam_sub, arr_constraint, nparams[order - 1]));
                                            } else {
                                                const_rotation_cross[order].insert(Constraint(nparam_sub, arr_constraint)); 
                                            }

                                            for(j = 0; j < nparam_sub; ++j){
                                                if(std::abs(arr_constraint[j]) > eps10) {
                                                    ofs_constraint << "index = " << std::setw(5) <<  j << ": " 
                                                        << std::setw(5) << arr_constraint[j] << std::endl;
                                                }
                                            }
                                        }

                                    } // nu
                                } // mu

                            } // ixyz

                        } while(g.next());

                    } // direction
                } // icrd
            }

            // additional constraint for the last order

            if (order == maxorder - 1) {

                nxyz2 = static_cast<int>(pow(static_cast<double>(3), order + 1));
                memory->allocate(xyzcomponent2, nxyz2, order + 1);
                fcs->get_xyzcomponent(order + 1, xyzcomponent2);

                CombinationWithRepetition<int> g_now(interaction_list_now.begin(), interaction_list_now.end(), order + 1);
                do {

                    std::vector<int> data = g_now.now();

                    for (unsigned int idata = 0; idata < data.size(); ++idata)  interaction_atom[idata + 1] = data[idata];
                    
                    for (ixyz = 0; ixyz < nxyz2; ++ixyz){

                        for (j = 0; j < order + 1; ++j) interaction_index[j + 1] = 3 * interaction_atom[j + 1] + xyzcomponent2[ixyz][j];
                        for (j = 0; j < order + 2; ++j){
                            ofs_constraint << "#$$ i_1 .. i_N = " << std::setw(5) << fcs->easyvizint(interaction_index[j]);
                        }
                        ofs_constraint << std::endl;

                        for (mu = 0; mu < 3; ++mu){
                            for (nu = 0; nu < 3; ++nu){

                                ofs_constraint << "#$$ mu = " << mu << " nu = " << nu << std::endl;

                                if (mu == nu) continue;

                                for (j = 0; j < nparams[order]; ++j) arr_constraint_self[j] = 0.0;

                                for (lambda = 0; lambda < order + 2; ++lambda){

                                    mu_lambda = interaction_index[lambda] % 3;

                                    for (jcrd = 0; jcrd < 3; ++jcrd){

                                        for (j = 0; j < order + 2; ++j) interaction_tmp[j] = interaction_index[j];

                                        interaction_tmp[lambda] = 3 * (interaction_index[lambda] / 3) + jcrd;

                                        levi_factor = 0;

                                        for (j = 0; j < 3; ++j){
                                            levi_factor += levi_civita(j, mu, nu)*levi_civita(j, mu_lambda, jcrd); 
                                        }

                                        if(levi_factor == 0) continue;

                                        fcs->sort_tail(order + 2, interaction_tmp);

                                        iter_found = list_found.find(FcProperty(order + 2, 1.0, interaction_tmp, 1));
                                        if(iter_found != list_found.end()){
                                            FcProperty arrtmp = *iter_found;
                                            arr_constraint_self[arrtmp.mother] += arrtmp.coef * static_cast<double>(levi_factor);                                
                                        }
                                    } // jcrd
                                } // lambda

                                if(!is_allzero(nparams[order], arr_constraint_self)){
                                    const_rotation_self[order].insert(Constraint(nparams[order], arr_constraint_self));         
                                }

                            } // nu
                        } // mu

                    } // ixyz


                } while(g_now.next());

                memory->deallocate(xyzcomponent2);
            }

        }

        std::cout << " done" << std::endl;
        memory->deallocate(arr_constraint);
        memory->deallocate(arr_constraint_self);
        memory->deallocate(interaction_tmp);
        memory->deallocate(interaction_index);
        memory->deallocate(interaction_atom);
    }

    for (order = 0; order < maxorder; ++order) {
        remove_redundant_rows(nparam_sub, const_rotation_cross[order]);
        remove_redundant_rows(nparams[order], const_rotation_self[order]);

     /*   if (order == 0) {
            nparam_sub = nparams[order];
        } else {
            nparam_sub = nparams[order] + nparams[order - 1];
        }
        std::cout << const_rotation_cross[order].size() << std::endl;
        for(std::set<Constraint>::iterator p = const_rotation_cross[order].begin(); p != const_rotation_cross[order].end(); ++p){
            Constraint const_tmp = *p;
            for (j = 0; j < nparam_sub; ++j){
                std::cout << std::setw(15) << std::scientific << const_tmp.w_const[j];
            }
            std::cout << std::endl;
        }

        std::cout << const_rotation_self[order].size() << std::endl;
        for(std::set<Constraint>::iterator p = const_rotation_self[order].begin(); p != const_rotation_self[order].end(); ++p){
            Constraint const_tmp = *p;
            for (j = 0; j < nparams[order]; ++j){
                std::cout << std::setw(15) << std::scientific << const_tmp.w_const[j];
            }
            std::cout << std::endl;
        }
        */
    }
    memory->deallocate(ind);
    memory->deallocate(nparams);
}

void Fitting::remove_redundant_rows(const int n, std::set<Constraint> &Constraint_Set)
{
    using namespace Eigen;

    int nrow = n;
    int ncol = Constraint_Set.size();
    double *arr_tmp;

    if(ncol > 0) {
        memory->allocate(arr_tmp, nrow);
        MatrixXd mat_tmp(nrow, ncol);

        int icol = 0;

        for (std::set<Constraint>::iterator p = Constraint_Set.begin(); p != Constraint_Set.end(); ++p){
            Constraint const_now = *p;
            for (int i = 0; i < nrow; ++i){
                mat_tmp(i, icol) = const_now.w_const[i];
            }
            ++icol;
        }

        FullPivLU<MatrixXd> lu_decomp(mat_tmp);
        int nrank = lu_decomp.rank();
        MatrixXd c_reduced = lu_decomp.image(mat_tmp);

        Constraint_Set.clear();

        for(icol = 0; icol < nrank; ++icol){
            for(int irow = 0; irow < nrow; ++irow){
                arr_tmp[irow] = c_reduced(irow, icol);
            }

            Constraint_Set.insert(Constraint(nrow, arr_tmp));
        }

        memory->deallocate(arr_tmp);
    }
}

bool Fitting::is_allzero(const int n, const double *arr, const int nshift){

    for(int i = nshift; i < n; ++i){
        if(std::abs(arr[i]) > eps10) {
            return false;
        }
    }
    return true;
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

int Fitting::levi_civita(const int i, const int j, const int k)
{
    int epsilon = (j - i) * (k - i) * (k - j) / 2;
    return epsilon;
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

