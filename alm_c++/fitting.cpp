#include "fitting.h"
#include "files.h"
#include "error.h"
#include "memory.h"
#include "symmetry.h"
#include "system.h"
#include "fcs.h"
#include "interaction.h"
#include <iostream>
#include <string>
#include <vector>
#include <boost/lexical_cast.hpp>
#include "timer.h"

using namespace ALM_NS;

Fitting::Fitting(ALM *alm): Pointers(alm){}

Fitting::~Fitting() {}

void Fitting::fitmain()
{
    files->ifs_disp_sym.open(files->file_disp_sym, std::ios::in | std::ios::binary);
    if(!files->ifs_disp_sym) error->exit("fitmain", "cannot open file disp_sym");

    files->ifs_force_sym.open(files->file_force_sym, std::ios::in | std::ios::binary);
    if(!files->ifs_force_sym) error->exit("fitmain", "cannot open file force_sym");


    int ntran = symmetry->ntran;
    int nat = system->nat;
    int natmin = symmetry->natmin;

    int i;
    int ndata = system->ndata;

    int N, M;
    int maxorder = interaction->maxorder;

    N = 0;
    for(i = 0; i < maxorder; ++i){
        N += fcs->ndup[i].size();
    }
    std::cout << "Number of Parameters : " << N << std::endl;
    M = 3 * natmin * ndata * ntran;

    memory->allocate(amat, M, N);
    fsum = new double [M];

    // calculate matrix elements for fitting
    calc_matrix_elements(M, N, nat, natmin, ntran, ndata, maxorder);
    timer->print_elapsed();

    // fitting with singular value decomposition
    if(constraint == 0) {
        fit_without_constraints(M, N);
    }

    // write force constants to file

    wrtfcs(fsum);

    memory->deallocate(amat);
    delete [] fsum;

    timer->print_elapsed();
}

void Fitting::fit_without_constraints(int M, int N)
{
    int i, j, k;
    int nrhs = 1, nrank, INFO, LWORK;
    double *WORK, *S, *amat_mod;
    double rcond = -1.0;
    double f_square = 0.0;

    std::cout << "Entering Fitting Routine: SVD without constraints" << std::endl << std::endl;

    LWORK = 3 * std::min<int>(M, N) + std::max<int>(2*std::min<int>(M, N), std::max<int>(M, N));
    LWORK = 2 * LWORK;

    WORK = new double [LWORK];
    S = new double [N];

    // transpose matrix A
    amat_mod = new double [M * N];
    k = 0;
    for (j = 0; j < N; ++j){
        for(i = 0; i < M; ++i){
            amat_mod[k++] = amat[i][j];
        }
    }
    for (i = 0; i < M; ++i){
        f_square += std::pow(fsum[i], 2);
    }

    // fitting with singular value decomposition
    dgelss_(&M, &N, &nrhs, amat_mod, &M, fsum, &M, S, &rcond, &nrank, WORK, &LWORK, &INFO);

    std::cout << "Finished !" << std::endl << std::endl;

    std::cout << "RANK of the MATRIX: " << nrank << std::endl;
    if(nrank < N) error->warn("fit_without_constraints", 
        "Matrix is rank-deficient. Force constants could not be determined uniquely :(");

    if(nrank == N) {
        double f_residual = 0.0;
        for (i = N; i < M; ++i){
            f_residual += std::pow(fsum[i], 2);
        }
        std::cout << std::endl << "Residual sum of squares for the solution: " << sqrt(f_residual) << std::endl;
        std::cout << "Fitting Error (%) : "<< sqrt(f_residual/f_square) * 100.0 << std::endl;
    }

    delete [] amat_mod;
    delete [] WORK;
    delete [] S;
}

void Fitting::calc_matrix_elements(const int M, const int N, const int nat, const int natmin,
    const int ntran, const int ndata, const int maxorder)
{
    int i, j, k, order;
    int itran, iat;
    double **u;
    double **f;
    int *ind;

    std::cout << "Calculation of Matrix Elements for Direct Fitting Started" << std::endl << std::endl;

    ind = new int [maxorder + 1];

    for (i = 0; i < M; ++i){
        for (j = 0; j < N; ++j){
            amat[i][j] = 0.0;
        }
        fsum[j] = 0.0;
    }

    memory->allocate(u, ntran, 3 * nat);
    memory->allocate(f, ntran, 3 * nat);

    int im = 0;
    int idata = 0;
    double amat_tmp;

    int iparam = 0;

    for(int data = 0; data < ndata; ++data){
        for (itran = 0; itran < ntran; ++itran){

            // read displacement and force

            for(i = 0; i < nat; ++i){
                for(j = 0; j < 3; ++j){
                    files->ifs_disp_sym.read((char *) &u[itran][3 * i + j], sizeof(double));
                    files->ifs_force_sym.read((char *) &f[itran][3 * i + j], sizeof(double));
                }
            }

            // construct r.h.s. vector B

            for (i = 0; i < natmin; ++i){
                iat = symmetry->map_p2s[i][0];
                for (j = 0; j < 3; ++j){
                    fsum[im++] = f[itran][3 * iat + j];
                }
            }

            //construct l.h.s. matrix A

            iparam = 0;

            for(order = 0; order < maxorder; ++order){

                int mm = 0;

                for(std::vector<int>::iterator iter = fcs->ndup[order].begin(); iter != fcs->ndup[order].end(); ++iter){
                    for (i = 0; i < *iter; ++i){
                        ind[0] = fcs->fc_set[order][mm].elems[0];
                        k = idata + inprim_index(fcs->fc_set[order][mm].elems[0]);
                        amat_tmp = 1.0;
                        for(j = 1; j < order + 2; ++j){
                            ind[j] = fcs->fc_set[order][mm].elems[j];
                            amat_tmp *= u[itran][fcs->fc_set[order][mm].elems[j]];
                        }
                        /*for(j = 0; j < order + 2; ++j){
                        std::cout << std::setw(5) << ind[j];
                        }
                        std::cout << std::setw(6) << gamma(order + 2, ind) << std::endl;*/
                        amat[k][iparam] -= gamma(order + 2, ind) * fcs->fc_set[order][mm].coef * amat_tmp;
                        ++mm;
                    }
                    ++iparam;
                }
            }

            idata += 3 * natmin;
        }
    }
    memory->deallocate(u);
    memory->deallocate(f);
    delete [] ind;

    std::cout << "Finished !" << std::endl;
}

int Fitting::inprim_index(const int n)
{

    int atmn = n / 3;
    int crdn = n % 3;

    for (int i = 0; i < symmetry->natmin; ++i){
        if(symmetry->map_p2s[i][0] == atmn){
            return 3 * i + crdn;
        }
    }
    error->exit("inprim_index", "This cannot happen");
}

void Fitting::wrtfcs(const double *params)
{
    int i, j, k, l, m;

    int maxorder = interaction->maxorder;
    std::string *str_fcs; 
    str_fcs = new std::string [maxorder];

    std::string str_tmp;

    for (i = 0; i < maxorder; ++i){
        str_fcs[i] = "*FC" + boost::lexical_cast<std::string>(i + 2);
    }

    files->ofs_fcs <<  "********************Force Constants (FCs)********************" << std::endl;
    files->ofs_fcs <<  "!     Force Constants will be printed in atomic unit        !" << std::endl;
    files->ofs_fcs <<  "!     FC2: Ry/a0^2     FC3: Ry/a0^3     FC4: Ry/a0^4   etc. !" << std::endl;
    files->ofs_fcs <<  "!     FC?: Ry/a0^?                                          !" << std::endl;
    files->ofs_fcs <<  "!     a0= Bohr radius                                       !" << std::endl;
    files->ofs_fcs << "*************************************************************"  << std::endl << std::endl;
    files->ofs_fcs << "---------------Symmetrically Independent FCs---------------" << std::endl;
    files->ofs_fcs << " Global No." << "  Local No." << "            FCs" << "            Pairs" << std::endl;

    k = 0;


    files->ofs_fcs.setf(std::ios::scientific);

    for (i = 0; i < maxorder; ++i){

        m = 0;
        
        if(fcs->ndup[i].size() > 0) {

            files->ofs_fcs << std::endl << std::setw(6) << str_fcs[i] << std::endl;

            for (j = 0; j < fcs->ndup[i].size(); ++j){

                files->ofs_fcs << std::setw(6) << k + 1 << std::setw(6) << j + 1 << std::setw(16) <<  params[k];
                for (l = 0; l < i + 2; ++l){
                    files->ofs_fcs << std::setw(7) << fcs->easyvizint(fcs->fc_set[i][m].elems[l]);    
                }
                files->ofs_fcs << std::endl;
                m += fcs->ndup[i][j];
                ++k;
            }
        }
    }

    for (i = 0; i < maxorder; ++i){
        str_fcs[i] = "**FC" + boost::lexical_cast<std::string>(i + 2);
    }

    files->ofs_fcs << std::endl << std::endl
        << "---------------All FCs below---------------" << std::endl;

    int ip = 0;
    int id;

    for (i = 0; i < maxorder; ++i){

        id = 0;

        if(fcs->ndup[i].size() > 0){
            files->ofs_fcs << std::endl << std::setw(6) << str_fcs[i] << std::endl;

            for (unsigned int iuniq = 0; iuniq < fcs->ndup[i].size(); ++iuniq){

                str_tmp = "# FC" + boost::lexical_cast<std::string>(i + 2) + "_";
                str_tmp += boost::lexical_cast<std::string>(iuniq + 1);

                files->ofs_fcs << str_tmp << std::setw(6) << fcs->ndup[i][iuniq] << std::setw(16) << params[ip] << std::endl;

                for (j = 0; j < fcs->ndup[i][iuniq]; ++j){
                    files->ofs_fcs << std::setw(5) << j + 1 << std::setw(16) << fcs->fc_set[i][id].coef;
                    for (k = 0; k < i + 2; ++k){
                        files->ofs_fcs << std::setw(6) << fcs->easyvizint(fcs->fc_set[i][id].elems[k]);
                    }
                    files->ofs_fcs << std::endl;
                    ++id;
                }
                files->ofs_fcs << std::endl;
                ++ip;
            }
            
        }
    }

    delete [] str_fcs;

    std::cout << "Force Constants are written to file: " << files->file_fcs << std::endl;
}

double Fitting::gamma(const int n, const int *arr)
{
    int *arr_tmp, *nsame;
    int i;
    int ind_front, nsame_to_front;

    arr_tmp = new int [n];
    nsame = new int [n];

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

    delete [] arr_tmp;
    delete [] nsame;

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