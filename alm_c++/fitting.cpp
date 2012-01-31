#include <iostream>
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
    std::cout << "Total Number of Parameters : " << N << std::endl;
    M = 3 * natmin * ndata * ntran;

    memory->allocate(amat, M, N);
    memory->allocate(fsum, M);

    // calculate matrix elements for fitting
    calc_matrix_elements(M, N, nat, natmin, ntran, ndata, maxorder);
    timer->print_elapsed();

    // fitting with singular value decomposition
    if(constraint == 0) {
        fit_without_constraints(M, N);
    } else {
        fit_with_constraints(M, N);
    }

    // write force constants to file

    wrtfcs(fsum);

    memory->deallocate(amat);
    memory->deallocate(fsum);

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

    memory->allocate(WORK, LWORK);
    memory->allocate(S, N);

    // transpose matrix A
    memory->allocate(amat_mod, M * N);
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

    memory->deallocate(amat_mod);
    memory->deallocate(WORK);
    memory->deallocate(S);
}

void Fitting::fit_with_constraints(int M, int N)
{
    std::cout << "Entering Fitting Routine: QRD with constraints" << std::endl << std::endl;
    translational_invariance();
}


void Fitting::calc_matrix_elements(const int M, const int N, const int nat, const int natmin,
    const int ntran, const int ndata, const int maxorder)
{
    int i, j, k, order;
    int itran, iat;
    double **u;
    double **f;
    int *ind;

    std::cout << "Calculation of Matrix Elements for Direct Fitting Started ...";

    memory->allocate(ind, maxorder + 1);

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
                        amat[k][iparam] -= gamma(order + 2, ind) * fcs->fc_set[order][mm].coef * amat_tmp;
                        ++mm;
                    }
                    ++iparam;
                }
            }

            idata += 3 * natmin;
        }
    }
    memory->deallocate(ind);
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

    //FcProperty list_tmp;

    memory->allocate(ind, maxorder + 1);
    memory->allocate(const_translation, maxorder);
    // const_translation = new std::set<Constraint> [maxorder];

    for (order = 0; order < maxorder; ++order){

        std::cout << std::setw(8) << interaction->str_order[order] << " ...";

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

        // number of parameters

        int nparams = fcs->ndup[order].size();

        const_translation[order].clear();

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

        std::cout << " done." << std::endl;
    }
    memory->deallocate(ind);

    std::cout << "Finished !" << std::endl << std::endl;
    for(order = 0;  order < maxorder; ++order){
        std::cout << "Number of Constraints for" << std::setw(9) << interaction->str_order[order] << " : " << const_translation[order].size() << std::endl;
    }
    std::cout << std::endl;
}


bool Fitting::is_allzero(const int n, const double *arr){

    for(int i = 0; i < n; ++i){
        if(std::abs(arr[i]) > eps10) {
            return false;
        }
    }
    return true;
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

    memory->allocate(str_fcs, maxorder);

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

    memory->deallocate(str_fcs);

    std::cout << "Force Constants are written to file: " << files->file_fcs << std::endl;
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