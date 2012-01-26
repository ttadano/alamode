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

using namespace ALM_NS;

Fitting::Fitting(ALM *alm): Pointers(alm){}

Fitting::~Fitting() {}

void Fitting::fitmain()
{
    files->ifs_disp_sym.open(files->file_disp_sym, std::ios::in | std::ios::binary);
    if(!files->ifs_disp_sym) error->exit("fitmain", "cannot open file disp_sym");

    files->ifs_force_sym.open(files->file_force_sym, std::ios::in | std::ios::binary);
    if(!files->ifs_force_sym) error->exit("fitmain", "cannot open file force_sym");

    double **u;
    double **f;
    int ntran = symmetry->ntran;
    int nat = system->nat;
    int natmin = symmetry->natmin;

    int i, j, k, itran;
    int ndata = system->ndata;

    int N, M;
    int maxorder = interaction->maxorder;

    int iat, order;

    N = 0;
    for(i = 0; i < maxorder; ++i){
        N += fcs->ndup[i].size();
    }
    std::cout << "Number of Parameters : " << N << std::endl;
    M = 3 * natmin * ndata * ntran;

    memory->allocate(amat, M, N);
    fsum = new double [M];

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
                    fsum[im++] = f[itran][3 * i + j];
                }
            }

            //construct l.h.s. matrix A

            iparam = 0;

            for(order = 0; order < maxorder; ++order){

                int mm = 0;

                for(std::vector<int>::iterator iter = fcs->ndup[order].begin(); iter != fcs->ndup[order].end(); ++iter){
                    for (i = 0; i < *iter; ++i){
                        /*   p_tmp = *p;
                        k = idata + inprim_index(p_tmp.elems[0]);
                        amat_tmp = 1.0;
                        for(j = 1; j < order + 2; ++j){
                        amat_tmp *= u[itran][p_tmp.elems[j]];
                        }
                        amat[k][iparam] -= p_tmp.coef * amat_tmp;
                        //   ++p;
                        */
                        k = idata + inprim_index(fcs->fc_set[order][mm].elems[0]);
                        amat_tmp = 1.0;
                        for(j = 1; j < order + 2; ++j){
                            amat_tmp *= u[itran][fcs->fc_set[order][mm].elems[j]];
                        }
                        amat[k][iparam] -= fcs->fc_set[order][mm].coef * amat_tmp;
                        ++mm;
                    }
                    ++iparam;
                }
            }

            idata += 3 * natmin;
        }
    }
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