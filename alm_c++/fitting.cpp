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

    int i, j, itran;
    int ndata = system->ndata;

    int N, M;
    int maxorder = interaction->maxorder;

    N = 0;
    for(i = 0; i < maxorder; ++i){
        N += fcs->ndup[i].size();
    }

    M = 3 * natmin * ndata * ntran;

    memory->allocate(amat, M, N);
    fsum = new double [M];

    memory->allocate(u, ntran, 3 * nat);
    memory->allocate(f, ntran, 3 * nat);

    while(!files->ifs_disp_sym.eof()){
        for (itran = 0; itran < ntran; ++itran){
            for(i = 0; i < nat; ++i){
                for(j = 0; j < 3; ++j){
                    files->ifs_disp_sym.read((char *) &u[itran][3 * i + j], sizeof(double));
                    files->ifs_force_sym.read((char *) &f[itran][3 * i + j], sizeof(double));
                }
            }
        }


    }
}