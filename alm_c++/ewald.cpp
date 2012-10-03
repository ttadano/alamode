#include "ewald.h"
#include "memory.h"
#include "error.h"
#include "system.h"
#include <fstream>

using namespace ALM_NS;

Ewald::Ewald(ALM *alm): Pointers(alm) {}

Ewald::~Ewald() {
    if (is_longrange) memory->deallocate(Q);
}

void Ewald::init()
{
    if (is_longrange) {
        memory->allocate(Q, system->nat);
        load_charge();
    }
}

void Ewald::load_charge()
{
    unsigned int i;
    std::ifstream ifs_long;

    ifs_long.open(file_longrange.c_str(), std::ios::in);
    if (!ifs_long) error->exit("load_charge", "cannot open file_longrange");

    ifs_long >> alpha;

    for (i = 0; i < system->nat; ++i){
        ifs_long >> Q[i];
    }

    ifs_long.close();
}

void Ewald::ewald_force(const int N, double **r, double **f)
{
    int i, j;

    for (i = 0; i < N; ++i) {
        for (j = 0; j < 3; ++j) {
        f[i][j] = 0.0;
        }
    }

    double **f_tmp;

    memory->allocate(f_tmp, N, 3);


    ewald_force_short(N, r, f_tmp);
}

void Ewald::ewald_force_short(const int N, double **r, double **f)
{
    int i, j;
    int iat, jat;
    double r_tmp[3];

    for (i = 0; i < N; ++i){
    for (j = 0; j < 3; ++j){
    f[i][j] = 0.0;
    }
    }

    for (iat = 0; iat < N; ++iat){

        for (jat = 0; jat < N; ++jat){
        
        
        }



        for (j = 0; j < 3; ++j){
       
        }
    }
}