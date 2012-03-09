#include "phonon_dos.h"
#include "kpoint.h"
#include "../alm_c++/constants.h"
#include "error.h"
#include "memory.h"
#include "dynamical.h"
#include <algorithm>
#include <vector>

using namespace PHON_NS;

Dos::Dos(PHON *phon): Pointers(phon){}

Dos::~Dos(){}

void Dos::setup()
{
    if(kpoint->kpoint_mode == 2) {
        flag_dos = true;
    } else {
        flag_dos = false;
    }

    if(flag_dos && delta_e < eps12) error->exit("dos_setup()", "invalid delta_e");
}

void Dos::calc_dos()
{
    int i;

    int n_energy;
    n_energy = static_cast<int>((emax - emin) / delta_e);

    memory->allocate(energy_dos, n_energy);
    memory->allocate(dos_phonon, n_energy);

    for (i = 0; i < n_energy; ++i){
        energy_dos[i] = emin + delta_e * static_cast<double>(i);
    }

    unsigned int nk = kpoint->nk;
    unsigned int nkx = kpoint->nkx;
    unsigned int nky = kpoint->nky;
    unsigned int nkz = kpoint->nkz;

    unsigned int ntetra = 6 * nk;
    memory->allocate(tetras, ntetra, 4);
    prepare_tetrahedron(nkx, nky, nkz);

    unsigned int neval = dynamical->neval;
    double **eval = dynamical->eval_phonon;

    for (i = 0; i < n_energy; ++i){
        dos_phonon[i] = dos_integration(neval, eval);
    }
}

double Dos::dos_integration(const unsigned int n, double **eval)
{
    double dos_ret = 0.0;

    unsigned int i, j, k;

    std::vector<double> e_tetra;

    for (i = 0; i < ntetra; ++i){
        for (j = 0; j < n; ++j){

            for (k = 0; k < 4; ++k){
            e_tetra.push_back(eval[tetras[i][k]][j]);
            }

            std::sort(e_tetra.begin(), e_tetra.end());
        
        }
    }
    return dos_ret;
}
void Dos::prepare_tetrahedron(const int nk1, const int nk2, const int nk3)
{
    int i, j, k, ii, jj, kk;
    int n1, n2, n3, n4, n5, n6, n7, n8;
    int m;

    int nk23 = nk2 * nk3;

    for (i = 0; i < nk1; ++i){
        for (j = 0; j < nk2; ++j){
            for (k = 0; k < nk3; ++k){

            ii = (i + 1) % nk1 + 1;
            jj = (j + 1) % nk2 + 1;
            kk = (k + 1) % nk3 + 1;

            n1 =  k +  j*nk3 +  i*nk23;
            n2 =  k +  j*nk3 + ii*nk23;
            n3 =  k + jj*nk3 +  i*nk23;
            n4 =  k + jj*nk3 + ii*nk23;
            n5 = kk +  j*nk3 +  i*nk23;
            n6 = kk +  j*nk3 + ii*nk23;
            n7 = kk + jj*nk3 +  i*nk23;
            n8 = kk + jj*nk3 + ii*nk23;

            m = 6 * (k + j*nk3 + i*nk23);

            tetras[m][0] = n1;
            tetras[m][1] = n2;
            tetras[m][2] = n3;
            tetras[m][3] = n6;

            ++m;

            tetras[m][0] = n2;
            tetras[m][1] = n3;
            tetras[m][2] = n4;
            tetras[m][3] = n6;

            ++m;

            tetras[m][0] = n1;
            tetras[m][1] = n3;
            tetras[m][2] = n5;
            tetras[m][3] = n6;

            ++m;

            tetras[m][0] = n3;
            tetras[m][1] = n4;
            tetras[m][2] = n6;
            tetras[m][3] = n8;

            ++m;

            tetras[m][0] = n3;
            tetras[m][1] = n6;
            tetras[m][2] = n7;
            tetras[m][3] = n8;

            ++m;

            tetras[m][0] = n3;
            tetras[m][1] = n5;
            tetras[m][2] = n6;
            tetras[m][3] = n7;
            }
        }
    }
}