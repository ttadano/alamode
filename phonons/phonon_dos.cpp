#include "phonon_dos.h"
#include "kpoint.h"
#include "../alm_c++/constants.h"
#include "error.h"
#include "memory.h"
#include "dynamical.h"
#include <algorithm>
#include <vector>

using namespace PHON_NS;

Dos::Dos(PHON *phon): Pointers(phon){
    Kayser2_to_Ry = amu_ry * std::pow(time_ry, 2) / std::pow(Hz_to_kayser, 2);
}

Dos::~Dos(){
    if(flag_dos) {
        memory->deallocate(energy_dos);
        memory->deallocate(dos_phonon);
    }
}

void Dos::setup()
{
    if(kpoint->kpoint_mode == 2) {
        flag_dos = true;
    } else {
        flag_dos = false;
    }

    if(flag_dos && delta_e < eps12) error->exit("dos_setup()", "invalid delta_e");

    if(flag_dos)
    {
        n_energy = static_cast<int>((emax - emin) / delta_e);
        memory->allocate(energy_dos, n_energy);
        memory->allocate(dos_phonon, n_energy);
    }
}

void Dos::calc_dos()
{
    int i;

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
        dos_phonon[i] = dos_integration(neval, ntetra, eval, std::pow(energy_dos[i],2)*Kayser2_to_Ry);
    }

    memory->deallocate(tetras);
}

double Dos::dos_integration(const unsigned int neval, const unsigned int ntetra, double **eval, const double e)
{
    double dos_ret = 0.0;

    unsigned int i, j, k;

    std::vector<double> e_tetra;
    double e1, e2, e3, e4;

    double f_ntetra = 1.0 / static_cast<double>(ntetra);

    for (i = 0; i < ntetra; ++i){
        for (j = 0; j < neval; ++j){

            e_tetra.clear();
            for (k = 0; k < 4; ++k){
                e_tetra.push_back(eval[tetras[i][k]][j]);
            }

            std::sort(e_tetra.begin(), e_tetra.end());
            e1 = e_tetra[0];
            e2 = e_tetra[1];
            e3 = e_tetra[2];
            e4 = e_tetra[3];

           /* std::cout << e1 << " " << e2 << " " << e3 << " " << e4 << std::endl;
            std::cout << "e = " << e << std::endl;*/

            if (e3 <= e && e < e4){
                dos_ret += f_ntetra*(3.0*std::pow((e4 - e), 2) / ((e4 - e1)*(e4 - e2)*(e4 - e3)));
            } else if (e2 <= e && e < e3) {
                dos_ret += f_ntetra*(3.0*(e2 - e1) + 6.0*(e - e2) - 3.0*(e4 + e3 - e2 - e1)*std::pow((e - e2), 2) / ((e3 - e2)*(e4 - e2))) / ((e3 - e1)*(e4 - e1));
            } else if (e1 <= e && e < e2) {
                dos_ret += f_ntetra*3.0*std::pow((e - e1), 2) / ((e2 - e1)*(e3 - e1)*(e4 - e1));
            }

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

                ii = (i + 1) % nk1;
                jj = (j + 1) % nk2;
                kk = (k + 1) % nk3;

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