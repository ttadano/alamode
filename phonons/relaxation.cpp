#include "relaxation.h"
#include "memory.h"
#include "kpoint.h"
#include "system.h"
#include "dynamical.h"
#include "fcs_phonon.h"
#include "../alm_c++/constants.h"

using namespace PHON_NS;

Relaxation::Relaxation(PHON *phon): Pointers(phon) {}
Relaxation::~Relaxation(){};


void Relaxation::setup_relaxation()
{
    unsigned int nk = kpoint->nk;
    unsigned int nband = dynamical->neval;

    memory->allocate(V3, nk*nband, nk*nband, nk*nband);
}

void Relaxation::calc_V3()
{
    unsigned int nband = dynamical->neval;
    unsigned int natmin = system->natmin;
    unsigned int nat = system->nat;
    unsigned int nk = kpoint->nk;

    unsigned int i, j, k, m;
    unsigned int icrd, jcrd, kcrd;
    unsigned int *atmn, *ind, *ind_mod;
    unsigned int *crdn;

    unsigned int k1, k2, k3;
    unsigned int b1, b2, b3;

    double mass_prod;

    std::complex<double> prod, prod_tmp;

    memory->allocate(atmn, 3);
    memory->allocate(crdn, 3);
    memory->allocate(ind, 3);
    memory->allocate(ind_mod, 3);

    std::set<FcsClass>::iterator location;
    std::vector<FcsClass>::iterator it;



    for (k1 = 0; k1 < nk; ++k1){
        for (b1 = 0; b1 < nband; ++b1){
            for (k2 = k1; k2 < nk; ++k2){
                for (b2 = 0; b2 < nband; ++b2){
                    for (k3 = k2; k3 < nk; ++k3){
                        for (b3 = 0; b3 < nband; ++b3){



                            /*
                            prod = (0.0, 0.0);

                            for (i = 0; i < natmin; ++i){

                            atmn[0] = system->map_p2s[i][0];

                            for (j = 0; j < nat; ++j){
                            for (k = 0; k < nat; ++k){

                            mass_prod = 1.0 / std::sqrt(system->mass[atmn[0]] * system->mass[j] * system->mass[k]);

                            prod_tmp = (0.0, 0.0);
                            for (icrd = 0; icrd < 3; ++icrd){
                            for (jcrd = 0; jcrd < 3; ++jcrd){
                            for (kcrd = 0; kcrd < 3; ++kcrd){

                            ind[0] = 3 * atmn[0] + icrd;
                            ind[1] = 3 * j + jcrd;
                            ind[2] = 3 * k + kcrd;

                            for (m = 0; m < 3; ++m) ind_mod[m] = ind[m];

                            if (ind[1] > ind[2]) std::swap(ind[1], ind[2]);

                            location = fcs_phonon->fcs_set[1].find(FcsClass(3, 1.0, ind));

                            if (location == fcs_phonon->fcs_set[1].end()) continue;

                            prod_tmp += (*location).fcs_val
                            * dynamical->evec_phonon[k1][b1][icrd]
                            * dynamical->evec_phonon[k2][b2][3 * system->map_s2p[j].atom_num + jcrd]
                            * dynamical->evec_phonon[k3][b3][3 * system->map_s2p[k].atom_num + kcrd];

                            }
                            }
                            }


                            if (std::norm(prod_tmp) > eps12) {
                            std::cout << "prod_tmp = " << prod_tmp << std::endl;
                            }
                            prod += mass_prod * prod_tmp;

                            }

                            }
                            }

                            if (std::norm(prod) > eps12){
                            std::cout << "prod = " << prod << std::endl;
                            }
                            */



                            prod = (0.0, 0.0);

                            for (it = fcs_phonon->force_constant[1].begin(); it != fcs_phonon->force_constant[1].end(); ++it){

                                for (m = 0; m < 3; ++m) {
                                    ind[m] = (*it).elems[m];
                                    atmn[m] = ind[m] / 3;
                                    crdn[m] = ind[m] % 3;
                                }

                                mass_prod = 1.0 / std::sqrt(system->mass[atmn[0]] * system->mass[atmn[1]] * system->mass[atmn[2]]);

                                prod += mass_prod * (*it).fcs_val 
                                    * dynamical->evec_phonon[k1][b1][3 * system->map_s2p[atmn[0]].atom_num + crdn[0]]
                                * dynamical->evec_phonon[k2][b2][3 * system->map_s2p[atmn[1]].atom_num + crdn[1]]
                                * dynamical->evec_phonon[k3][b3][3 * system->map_s2p[atmn[2]].atom_num + crdn[2]];
                            }

                            if (std::norm(prod) > eps12){

                                std::cout << "k = " << k1 << " " << k2 << " " << k3 << std::endl;
                                std::cout << "band = " << b1 << " " << b2 << " " << b3 << std::endl;
                                std::cout << "prod = " << prod << std::endl;
                            }


                        }
                    }
                }
            }
        }
    }
}