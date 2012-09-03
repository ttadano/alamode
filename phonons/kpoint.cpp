#include "kpoint.h"
#include "memory.h"
#include "error.h"
#include "system.h"
#include "phonon_dos.h"
#include "symmetry_core.h"
#include "dynamical.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <set>
#include <numeric>

using namespace PHON_NS;

Kpoint::Kpoint(PHON *phon): Pointers(phon) {
    npath = 0;
}

Kpoint::~Kpoint() {
    memory->deallocate(xk);
}

void Kpoint::kpoint_setups()
{
    unsigned int i, j;

    if (phon->mode == "boltzmann") kpoint_mode = 3;

    std::cout << "**k-points**" << std::endl;

    switch (kpoint_mode){

    case 0:
        std::cout << " kpoint_mode = 0: calculation on given k-points" << std::endl;
        std::cin >> nk;
        memory->allocate(xk, nk, 3);

        for(i = 0; i < nk; ++i){
            std::cin >> xk[i][0] >> xk[i][1] >> xk[i][2];
        };
        break;

    case 1:
        std::cout << " kpoint_mode = 1: band structure calculation" << std::endl;
        std::cin >> npath;
        std::cout << " Number of paths: " << npath << std::endl;
        memory->allocate(kp_symbol, npath, 2);
        memory->allocate(kp_bound, npath, 2, 3);
        memory->allocate(nkp, npath);

        nk = 0;
        for (i = 0; i < npath; ++i){
            std::cin >> kp_symbol[i][0] >> kp_bound[i][0][0] >> kp_bound[i][0][1] >> kp_bound[i][0][2]
            >> kp_symbol[i][1] >> kp_bound[i][1][0] >> kp_bound[i][1][1] >> kp_bound[i][1][2] >> nkp[i];
            nk += nkp[i];
        };
        memory->allocate(xk, nk, 3);
        memory->allocate(kpoint_direction, nk, 3);

        gen_kpoints_band();
        memory->deallocate(kp_bound);

        double **xk_c, tmp[3];
        memory->allocate(kaxis, nk);
        memory->allocate(xk_c, nk, 3);

        for (i = 0; i < nk; ++i){
            system->rotvec(xk_c[i], xk[i], system->rlavec_p, 'T');
        }

        unsigned int j, k;
        unsigned int ik;

        ik = 0;
        std::cout << std::endl << " ---------- kpval at the edges ----------" << std::endl;
        std::cout << " kpval";
        std::cout.unsetf(std::ios::scientific);

        for (i = 0; i < npath; ++i){
            for (j = 0; j < nkp[i]; ++j){
                if (j == 0){
                    if(ik == 0) {
                        kaxis[ik] = 0.0;
                    } else {
                        kaxis[ik] = kaxis[ik - 1];
                    }
                    std::cout << std::setw(10) << kaxis[ik];
                } else {
                    for (k = 0; k < 3; ++k){
                        tmp[k] = xk_c[ik][k] - xk_c[ik - 1][k];
                    }
                    kaxis[ik] = kaxis[ik - 1] + std::sqrt(tmp[0] * tmp[0] + tmp[1] * tmp[1] + tmp[2] * tmp[2]);
                }
                ++ik;
            }
        }

        std::cout << std::setw(10) << kaxis[ik - 1];
        std::cout << std::endl;
        std::cout << " ----------------------------------------" << std::endl;

        memory->deallocate(xk_c);
        memory->deallocate(nkp);

        std::cout.setf(std::ios::scientific);

        break;
    case 2:
        std::cout << " kpoint_mode = 2: DOS calculation" << std::endl;
        std::cin >> nkx >> nky >> nkz;
        nk = nkx * nky * nkz;
        memory->allocate(xk, nk, 3);

        double emin, emax, delta_e;
        std::cin >> emin >> emax >> delta_e;
        dos->emin = emin;
        dos->emax = emax;
        dos->delta_e = delta_e;

        gen_kmesh();
        std::cout << " nkx: " << std::setw(6) << nkx;
        std::cout << " nky: " << std::setw(6) << nky;
        std::cout << " nkz: " << std::setw(6) << nkz;
        std::cout << std::endl;
        std::cout << " Number of k-points: " << nk << std::endl << std::endl;
        break;
    case 3:
        std::cout << " kpoint_mode = 3: Uniform k-point grid for Boltzmann" << std::endl;
        std::cin >> nkx >> nky >> nkz;
        nk = nkx * nky * nkz;
        memory->allocate(xk, nk, 3);
        gen_kmesh();
        memory->allocate(knum_minus, nk);
        gen_nkminus();
        reduce_kpoints();

        std::cout << " nkx: " << std::setw(6) << nkx;
        std::cout << " nky: " << std::setw(6) << nky;
        std::cout << " nkz: " << std::setw(6) << nkz;
        std::cout << std::endl;
        std::cout << " Number of irreducible k-points: " << nk_equiv.size() << std::endl << std::endl;

        std::cout << "#k-points (in unit of reciprocal vector), weights" << std::endl;
        
        ik = 0;
        for (i = 0; i < nk_equiv.size(); ++i){
            std::cout << std::setw(8) << i + 1 << ":";
            for (j = 0; j < 3; ++j){
                std::cout << std::setw(15) << kpIBZ[ik].kval[j];
            }
            std::cout << std::setw(15) << weight_k[i] << std::endl;
            ik += nk_equiv[i];
        }
        std::cout << std::endl;

        break;
    default:
        error->exit("read_kpoints", "invalid kpoint_mode = ", kpoint_mode);
    }


}

void Kpoint::gen_kpoints_band()
{
    unsigned int i, j, k;
    double xk_s[3], xk_e[3];
    double xk_direction[3], norm;
    unsigned int ik = 0;

    for (i = 0; i < npath; ++i){
        for (j = 0; j < 3; ++j){
            xk_s[j] = kp_bound[i][0][j];
            xk_e[j] = kp_bound[i][1][j];
            xk_direction[j] = xk_e[j] - xk_s[j];
        }

        system->rotvec(xk_direction, xk_direction, system->rlavec_p, 'T');
        norm = std::pow(xk_direction[0], 2) + std::pow(xk_direction[1], 2) + std::pow(xk_direction[2], 2);
        norm = std::sqrt(norm);

        for (j = 0; j < 3; ++j) xk_direction[j] /= norm;

        for(j = 0; j < nkp[i]; ++j){
            for(k = 0; k < 3; ++k){
                xk[ik][k] = xk_s[k] + (xk_e[k] - xk_s[k]) * static_cast<double>(j) / static_cast<double>(nkp[i] - 1);
                kpoint_direction[ik][k] = xk_direction[k];
            }
            ++ik;
        }
    }
}

void Kpoint::gen_kmesh()
{
    unsigned int ix, iy, iz, ik;
    unsigned int i;
    double **xkr, xk_tmp;

    memory->allocate(xkr, nk, 3);

    for (ix = 0; ix < nkx; ++ix){
        for (iy = 0; iy < nky; ++iy){
            for (iz = 0; iz < nkz; ++iz){
                ik = iz + iy * nkz + ix * nkz * nky;
                xkr[ik][0] = static_cast<double>(ix) / static_cast<double>(nkx);
                xkr[ik][1] = static_cast<double>(iy) / static_cast<double>(nky);
                xkr[ik][2] = static_cast<double>(iz) / static_cast<double>(nkz);
            }
        }
    }

    for (ik = 0; ik < nk; ++ik){
        for (i = 0; i < 3; ++i){
            if (xkr[ik][i] >= 0.5){
                xk_tmp = 1.0;
            } else {
                xk_tmp = 0.0;
            }
            xk[ik][i] =  xk_tmp - xkr[ik][i]; 
        }
    }
    memory->deallocate(xkr);
}

void Kpoint::gen_nkminus()
{
    unsigned int i;
    unsigned int ik, jk;
    double xk_tmp, norm;
    std::vector<KpointList> ksets;
    std::set<KpointList>::iterator it;
    std::vector<double> ktmp;
    bool found_same;

    ksets.clear();

    for (ik = 0; ik < nk; ++ik){

        ktmp.clear();

        ktmp.push_back(xk[ik][0]);
        ktmp.push_back(xk[ik][1]);
        ktmp.push_back(xk[ik][2]);

        ksets.push_back(KpointList(ik, ktmp));
    }

    for (ik = 0; ik < nk; ++ik){

        found_same = false;

        ktmp.clear();
        ktmp.push_back(-xk[ik][0]);
        ktmp.push_back(-xk[ik][1]);
        ktmp.push_back(-xk[ik][2]);

        for (std::vector<KpointList>::const_iterator it = ksets.begin(); it != ksets.end(); ++it){
            norm = std::pow(std::fmod(ktmp[0]-(*it).kval[0], 1.0), 2) 
                + std::pow(std::fmod(ktmp[1]-(*it).kval[1], 1.0), 2) 
                + std::pow(std::fmod(ktmp[2]-(*it).kval[2], 1.0), 2);

            if (std::sqrt(norm) < eps12) {
                knum_minus[ik] = (*it).knum;
                found_same = true;
                break;
            } 
        }

        if (!found_same) {
            error->exit("gen_nkminus", "k-point of the inverse point is not found");
        }
    }

    ksets.clear();
}

void Kpoint::reduce_kpoints()
{
    unsigned int ik;
    unsigned int i, j;
    std::set<KpointList> ksets;
    std::set<KpointList>::iterator it;
    std::vector<double> ktmp;

    unsigned int *kequiv;
    unsigned int nsame, knum_found;
    double xk_sym[3], srot[3][3];

    kpIBZ.clear();
    nk_equiv.clear();
    weight_k.clear();

    memory->allocate(kequiv, nk);

    ksets.clear();

    for (ik = 0; ik < nk; ++ik){

        ktmp.clear();

        ktmp.push_back(xk[ik][0]);
        ktmp.push_back(xk[ik][1]);
        ktmp.push_back(xk[ik][2]);

        ksets.insert(KpointList(ik, ktmp));

        kequiv[ik] = ik;
    }

    for (ik = 0; ik < nk; ++ik){

        if(kequiv[ik] != ik) continue;

        nsame = 1;

        ktmp.clear();
        ktmp.push_back(xk[ik][0]);
        ktmp.push_back(xk[ik][1]);
        ktmp.push_back(xk[ik][2]);

        kpIBZ.push_back(KpointList(ik, ktmp));

        for (std::vector<SymmetryOperation>::iterator isym = symmetry->SymmList.begin(); isym != symmetry->SymmList.end(); ++isym){
            for (i = 0; i < 3; ++i){
                for (j = 0; j < 3; ++j){
                    srot[i][j] = static_cast<double>((*isym).symop[3 * i + j]);
                }
            }
            system->rotvec(xk_sym, xk[ik], srot, 'T');

            for (i = 0; i < 3; ++i){
                xk_sym[i] = dynamical->fold(xk_sym[i]);
            }    

            ktmp.clear();  
            ktmp.push_back(xk_sym[0]);
            ktmp.push_back(xk_sym[1]);
            ktmp.push_back(xk_sym[2]);

            it  = ksets.find(KpointList(ik, ktmp));

            if (it != ksets.end()){
                knum_found = (*it).knum;

                if(knum_found > ik && kequiv[knum_found] == knum_found) {
                    kequiv[knum_found] = ik;
                    ++nsame;
                    kpIBZ.push_back(KpointList(knum_found, (*it).kval));
                }
            }
        }
        if (nsame > 0) {
            nk_equiv.push_back(nsame);
        }
    }

    for (std::vector<unsigned int>::iterator p = nk_equiv.begin(); p != nk_equiv.end(); ++p){
        weight_k.push_back(static_cast<double>(*p)/static_cast<double>(nk));
    }
    std::cout << std::endl;
}
