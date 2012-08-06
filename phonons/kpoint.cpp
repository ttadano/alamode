#include "kpoint.h"
#include "memory.h"
#include "error.h"
#include "system.h"
#include "phonon_dos.h"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace PHON_NS;

Kpoint::Kpoint(PHON *phon): Pointers(phon) {
    npath = 0;
}

Kpoint::~Kpoint() {
    memory->deallocate(xk);
}

void Kpoint::kpoint_setups()
{
    unsigned int i;

    if (phon->mode == "boltzmann") kpoint_mode = 3;

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
            system->rotvec(xk_c[i], xk[i], system->rlavec_p);
        }

        unsigned int j, k;
        unsigned int ik;

        ik = 0;

        for (i = 0; i < npath; ++i){
            for (j = 0; j < nkp[i]; ++j){
                if(j == 0){
                    if(ik == 0) {
                        kaxis[ik] = 0.0;
                    } else {
                        kaxis[ik] = kaxis[ik - 1];
                    }
                } else {
                    for (k = 0; k < 3; ++k){
                        tmp[k] = xk_c[ik][k] - xk_c[ik - 1][k];
                    }
                    kaxis[ik] = kaxis[ik - 1] + std::sqrt(tmp[0] * tmp[0] + tmp[1] * tmp[1] + tmp[2] * tmp[2]);
                }                
                ++ik;
            }
        }

        memory->deallocate(xk_c);
        memory->deallocate(nkp);
        
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
        break;
    case 3:
        std::cout << " kpoint_mode = 3: Uniform k-point grid for Boltzmann" << std::endl;
        std::cin >> nkx >> nky >> nkz;
        nk = nkx * nky * nkz;
        memory->allocate(xk, nk, 3);
        gen_kmesh();
        break;
    default:
        error->exit("read_kpoints", "invalid kpoint_mode = ", kpoint_mode);
    }
    std::cout << " Number of k-points: " << nk << std::endl << std::endl;

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

        system->rotvec(xk_direction, xk_direction, system->rlavec_p);
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

    for(ix = 0; ix < nkx; ++ix){
        for(iy = 0; iy < nky; ++iy){
            for(iz = 0; iz < nkz; ++iz){
                ik = iz + iy * nkz + ix * nkz * nky;
                xkr[ik][0] = static_cast<double>(ix) / static_cast<double>(nkx);
                xkr[ik][1] = static_cast<double>(iy) / static_cast<double>(nky);
                xkr[ik][2] = static_cast<double>(iz) / static_cast<double>(nkz);
            }
        }
    }

    for(ik = 0; ik < nk; ++ik){
        for(i = 0; i < 3; ++i){
            if(xkr[ik][i] >= 0.5){
                xk_tmp = 1.0;
            } else {
                xk_tmp = 0.0;
            }
            xk[ik][i] =  xk_tmp - xkr[ik][i]; 
        }
    }
    memory->deallocate(xkr);
}
