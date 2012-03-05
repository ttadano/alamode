#include "kpoint.h"
#include "memory.h"
#include "error.h"
#include <iostream>

using namespace PHON_NS;

Kpoint::Kpoint(PHON *phon): Pointers(phon) {
    npath = 0;
}

Kpoint::~Kpoint() {}

void Kpoint::read_kpoints()
{
    unsigned int i;

    switch (kpoint_mode){
    case 0:
        std::cout << "kpoint_mode = 0: calculation on given k-points" << std::endl;
        std::cin >> nk;
        memory->allocate(xk, nk, 3);
        for(i = 0; i < nk; ++i){
            std::cin >> xk[i][0] >> xk[i][1] >> xk[i][2];
        };
        break;
    case 1:
        std::cout << "kpoint_mode = 1: band structure calculation" << std::endl;
        std::cin >> npath;
        std::cout << "Number of paths: " << npath << std::endl;
        memory->allocate(kp_symbol, npath, 2);
        memory->allocate(kp_bound, npath, 2, 3);
        memory->allocate(nkp, npath);

        for (i = 0; i < npath; ++i){
            std::cin >> kp_symbol[i][0] >> kp_bound[i][0][0] >> kp_bound[i][0][1] >> kp_bound[i][0][2]
            >> kp_symbol[i][1] >> kp_bound[i][1][0] >> kp_bound[i][1][1] >> kp_bound[i][1][2] >> nkp[i];
        };
        break;
    case 2:
        std::cout << "kpoint_mode = 2: DOS calculation" << std::endl;
        std::cin >> nkx >> nky >> nkz;
        double emin, emax, delta_e;
        std::cin >> emin >> emax >> delta_e;
        break;
    default:
        error->exit("read_kpoints", "invalid kpoint_mode = ", kpoint_mode);
    }
}