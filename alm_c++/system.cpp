#include "system.h"
#include "constants.h"
#include <iostream>

using namespace ALM_NS;

System::System(ALM *alm): Pointers(alm) {}

System::~System() {}

void System::init(){

    recips(lavec, rlavec);


    std::cout.setf(std::ios::scientific);

    std::cout << "Lattice Vector" << std::endl;
    for (int i = 0; i < 3; i++){
        for (int j = 0; j < 3; j++){
            std::cout <<  " " << lavec[i][j] ;
        }
        std::cout << std::endl;
    }

    std::cout << std::endl << "Reciprocal Lattice Vector" << std::endl;    
    for (int i = 0; i < 3; i++){
        for (int j = 0; j < 3; j++){
            std::cout <<  " " << rlavec[i][j] ;
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

void System::recips(double vec[3][3], double inverse[3][3])
{
    double det;
    det = vec[0][0] * vec[1][1] * vec[2][2] 
    + vec[1][0] * vec[2][1] * vec[0][2] 
    + vec[2][0] * vec[0][1] * vec[1][2]
    - vec[0][0] * vec[2][1] * vec[1][2] 
    - vec[2][0] * vec[1][1] * vec[0][2]
    - vec[1][0] * vec[0][1] * vec[2][2];

    if(det < 1.0e-12) {};
    double factor = 2.0 * pi / det;
    inverse[0][0] = (vec[1][1] * vec[2][2] - vec[1][2] * vec[2][1]) * factor;
    inverse[0][1] = (vec[0][2] * vec[2][1] - vec[0][1] * vec[2][2]) * factor;
    inverse[0][2] = (vec[0][1] * vec[1][2] - vec[0][2] * vec[1][1]) * factor;

    inverse[1][0] = (vec[1][2] * vec[2][0] - vec[1][0] * vec[2][2]) * factor;
    inverse[1][1] = (vec[0][0] * vec[2][2] - vec[0][2] * vec[2][0]) * factor;
    inverse[1][2] = (vec[0][2] * vec[1][0] - vec[0][0] * vec[1][2]) * factor;

    inverse[2][0] = (vec[1][0] * vec[2][1] - vec[1][1] * vec[2][0]) * factor;
    inverse[2][1] = (vec[0][1] * vec[2][0] - vec[0][0] * vec[2][1]) * factor;
    inverse[2][2] = (vec[0][0] * vec[1][1] - vec[0][1] * vec[1][0]) * factor;
}