#include <iostream>
#include <iomanip>
#include "system.h"
#include "constants.h"
#include "timer.h"
#include "memory.h"
#include "error.h"

using namespace ALM_NS;

System::System(ALM *alm): Pointers(alm) {}

System::~System() {}

void System::init(){

    int i, j;

    recips(lavec, rlavec);

    std::cout.setf(std::ios::scientific);

    std::cout << "Lattice Vector" << std::endl;
    for (i = 0; i < 3; ++i){
        for (j = 0; j < 3; ++j){
            std::cout <<  " " << lavec[i][j] ;
        }
        std::cout << std::endl;
    }

    std::cout << std::endl << "Reciprocal Lattice Vector" << std::endl;    
    for (i = 0; i < 3; ++i){
        for (j = 0; j < 3; ++j){
            std::cout <<  " " << rlavec[i][j] ;
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;

    std::cout << "Atomic positions in fractional coordinate and atomic species" << std::endl;
    for (i = 0; i < nat; ++i) {
        std::cout << std::setw(5) << i + 1;
        std::cout << " " << xcoord[i][0];
        std::cout << " " << xcoord[i][1];
        std::cout << " " << xcoord[i][2];
        std::cout << " " << kd[i] << std::endl;
    }
    std::cout << std::endl;
    std::cout << "Number of input data: " << ndata << std::endl;

    std::cout.unsetf(std::ios::scientific);
    timer->print_elapsed();

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

    if(det < eps12) {
        error->exit("recips", "Lattice Vector is singular");
    }

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

void System::frac2cart(double **xf)
{
    int i, j;

    double **x_tmp;
    memory->allocate(x_tmp, nat, 3);

    for (i = 0; i < nat; ++i){
        for (j = 0; j < 3; ++j){
            x_tmp[i][j] = lavec[j][0] * xf[i][0] + lavec[j][1] * xf[i][1] + lavec[j][2] * xf[i][2];
        }
    }
    for (i = 0; i < nat; ++i){
        for (j = 0; j < 3; ++j){
            xf[i][j] = x_tmp[i][j];   
        }
    }
    memory->deallocate(x_tmp);

}