#include "system.h"
#include "fcs_phonon.h"
#include "error.h"
#include "memory.h"
#include "../alm_c++/constants.h"
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>

using namespace PHON_NS;

System::System(PHON *phon): Pointers(phon) {}

System::~System() {}

void System::setup()
{
    unsigned int i, j;

    load_system_info();

    recips(lavec_s, rlavec_s);
    recips(lavec_p, rlavec_p);

    memory->allocate(xr_p, nat, 3);
    memory->allocate(xc, nat, 3);

    for (i = 0; i < nat; ++i){
        rotvec(xc[i], xr_s[i], lavec_s);
    }

    for (i = 0; i < nat; ++i){
        rotvec(xr_p[i], xc[i], rlavec_s);
        for(j = 0; j < 3; ++j){
            xr_p[i][j] /=  2.0 * pi;
        }
    }

    std::cout << "**Lattice Vectors**" << std::endl;
    std::cout.setf(std::ios::scientific);

    std::cout << " *Super Cell* " << std::endl << std::endl;
    std::cout << " " << lavec_s[0][0] << " " << lavec_s[1][0] << " " << lavec_s[2][0] << " : a1" << std::endl;
    std::cout << " " << lavec_s[0][1] << " " << lavec_s[1][1] << " " << lavec_s[2][1] << " : a2" << std::endl;
    std::cout << " " << lavec_s[0][2] << " " << lavec_s[1][2] << " " << lavec_s[2][2] << " : a3" << std::endl;
    std::cout << std::endl;

    std::cout << " " << rlavec_s[0][0] << " " << rlavec_s[0][1] << " " << rlavec_s[0][2] << " : b1" << std::endl;
    std::cout << " " << rlavec_s[1][0] << " " << rlavec_s[1][1] << " " << rlavec_s[1][2] << " : b2" << std::endl;
    std::cout << " " << rlavec_s[2][0] << " " << rlavec_s[2][1] << " " << rlavec_s[2][2] << " : b3" << std::endl;
    std::cout << std::endl << std::endl;

    std::cout << " *Primitive Cell* " << std::endl << std::endl;

    std::cout << " *Super Cell* " << std::endl << std::endl;
    std::cout << " " << lavec_p[0][0] << " " << lavec_p[1][0] << " " << lavec_p[2][0] << " : a1" << std::endl;
    std::cout << " " << lavec_p[0][1] << " " << lavec_p[1][1] << " " << lavec_p[2][1] << " : a2" << std::endl;
    std::cout << " " << lavec_p[0][2] << " " << lavec_p[1][2] << " " << lavec_p[2][2] << " : a3" << std::endl;
    std::cout << std::endl;

    std::cout << " " << rlavec_p[0][0] << " " << rlavec_p[0][1] << " " << rlavec_p[0][2] << " : b1" << std::endl;
    std::cout << " " << rlavec_p[1][0] << " " << rlavec_p[1][1] << " " << rlavec_p[1][2] << " : b2" << std::endl;
    std::cout << " " << rlavec_p[2][0] << " " << rlavec_p[2][1] << " " << rlavec_p[2][2] << " : b3" << std::endl;
    std::cout << std::endl << std::endl;

    std::cout << " Number of Atoms: " << nat << std::endl << std::endl;

    // Atomic masses in Rydberg unit

    memory->allocate(mass, nat);
    for (i = 0; i < nat; ++i){
        mass[i] = mass_kd[kd[i]]*amu_ry;
    }
}

void System::load_system_info()
{
    unsigned int i;

    std::string file_fcs = fcs_phonon->file_fcs;
    std::ifstream ifs_fcs;

    ifs_fcs.open(file_fcs.c_str(), std::ios::in);
    if(!ifs_fcs) error->exit("load_system_info", "cannot open file file_fcs");

    std::string str_tmp;

    while(!ifs_fcs.eof())
    {
        std::getline(ifs_fcs, str_tmp);
        if(str_tmp == "##SYSTEM INFO"){

            std::getline(ifs_fcs, str_tmp);

            for (i = 0; i < 3; ++i){
                ifs_fcs >> lavec_s[0][i] >> lavec_s[1][i] >> lavec_s[2][i];
            }
            ifs_fcs.ignore();
            std::getline(ifs_fcs, str_tmp);
            ifs_fcs >> nkd;
            memory->allocate(symbol_kd, nkd);
            memory->allocate(mass_kd, nkd);

            for (i = 0; i < nkd; ++i){
                ifs_fcs >> str_tmp >> symbol_kd[i] >> mass_kd[i];
            }
            ifs_fcs.ignore();
            std::getline(ifs_fcs, str_tmp);
            ifs_fcs >> nat >> natmin >> ntran;

            memory->allocate(xr_s, nat, 3);
            memory->allocate(kd, nat);
            memory->allocate(map_p2s, natmin, ntran);
            memory->allocate(map_s2p, nat);

            unsigned int ikd, itran, icell;
            std::getline(ifs_fcs, str_tmp);
            std::getline(ifs_fcs, str_tmp);
            for (i = 0; i < nat; ++i){
                ifs_fcs >> str_tmp >> ikd >> xr_s[i][0] >> xr_s[i][1] >> xr_s[i][2] >> itran >> icell;
                kd[i] = ikd - 1;
                map_p2s[icell - 1][itran - 1] = i;
                map_s2p[i].atom_num = icell - 1;
                map_s2p[i].tran_num = itran - 1;
            }
        }
    }
    ifs_fcs.close();
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

void System::rotvec(double vec_out[3], double vec_in[3], double mat[3][3], char mode)
{
    // Perform matrix x vector multiplication. 
    //
    // vec_out = mat      * vec_in   (mode = 'N')
    //          (mat)^{t} * vec_in   (mode = 'T')
    //

    unsigned int i;
    double vec_tmp[3];

    for (i = 0; i < 3; ++i){
        vec_tmp[i] = vec_in[i];
    }

    if (mode == 'N') {
        for (i = 0; i < 3; ++i){
            vec_out[i] = mat[i][0] * vec_tmp[0] + mat[i][1] * vec_tmp[1] + mat[i][2] * vec_tmp[2];
        }
    } else if (mode == 'T'){
        for (i = 0; i < 3; ++i){
            vec_out[i] = mat[0][i] * vec_tmp[0] + mat[1][i] * vec_tmp[1] + mat[2][i] * vec_tmp[2];
        }
    } else {
        error->exit("rotvec", "invalid input valiable for mode");
    }
}