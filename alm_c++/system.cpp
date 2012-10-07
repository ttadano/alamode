#include <iostream>
#include <iomanip>
#include <fstream>
#include "system.h"
#include "constants.h"
#include "timer.h"
#include "memory.h"
#include "error.h"
#include "constraint.h"
#include "fcs.h"
#include "symmetry.h"
#include "fitting.h"

using namespace ALM_NS;

System::System(ALM *alm): Pointers(alm) {}

System::~System() {}

void System::init(){

    int i, j;

    recips(lavec, rlavec);

    std::cout.setf(std::ios::scientific);

    std::cout << "Lattice Vector" << std::endl;
    std::cout << " " << lavec[0][0] << " " << lavec[1][0] << " " << lavec[2][0] << " : a1" << std::endl;
    std::cout << " " << lavec[0][1] << " " << lavec[1][1] << " " << lavec[2][1] << " : a2" << std::endl;
    std::cout << " " << lavec[0][2] << " " << lavec[1][2] << " " << lavec[2][2] << " : a3" << std::endl;
    std::cout << std::endl;

    std::cout << std::endl << "Reciprocal Lattice Vector" << std::endl;
    std::cout << " " << rlavec[0][0] << " " << rlavec[0][1] << " " << rlavec[0][2] << " : b1" << std::endl;
    std::cout << " " << rlavec[1][0] << " " << rlavec[1][1] << " " << rlavec[1][2] << " : b2" << std::endl;
    std::cout << " " << rlavec[2][0] << " " << rlavec[2][1] << " " << rlavec[2][2] << " : b3" << std::endl;
    std::cout << std::endl;

    double vec_tmp[3][3];
     for (i = 0; i < 3; ++i){
            for (j = 0; j < 3; ++j){
                vec_tmp[i][j] = lavec[j][i];
            }
        }

    cell_volume = volume(vec_tmp[0], vec_tmp[1], vec_tmp[2]);
    std::cout << " Cell volume = " << cell_volume << " (a.u)^3" << std::endl;

    std::cout << "Atomic positions in fractional coordinate and atomic species" << std::endl;
    for (i = 0; i < nat; ++i) {
        std::cout << std::setw(5) << i + 1;
        std::cout << std::setw(15) << xcoord[i][0];
        std::cout << std::setw(15) << xcoord[i][1];
        std::cout << std::setw(15) << xcoord[i][2];
        std::cout << std::setw(5) << kd[i] << std::endl;
    }
    std::cout << std::endl;
    std::cout << "Number of input data: " << ndata << std::endl;
    std::cout << "nstart = " << nstart << ", nend = " << nend << ", nskip = " << nskip << std::endl;

    std::cout.unsetf(std::ios::scientific);

    // generate cartesian coordinate

    memory->allocate(x_cartesian, nat, 3);

    for (i = 0; i < nat; ++i){
        for (j = 0; j < 3; ++j){
            x_cartesian[i][j] = xcoord[i][j];
        }
    }
    frac2cart(x_cartesian);

    timer->print_elapsed();
}

void System::recips(double aa[3][3], double bb[3][3])
{
    /*
    Calculate Reciprocal Lattice Vectors

    Here, BB is just the inverse matrix of AA (multiplied by factor 2 Pi)

    BB = 2 Pi AA^{-1},
       = t(b1, b2, b3)
                     
         (b11 b12 b13)
       = (b21 b22 b23)
         (b31 b32 b33),

    b1 = t(b11, b12, b13) etc.
    */

    double det;
    det = aa[0][0] * aa[1][1] * aa[2][2] 
    + aa[1][0] * aa[2][1] * aa[0][2] 
    + aa[2][0] * aa[0][1] * aa[1][2]
    - aa[0][0] * aa[2][1] * aa[1][2] 
    - aa[2][0] * aa[1][1] * aa[0][2]
    - aa[1][0] * aa[0][1] * aa[2][2];

    if(det < eps12) {
        error->exit("recips", "Lattice Vector is singular");
    }

    double factor = 2.0 * pi / det;

    bb[0][0] = (aa[1][1] * aa[2][2] - aa[1][2] * aa[2][1]) * factor;
    bb[0][1] = (aa[0][2] * aa[2][1] - aa[0][1] * aa[2][2]) * factor;
    bb[0][2] = (aa[0][1] * aa[1][2] - aa[0][2] * aa[1][1]) * factor;

    bb[1][0] = (aa[1][2] * aa[2][0] - aa[1][0] * aa[2][2]) * factor;
    bb[1][1] = (aa[0][0] * aa[2][2] - aa[0][2] * aa[2][0]) * factor;
    bb[1][2] = (aa[0][2] * aa[1][0] - aa[0][0] * aa[1][2]) * factor;

    bb[2][0] = (aa[1][0] * aa[2][1] - aa[1][1] * aa[2][0]) * factor;
    bb[2][1] = (aa[0][1] * aa[2][0] - aa[0][0] * aa[2][1]) * factor;
    bb[2][2] = (aa[0][0] * aa[1][1] - aa[0][1] * aa[1][0]) * factor;
}

void System::frac2cart(double **xf)
{
    // x_cartesian = A x_fractional

    int i, j;

    double *x_tmp;
    memory->allocate(x_tmp, 3);

    for (i = 0; i < nat; ++i){

        rotvec(x_tmp, xf[i], lavec);
        
        for (j = 0; j < 3; ++j){
            xf[i][j] = x_tmp[j];
        }
    }
    memory->deallocate(x_tmp);
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

void System::load_reference_system()
{
    int i, j;
    int iat, jat;
    int icrd, jcrd;

    std::ifstream ifs_fc2;

    ifs_fc2.open(constraint->fc2_file.c_str(), std::ios::in);
    if(!ifs_fc2) error->exit("calc_constraint_matrix", "cannot open file fc2_file");

    bool is_found_system = false;

    int nparam_harmonic_ref;
    int nparam_harmonic = fcs->ndup[0].size();

    std::string str_tmp;

    while(!ifs_fc2.eof() && !is_found_system)
    {
        std::getline(ifs_fc2, str_tmp);
        if(str_tmp == "##SYSTEM INFO"){

            is_found_system = true;

            std::getline(ifs_fc2, str_tmp);
            for (i = 0; i < 3; ++i){
                ifs_fc2 >> lavec_s[0][i] >> lavec_s[1][i] >> lavec_s[2][i];
            }
            ifs_fc2.ignore();
            std::getline(ifs_fc2, str_tmp);
            ifs_fc2 >> nkd_s;
            ifs_fc2.ignore();
            for (i = 0; i < nkd_s; ++i){
                std::getline(ifs_fc2, str_tmp);
            }
            std::getline(ifs_fc2, str_tmp);
            ifs_fc2.ignore();

            ifs_fc2 >> nat_s >> symmetry->natmin_s >> symmetry->ntran_s;

            if (symmetry->natmin_s != symmetry->natmin) {
                error->exit("load_reference_system", "The number of atoms in the primitive cell is not consistent");
            }

            if (nat_s != nat) {
                std::cout << "The number of atoms in the reference system differs from input." << std::endl;
                std::cout << "Trying to map the related force constants (^o^)" << std::endl << std::endl;
            }

            memory->allocate(xcoord_s, nat_s, 3);
            memory->allocate(kd_s, nat_s);
            memory->allocate(symmetry->map_p2s_s, symmetry->natmin_s, symmetry->ntran_s);
            memory->allocate(symmetry->map_s2p_s, nat_s);

            unsigned int ikd, itran, icell;
            std::getline(ifs_fc2, str_tmp);
            std::getline(ifs_fc2, str_tmp);
            for (i = 0; i < nat_s; ++i){
                ifs_fc2 >> str_tmp >> ikd >> xcoord_s[i][0] >> xcoord_s[i][1] >> xcoord_s[i][2] >> itran >> icell;
                kd_s[i] = ikd;
                symmetry->map_p2s_s[icell - 1][itran - 1] = i;
                symmetry->map_s2p_s[i].atom_num = icell - 1;
                symmetry->map_s2p_s[i].tran_num = itran - 1;
            }
        }
    }
    if(!is_found_system) error->exit("load_reference_system", "SYSTEM INFO flag not found in the fc2_file");

    //
    // Generate Mapping Information (big supercell -> small supercell)
    //

    double *xtmp;
    double *xdiff;
    int **intpair_tmp;

    memory->allocate(xtmp, 3);
    memory->allocate(xdiff, 3);
    memory->allocate(map_ref, nat_s);

    bool map_found;
    double dist;

    for (iat = 0; iat < nat_s; ++iat){
        map_found = false;

        rotvec(xtmp, xcoord_s[iat], lavec_s);
        rotvec(xtmp, xtmp, rlavec);

        for (icrd = 0; icrd < 3; ++icrd) xtmp[icrd] /= 2.0 * pi;

        for (jat = 0; jat < nat; ++jat){
            for (icrd = 0; icrd < 3; ++icrd){
                xdiff[icrd] = xtmp[icrd] - xcoord[jat][icrd];
                xdiff[icrd] = std::fmod(xdiff[icrd], 1.0);
            }
            dist = xdiff[0] * xdiff[0] + xdiff[1] * xdiff[1] + xdiff[2] * xdiff[2];

            if (dist < eps12 && kd_s[iat] == kd[jat]) {
                map_ref[iat] = jat;
                map_found = true;
                break;
            }
        }
        if(!map_found) error->exit("load_reference_system", "Could not find an equivalent atom for atom ", iat + 1);
    }

    memory->deallocate(xtmp);
    memory->deallocate(xdiff);

    ifs_fc2.clear();
    ifs_fc2.seekg(0, std::ios_base::beg);

    bool is_found_fc2 = false;

    while(!ifs_fc2.eof() && !is_found_fc2)
    {
        std::getline(ifs_fc2, str_tmp);
        if(str_tmp == "##HARMONIC FORCE CONSTANTS")
        {
            ifs_fc2 >> nparam_harmonic_ref;
            if(nparam_harmonic_ref < nparam_harmonic) {
                error->exit("load_reference_system", "Reference file doesn't contain necessary fc2. (too few)");
            } else if (nparam_harmonic_ref > nparam_harmonic){
                error->exit("load_reference_system","Reference file contains extra force constants." );
            }

            is_found_fc2 = true;

            memory->allocate(fitting->fc2_ref, nparam_harmonic);
            memory->allocate(intpair_tmp, nparam_harmonic, 2);

            for (i = 0; i < nparam_harmonic; ++i){
                ifs_fc2 >> fitting->fc2_ref[i] >> intpair_tmp[i][0] >> intpair_tmp[i][1];
            }

            std::set<FcProperty> list_found;
            std::set<FcProperty>::iterator iter_found;
            int *ind;
            memory->allocate(ind, 2);

            list_found.clear();
            for(std::vector<FcProperty>::iterator p = fcs->fc_set[0].begin(); p != fcs->fc_set[0].end(); ++p){
                FcProperty list_tmp = *p; // Using copy constructor
                for (i = 0; i < 2; ++i){
                    ind[i] = list_tmp.elems[i];
                }
                list_found.insert(FcProperty(2, list_tmp.coef, ind, list_tmp.mother));
            }

            for (i = 0; i < nparam_harmonic; ++i){
                constraint->const_mat[i][i] = 1.0;
            }

            for (i = 0; i < nparam_harmonic; ++i){

                iter_found = list_found.find(FcProperty(2, 1.0, intpair_tmp[i], 1));
                if(iter_found == list_found.end()) {
                    error->exit("load_reference_system", "Cannot find equivalent force constant, number: ", i + 1);
                }
                FcProperty arrtmp = *iter_found;
                constraint->const_rhs[arrtmp.mother] = fitting->fc2_ref[i];
            }

            memory->deallocate(intpair_tmp);
            memory->deallocate(ind);
            list_found.clear();
        }
    }

    if(!is_found_fc2) error->exit("load_reference_system", "HARMONIC FORCE CONSTANTS flag not found in the fc2_file");
    ifs_fc2.close();
}

double System::volume(double vec1[3], double vec2[3], double vec3[3])
{
    double vol;

    vol = std::abs(vec1[0]*(vec2[1]*vec3[2] - vec2[2]*vec3[1]) 
        + vec1[1]*(vec2[2]*vec3[0] - vec2[0]*vec3[2]) 
        + vec1[2]*(vec2[0]*vec3[1] - vec2[1]*vec3[0]));

    return vol;
}
