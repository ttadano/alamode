#include "symmetry_core.h"
#include "memory.h"
#include "../alm_c++/constants.h"
#include "error.h"
#include "system.h"
#include <Eigen/Core>
#include <iomanip>
#include <fstream>

using namespace PHON_NS;

Symmetry::Symmetry(PHON *phon): Pointers(phon){
    file_sym = "SYMM_INFO_PRIM";
//    time_reversal_sym = true;
    time_reversal_sym = false;
}
Symmetry::~Symmetry(){}

void Symmetry::setup_symmetry()
{
    unsigned int natmin = system->natmin;
    double **xtmp;
    unsigned int *kdtmp;

    memory->allocate(xtmp, natmin, 3);
    memory->allocate(kdtmp, natmin);

    unsigned int i, j;

    for (i = 0; i < natmin; ++i){
        system->rotvec(xtmp[i], system->xr_s[system->map_p2s[i][0]], system->lavec_s);
        system->rotvec(xtmp[i], xtmp[i], system->rlavec_p);
        
        for (j = 0; j < 3; ++j) xtmp[i][j] /= 2.0 * pi;
        
        kdtmp[i] = system->kd[system->map_p2s[i][0]];
    }

    std::cout << "Atoms in the primitive cell" << std::endl;
    for (i = 0; i < natmin; ++i){
    std::cout << std::setw(6) << i + 1 << ":";
    for (j = 0; j < 3; ++j) {
        std::cout << std::setw(15) << xtmp[i][j];
    }
    std::cout << std::setw(5) << kdtmp[i] << std::endl;
    }
    std::cout << std::endl;

    SymmList.clear();
    gensym(natmin, nsym, nnp, system->lavec_p, system->rlavec_p, xtmp, kdtmp);

    std::cout << "**Symmetry" << std::endl;
    std::cout << "Number of Symmetry Operations = " << nsym << std::endl << std::endl;
}

void Symmetry::gensym(unsigned int nat, unsigned int &nsym, unsigned int nnp, double aa[3][3], double bb[3][3], double **x, unsigned int *kd)
{
    unsigned int i, j;

    std::ofstream ofs_sym;    
    std::ifstream ifs_sym;

    if(nsym == 0) {

        // Automatically find symmetries.

        std::cout << "Generating Symmetry Operations: This can take a while." << std::endl << std::endl;

        // findsym(nat, nnp, kd, aa, bb, x, nsym, rot, tran_int);
        findsym(nat, nnp, kd, aa, bb, x);
        nsym = SymmList.size();

        ofs_sym.open(file_sym.c_str(), std::ios::out);
        if (!ofs_sym) error->exit("gensym", "cannot open file_sym");

        ofs_sym << nsym << std::endl;
        ofs_sym << nnp << std::endl;

        for(std::vector<SymmetryOperation>::iterator p = SymmList.begin(); p != SymmList.end(); ++p){
            SymmetryOperation sym_tmp = *p;
            for (i = 0; i < 9; ++i){
                ofs_sym << std::setw(4) << sym_tmp.symop[i];
            }
            ofs_sym << std::setw(7) << " ";
            for(i = 9; i < 12; ++i){
                ofs_sym << sym_tmp.symop[i] << std::setw(4);
            }
            ofs_sym << std::endl;
        }

        ofs_sym.close();
    } 
    else if(nsym == 1) {

        // Identity operation only !

        int rot_tmp[3][3], tran_tmp[3];

        for (i = 0; i < 3; ++i){
            for (j = 0; j < 3; ++j){
                if(i == j) {
                    rot_tmp[i][j] = 1;
                } else {
                    rot_tmp[i][j] = 0;
                }
            }
            tran_tmp[i] = 0;
        }

        SymmList.push_back(SymmetryOperation(rot_tmp, tran_tmp));
    } 
    else {
        unsigned int nsym2;
        int rot_tmp[3][3], tran_tmp[3];

        ifs_sym.open(file_sym.c_str(), std::ios::in);
        ifs_sym >> nsym2 >> nnp;

        if(nsym != nsym2) error->exit("gensym", "nsym in the given file and the input file are not consistent.");

        for (i = 0; i < nsym; ++i) {
            ifs_sym >> rot_tmp[0][0] >> rot_tmp[0][1] >> rot_tmp[0][2]
                    >> rot_tmp[1][0] >> rot_tmp[1][1] >> rot_tmp[1][2] 
                    >> rot_tmp[2][0] >> rot_tmp[2][1] >> rot_tmp[2][2]
                    >> tran_tmp[0] >> tran_tmp[1] >> tran_tmp[2];
  
            SymmList.push_back(SymmetryOperation(rot_tmp, tran_tmp));
        }
        ifs_sym.close();
    }
}

void Symmetry::findsym(unsigned int nat, unsigned int nnp, unsigned int *kd, double aa[3][3], double bb[3][3], double **x)
{
    int i, j;

    int m11, m12, m13, m21, m22, m23, m31, m32, m33;
    unsigned int  np1, np2, np3;
    int det;

    Eigen::Matrix3d amat, bmat;
    Eigen::Matrix3d rot2;

    int rot_tmp[3][3], rot_reciprocal[3][3];
    int tran_tmp[3];

    for (i = 0; i < 3; ++i){
        for (j = 0; j < 3; ++j){

            amat(i,j) = aa[i][j];
            bmat(i,j) = bb[i][j];

            if(i == j) {
                rot_tmp[i][j] = 1;
            } else {
                rot_tmp[i][j] = 0;
            }
        }
        tran_tmp[i] = 0;
    }

    // Add the identity operation to the list

    SymmList.push_back(SymmetryOperation(rot_tmp, tran_tmp));

    int nnps = nnp * nnp * nnp;
    int **arr_trans;

    memory->allocate(arr_trans, nnps, 3);

    int itran = 0;

    for (np1 = 0; np1 < nnp; ++np1){
        for (np2 = 0; np2 < nnp; ++np2){
            for (np3 = 0; np3 < nnp; ++np3){
                arr_trans[itran][0] = np1;
                arr_trans[itran][1] = np2;
                arr_trans[itran][2] = np3;
                ++itran;
            }
        }
    }   

    for (m11 = -1; m11 <= 1; ++m11){
        for (m12 = -1; m12 <= 1; ++m12) {
            for (m13 = -1; m13 <= 1; ++m13){
                for (m21 = -1; m21 <= 1; ++m21){
                    for (m22 = -1; m22 <= 1; ++m22){
                        for (m23 = -1; m23 <= 1; ++m23){
                            for (m31 = -1; m31 <= 1; ++m31){
                                for (m32 = -1; m32 <= 1; ++m32){
                                    for (m33 = -1; m33 <= 1; ++m33){

                                        det = m11 * (m22 * m33 - m32 * m23)
                                            - m21 * (m12 * m33 - m32 * m13)
                                            + m31 * (m12 * m23 - m22 * m13);

                                        if (det != 1 && det != -1) continue;

                                        rot_tmp[0][0] = m11;
                                        rot_tmp[0][1] = m12;
                                        rot_tmp[0][2] = m13;
                                        rot_tmp[1][0] = m21;
                                        rot_tmp[1][1] = m22;
                                        rot_tmp[1][2] = m23;
                                        rot_tmp[2][0] = m31;
                                        rot_tmp[2][1] = m32;
                                        rot_tmp[2][2] = m33;

                                        det = 1 / det;

                                        rot_reciprocal[0][0] = (m22 * m33 - m23 * m32) * det ;
                                        rot_reciprocal[0][1] = (m23 * m31 - m21 * m33) * det ;
                                        rot_reciprocal[0][2] = (m21 * m32 - m22 * m31) * det ;
                                        rot_reciprocal[1][0] = (m32 * m13 - m33 * m12) * det ;
                                        rot_reciprocal[1][1] = (m33 * m11 - m31 * m13) * det ;
                                        rot_reciprocal[1][2] = (m31 * m12 - m32 * m11) * det ;
                                        rot_reciprocal[2][0] = (m12 * m23 - m13 * m22) * det ;
                                        rot_reciprocal[2][1] = (m13 * m21 - m11 * m23) * det ;
                                        rot_reciprocal[2][2] = (m11 * m22 - m12 * m21) * det ;

                                        for (i = 0; i < 3; ++i) {
                                            for (j = 0; j < 3; ++j){
                                                rot2(i,j) = static_cast<double>(rot_reciprocal[i][j]);
                                            }
                                        }
                                        
                                        if(!is_ortho(rot2, amat, bmat)) continue;

#pragma omp parallel for private(np1, np2, np3)
                                        for (itran = 0; itran < nnps; ++itran){

                                            np1 = arr_trans[itran][0];
                                            np2 = arr_trans[itran][1];
                                            np3 = arr_trans[itran][2];

                                            if(m11 == 1 && m12 == 0 && m13 ==0 &&
                                                m21 == 0 && m22 == 1 && m23 == 0 &&
                                                m31 == 0 && m32 == 0 && m33 == 1 &&
                                                np1 == 0 && np2 == 0 && np3 == 0) continue;

                                            if(!is_invariant(rot2, nat, kd, x, arr_trans[itran], nnp)) continue;

                                            // STL containers are not thread-safe
#pragma omp critical
                                            SymmList.push_back(SymmetryOperation(rot_tmp, arr_trans[itran]));
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    std::sort(SymmList.begin() + 1, SymmList.end());
    memory->deallocate(arr_trans);
}

bool Symmetry::is_ortho(Eigen::Matrix3d rot, Eigen::Matrix3d amat, Eigen::Matrix3d bmat)
{
    double pi2 = 2.0 * pi;

    Eigen::Matrix3d sat, unit;

    double tmp;

    sat = rot * amat.transpose();
    unit = (sat.transpose() * (bmat * (bmat.transpose() * sat)));
    unit /= pow(pi2, 2);

    tmp = pow((unit(0,0) - 1.0), 2) + pow((unit(1,1) - 1.0), 2) + pow((unit(2,2) - 1.0), 2)
        + pow(unit(0,1), 2) + pow(unit(0,2), 2)
        + pow(unit(1,0), 2) + pow(unit(1,2), 2)
        + pow(unit(2,0), 2) + pow(unit(2,1), 2);

    if(tmp > eps) {
        return false;
    } else {
        return true;
    }
}

bool Symmetry::is_invariant(Eigen::Matrix3d rot, unsigned int nat, unsigned int *kd, double **x, int tran[3], unsigned int nnp)
{

    int i, j, k, l;
    Eigen::Vector3d wsi, usi, vsi, tmp;

    bool value = true;

    for (i = 0; i < nat; ++i){

        for (j = 0; j < 3; ++j){   
            wsi(j) = x[i][j] - static_cast<double>(tran[j]) / static_cast<double>(nnp);
        }

        usi = rot.transpose() * wsi;

        l = -1;

        for (j = 0; j < nat; ++j){

            if(kd[j] == kd[i]) {

                for (k = 0; k < 3; ++k) { 
                    vsi(k) = x[j][k];
                    tmp(k) = fmod(std::abs(usi(k) - vsi(k)), 1.0); 
                    // Need "std" to specify floating point operation
                    // especially for intel compiler (there was no problem in MSVC)
                    tmp(k) = std::min<double>(tmp(k), 1.0 - tmp(k)) ;
                }
                double diff = tmp.dot(tmp);
                if (diff < eps12) l = j;
            }
        }

        if(l == -1) value = false;

    }
    return value;
}
