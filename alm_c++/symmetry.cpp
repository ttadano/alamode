#include <cmath>
#include <iostream>
#include <iomanip>
#include "symmetry.h"
#include "system.h"
#include "memory.h"
#include "constants.h"
#include "timer.h"
#include "error.h"
#include "files.h"
#include <Eigen/Core>

using namespace ALM_NS;

Symmetry::Symmetry(ALM *alm) : Pointers(alm) 
{
    file_sym = "SYMM_INFO";
    maxsym = 50000;
}

Symmetry::~Symmetry() {
    memory->deallocate(tnons);
    memory->deallocate(symrel_int);
    memory->deallocate(symrel);
    memory->deallocate(map_sym);
    memory->deallocate(map_p2s);
    memory->deallocate(map_s2p);
}

void Symmetry::init()
{
    int nat = system->nat;

    memory->allocate(tnons, maxsym, 3);
    memory->allocate(symrel_int, maxsym, 3, 3);

    gensym(nat, nsym, nnp, system->lavec, system->rlavec, system->xcoord,
        system->kd, symrel_int, tnons);

    std::cout << std::endl << "Number of symmetry operations = " << nsym << std::endl;
    timer->print_elapsed();

    memory->allocate(symrel, nsym, 3, 3);
    symop_in_cart(system->lavec, system->rlavec);

    pure_translations();

    memory->allocate(map_sym, nat, nsym);
    memory->allocate(map_p2s, natmin, ntran);
    memory->allocate(map_s2p, nat);

    genmaps(nat, system->xcoord, map_sym, map_p2s, map_s2p);

    std::cout << std::endl;
    std::cout << "**Cell-Atom Correspondens Below**" << std::endl;
    std::cout << std::setw(6) << "CELL" << " | " << std::setw(5) << "ATOM" << std::endl;

    for (int i = 0; i < ntran; ++i){
        std::cout << std::setw(6) << i + 1 << " | ";
        for (int j = 0; j < natmin; ++j)  {
            std::cout << std::setw(5) << map_p2s[j][i];
            if((j + 1)%5 == 0) {
                std::cout << std::endl << "       | ";
            }
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;

    if(multiply_data) data_multiplier(nat, system->ndata);

    timer->print_elapsed();
}

void Symmetry::gensym(int nat, int &nsym, int nnp,
    double aa[3][3], double bb[3][3], double **x, int *kd, int ***rot, double **tran)
{
    int i, j;

    if(nsym == 0) {
        // Automatically find symmetries.

        int **tran_int;
        memory->allocate(tran_int, maxsym, 3);
        findsym(nat, nnp, kd, aa, bb, x, nsym, rot, tran_int);

        ofs_sym.open(file_sym.c_str(), std::ios::out);
        ofs_sym << nsym << std::endl;
        ofs_sym << nnp << std::endl;

        for (i = 0; i < nsym; ++i){
            ofs_sym << std::setw(4) << rot[i][0][0] << std::setw(4) <<rot[i][0][1] << std::setw(4) << rot[i][0][2]
            << std::setw(4) << rot[i][1][0] << std::setw(4) << rot[i][1][1] << std::setw(4) << rot[i][1][2] 
            << std::setw(4) << rot[i][2][0] << std::setw(4) << rot[i][2][1] << std::setw(4) << rot[i][2][2] 
            << std::setw(7) << tran_int[i][0] << std::setw(4) << tran_int[i][1] << std::setw(4) << tran_int[i][2] << std::endl;
        }

        for(i = 0; i < nsym; ++i){
            for (j = 0; j <3; ++j){
                tran[i][j] = static_cast<double>(tran_int[i][j]) / static_cast<double>(nnp);
            }
        }
        memory->deallocate(tran_int);
        ofs_sym.close();
    } 
    else if(nsym == 1) {
        // Just identity operation.

        for (i = 0; i < 3; ++i){
            for (j = 0; j < 3; ++j){
                if(i == j) {
                    rot[0][i][j] = 1;
                } else {
                    rot[0][i][j] = 0;
                }
            }
            tran[0][i] = 0.0;
        }
    } 
    else {
        int nsym2;
        int **tran_int;
        ifs_sym.open(file_sym.c_str(), std::ios::in);
        ifs_sym >> nsym2 >> nnp;

        if(nsym != nsym2) error->exit("gensym", "nsym in the given file and the input file are not consistent.");
        memory->allocate(tran_int, nsym, 3);

        for (i = 0; i < nsym; ++i) {
            ifs_sym >> rot[i][0][0] >> rot[i][0][1] >> rot[i][0][2]
            >> rot[i][1][0] >> rot[i][1][1] >> rot[i][1][2] 
            >> rot[i][2][0] >> rot[i][2][1] >> rot[i][2][2]
            >> tran_int[i][0] >> tran_int[i][1] >> tran_int[i][2];

            for (j = 0; j < 3; ++j){
                tran[i][j] = static_cast<double>(tran_int[i][j])/static_cast<double>(nnp);
            }
        }

        memory->deallocate(tran_int);
        ifs_sym.close();
    }
}

void Symmetry::findsym(int nat, int nnp, int *kd, double aa[3][3], double bb[3][3],
    double **x, int &nsym, int ***rot, int **tran_int)
{
    // Symmetry Finder (originally from TAPP code)

    int i, j;

    int m11, m12, m13, m21, m22, m23, m31, m32, m33;
    int det, np1, np2, np3;

    Eigen::Matrix3d amat, bmat;
    Eigen::Matrix3d rot2;

    for (i = 0; i < 3; ++i){
        for (j = 0; j < 3; ++j){

            amat(i,j) = aa[i][j];
            bmat(i,j) = bb[i][j];

            if(i == j) {
                rot[0][i][j] = 1;
            } else {
                rot[0][i][j] = 0;
            }
        }
        tran_int[0][i] = 0;
    }

    nsym = 1;

    int rot_tmp[3][3], rot_reciprocal[3][3];
    int tran_tmp[3];

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

                                        rot_reciprocal[0][0] = (m22 * m33 - m23 * m32) / det ;
                                        rot_reciprocal[0][1] = (m23 * m31 - m21 * m33) / det ;
                                        rot_reciprocal[0][2] = (m21 * m32 - m22 * m31) / det ;
                                        rot_reciprocal[1][0] = (m32 * m13 - m33 * m12) / det ;
                                        rot_reciprocal[1][1] = (m33 * m11 - m31 * m13) / det ;
                                        rot_reciprocal[1][2] = (m31 * m12 - m32 * m11) / det ;
                                        rot_reciprocal[2][0] = (m12 * m23 - m13 * m22) / det ;
                                        rot_reciprocal[2][1] = (m13 * m21 - m11 * m23) / det ;
                                        rot_reciprocal[2][2] = (m11 * m22 - m12 * m21) / det ;

                                        for (i = 0; i < 3; ++i) {
                                            for (j = 0; j < 3; ++j){
                                                rot2(i,j) = static_cast<double>(rot_reciprocal[i][j]);
                                            }
                                        }

                                        if(!is_ortho(rot2, amat, bmat)) continue;

                                        for (np1 = 0; np1 < nnp; ++np1){
                                            for (np2 = 0; np2 < nnp; ++np2){
                                                for (np3 = 0; np3 < nnp; ++np3){

                                                    if(m11 == 1 && m12 == 0 && m13 ==0 &&
                                                        m21 == 0 && m22 == 1 && m23 == 0 &&
                                                        m31 == 0 && m32 == 0 && m33 == 1 &&
                                                        np1 == 0 && np2 == 0 && np3 == 0) continue;

                                                    tran_tmp[0] = np1;
                                                    tran_tmp[1] = np2;
                                                    tran_tmp[2] = np3;

                                                    if(!is_invariant(rot2, nat, kd, x, tran_tmp, nnp)) continue;

                                                    ++nsym;

                                                    for (i = 0; i < 3; ++i){
                                                        for (j = 0; j < 3; ++j){
                                                            rot[nsym - 1][i][j] = rot_tmp[i][j];
                                                        }
                                                        tran_int[nsym - 1][i] = tran_tmp[i];
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
            }
        }
    }
}

bool Symmetry::is_ortho(Eigen::Matrix3d rot, Eigen::Matrix3d amat, Eigen::Matrix3d bmat)
{
    double pi2 = 2.0 * pi;

    Eigen::Matrix3d sat, unit;

    double tmp;

    sat = rot * amat.transpose();
    unit = (sat.transpose() * (bmat.transpose() * (bmat * sat)));
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

bool Symmetry::is_invariant(Eigen::Matrix3d rot, int nat, int *kd, double **x, int tran[3], int nnp)
{

    int i, j, k, l;
    Eigen::Vector3d wsi, usi, vsi, tmp;

    bool value = true;

    for (i = 0; i < nat; ++i){

        for (j = 0; j < 3; ++j){   
            wsi(j) = x[i][j] - static_cast<double>(tran[j]) / nnp;
        }

        usi = rot * wsi;

        l = -1;

        for (j = 0; j < nat; ++j){

            if(kd[j] == kd[i]) {

                for (k = 0; k < 3; ++k) { 
                    vsi(k) = x[j][k]; 
                    tmp(k) = fmod(abs(usi(k) - vsi(k)), 1.0);
                    tmp(k) = std::min<double>(tmp(k), 1.0 - tmp(k)) ;
                }

                if (tmp.dot(tmp) < eps12) l = j;
            }
        }

        if(l == -1) value = false;

    }

    return value;
}

void Symmetry::symop_in_cart(double lavec[3][3], double rlavec[3][3])
{
    int i, j;
    Eigen::Matrix3d aa, bb, sym_tmp;
    Eigen::Matrix3d sym_crt;

    for (i = 0; i < 3; ++i){
        for (j = 0; j < 3; ++j){
            aa(i,j) = lavec[i][j];
            bb(i,j) = rlavec[i][j];
        }
    }

    for (int isym = 0; isym < nsym; ++isym) {

        for (i = 0; i < 3; ++i){
            for (j = 0; j < 3; ++j){
                sym_tmp(i,j) = static_cast<double>(symrel_int[isym][i][j]);
            }
        }
        sym_crt = (aa * (sym_tmp * bb.transpose())) / (2.0 * pi);

        for (i = 0; i < 3; ++i){
            for (j = 0; j < 3; ++j){
                symrel[isym][i][j] = sym_crt(i,j);
            }
        }
    }

#ifdef _DEBUG

    std::cout << "Symmetry Operations in Cartesian Coordinate" << std::endl;
    for (int isym = 0; isym < nsym; ++isym){
        for (i = 0; i < 3; ++i){
            for (j = 0; j < 3; ++j){
                std::cout << std::setw(8) << symrel[isym][i][j];    
            }
        }
        std::cout << std::endl;
    }
#endif
}

void Symmetry::pure_translations()
{
    int i;

    ntran = 0;
    for(i = 0; i < nsym; ++i){
        if(symrel_int[i][0][0] == 1 && symrel_int[i][1][1] == 1 && symrel_int[i][2][2] == 1) {
            ++ntran;
        }
    }

    natmin = system->nat / ntran;

    if(ntran > 1) {
        std::cout << "Given system is not primitive cell;" << std::endl;
        std::cout << std::setw(8) <<  ntran << " translation operations exist." << std::endl;
    } else {
        std::cout << "Given system is a primitive cell." << std::endl;
    }
    std::cout << "Each cell contains " << natmin << " atoms" << std::endl;

    symnum_tran = new int[ntran];
    int isym = 0;

    for (i = 0; i < nsym; ++i){
        if(symrel_int[i][0][0] == 1 && symrel_int[i][1][1] == 1 && symrel_int[i][2][2] == 1) {
            symnum_tran[isym++] = i; 
        }
    }
}

void Symmetry::genmaps(int nat, double **x, int **map_sym, int **map_p2s, Maps *map_s2p)
{
    double **xnew;
    int isym, iat, jat;
    int i;
    double tmp[3], dist; 

    memory->allocate(xnew, nat, 3);

    for(iat = 0; iat < nat; ++iat){
        for(isym = 0; isym < nsym; isym++){
            map_sym[iat][isym] = -1;
        }
    }

    for (isym = 0; isym < nsym; ++isym){

        for (iat = 0; iat < nat; ++iat){

            for (i = 0; i < 3; ++i){
                xnew[iat][i] = static_cast<double>(symrel_int[isym][i][0]) * x[iat][0] 
                + static_cast<double>(symrel_int[isym][i][1]) * x[iat][1] 
                + static_cast<double>(symrel_int[isym][i][2]) * x[iat][2] 
                + tnons[isym][i];
            }

            for (jat = 0; jat < nat; ++jat){

                for (i = 0; i < 3; ++i){
                    tmp[i] = fmod(abs(xnew[iat][i] - x[jat][i]), 1.0);
                    tmp[i] = std::min<double>(tmp[i], 1.0 - tmp[i]);
                }

                dist = tmp[0] * tmp[0] + tmp[1] * tmp[1] + tmp[2] * tmp[2];
                if(dist < eps10) {
                    map_sym[iat][isym] = jat;
                    break;
                }
            }
            if (map_sym[iat][isym] == -1) error->exit("genmaps", "cannot find symmetry for operation # ", isym + 1);
        }
    }
    memory->deallocate(xnew);    

    bool *is_checked;
    memory->allocate(is_checked, nat);

    for (i = 0; i < nat; ++i) is_checked[i] = false;

    jat = 0;
    int atomnum_translated;
    for (iat = 0; iat < nat; ++iat){

        if(is_checked[iat]) continue;
        for (i = 0; i < ntran; ++i){
            atomnum_translated = map_sym[iat][symnum_tran[i]];
            map_p2s[jat][i] = atomnum_translated;
            is_checked[atomnum_translated] = true;
        }
        ++jat;
    }

    memory->deallocate(is_checked);

    for (iat = 0; iat < natmin; ++iat){
        for (i  = 0; i < ntran; ++i){
            atomnum_translated = map_p2s[iat][i];
            map_s2p[atomnum_translated].atom_num = iat;
            map_s2p[atomnum_translated].tran_num = i;
        }
    }
}

void Symmetry::data_multiplier(int nat, int ndata)
{
    int i, j, k, itran;
    double **u, **f;
    double ***u_sym, ***f_sym;

    memory->allocate(u, nat, 3);
    memory->allocate(f, nat, 3);

    memory->allocate(u_sym, ntran, nat, 3);
    memory->allocate(f_sym, ntran, nat, 3);

    files->ofs_disp_sym.open(files->file_disp_sym.c_str(), std::ios::out | std::ios::binary);
    files->ofs_force_sym.open(files->file_force_sym.c_str(), std::ios::out | std::ios::binary);

    if(!files->ofs_disp_sym) error->exit("data_multiplier", "cannot open file_disp"); 
    if(!files->ofs_force_sym) error->exit("data_multiplier", "cannot open file_force");

    for (i = 0; i < ndata; ++i){
        for (j = 0; j < nat; ++j){
            files->ifs_disp >> u[j][0] >> u[j][1] >> u[j][2];
            files->ifs_force >> f[j][0] >> f[j][1] >> f[j][2];
        }


        for (itran = 0; itran < ntran; ++itran){
            for (j = 0; j < nat; ++j){
                for (k = 0; k < 3; ++k){
                    u_sym[itran][map_sym[j][symnum_tran[itran]]][k] = u[j][k];
                    f_sym[itran][map_sym[j][symnum_tran[itran]]][k] = f[j][k];
                }
            }
        }

        for (itran = 0; itran < ntran; ++itran){
            for (j = 0; j < nat; ++j){
                files->ofs_disp_sym.write((char *) &u_sym[itran][j][0], sizeof(double));
                files->ofs_disp_sym.write((char *) &u_sym[itran][j][1], sizeof(double));
                files->ofs_disp_sym.write((char *) &u_sym[itran][j][2], sizeof(double));
                files->ofs_force_sym.write((char *) &f_sym[itran][j][0], sizeof(double));
                files->ofs_force_sym.write((char *) &f_sym[itran][j][1], sizeof(double));
                files->ofs_force_sym.write((char *) &f_sym[itran][j][2], sizeof(double));
            }
        }
    }

    files->ofs_disp_sym.close();
    files->ofs_force_sym.close();
    std::cout << "Symmetrycally equivalent displacements and forces data are" << std::endl;
    std::cout << "stored in files: " << files->file_disp_sym << " " << files->file_force_sym << std::endl;

    memory->deallocate(u);
    memory->deallocate(f);
    memory->deallocate(u_sym);
    memory->deallocate(f_sym);
}
