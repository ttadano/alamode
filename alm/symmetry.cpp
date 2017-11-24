/*
 symmetry.cpp

 Copyright (c) 2014, 2015, 2016 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#include <cmath>
#include <iostream>
#include <iomanip>
#include "mathfunctions.h"
#include "symmetry.h"
#include "system.h"
#include "memory.h"
#include "constants.h"
#include "timer.h"
#include "error.h"
#include "interaction.h"
#include "files.h"
#include <vector>
#include <algorithm>

#ifdef _USE_EIGEN
#include <Eigen/Core>
#endif

using namespace ALM_NS;

Symmetry::Symmetry(ALM *alm) : Pointers(alm)
{
}

Symmetry::~Symmetry()
{
    memory->deallocate(map_sym);
    memory->deallocate(map_p2s);
    memory->deallocate(map_s2p);
    memory->deallocate(symnum_tran);
    SymmData.clear();
}

void Symmetry::init()
{
    int i, j;
    int nat = system->nat;

    std::cout << " SYMMETRY" << std::endl;
    std::cout << " ========" << std::endl << std::endl;

    setup_symmetry_operation(nat, nsym,
                             system->lavec, system->rlavec,
                             system->xcoord, system->kd);

    std::cout << "  Number of symmetry operations = " << SymmData.size() << std::endl;

    // int nsym_fc;

    //    if (nsym_fc == nsym) {
    //        std::cout << "  All symmetry operations will be used to" << std::endl;
    //        std::cout << "  reduce the number of force constants." << std::endl;
    //    } else {
    //        std::cout << "  " << nsym_fc << " symmetry operations out of "
    //            << nsym << " will be used to reduce the number of parameters." << std::endl;
    //        std::cout << "  Other " << nsym - nsym_fc
    //            << " symmetry operations will be imposed as constraints." << std::endl;
    //    }
    //    std::cout << std::endl;

    //    int counter = 0;
    //    for (auto it = SymmData.begin(); it != SymmData.end(); ++it) {
    //        std::cout << "Symm. No. : " << std::setw(4) << counter + 1;
    //        std::cout << "( " << (*it).compatible_with_lattice << " " << (*it).compatible_with_cartesian << ")" << std::endl;
    //        for (i = 0; i < 3; ++i) {
    //            for (j = 0; j < 3; ++j) {
    //                std::cout << std::setw(15) << (*it).rotation[i][j];
    //            }
    //            std::cout << "           ";
    //            for (j = 0; j < 3; ++j) {
    //                std::cout << std::setw(15) << (*it).rotation_cart[i][j];
    //            }
    //
    //            std::cout << std::endl;
    //        }
    //        std::cout << std::endl;
    //        ++counter;
    //    }

    pure_translations();

    memory->allocate(map_sym, nat, nsym);
    memory->allocate(map_p2s, nat_prim, ntran);
    memory->allocate(map_s2p, nat);

    genmaps(nat, system->xcoord, map_sym, map_p2s, map_s2p);

    std::cout << std::endl;
    std::cout << "  **Cell-Atom Correspondens Below**" << std::endl;
    std::cout << std::setw(6) << " CELL" << " | " << std::setw(5) << "ATOM" << std::endl;

    for (int i = 0; i < ntran; ++i) {
        std::cout << std::setw(6) << i + 1 << " | ";
        for (int j = 0; j < nat_prim; ++j) {
            std::cout << std::setw(5) << map_p2s[j][i] + 1;
            if ((j + 1) % 5 == 0) {
                std::cout << std::endl << "       | ";
            }
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;

    timer->print_elapsed();
    std::cout << " -------------------------------------------------------------------" << std::endl;
    std::cout << std::endl;
}

void Symmetry::setup_symmetry_operation(int nat,
                                        unsigned int &nsym,
                                        double aa[3][3],
                                        double bb[3][3],
                                        double **x,
                                        int *kd)
{
    int i, j;

    SymmData.clear();

    if (nsym == 0) {

        // Automatically find symmetries.

        std::cout << "  NSYM = 0 : Trying to find symmetry operations." << std::endl;
        std::cout << "             Please be patient. " << std::endl;
        std::cout << "             This can take a while for a large supercell." << std::endl << std::endl;

        findsym(nat, aa, x, SymmData);
        // The order in SymmData changes for each run because it was generated
        // with OpenMP. Therefore, we sort the list here to have the same result. 
        std::sort(SymmData.begin() + 1, SymmData.end());
        nsym = SymmData.size();

        if (is_printsymmetry) {
            std::ofstream ofs_sym;
            std::cout << "  PRINTSYM = 1: Symmetry information will be stored in SYMM_INFO file."
                << std::endl << std::endl;
            ofs_sym.open(file_sym.c_str(), std::ios::out);
            ofs_sym << nsym << std::endl;

            for (auto p = SymmData.begin(); p != SymmData.end(); ++p) {
                for (i = 0; i < 3; ++i) {
                    for (j = 0; j < 3; ++j) {
                        ofs_sym << std::setw(4) << (*p).rotation[i][j];
                    }
                }
                ofs_sym << "  ";
                for (i = 0; i < 3; ++i) {
                    ofs_sym << std::setprecision(15) << std::setw(21) << (*p).tran[i];
                }
                ofs_sym << std::endl;
            }

            ofs_sym.close();
        }

    } else if (nsym == 1) {

        // Identity operation only !

        std::cout << "  NSYM = 1 : Only the identity matrix will be considered."
            << std::endl << std::endl;

        int rot_tmp[3][3];
        double rot_cart_tmp[3][3];
        double tran_tmp[3];

        for (i = 0; i < 3; ++i) {
            for (j = 0; j < 3; ++j) {
                if (i == j) {
                    rot_tmp[i][j] = 1;
                    rot_cart_tmp[i][j] = 1.0;
                } else {
                    rot_tmp[i][j] = 0;
                    rot_cart_tmp[i][j] = 0.0;
                }
            }
            tran_tmp[i] = 0.0;
        }

        SymmData.push_back(SymmetryOperation(rot_tmp,
                                             tran_tmp,
                                             rot_cart_tmp,
                                             true,
                                             true,
                                             true));

    } else {

        std::cout << "  NSYM > 1 : Symmetry operations will be read from SYMM_INFO file"
            << std::endl << std::endl;

        int nsym2;
        int rot_tmp[3][3];
        double rot_cart_tmp[3][3];
        double tran_tmp[3];
        std::ifstream ifs_sym;

        ifs_sym.open(file_sym.c_str(), std::ios::in);
        ifs_sym >> nsym2;

        if (nsym != nsym2)
            error->exit("setup_symmetry_operations",
                        "nsym in the given file and the input file are not consistent.");

        for (i = 0; i < nsym; ++i) {
            ifs_sym
                >> rot_tmp[0][0] >> rot_tmp[0][1] >> rot_tmp[0][2]
                >> rot_tmp[1][0] >> rot_tmp[1][1] >> rot_tmp[1][2]
                >> rot_tmp[2][0] >> rot_tmp[2][1] >> rot_tmp[2][2]
                >> tran_tmp[0] >> tran_tmp[1] >> tran_tmp[2];

            symop_in_cart(rot_cart_tmp, rot_tmp, system->lavec, system->rlavec);
            SymmData.push_back(SymmetryOperation(rot_tmp,
                                                 tran_tmp,
                                                 rot_cart_tmp,
                                                 is_compatible(rot_tmp),
                                                 is_compatible(rot_cart_tmp),
                                                 is_translation(rot_tmp)));
        }
        ifs_sym.close();
    }

#ifdef _DEBUG
    //  print_symmetrized_coordinate(x);
#endif
}

void Symmetry::findsym(int nat,
                       double aa[3][3],
                       double **x,
                       std::vector<SymmetryOperation> &symop_all)
{
    std::vector<RotationMatrix> LatticeSymmList;

    // Generate rotational matrices that don't change the metric tensor
    LatticeSymmList.clear();
    find_lattice_symmetry(aa, LatticeSymmList);

    // Generate all the space group operations with translational vectors
    symop_all.clear();
    find_crystal_symmetry(nat, system->nclassatom, system->atomlist_class, x,
                          LatticeSymmList, symop_all);

    LatticeSymmList.clear();
}

void Symmetry::find_lattice_symmetry(double aa[3][3],
                                     std::vector<RotationMatrix> &LatticeSymmList)
{
    /*
    Find the rotational matrices that leave the metric tensor invariant.

    Metric tensor G = (g)_{ij} = a_{i} * a_{j} is invariant under crystal symmetry operations T,
    i.e. T^{t}GT = G. Since G can be written as G = A^{t}A, the invariance condition is given by
    (AT)^{t}(AT) = G0 (original).
    */

    int i, j, k;
    int m11, m12, m13, m21, m22, m23, m31, m32, m33;

    int nsym_tmp = 0;
    int mat_tmp[3][3];
    double det, res;
    double rot_tmp[3][3];
    double aa_rot[3][3];

    double metric_tensor[3][3];
    double metric_tensor_rot[3][3];


    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            metric_tensor[i][j] = 0.0;
            for (k = 0; k < 3; ++k) {
                metric_tensor[i][j] += aa[k][i] * aa[k][j];
            }
        }
    }

    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            if (i == j) {
                mat_tmp[i][i] = 1;
            } else {
                mat_tmp[i][j] = 0;
            }
        }
    }

    // Identity matrix should be the first entry.
    LatticeSymmList.push_back(mat_tmp);

    for (m11 = -1; m11 <= 1; ++m11) {
        for (m12 = -1; m12 <= 1; ++m12) {
            for (m13 = -1; m13 <= 1; ++m13) {
                for (m21 = -1; m21 <= 1; ++m21) {
                    for (m22 = -1; m22 <= 1; ++m22) {
                        for (m23 = -1; m23 <= 1; ++m23) {
                            for (m31 = -1; m31 <= 1; ++m31) {
                                for (m32 = -1; m32 <= 1; ++m32) {
                                    for (m33 = -1; m33 <= 1; ++m33) {

                                        if (m11 == 1 && m12 == 0 && m13 == 0 &&
                                            m21 == 0 && m22 == 1 && m23 == 0 &&
                                            m31 == 0 && m32 == 0 && m33 == 1)
                                            continue;

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

                                        // Here, aa_rot = aa * rot_tmp is correct.
                                        matmul3(aa_rot, aa, rot_tmp);

                                        for (i = 0; i < 3; ++i) {
                                            for (j = 0; j < 3; ++j) {
                                                metric_tensor_rot[i][j] = 0.0;
                                                for (k = 0; k < 3; ++k) {
                                                    metric_tensor_rot[i][j] += aa_rot[k][i] * aa_rot[k][j];
                                                }
                                            }
                                        }

                                        res = 0.0;
                                        for (i = 0; i < 3; ++i) {
                                            for (j = 0; j < 3; ++j) {
                                                res += std::pow(metric_tensor[i][j] - metric_tensor_rot[i][j], 2.0);
                                            }
                                        }

                                        // Metric tensor is invariant under symmetry operations.

                                        if (res < tolerance * tolerance) {
                                            ++nsym_tmp;
                                            for (i = 0; i < 3; ++i) {
                                                for (j = 0; j < 3; ++j) {
                                                    mat_tmp[i][j] = static_cast<int>(rot_tmp[i][j]);
                                                }
                                            }
                                            LatticeSymmList.push_back(mat_tmp);
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

    if (LatticeSymmList.size() > 48) {
        error->exit("find_lattice_symmetry", "Number of lattice symmetry is larger than 48.");
    }
}

void Symmetry::find_crystal_symmetry(int nat,
                                     int nclass,
                                     std::vector<unsigned int> *atomclass,
                                     double **x,
                                     std::vector<RotationMatrix> LatticeSymmList,
                                     std::vector<SymmetryOperation> &CrystalSymmList)
{
    unsigned int i, j;
    unsigned int iat, jat, kat, lat;
    double x_rot[3];
    double rot[3][3], rot_tmp[3][3], rot_cart[3][3];
    double mag[3], mag_rot[3];
    double tran[3];
    double x_rot_tmp[3];
    double tmp[3];
    double diff;

    int rot_int[3][3];

    int ii, jj, kk;
    unsigned int itype;

    bool is_found;
    bool isok;
    bool mag_sym1, mag_sym2;

    bool is_identity_matrix;


    // Add identity matrix first.
    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            if (i == j) {
                rot_int[i][j] = 1;
                rot_cart[i][j] = 1.0;
            } else {
                rot_int[i][j] = 0;
                rot_cart[i][j] = 0.0;
            }
        }
        tran[i] = 0.0;
    }

    CrystalSymmList.push_back(SymmetryOperation(rot_int,
                                                tran,
                                                rot_cart,
                                                is_compatible(rot_int),
                                                is_compatible(rot_cart),
                                                is_translation(rot_int)));


    for (auto it_latsym = LatticeSymmList.begin(); it_latsym != LatticeSymmList.end(); ++it_latsym) {

        iat = atomclass[0][0];

        for (i = 0; i < 3; ++i) {
            for (j = 0; j < 3; ++j) {
                rot[i][j] = static_cast<double>((*it_latsym).mat[i][j]);
            }
        }

        rotvec(x_rot, x[iat], rot);

#ifdef _OPENMP
#pragma omp parallel for private(jat, tran, isok, kat, x_rot_tmp, is_found, lat, tmp, diff, \
    i, j, itype, jj, kk, is_identity_matrix, mag, mag_rot, rot_tmp, rot_cart, mag_sym1, mag_sym2)
#endif
        for (ii = 0; ii < atomclass[0].size(); ++ii) {
            jat = atomclass[0][ii];

            for (i = 0; i < 3; ++i) {
                tran[i] = x[jat][i] - x_rot[i];
                tran[i] = tran[i] - nint(tran[i]);
            }

            if ((std::abs(tran[0]) > eps12 && !interaction->is_periodic[0]) ||
                (std::abs(tran[1]) > eps12 && !interaction->is_periodic[1]) ||
                (std::abs(tran[2]) > eps12 && !interaction->is_periodic[2]))
                continue;

            is_identity_matrix =
            (std::pow(rot[0][0] - 1.0, 2) + std::pow(rot[0][1], 2) + std::pow(rot[0][2], 2)
                + std::pow(rot[1][0], 2) + std::pow(rot[1][1] - 1.0, 2) + std::pow(rot[1][2], 2)
                + std::pow(rot[2][0], 2) + std::pow(rot[2][1], 2) + std::pow(rot[2][2] - 1.0, 2)
                + std::pow(tran[0], 2) + std::pow(tran[1], 2) + std::pow(tran[2], 2)) < eps12;
            if (is_identity_matrix) continue;

            isok = true;

            for (itype = 0; itype < nclass; ++itype) {

                for (jj = 0; jj < atomclass[itype].size(); ++jj) {

                    kat = atomclass[itype][jj];

                    rotvec(x_rot_tmp, x[kat], rot);

                    for (i = 0; i < 3; ++i) {
                        x_rot_tmp[i] += tran[i];
                    }

                    is_found = false;

                    for (kk = 0; kk < atomclass[itype].size(); ++kk) {

                        lat = atomclass[itype][kk];

                        for (i = 0; i < 3; ++i) {
                            tmp[i] = std::fmod(std::abs(x[lat][i] - x_rot_tmp[i]), 1.0);
                            tmp[i] = std::min<double>(tmp[i], 1.0 - tmp[i]);
                        }
                        diff = tmp[0] * tmp[0] + tmp[1] * tmp[1] + tmp[2] * tmp[2];
                        if (diff < tolerance * tolerance) {
                            is_found = true;
                            break;
                        }
                    }

                    if (!is_found) isok = false;
                }
            }

            if (isok) {
                matmul3(rot_tmp, rot, system->rlavec);
                matmul3(rot_cart, system->lavec, rot_tmp);

                for (i = 0; i < 3; ++i) {
                    for (j = 0; j < 3; ++j) {
                        rot_cart[i][j] /= (2.0 * pi);
                    }
                }

                if (system->lspin && system->noncollinear) {
                    for (i = 0; i < 3; ++i) {
                        mag[i] = system->magmom[jat][i];
                        mag_rot[i] = system->magmom[iat][i];
                    }


                    rotvec(mag_rot, mag_rot, rot_cart);

                    // In the case of improper rotation, the factor -1 should be multiplied
                    // because the inversion operation doesn't flip the spin.
                    if (!is_proper(rot_cart)) {
                        for (i = 0; i < 3; ++i) {
                            mag_rot[i] = -mag_rot[i];
                        }
                    }

                    mag_sym1 = (std::pow(mag[0] - mag_rot[0], 2.0)
                        + std::pow(mag[1] - mag_rot[1], 2.0)
                        + std::pow(mag[2] - mag_rot[2], 2.0)) < eps6;

                    mag_sym2 = (std::pow(mag[0] + mag_rot[0], 2.0)
                        + std::pow(mag[1] + mag_rot[1], 2.0)
                        + std::pow(mag[2] + mag_rot[2], 2.0)) < eps6;

                    if (!mag_sym1 && !mag_sym2) {
                        isok = false;
                    } else if (!mag_sym1 && mag_sym2 && !trev_sym_mag) {
                        isok = false;
                    }
                }
            }


            if (isok) {
#ifdef _OPENMP
#pragma omp critical
#endif
                CrystalSymmList.push_back(SymmetryOperation((*it_latsym).mat,
                                                            tran,
                                                            rot_cart,
                                                            is_compatible((*it_latsym).mat),
                                                            is_compatible(rot_cart),
                                                            is_translation((*it_latsym).mat)));
            }
        }

    }
}


void Symmetry::symop_in_cart(double rot_cart[3][3],
                             const int rot_lattice[3][3],
                             const double lavec[3][3],
                             const double rlavec[3][3])
{
    int i, j;
    double sym_tmp[3][3];
    double tmp[3][3];

    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            sym_tmp[i][j] = static_cast<double>(rot_lattice[i][j]);
        }
    }

    matmul3(tmp, sym_tmp, rlavec);
    matmul3(rot_cart, lavec, tmp);

    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            rot_cart[i][j] = rot_cart[i][j] / (2.0 * pi);
        }
    }
}


void Symmetry::pure_translations()
{
    int i;

    ntran = 0;
    for (i = 0; i < nsym; ++i) {
        if (SymmData[i].is_translation) ++ntran;
    }

    nat_prim = system->nat / ntran;

    if (ntran > 1) {
        std::cout << "  Given system is not primitive cell." << std::endl;
        std::cout << "  There are " << std::setw(5)
            << ntran << " translation operations." << std::endl;
    } else {
        std::cout << "  Given system is a primitive cell." << std::endl;
    }
    std::cout << "  Primitive cell contains " << nat_prim << " atoms" << std::endl;

    if (system->nat % ntran) {
        error->exit("pure_translations",
                    "nat != nat_prim * ntran. Something is wrong in the structure.");
    }

    memory->allocate(symnum_tran, ntran);

    int isym = 0;

    for (i = 0; i < nsym; ++i) {
        if (SymmData[i].is_translation) symnum_tran[isym++] = i;
    }
}

void Symmetry::genmaps(int nat,
                       double **x,
                       int **map_sym,
                       int **map_p2s,
                       Maps *map_s2p)
{
    int isym, iat, jat;
    int i, j;
    int itype;
    int ii, jj;
    double xnew[3];
    double tmp[3], diff;
    double rot_double[3][3];

    for (iat = 0; iat < nat; ++iat) {
        for (isym = 0; isym < nsym; ++isym) {
            map_sym[iat][isym] = -1;
        }
    }

#ifdef _OPENMP
#pragma omp parallel for private(i, j, rot_double, itype, ii, iat, xnew, jj, jat, tmp, diff, isym)
#endif
    for (isym = 0; isym < nsym; ++isym) {

        for (i = 0; i < 3; ++i) {
            for (j = 0; j < 3; ++j) {
                rot_double[i][j] = static_cast<double>(SymmData[isym].rotation[i][j]);
            }
        }

        for (itype = 0; itype < system->nclassatom; ++itype) {

            for (ii = 0; ii < system->atomlist_class[itype].size(); ++ii) {

                iat = system->atomlist_class[itype][ii];

                rotvec(xnew, x[iat], rot_double);

                for (i = 0; i < 3; ++i) xnew[i] += SymmData[isym].tran[i];

                for (jj = 0; jj < system->atomlist_class[itype].size(); ++jj) {

                    jat = system->atomlist_class[itype][jj];

                    for (i = 0; i < 3; ++i) {
                        tmp[i] = std::fmod(std::abs(x[jat][i] - xnew[i]), 1.0);
                        tmp[i] = std::min<double>(tmp[i], 1.0 - tmp[i]);
                    }
                    diff = tmp[0] * tmp[0] + tmp[1] * tmp[1] + tmp[2] * tmp[2];
                    if (diff < tolerance * tolerance) {
                        map_sym[iat][isym] = jat;
                        break;
                    }
                }
                if (map_sym[iat][isym] == -1) {
                    error->exit("genmaps",
                                "cannot find symmetry for operation # ",
                                isym + 1);
                }
            }
        }
    }

    bool *is_checked;
    memory->allocate(is_checked, nat);

    for (i = 0; i < nat; ++i) is_checked[i] = false;

    jat = 0;
    int atomnum_translated;
    for (iat = 0; iat < nat; ++iat) {

        if (is_checked[iat]) continue;
        for (i = 0; i < ntran; ++i) {
            atomnum_translated = map_sym[iat][symnum_tran[i]];
            map_p2s[jat][i] = atomnum_translated;
            is_checked[atomnum_translated] = true;
        }
        ++jat;
    }

    memory->deallocate(is_checked);

    for (iat = 0; iat < nat_prim; ++iat) {
        for (i = 0; i < ntran; ++i) {
            atomnum_translated = map_p2s[iat][i];
            map_s2p[atomnum_translated].atom_num = iat;
            map_s2p[atomnum_translated].tran_num = i;
        }
    }
}

bool Symmetry::is_translation(const int rot[3][3])
{
    bool ret;

    ret =
        rot[0][0] == 1 && rot[0][1] == 0 && rot[0][2] == 0 &&
        rot[1][0] == 0 && rot[1][1] == 1 && rot[1][2] == 0 &&
        rot[2][0] == 0 && rot[2][1] == 0 && rot[2][2] == 1;

    return ret;
}


template <typename T>
bool Symmetry::is_compatible(const T rot[3][3],
                             const double tolerance_zero)
{
    int i, j;
    int nfinite;
    double rot_double[3][3];

    nfinite = 0;
    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            rot_double[i][j] = static_cast<double>(rot[i][j]);
            if (std::abs(rot_double[i][j]) > tolerance_zero) ++nfinite;
        }
    }

    if (nfinite == 3) return true;

    return false;
}

void Symmetry::print_symmetrized_coordinate(double **x)
{
    int i, j, k, l;
    int isym = 0;
    int nat = system->nat;
    int m11, m12, m13, m21, m22, m23, m31, m32, m33;
    int det;
    double tran[3];
    double **x_symm, **x_avg;
#ifdef _USE_EIGEN
    Eigen::Matrix3d rot;
    Eigen::Vector3d wsi, usi, vsi, tmp;
#else 
    double rot[3][3];
    double wsi[3], usi[3], vsi[3], tmp[3];
#endif

    memory->allocate(x_symm, nat, 3);
    memory->allocate(x_avg, nat, 3);

    for (i = 0; i < nat; ++i) {
        for (j = 0; j < 3; ++j) {
            x_avg[i][j] = 0.0;
        }
    }

    for (std::vector<SymmetryOperation>::iterator it = SymmData.begin();
         it != SymmData.end(); ++it) {

        ++isym;
        std::cout << "Symmetry No. : " << std::setw(5) << isym << std::endl;

        m11 = (*it).rotation[0][0];
        m12 = (*it).rotation[0][1];
        m13 = (*it).rotation[0][2];
        m21 = (*it).rotation[1][0];
        m22 = (*it).rotation[1][1];
        m23 = (*it).rotation[1][2];
        m31 = (*it).rotation[2][0];
        m32 = (*it).rotation[2][1];
        m33 = (*it).rotation[2][2];

        for (i = 0; i < 3; ++i) tran[i] = (*it).tran[i];

        det = m11 * (m22 * m33 - m32 * m23)
            - m21 * (m12 * m33 - m32 * m13)
            + m31 * (m12 * m23 - m22 * m13);

#ifdef _USE_EIGEN
        rot(0, 0) = static_cast<double>((m22 * m33 - m23 * m32) * det);
        rot(0, 1) = static_cast<double>((m23 * m31 - m21 * m33) * det);
        rot(0, 2) = static_cast<double>((m21 * m32 - m22 * m31) * det);
        rot(1, 0) = static_cast<double>((m32 * m13 - m33 * m12) * det);
        rot(1, 1) = static_cast<double>((m33 * m11 - m31 * m13) * det);
        rot(1, 2) = static_cast<double>((m31 * m12 - m32 * m11) * det);
        rot(2, 0) = static_cast<double>((m12 * m23 - m13 * m22) * det);
        rot(2, 1) = static_cast<double>((m13 * m21 - m11 * m23) * det);
        rot(2, 2) = static_cast<double>((m11 * m22 - m12 * m21) * det);
#else
        rot[0][0] = static_cast<double>((m22 * m33 - m23 * m32) * det);
        rot[0][1] = static_cast<double>((m23 * m31 - m21 * m33) * det);
        rot[0][2] = static_cast<double>((m21 * m32 - m22 * m31) * det);
        rot[1][0] = static_cast<double>((m32 * m13 - m33 * m12) * det);
        rot[1][1] = static_cast<double>((m33 * m11 - m31 * m13) * det);
        rot[1][2] = static_cast<double>((m31 * m12 - m32 * m11) * det);
        rot[2][0] = static_cast<double>((m12 * m23 - m13 * m22) * det);
        rot[2][1] = static_cast<double>((m13 * m21 - m11 * m23) * det);
        rot[2][2] = static_cast<double>((m11 * m22 - m12 * m21) * det);
#endif

        for (i = 0; i < nat; ++i) {
            for (j = 0; j < 3; ++j) {
#ifdef _USE_EIGEN
                wsi(j) = x[i][j] - tran[j];
#else 
                wsi[j] = x[i][j] - tran[j];
#endif
            }

#ifdef _USE_EIGEN
            usi = rot * wsi;
#else
            rotvec(usi, wsi, rot);
#endif

            l = -1;

            for (j = 0; j < nat; ++j) {
                for (k = 0; k < 3; ++k) {
#ifdef _USE_EIGEN
                    vsi(k) = x[j][k];
                    tmp(k) = std::fmod(std::abs(usi(k) - vsi(k)), 1.0);
                    // need "std" to specify floating point operation
                    // especially for Intel compiler (there was no problem in MSVC)
                    tmp(k) = std::min<double>(tmp(k), 1.0 - tmp(k)) ;
#else
                    vsi[k] = x[j][k];
                    tmp[k] = std::fmod(std::abs(usi[k] - vsi[k]), 1.0);
                    tmp[k] = std::min<double>(tmp[k], 1.0 - tmp[k]);
#endif
                }
#ifdef _USE_EIGEN
                double diff = tmp.dot(tmp);
#else
                double diff = tmp[0] * tmp[0] + tmp[1] * tmp[1] + tmp[2] * tmp[2];
#endif
                if (diff < tolerance * tolerance) {
                    l = j;
                    break;
                }
            }
            if (l == -1)
                error->exit("print_symmetrized_coordinate",
                            "This cannot happen.");

            for (j = 0; j < 3; ++j) {
#ifdef _USE_EIGEN
                x_symm[l][j] = usi(j);
#else 
                x_symm[l][j] = usi[j];
#endif
                /*
                do {
                if (x_symm[l][j] < 0.0) {
                x_symm[l][j] += 1.0;
                } else if (x_symm[l][j] >= 1.0) {
                x_symm[l][j] -= 1.0;
                }
                } while(x_symm[l][j] < 0.0 || x_symm[l][j] >= 1.0);
                */
            }

        }

        for (i = 0; i < nat; ++i) {
            for (j = 0; j < 3; ++j) {
                std::cout << std::setw(20) << std::scientific << x_symm[i][j];
            }
            std::cout << " ( ";
            for (j = 0; j < 3; ++j) {
                std::cout << std::setw(20) << std::scientific << x_symm[i][j] - x[i][j];
            }
            std::cout << " )" << std::endl;

            for (j = 0; j < 3; ++j) {
                x_avg[i][j] += x_symm[i][j];
            }
        }

    }

    for (i = 0; i < nat; ++i) {
        for (j = 0; j < 3; ++j) {
            x_avg[i][j] /= static_cast<double>(SymmData.size());
        }
    }

    std::cout << "Symmetry Averaged Coordinate" << std::endl;
    for (i = 0; i < nat; ++i) {
        for (j = 0; j < 3; ++j) {
            std::cout << std::setw(20) << std::scientific
                << std::setprecision(9) << x_avg[i][j];
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;

    std::cout.setf(std::ios::floatfield);

    memory->deallocate(x_symm);
    memory->deallocate(x_avg);
}

bool Symmetry::is_proper(const double rot[3][3])
{
    double det;
    bool ret;

    det = rot[0][0] * (rot[1][1] * rot[2][2] - rot[2][1] * rot[1][2])
        - rot[1][0] * (rot[0][1] * rot[2][2] - rot[2][1] * rot[0][2])
        + rot[2][0] * (rot[0][1] * rot[1][2] - rot[1][1] * rot[0][2]);

    if (std::abs(det - 1.0) < eps12) {
        ret = true;
    } else if (std::abs(det + 1.0) < eps12) {
        ret = false;
    } else {
        error->exit("is_proper", "This cannot happen.");
    }

    return ret;
}
