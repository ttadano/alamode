/*
 symmetry_core.cpp

 Copyright (c) 2014, 2015, 2016 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#include "mpi_common.h"
#include "symmetry_core.h"
#include "constants.h"
#include "error.h"
#include "mathfunctions.h"
#include "memory.h"
#include "system.h"
#include <iomanip>
#include <fstream>
#include <algorithm>

using namespace PHON_NS;

Symmetry::Symmetry(PHON *phon): Pointers(phon)
{
    set_default_variables();
}

Symmetry::~Symmetry() {}


void Symmetry::set_default_variables()
{
    file_sym = "SYMM_INFO_PRIM";
    time_reversal_sym = false;
    nsym = 0;
    symmetry_flag = true;
    printsymmetry = false;
    tolerance = 1.0e-6;
}


void Symmetry::setup_symmetry()
{
    unsigned int natmin = system->natmin;
    double **xtmp;
    unsigned int *kdtmp;

    memory->allocate(xtmp, natmin, 3);
    memory->allocate(kdtmp, natmin);

    for (auto i = 0; i < natmin; ++i) {
        rotvec(xtmp[i], system->xr_s[system->map_p2s[i][0]], system->lavec_s);
        rotvec(xtmp[i], xtmp[i], system->rlavec_p);

        for (auto j = 0; j < 3; ++j) xtmp[i][j] /= 2.0 * pi;

        kdtmp[i] = system->kd[system->map_p2s[i][0]];
    }

    SymmList.clear();

    if (mympi->my_rank == 0) {
        std::cout << " Symmetry" << std::endl;
        std::cout << " ========" << std::endl << std::endl;
        setup_symmetry_operation(natmin,
                                 nsym,
                                 system->lavec_p,
                                 system->rlavec_p,
                                 xtmp,
                                 kdtmp);
    }

    MPI_Bcast(&nsym, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    broadcast_symmlist(SymmList);

    if (mympi->my_rank == 0) {
        std::cout << "  Number of symmetry operations : "
            << nsym << std::endl << std::endl;
        gensym_withmap(xtmp, kdtmp);
    }
    memory->deallocate(xtmp);
    memory->deallocate(kdtmp);
}

void Symmetry::setup_symmetry_operation(int N,
                                        unsigned int &nsym,
                                        double aa[3][3],
                                        double bb[3][3],
                                        double **x,
                                        unsigned int *kd)
{
    int i, j;
    std::ofstream ofs_sym;
    std::ifstream ifs_sym;
    SymmList.clear();

    if (nsym == 0) {

        // Automatically find symmetries.

        std::cout << "  NSYM = 0 is given: Trying to find symmetry operations." << std::endl;

        findsym(N, aa, x, SymmList);

        std::sort(SymmList.begin() + 1, SymmList.end());
        nsym = SymmList.size();

        if (printsymmetry) {
            std::cout
                << "  PRINTSYMM = 1: Symmetry information will be stored in SYMM_INFO_PRIM file."
                << std::endl << std::endl;
            ofs_sym.open(file_sym.c_str(), std::ios::out);
            ofs_sym << nsym << std::endl;

            for (auto p = SymmList.begin(); p != SymmList.end(); ++p) {
                for (i = 0; i < 3; ++i) {
                    for (j = 0; j < 3; ++j) {
                        ofs_sym << std::setw(4) << (*p).rot[i][j];
                    }
                }
                ofs_sym << "  ";
                for (i = 0; i < 3; ++i) {
                    ofs_sym << std::setprecision(15) << std::setw(20) << (*p).tran[i];
                }
                ofs_sym << std::endl;
            }

            ofs_sym.close();
        }

    } else if (nsym == 1) {

        // Identity operation only !

        std::cout << "  NSYM = 1 is given: Only the identity matrix will be considered."
            << std::endl << std::endl;

        int rot_tmp[3][3];
        double tran_tmp[3];

        for (i = 0; i < 3; ++i) {
            for (j = 0; j < 3; ++j) {
                if (i == j) {
                    rot_tmp[i][j] = 1;
                } else {
                    rot_tmp[i][j] = 0;
                }
            }
            tran_tmp[i] = 0.0;
        }

        SymmList.emplace_back(rot_tmp, tran_tmp);

    } else {

        std::cout
            << "  NSYM > 1 is given: Symmetry operations will be read from SYMM_INFO_PRIM file"
            << std::endl << std::endl;

        int nsym2;
        int rot_tmp[3][3];
        double tran_tmp[3];

        ifs_sym.open(file_sym.c_str(), std::ios::in);
        ifs_sym >> nsym2;

        if (nsym != nsym2)
            error->exit("setup_symmetry_operation",
                        "nsym in the given file and the input file are not consistent.");

        for (i = 0; i < nsym; ++i) {
            ifs_sym
                >> rot_tmp[0][0] >> rot_tmp[0][1] >> rot_tmp[0][2]
                >> rot_tmp[1][0] >> rot_tmp[1][1] >> rot_tmp[1][2]
                >> rot_tmp[2][0] >> rot_tmp[2][1] >> rot_tmp[2][2]
                >> tran_tmp[0] >> tran_tmp[1] >> tran_tmp[2];

            SymmList.emplace_back(rot_tmp, tran_tmp);
        }
        ifs_sym.close();
    }
}


void Symmetry::findsym(int N,
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
    find_crystal_symmetry(N, system->nclassatom, system->atomlist_class, x,
                          LatticeSymmList, symop_all);

    LatticeSymmList.clear();
}

void Symmetry::find_lattice_symmetry(double aa[3][3],
                                     std::vector<RotationMatrix> &LatticeSymmList) const
{
    /*
    Find the rotational matrices that leave the metric tensor invariant.

    Metric tensor G = (g)_{ij} = a_{i} * a_{j} is invariant under crystal symmetry operations T,
    i.e. T^{t}GT = G. Since G can be written as G = A^{t}A, the invariance condition is given by
    (AT)^{t}(AT) = G0 (original).
    */

    int i, j, k;

    int nsym_tmp = 0;
    int mat_tmp[3][3];
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
    LatticeSymmList.emplace_back(mat_tmp);

    for (int m11 = -1; m11 <= 1; ++m11) {
        for (int m12 = -1; m12 <= 1; ++m12) {
            for (int m13 = -1; m13 <= 1; ++m13) {
                for (int m21 = -1; m21 <= 1; ++m21) {
                    for (int m22 = -1; m22 <= 1; ++m22) {
                        for (int m23 = -1; m23 <= 1; ++m23) {
                            for (int m31 = -1; m31 <= 1; ++m31) {
                                for (int m32 = -1; m32 <= 1; ++m32) {
                                    for (int m33 = -1; m33 <= 1; ++m33) {

                                        if (m11 == 1 && m12 == 0 && m13 == 0 &&
                                            m21 == 0 && m22 == 1 && m23 == 0 &&
                                            m31 == 0 && m32 == 0 && m33 == 1)
                                            continue;

                                        double det = m11 * (m22 * m33 - m32 * m23)
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

                                        double res = 0.0;
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
                                            LatticeSymmList.emplace_back(mat_tmp);
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

void Symmetry::find_crystal_symmetry(int N,
                                     int nclass,
                                     std::vector<unsigned int> *atomclass,
                                     double **x,
                                     const std::vector<RotationMatrix> &LatticeSymmList,
                                     std::vector<SymmetryOperation> &CrystalSymmList) const
{
    unsigned int i, j;
    unsigned int jat, kat, lat;
    double x_rot[3];
    double rot[3][3], rot_tmp[3][3], rot_cart[3][3];
    double tran[3];
    double x_rot_tmp[3];
    double tmp[3];
    double mag[3], mag_rot[3];
    double diff;

    int rot_int[3][3];

    bool is_found;
    bool isok;
    bool mag_sym1, mag_sym2;
    bool is_identity_matrix;


    // Add identity matrix first.
    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            if (i == j) {
                rot_int[i][j] = 1;
            } else {
                rot_int[i][j] = 0;
            }
        }
        tran[i] = 0.0;
    }

    CrystalSymmList.emplace_back(rot_int, tran);

    for (const auto &it_latsym : LatticeSymmList) {

        unsigned int iat = atomclass[0][0];

        for (i = 0; i < 3; ++i) {
            for (j = 0; j < 3; ++j) {
                rot[i][j] = static_cast<double>(it_latsym.mat[i][j]);
            }
        }

        rotvec(x_rot, x[iat], rot);

#ifdef _OPENMP
#pragma omp parallel for private(jat, tran, isok, kat, x_rot_tmp, is_found, lat, tmp, diff, \
    i, j, is_identity_matrix, mag, mag_rot, rot_tmp, rot_cart, mag_sym1, mag_sym2)
#endif
        for (int ii = 0; ii < atomclass[0].size(); ++ii) {
            jat = atomclass[0][ii];

            for (i = 0; i < 3; ++i) {
                tran[i] = x[jat][i] - x_rot[i];
                tran[i] = tran[i] - nint(tran[i]);
            }

            isok = true;

            is_identity_matrix =
                std::pow(rot[0][0] - 1.0, 2) + std::pow(rot[0][1], 2) + std::pow(rot[0][2], 2)
                + std::pow(rot[1][0], 2) + std::pow(rot[1][1] - 1.0, 2) + std::pow(rot[1][2], 2)
                + std::pow(rot[2][0], 2) + std::pow(rot[2][1], 2) + std::pow(rot[2][2] - 1.0, 2)
                + std::pow(tran[0], 2) + std::pow(tran[1], 2) + std::pow(tran[2], 2) < eps12;

            if (is_identity_matrix) continue;

            for (unsigned int itype = 0; itype < nclass; ++itype) {

                for (int jj = 0; jj < atomclass[itype].size(); ++jj) {

                    kat = atomclass[itype][jj];

                    rotvec(x_rot_tmp, x[kat], rot);

                    for (i = 0; i < 3; ++i) {
                        x_rot_tmp[i] += tran[i];
                    }

                    is_found = false;

                    for (int kk = 0; kk < atomclass[itype].size(); ++kk) {

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

            if (isok && system->lspin && system->noncollinear) {
                for (i = 0; i < 3; ++i) {
                    mag[i] = system->magmom[jat][i];
                    mag_rot[i] = system->magmom[iat][i];
                }

                matmul3(rot_tmp, rot, system->rlavec_p);
                matmul3(rot_cart, system->lavec_p, rot_tmp);

                for (i = 0; i < 3; ++i) {
                    for (j = 0; j < 3; ++j) {
                        rot_cart[i][j] /= 2.0 * pi;
                    }
                }
                rotvec(mag_rot, mag_rot, rot_cart);

                // In the case of improper rotation, the factor -1 should be multiplied
                // because the inversion operation doesn't flip the spin.
                if (!is_proper(rot_cart)) {
                    for (i = 0; i < 3; ++i) {
                        mag_rot[i] = -mag_rot[i];
                    }
                }

                mag_sym1 = std::pow(mag[0] - mag_rot[0], 2.0)
                    + std::pow(mag[1] - mag_rot[1], 2.0)
                    + std::pow(mag[2] - mag_rot[2], 2.0) < eps6;

                mag_sym2 = std::pow(mag[0] + mag_rot[0], 2.0)
                    + std::pow(mag[1] + mag_rot[1], 2.0)
                    + std::pow(mag[2] + mag_rot[2], 2.0) < eps6;

                if (!mag_sym1 && !mag_sym2) {
                    isok = false;
                } else if (!mag_sym1 && mag_sym2 && !trev_sym_mag) {
                    isok = false;
                }
            }


            if (isok) {
#ifdef _OPENMP
#pragma omp critical
#endif
                CrystalSymmList.emplace_back(it_latsym.mat, tran);
            }
        }

    }
}


void Symmetry::gensym_withmap(double **x,
                              const unsigned int *kd)
{
    // Generate symmetry operations in Cartesian coordinate with the atom-mapping information.

    double S[3][3], T[3][3], S_recip[3][3], mat_tmp[3][3];
    double shift[3], x_mod[3], tmp[3];
    unsigned int *map_tmp;
    int i, j;
    unsigned int natmin = system->natmin;

    SymmListWithMap.clear();

    memory->allocate(map_tmp, natmin);

    for (const auto &isym : SymmList) {

        for (i = 0; i < 3; ++i) {
            for (j = 0; j < 3; ++j) {
                T[i][j] = static_cast<double>(isym.rot[i][j]);
            }
        }

        for (i = 0; i < 3; ++i) {
            shift[i] = isym.tran[i];
        }

        invmat3(mat_tmp, T);
        for (i = 0; i < 3; ++i) {
            for (j = 0; j < 3; ++j) {
                S_recip[i][j] = mat_tmp[j][i];
            }
        }

        // Convert to Cartesian coordinate
        matmul3(mat_tmp, T, system->rlavec_p);
        matmul3(S, system->lavec_p, mat_tmp);
        for (i = 0; i < 3; ++i) {
            for (j = 0; j < 3; ++j) {
                S[i][j] /= 2.0 * pi;
            }
        }

        // Generate mapping information

        for (i = 0; i < natmin; ++i) {

            rotvec(x_mod, x[i], T);

            for (j = 0; j < 3; ++j) {
                x_mod[j] += shift[j];
            }

            int num_mapped = -1;

            for (j = 0; j < natmin; ++j) {

                if (kd[j] == kd[i]) {

                    for (int k = 0; k < 3; ++k) {
                        tmp[k] = std::fmod(std::abs(x_mod[k] - x[j][k]), 1.0);
                        tmp[k] = std::min<double>(tmp[k], 1.0 - tmp[k]);
                    }
                    double diff = tmp[0] * tmp[0] + tmp[1] * tmp[1] + tmp[2] * tmp[2];
                    if (diff < tolerance * tolerance) {
                        num_mapped = j;
                        break;
                    }
                }
            }

            if (num_mapped == -1) {
                error->exit("gensym_withmap", "cannot find a equivalent atom");
            }
            map_tmp[i] = num_mapped;
        }

        // Add to vector

        SymmListWithMap.emplace_back(S,
                                     T,
                                     S_recip,
                                     map_tmp,
                                     natmin,
                                     shift);
    }
}


void Symmetry::broadcast_symmlist(std::vector<SymmetryOperation> &sym) const
{
    int i, j, k;
    int n;
    std::vector<int> sym_entry;
    int ***rot_tmp, rot[3][3];
    double **tran_tmp, tran[3];

    if (mympi->my_rank == 0) n = sym.size();
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

    memory->allocate(rot_tmp, n, 3, 3);
    memory->allocate(tran_tmp, n, 3);

    if (mympi->my_rank == 0) {
        for (i = 0; i < n; ++i) {
            for (j = 0; j < 3; ++j) {
                for (k = 0; k < 3; ++k) {
                    rot_tmp[i][j][k] = sym[i].rot[j][k];
                }
                tran_tmp[i][j] = sym[i].tran[j];
            }
        }
    }
    MPI_Bcast(&rot_tmp[0][0][0], 9 * n, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&tran_tmp[0][0], 3 * n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (mympi->my_rank > 0) {
        for (i = 0; i < n; ++i) {
            for (j = 0; j < 3; ++j) {
                for (k = 0; k < 3; ++k) {
                    rot[j][k] = rot_tmp[i][j][k];
                }
                tran[j] = tran_tmp[i][j];
            }
            sym.emplace_back(rot, tran);
        }
    }

    memory->deallocate(rot_tmp);
    memory->deallocate(tran_tmp);
}

bool Symmetry::is_proper(double rot[3][3]) const
{
    bool ret;

    double det = rot[0][0] * (rot[1][1] * rot[2][2] - rot[2][1] * rot[1][2])
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
