#include "mpi_common.h"
#include "symmetry_core.h"
#include "memory.h"
#include "constants.h"
#include "error.h"
#include "system.h"
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <set>
#include "mathfunctions.h"

#ifdef _USE_EIGEN
#include <Eigen/Core>
#endif

using namespace PHON_NS;

Symmetry::Symmetry(PHON *phon): Pointers(phon){
    file_sym = "SYMM_INFO_PRIM";
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
        rotvec(xtmp[i], system->xr_s[system->map_p2s[i][0]], system->lavec_s);
        rotvec(xtmp[i], xtmp[i], system->rlavec_p);

        for (j = 0; j < 3; ++j) xtmp[i][j] /= 2.0 * pi;

        kdtmp[i] = system->kd[system->map_p2s[i][0]];
    }

    SymmList.clear();

    if (mympi->my_rank == 0) {
        setup_symmetry_operation(natmin, nsym, nnp, system->lavec_p, system->rlavec_p, xtmp, kdtmp);
    }

    MPI_Bcast(&nsym, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    broadcast_symmlist(SymmList);

    if (mympi->my_rank == 0) {
        std::cout << "**Symmetry" << std::endl;
        std::cout << "Number of Symmetry Operations = " << nsym << std::endl << std::endl;
    }
    if (mympi->my_rank == 0) {
        gensym_withmap(xtmp, kdtmp);
    }
}

void Symmetry::setup_symmetry_operation(int N, unsigned int &nsym, unsigned int &nnp, double aa[3][3], double bb[3][3], double **x, unsigned int *kd)
{
    int i, j;

    SymmList.clear();

    if(nsym == 0) {

        // Automatically find symmetries.

        std::cout << "NSYM = 0 is given: Trying to find symmetry operations." << std::endl;

        findsym(N, aa, x, SymmList);

        std::sort(SymmList.begin() + 1, SymmList.end());
        nsym = SymmList.size();

        if (printsymmetry) {
            std::cout << "PRINTSYMM = 1: Symmetry information will be stored in SYMM_INFO_PRIM file." << std::endl << std::endl;
            ofs_sym.open(file_sym.c_str(), std::ios::out);
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
    } 
    else if(nsym == 1) {

        // Identity operation only !

        std::cout << "NSYM = 1 is given: Only the identity matrix will be considered." << std::endl << std::endl;

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
        nnp = 1;
    } 
    else {

        std::cout << "NSYM > 1 is given: Symmetry operations will be read from SYMM_INFO file" << std::endl << std::endl;

        int nsym2;
        int rot_tmp[3][3], tran_tmp[3];

        ifs_sym.open(file_sym.c_str(), std::ios::in);
        ifs_sym >> nsym2 >> nnp;

        if(nsym != nsym2) error->exit("setup_symmetry_operation", "nsym in the given file and the input file are not consistent.");

        for (i = 0; i < nsym; ++i) {
            ifs_sym >> rot_tmp[0][0] >> rot_tmp[0][1] >> rot_tmp[0][2]
            >> rot_tmp[1][0] >> rot_tmp[1][1] >> rot_tmp[1][2] 
            >> rot_tmp[2][0] >> rot_tmp[2][1] >> rot_tmp[2][2]
            >> tran_tmp[0] >> tran_tmp[1] >> tran_tmp[2];

            SymmList.push_back(SymmetryOperation(rot_tmp, tran_tmp));
        }
        ifs_sym.close();
    }

#ifdef _DEBUG
    // print_symmetrized_coordinate(x);
#endif
}


void Symmetry::findsym(int N, double aa[3][3], double **x, std::vector<SymmetryOperation> &symop) {

    unsigned int i;
    int tran_int[3];

    std::vector<RotationMatrix> LatticeSymmList;
    std::vector<SymmetryOperationTransFloat> CrystalSymmList;

    LatticeSymmList.clear();
    find_lattice_symmetry(aa, LatticeSymmList);
    CrystalSymmList.clear();
    find_crystal_symmetry(N, system->nclassatom, system->atomlist_class, x, LatticeSymmList, CrystalSymmList);

    find_nnp_for_translation(nnp, CrystalSymmList);

    symop.clear();
    for(std::vector<SymmetryOperationTransFloat>::iterator p = CrystalSymmList.begin(); p != CrystalSymmList.end(); ++p){
        SymmetryOperationTransFloat sym_tmp = *p;

        for(i = 0; i < 3; ++i){
            tran_int[i] = nint(static_cast<double>(nnp) * sym_tmp.tran[i]);
            if (tran_int[i] < 0) tran_int[i] += nnp;
        }

        symop.push_back(SymmetryOperation(sym_tmp.rot, tran_int));
    }
    LatticeSymmList.clear();
    CrystalSymmList.clear();
}

void Symmetry::find_lattice_symmetry(double aa[3][3], std::vector<RotationMatrix> &LatticeSymmList) {

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

    for (m11 = -1; m11 <= 1; ++m11){
        for (m12 = -1; m12 <= 1; ++m12) {
            for (m13 = -1; m13 <= 1; ++m13){
                for (m21 = -1; m21 <= 1; ++m21){
                    for (m22 = -1; m22 <= 1; ++m22){
                        for (m23 = -1; m23 <= 1; ++m23){
                            for (m31 = -1; m31 <= 1; ++m31){
                                for (m32 = -1; m32 <= 1; ++m32){
                                    for (m33 = -1; m33 <= 1; ++m33){

                                        if (m11 == 1 && m12 == 0 && m13 == 0 &&
                                            m21 == 0 && m22 == 1 && m23 == 0 &&
                                            m31 == 0 && m32 == 0 && m33 == 1) continue;

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

void Symmetry::find_crystal_symmetry(int N, int nclass, std::vector<unsigned int> *atomclass, double **x, std::vector<RotationMatrix> LatticeSymmList, std::vector<SymmetryOperationTransFloat> &CrystalSymmList){

    unsigned int i, j;
    unsigned int iat, jat, kat, lat;
    double x_rot[3];
    double rot[3][3];
    double tran[3];
    double x_rot_tmp[3];
    double tmp[3];
    double diff;

    int rot_int[3][3];

    int ii, jj, kk;
    unsigned int itype;

    bool is_found;
    bool isok;

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

    CrystalSymmList.push_back(SymmetryOperationTransFloat(rot_int, tran));


    for (std::vector<RotationMatrix>::iterator it_latsym = LatticeSymmList.begin(); it_latsym != LatticeSymmList.end(); ++it_latsym) {

        iat = atomclass[0][0];

        for (i = 0; i < 3; ++i) {
            for (j = 0; j < 3; ++j) {
                rot[i][j] = static_cast<double>((*it_latsym).mat[i][j]);
            }
        }

        rotvec(x_rot, x[iat], rot);

#ifdef _OPENMP
#pragma omp parallel for private(jat, tran, isok, kat, x_rot_tmp, is_found, lat, tmp, diff, i, itype, jj, kk, is_identity_matrix)
#endif
        for (ii = 0; ii < atomclass[0].size(); ++ii) {
            jat = atomclass[0][ii];

            for (i = 0; i < 3; ++i) {
                tran[i] = x[jat][i] - x_rot[i];
                tran[i] = tran[i] - nint(tran[i]);
            }

            isok = true;

            is_identity_matrix = 
                ( std::pow(rot[0][0] - 1.0, 2) + std::pow(rot[0][1], 2) + std::pow(rot[0][2], 2) 
                + std::pow(rot[1][0], 2) + std::pow(rot[1][1] - 1.0, 2) + std::pow(rot[1][2], 2)
                + std::pow(rot[2][0], 2) + std::pow(rot[2][1], 2) + std::pow(rot[2][2] - 1.0, 2)
                + std::pow(tran[0], 2) + std::pow(tran[1], 2) + std::pow(tran[2], 2) ) < eps12;

            if (is_identity_matrix) continue;

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
#ifdef _OPENMP
#pragma omp critical
#endif
                CrystalSymmList.push_back(SymmetryOperationTransFloat((*it_latsym).mat, tran));
            }
        }

    }
}

void Symmetry::find_nnp_for_translation(unsigned int &ret, std::vector<SymmetryOperationTransFloat> symminfo) {

    int i;

    ret = 1;

    std::set<double> translation_set;
    double tran_tmp;
    double tran_numerator;
    bool is_integer;
    bool is_found;

    translation_set.clear();

    for (std::vector<SymmetryOperationTransFloat>::iterator it = symminfo.begin(); it != symminfo.end(); ++it) {

        for (i = 0; i < 3; ++i) {
            tran_tmp = std::abs((*it).tran[i]);

            if (translation_set.find(tran_tmp) == translation_set.end()) {
                translation_set.insert(tran_tmp);
            }
        }
    }

    is_found = false;

    while(1) {

        is_integer = true;

        for (std::set<double>::iterator it = translation_set.begin(); it != translation_set.end(); ++it) {

            tran_numerator = (*it) * static_cast<double>(ret);
            if (std::abs(tran_numerator - static_cast<double>(nint(tran_numerator))) > tolerance) {
                is_integer = false;
                break;
            }
        }

        if (is_integer) {
            is_found = true;
            break;
        }

        if (ret > 1000) break;
        ++ret;
    }


    if (!is_found) {
        // This should not happen.
        error->exit("find_nnp_for_translation", "Cannot find nnp.");
    }
}


void Symmetry::gensym_withmap(double **x, unsigned int *kd)
{
    // Generate symmetry operations in Cartesian coordinate with the atom-mapping information.

    double S[3][3], T[3][3], S_recip[3][3], mat_tmp[3][3];
    double shift[3], x_mod[3], S_double[3][3], tmp[3];
    double diff;
    unsigned int *map_tmp;
    int i, j, k;
    int num_mapped;
    unsigned int natmin = system->natmin;

    SymmListWithMap.clear();

    memory->allocate(map_tmp, natmin);

    for (std::vector<SymmetryOperation>::iterator isym = SymmList.begin(); isym != SymmList.end(); ++isym) {

        for (i = 0; i < 3; ++i) {
            for (j = 0; j < 3; ++j) {
                T[i][j] = static_cast<double>((*isym).symop[3 * i + j]);
            }
        }

        for (i = 0; i < 3; ++i) {
            shift[i] = static_cast<double>((*isym).symop[i + 9]) / static_cast<double>(nnp);
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

            num_mapped = -1;

            for (j = 0; j < natmin; ++j) {

                if (kd[j] == kd[i]) {

                    for (k = 0; k < 3; ++k) {
                        tmp[k] = std::fmod(std::abs(x_mod[k] - x[j][k]), 1.0);
                        tmp[k] = std::min<double>(tmp[k], 1.0 - tmp[k]);	
                    }
                    diff = tmp[0] * tmp[0] + tmp[1] * tmp[1] + tmp[2] * tmp[2];
                    if (diff < eps12) {
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

        SymmListWithMap.push_back(SymmetryOperationWithMapping(S, T, S_recip, map_tmp, natmin, shift));
    }
}


void Symmetry::broadcast_symmlist(std::vector<SymmetryOperation> &sym)
{
    int i, j;
    int n;
    std::vector<int> sym_entry;
    int **sym_tmp;

    if (mympi->my_rank == 0) n = sym.size();
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

    memory->allocate(sym_tmp, n, 9);

    if (mympi->my_rank == 0) {
        for (i = 0; i < n; ++i) {
            for (j = 0; j < 9; ++j) {
                sym_tmp[i][j] = sym[i].symop[j];
            }
        }
    }
    MPI_Bcast(&sym_tmp[0][0], 9*n, MPI_INT, 0, MPI_COMM_WORLD);

    if (mympi->my_rank > 0) {
        for (i = 0; i < n; ++i) {
            sym_entry.clear();
            for (j = 0; j < 9; ++j) {
                sym_entry.push_back(sym_tmp[i][j]);
            }
            sym.push_back(SymmetryOperation(sym_entry));
        }
    }

    memory->deallocate(sym_tmp);
}