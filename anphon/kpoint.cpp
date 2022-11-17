/*
kpoint.cpp

Copyright (c) 2014 Terumasa Tadano

This file is distributed under the terms of the MIT license.
Please see the file 'LICENCE.txt' in the root directory
or http://opensource.org/licenses/mit-license.php for information.
*/

#include "mpi_common.h"
#include "kpoint.h"
#include "constants.h"
#include "memory.h"
#include "error.h"
#include "phonon_dos.h"
#include "system.h"
#include "symmetry_core.h"
#include "timer.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <set>
#include <map>
#include <algorithm>
#include "parsephon.h"
#include "mathfunctions.h"

using namespace PHON_NS;

Kpoint::Kpoint(PHON *phon) : Pointers(phon)
{
    set_default_variables();
}

Kpoint::~Kpoint()
{
    deallocate_variables();
}

void Kpoint::set_default_variables()
{
    kp_planes = nullptr;
    kp_planes_tri = nullptr;
    kpoint_bs = nullptr;
    kpoint_general = nullptr;
}

void Kpoint::deallocate_variables()
{
    if (kp_planes) {
        deallocate(kp_planes);
    }
    if (kp_planes_tri) {
        deallocate(kp_planes_tri);
    }
    if (kpoint_bs) delete kpoint_bs;
    if (kpoint_general) delete kpoint_general;
}

void Kpoint::kpoint_setups(const std::string mode)
{
    unsigned int i, j;
    std::string str_tmp;

    MPI_Bcast(&kpoint_mode, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (mympi->my_rank == 0) {
        std::cout << " k points" << std::endl;
        std::cout << " ========" << std::endl << std::endl;
    }

    switch (kpoint_mode) {
        case 0:

            if (mympi->my_rank == 0) {
                std::cout << "  KPMODE = 0 : Calculation on given k points" << std::endl;
            }

            setup_kpoint_given(kpInp, system->rlavec_p);

            if (mympi->my_rank == 0) {
                std::cout << "  Number of k points : " << kpoint->kpoint_general->nk << std::endl << std::endl;
                std::cout << "  List of k points : " << std::endl;
                for (i = 0; i < kpoint->kpoint_general->nk; ++i) {
                    std::cout << std::setw(5) << i + 1 << ":";
                    for (j = 0; j < 3; ++j) {
                        std::cout << std::setw(15) << kpoint->kpoint_general->xk[i][j];
                    }
                    std::cout << std::endl;
                }
                std::cout << std::endl;
            }

            break;

        case 1:

            if (mympi->my_rank == 0) {
                std::cout << "  KPMODE = 1: Band structure calculation" << std::endl;
            }

            setup_kpoint_band(kpInp, system->rlavec_p);
            if (mympi->my_rank == 0) {
                std::cout << "  Number of paths : " << kpInp.size() << std::endl << std::endl;
                std::cout << "  List of k paths : " << std::endl;

                for (i = 0; i < kpInp.size(); ++i) {
                    std::cout << std::setw(4) << i + 1 << ":";
                    std::cout << std::setw(3) << kpInp[i].kpelem[0];
                    std::cout << " (";
                    for (int k = 0; k < 3; ++k) {
                        std::cout << std::setprecision(4) << std::setw(8)
                                  << std::atof(kpInp[i].kpelem[k + 1].c_str());
                    }
                    std::cout << ")";
                    std::cout << std::setw(3) << kpInp[i].kpelem[4];
                    std::cout << " (";
                    for (int k = 0; k < 3; ++k) {
                        std::cout << std::setprecision(4) << std::setw(8)
                                  << std::atof(kpInp[i].kpelem[k + 5].c_str());
                    }
                    std::cout << ")";
                    std::cout << std::setw(4) << kpInp[i].kpelem[8] << std::endl;
                }
                std::cout << std::endl;
                std::cout << "  Number of k points : " << kpoint_bs->nk << std::endl << std::endl;

            }

            break;

        case 2:

            if (mympi->my_rank == 0) {
                std::cout << "  KPMODE = 2: Uniform grid" << std::endl;
            }

            unsigned int nk_tmp[3];
            nk_tmp[0] = 0;
            nk_tmp[1] = 0;
            nk_tmp[2] = 0;
            if (mympi->my_rank == 0) {
                for (i = 0; i < 3; ++i) {
                    nk_tmp[i] = std::atoi(kpInp[0].kpelem[i].c_str());
                }
            }
            MPI_Bcast(&nk_tmp[0], 3, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
            dos->kmesh_dos = new KpointMeshUniform(nk_tmp);
            dos->kmesh_dos->setup(symmetry->SymmList,
                                  system->rlavec_p,
                                  symmetry->time_reversal_sym);

            if (mympi->my_rank == 0) {
                std::cout << "  Gamma-centered uniform grid with the following mesh density: " << std::endl;
                std::cout << "  nk1:" << std::setw(4) << dos->kmesh_dos->nk_i[0] << std::endl;
                std::cout << "  nk2:" << std::setw(4) << dos->kmesh_dos->nk_i[1] << std::endl;
                std::cout << "  nk3:" << std::setw(4) << dos->kmesh_dos->nk_i[2] << std::endl;
                std::cout << std::endl;
                std::cout << "  Number of k points : " << dos->kmesh_dos->nk << std::endl;
                std::cout << "  Number of irreducible k points : " << dos->kmesh_dos->nk_irred << std::endl <<
                          std::endl;
                std::cout << "  List of irreducible k points (reciprocal coordinate, weight) : " << std::endl;

                for (i = 0; i < dos->kmesh_dos->nk_irred; ++i) {
                    std::cout << "  " << std::setw(5) << i + 1 << ":";
                    for (j = 0; j < 3; ++j) {
                        std::cout << std::setprecision(5) << std::setw(14)
                                  << std::scientific << dos->kmesh_dos->kpoint_irred_all[i][0].kval[j];
                    }
                    std::cout << std::setprecision(6) << std::setw(11)
                              << std::fixed << dos->kmesh_dos->weight_k[i] << std::endl;
                }
                std::cout << std::endl;
            }

            break;

        default:
            exit("setup_kpoints", "This cannot happen.");
    }
}

void Kpoint::setup_kpoint_given(const std::vector<KpointInp> &kpinfo,
                                const double rlavec_p[3][3])
{
    int i;
    double **k, **kdirec;
    unsigned int n = kpinfo.size();

    MPI_Bcast(&n, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    allocate(k, n, 3);
    allocate(kdirec, n, 3);

    if (mympi->my_rank == 0) {
        int j = 0;
        for (const auto &it: kpinfo) {
            for (i = 0; i < 3; ++i) {
                k[j][i] = std::atof(it.kpelem[i].c_str());
            }

            rotvec(kdirec[j], k[j], rlavec_p, 'T');

            const auto norm = kdirec[j][0] * kdirec[j][0]
                              + kdirec[j][1] * kdirec[j][1]
                              + kdirec[j][2] * kdirec[j][2];

            if (norm > eps) {
                for (i = 0; i < 3; ++i) kdirec[j][i] /= std::sqrt(norm);
            }

            ++j;
        }
    }

    MPI_Bcast(&k[0][0], 3 * n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&kdirec[0][0], 3 * n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    kpoint_general = new KpointGeneral(n, k, kdirec);

    deallocate(k);
    deallocate(kdirec);
}

void Kpoint::setup_kpoint_band(const std::vector<KpointInp> &kpinfo,
                               const double rlavec_p[3][3])
{
    int j, k;

    double **xk_tmp;
    double **kdirec_tmp;
    double *axis_tmp;
    unsigned int n = 0;

    if (mympi->my_rank == 0) {

        std::string **kp_symbol;
        unsigned int *nk_path;
        double **k_start, **k_end;

        const auto npath = kpinfo.size();

        allocate(kp_symbol, npath, 2);
        allocate(k_start, npath, 3);
        allocate(k_end, npath, 3);
        allocate(nk_path, npath);

        n = 0;
        int i = 0;

        for (const auto &it: kpinfo) {
            kp_symbol[i][0] = it.kpelem[0];
            kp_symbol[i][1] = it.kpelem[4];

            for (j = 0; j < 3; ++j) {
                k_start[i][j] = std::atof(it.kpelem[j + 1].c_str());
                k_end[i][j] = std::atof(it.kpelem[j + 5].c_str());

            }
            nk_path[i] = std::atoi(it.kpelem[8].c_str());
            n += nk_path[i];
            ++i;
        }

        allocate(xk_tmp, n, 3);
        allocate(kdirec_tmp, n, 3);
        allocate(axis_tmp, n);

        unsigned int ik = 0;
        double direc_tmp[3], tmp[3];

        for (i = 0; i < npath; ++i) {
            for (j = 0; j < 3; ++j) {
                direc_tmp[j] = k_end[i][j] - k_start[i][j];
            }

            rotvec(direc_tmp, direc_tmp, rlavec_p, 'T');
            auto norm = std::pow(direc_tmp[0], 2)
                        + std::pow(direc_tmp[1], 2)
                        + std::pow(direc_tmp[2], 2);
            norm = std::sqrt(norm);

            if (norm > eps) {
                for (j = 0; j < 3; ++j) direc_tmp[j] /= norm;
            }

            for (j = 0; j < nk_path[i]; ++j) {
                for (k = 0; k < 3; ++k) {
                    xk_tmp[ik][k] = k_start[i][k]
                                    + (k_end[i][k] - k_start[i][k])
                                      * static_cast<double>(j) / static_cast<double>(nk_path[i] - 1);

                    kdirec_tmp[ik][k] = direc_tmp[k];
                }

                if (ik == 0) {
                    axis_tmp[ik] = 0.0;
                } else {
                    if (j == 0) {
                        axis_tmp[ik] = axis_tmp[ik - 1];
                    } else {
                        for (k = 0; k < 3; ++k) tmp[k] = xk_tmp[ik][k] - xk_tmp[ik - 1][k];
                        rotvec(tmp, tmp, rlavec_p, 'T');
                        axis_tmp[ik] = axis_tmp[ik - 1]
                                       + std::sqrt(tmp[0] * tmp[0]
                                                   + tmp[1] * tmp[1]
                                                   + tmp[2] * tmp[2]);
                    }
                }
                ++ik;
            }
        }
        deallocate(nk_path);
        deallocate(k_start);
        deallocate(k_end);
        deallocate(kp_symbol);
    }

    MPI_Bcast(&n, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

    if (mympi->my_rank > 0) {
        allocate(xk_tmp, n, 3);
        allocate(kdirec_tmp, n, 3);
        allocate(axis_tmp, n);
    }

    MPI_Bcast(&xk_tmp[0][0], 3 * n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&kdirec_tmp[0][0], 3 * n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&axis_tmp[0], n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    kpoint_bs = new KpointBandStructure(n, xk_tmp, kdirec_tmp, axis_tmp);

    deallocate(xk_tmp);
    deallocate(kdirec_tmp);
    deallocate(axis_tmp);
}

void KpointMeshUniform::setup(const std::vector<SymmetryOperation> &symmlist,
                              const double rlavec_p[3][3],
                              const bool time_reversal_symmetry)
{
    const bool usesym = true;

    gen_kmesh(symmlist, usesym, time_reversal_symmetry);

    for (auto i = 0; i < nk; ++i) {
        for (auto j = 0; j < 3; ++j) kvec_na[i][j] = xk[i][j];

        rotvec(&kvec_na[i][0], &kvec_na[i][0], rlavec_p, 'T');
        const auto norm = kvec_na[i][0] * kvec_na[i][0]
                          + kvec_na[i][1] * kvec_na[i][1]
                          + kvec_na[i][2] * kvec_na[i][2];

        if (norm > eps) {
            for (auto j = 0; j < 3; ++j) kvec_na[i][j] /= std::sqrt(norm);
        }
    }

    nk_irred = kpoint_irred_all.size();
    weight_k.resize(nk_irred);
    for (auto i = 0; i < nk_irred; ++i) {
        weight_k[i]
                = static_cast<double>(kpoint_irred_all[i].size())
                  / static_cast<double>(nk);
    }
    gen_nkminus();

    kmap_to_irreducible.resize(nk);
    for (auto i = 0; i < nk_irred; ++i) {
        for (auto j = 0; j < kpoint_irred_all[i].size(); ++j) {
            kmap_to_irreducible[kpoint_irred_all[i][j].knum] = i;
        }
    }
    // Compute small group of every irreducible k points for later use
    small_group_of_k.resize(nk_irred);
    set_small_groups_k_irred(usesym, symmlist);
}

void KpointMeshUniform::gen_kmesh(const std::vector<SymmetryOperation> &symmlist,
                                  const bool usesym,
                                  const bool time_reversal_symmetry)
{
    unsigned int ik;
    double **xkr;
    unsigned int nsym;

    allocate(xkr, nk, 3);

    for (unsigned int ix = 0; ix < nk_i[0]; ++ix) {
        for (unsigned int iy = 0; iy < nk_i[1]; ++iy) {
            for (unsigned int iz = 0; iz < nk_i[2]; ++iz) {
                ik = iz + iy * nk_i[2] + ix * nk_i[2] * nk_i[1];
                xkr[ik][0] = static_cast<double>(ix) / static_cast<double>(nk_i[0]);
                xkr[ik][1] = static_cast<double>(iy) / static_cast<double>(nk_i[1]);
                xkr[ik][2] = static_cast<double>(iz) / static_cast<double>(nk_i[2]);
            }
        }
    }
    if (usesym) {
        nsym = symmlist.size();
    } else {
        nsym = 1;
    }
    reduce_kpoints(nsym, symmlist, time_reversal_symmetry, xkr);

    for (ik = 0; ik < nk; ++ik) {
        for (unsigned int i = 0; i < 3; ++i) {
            xk[ik][i] = xkr[ik][i] - static_cast<double>(nint(xkr[ik][i]));
        }
    }

    deallocate(xkr);
}

void KpointMeshUniform::reduce_kpoints(const unsigned int nsym,
                                       const std::vector<SymmetryOperation> &symmlist,
                                       const bool time_reversal_symmetry,
                                       double **xkr)
{
    unsigned int ik;
    unsigned int i, j;
    int isym;

    bool *k_found;

    std::vector<KpointList> k_group;
    std::vector<double> ktmp;

    double srot[3][3];
    double xk_sym[3], xk_orig[3];
    double srot_inv[3][3], srot_inv_t[3][3];

    double ***symop_k;

    allocate(symop_k, nsym, 3, 3);

    for (isym = 0; isym < nsym; ++isym) {

        for (i = 0; i < 3; ++i) {
            for (j = 0; j < 3; ++j) {
                srot[i][j] = static_cast<double>(symmlist[isym].rot[i][j]);
            }
        }

        invmat3(srot_inv, srot);
        transpose3(srot_inv_t, srot_inv);

        for (i = 0; i < 3; ++i) {
            for (j = 0; j < 3; ++j) {
                symop_k[isym][i][j] = srot_inv_t[i][j];
            }
        }
    }

    kpoint_irred_all.clear();

    allocate(k_found, nk);

    for (ik = 0; ik < nk; ++ik) k_found[ik] = false;

    for (ik = 0; ik < nk; ++ik) {

        if (k_found[ik]) continue;

        k_group.clear();

        for (i = 0; i < 3; ++i) xk_orig[i] = xkr[ik][i];

        for (isym = 0; isym < nsym; ++isym) {

            rotvec(xk_sym, xk_orig, symop_k[isym]);

            for (i = 0; i < 3; ++i) xk_sym[i] = xk_sym[i] - nint(xk_sym[i]);

            int nloc = get_knum(xk_sym);

            if (nloc == -1) {

                //     exit("reduce_kpoints", "Cannot find the kpoint");

            } else {

                if (!k_found[nloc]) {
                    k_found[nloc] = true;
                    ktmp.clear();
                    ktmp.push_back(xk_sym[0]);
                    ktmp.push_back(xk_sym[1]);
                    ktmp.push_back(xk_sym[2]);

                    k_group.emplace_back(nloc, ktmp);
                }

            }

            // Time-reversal symmetry

            if (time_reversal_symmetry) {

                for (i = 0; i < 3; ++i) xk_sym[i] *= -1.0;

                nloc = get_knum(xk_sym);

                if (nloc == -1) {

                    //     exit("reduce_kpoints", "Cannot find the kpoint");

                } else {

                    if (!k_found[nloc]) {
                        k_found[nloc] = true;
                        ktmp.clear();
                        ktmp.push_back(xk_sym[0]);
                        ktmp.push_back(xk_sym[1]);
                        ktmp.push_back(xk_sym[2]);

                        k_group.emplace_back(nloc, ktmp);
                    }
                }
            }
        }
        kpoint_irred_all.push_back(k_group);
    }

    deallocate(k_found);
    deallocate(symop_k);
}

void KpointMeshUniform::gen_nkminus()
{
    kindex_minus_xk.resize(nk);
    double minus_xk[3];

    for (unsigned int ik = 0; ik < nk; ++ik) {

        for (auto i = 0; i < 3; ++i) minus_xk[i] = -xk[ik][i];

        const auto ik_minus = get_knum(minus_xk);

        if (ik_minus == -1) {
//            exit("gen_nkminus",
//                        "-xk doesn't exist on the mesh point.");
        }

        if (ik_minus < ik) continue;

        kindex_minus_xk[ik] = ik_minus;
        kindex_minus_xk[ik_minus] = ik;
    }
}

void KpointMeshUniform::set_small_groups_k_irred(const bool usesym,
                                                 const std::vector<SymmetryOperation> &symmlist)
{
    small_group_of_k.resize(nk_irred);
    for (auto ik = 0; ik < nk_irred; ++ik) {
        small_group_of_k[ik]
                = get_small_group_of_k(kpoint_irred_all[ik][0].knum,
                                       usesym,
                                       symmlist);
    }
}

std::vector<int> KpointMeshUniform::get_small_group_of_k(const unsigned int ik,
                                                         const bool usesym,
                                                         const std::vector<SymmetryOperation> &symmlist) const
{
    std::vector<int> small_group;
    small_group.clear();
    unsigned int nsym;
    if (usesym) {
        nsym = symmlist.size();
    } else {
        nsym = 1;
    }
    for (auto isym = 0; isym < nsym; ++isym) {
        const auto ksym = knum_sym(ik, symmlist[isym].rot);
        if (ksym == ik) {
            small_group.push_back(isym);
        }
    }
    return small_group;
}

int KpointMeshUniform::knum_sym(const unsigned int ik,
                                const int rot[3][3]) const
{
    // Returns kpoint index of S(symop_num)*xk[ik_in]
    // Works only for gamma-centered mesh calculations
    int i;

    double srot[3][3];
    double srot_inv[3][3], srot_inv_t[3][3];
    double xk_orig[3], xk_sym[3];

    for (i = 0; i < 3; ++i) {
        for (auto j = 0; j < 3; ++j) {
            srot[i][j] = static_cast<double>(rot[i][j]);
        }
    }

    invmat3(srot_inv, srot);
    transpose3(srot_inv_t, srot_inv);

    for (i = 0; i < 3; ++i) xk_orig[i] = xk[ik][i];

    rotvec(xk_sym, xk_orig, srot_inv_t);
    for (i = 0; i < 3; ++i) {
        xk_sym[i] = xk_sym[i] - nint(xk_sym[i]);
    }

    return get_knum(xk_sym);
}

int KpointMeshUniform::get_knum(const double xk[3]) const
{
    int i;
    double diff[3];
    double dnk[3];

    for (i = 0; i < 3; ++i) dnk[i] = static_cast<double>(nk_i[i]);
    for (i = 0; i < 3; ++i) diff[i] = static_cast<double>(nint(xk[i] * dnk[i])) - xk[i] * dnk[i];

    const auto norm = std::sqrt(diff[0] * diff[0]
                                + diff[1] * diff[1]
                                + diff[2] * diff[2]);

    if (norm >= eps12) return -1;

    const int iloc = nint(xk[0] * dnk[0] + 2.0 * dnk[0]) % nk_i[0];
    const int jloc = nint(xk[1] * dnk[1] + 2.0 * dnk[1]) % nk_i[1];
    const int kloc = nint(xk[2] * dnk[2] + 2.0 * dnk[2]) % nk_i[2];

    return kloc + nk_i[2] * jloc + nk_i[1] * nk_i[2] * iloc;
}

void Kpoint::mpi_broadcast_kplane_vector(const unsigned int nplane,
                                         std::vector<KpointPlane> *&kp_plane) const
{
    int j;
    int **naxis;
    double **xk_plane;

    for (int i = 0; i < nplane; ++i) {
        int nkp = kp_plane[i].size();

        MPI_Bcast(&nkp, 1, MPI_INT, 0, MPI_COMM_WORLD);

        allocate(naxis, nkp, 2);
        allocate(xk_plane, nkp, 3);

        if (mympi->my_rank == 0) {
            for (j = 0; j < nkp; ++j) {
                naxis[j][0] = kp_plane[i][j].n[0];
                naxis[j][1] = kp_plane[i][j].n[1];
                xk_plane[j][0] = kp_plane[i][j].k[0];
                xk_plane[j][1] = kp_plane[i][j].k[1];
                xk_plane[j][2] = kp_plane[i][j].k[2];
            }
        }

        MPI_Bcast(&naxis[0][0], 2 * nkp, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&xk_plane[0][0], 3 * nkp, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        if (mympi->my_rank > 0) {
            for (j = 0; j < nkp; ++j) {
                kp_plane[i].emplace_back(xk_plane[j], naxis[j]);
            }
        }
        deallocate(naxis);
        deallocate(xk_plane);
    }
}

void Kpoint::setup_kpoint_plane(const std::vector<KpointInp> &kpinfo,
                                unsigned int &nplane,
                                std::vector<KpointPlane> *&kp_plane)
{
    nplane = kpinfo.size();
    MPI_Bcast(&nplane, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

    allocate(kp_plane, nplane);
    allocate(kp_planes_tri, nplane);

    if (mympi->my_rank == 0) {
        gen_kpoints_plane(kpinfo, kp_plane, kp_planes_tri);
    }

    mpi_broadcast_kplane_vector(nplane, kp_plane);
}

void Kpoint::gen_kpoints_plane(const std::vector<KpointInp> &kplist,
                               std::vector<KpointPlane> *kpout,
                               std::vector<KpointPlaneTriangle> *kpout_tri)
{
    int j;
    const int nplane = kplist.size();

    int ik1, ik2;

    int **triangle;

    double xk_tmp[3];
    double xk0[3], xk1[3], xk2[3], xk3[3];
    double **xk;
    int n_in[2];

    for (int i = 0; i < nplane; ++i) {
        const auto N1 = std::atoi(kplist[i].kpelem[3].c_str());
        const auto N2 = std::atoi(kplist[i].kpelem[7].c_str());

        n_in[0] = N1;
        n_in[1] = N2;

        const auto frac1 = 1.0 / static_cast<double>(N1 - 1);
        const auto frac2 = 1.0 / static_cast<double>(N2 - 1);

        for (j = 0; j < 3; ++j) {
            xk0[j] = 0.0;
            xk1[j] = std::atof(kplist[i].kpelem[j].c_str());
            xk2[j] = std::atof(kplist[i].kpelem[4 + j].c_str());
        }

        auto dprod = 0.0;
        auto norm1 = 0.0;
        auto norm2 = 0.0;

        for (j = 0; j < 3; ++j) {
            dprod += xk1[j] * xk2[j];
            norm1 += xk1[j] * xk1[j];
            norm2 += xk2[j] * xk2[j];
        }
        const auto costheta = dprod / std::sqrt(norm1 * norm2);
        if (std::abs(std::abs(costheta) - 1.0) < eps12) {
            exit("gen_kpoints_plane",
                 "Two vectors have to be linearly independent with each other.");
        }

        kp_plane_geometry.emplace_back(xk0, xk1, xk2, n_in);

        for (ik1 = 0; ik1 < N1; ++ik1) {
            for (ik2 = 0; ik2 < N2; ++ik2) {

                for (j = 0; j < 3; ++j) {
                    xk_tmp[j] = static_cast<double>(ik1) * frac1 * xk1[j]
                                + static_cast<double>(ik2) * frac2 * xk2[j];
                }
                if (in_first_BZ(xk_tmp)) {
                    n_in[0] = ik1;
                    n_in[1] = ik2;
                    kpout[i].emplace_back(xk_tmp, n_in);
                }
            }
        }

        allocate(xk, N1 * N2, 3);

        int m = 0;
        for (ik1 = 0; ik1 < N1; ++ik1) {
            for (ik2 = 0; ik2 < N2; ++ik2) {

                for (j = 0; j < 3; ++j) {
                    xk[m][j] = static_cast<double>(ik1) * frac1 * xk1[j]
                               + static_cast<double>(ik2) * frac2 * xk2[j];
                }
                ++m;
            }
        }

        const auto number_of_tiles = (N1 - 1) * (N2 - 1);
        const auto number_of_triangle_tiles = 2 * number_of_tiles;

        allocate(triangle, number_of_triangle_tiles, 3);

        for (ik1 = 0; ik1 < N1 - 1; ++ik1) {
            for (ik2 = 0; ik2 < N2 - 1; ++ik2) {

                const auto n1 = ik2 + ik1 * N2;
                const auto n2 = ik2 + (ik1 + 1) * N2;
                const auto n3 = ik2 + 1 + ik1 * N2;
                const auto n4 = ik2 + 1 + (ik1 + 1) * N2;

                m = 2 * (ik2 + ik1 * (N2 - 1));

                triangle[m][0] = n1;
                triangle[m][1] = n2;
                triangle[m][2] = n4;

                ++m;

                triangle[m][0] = n1;
                triangle[m][1] = n3;
                triangle[m][2] = n4;

            }
        }

        for (int itri = 0; itri < number_of_triangle_tiles; ++itri) {

            for (j = 0; j < 3; ++j) {
                xk1[j] = xk[triangle[itri][0]][j];
                xk2[j] = xk[triangle[itri][1]][j];
                xk3[j] = xk[triangle[itri][2]][j];
            }

            const auto is_inside_FBZ = in_first_BZ(xk1) || in_first_BZ(xk2) || in_first_BZ(xk3);
            if (is_inside_FBZ) {
                kpout_tri[i].emplace_back(itri, triangle[itri]);
            }
        }

        deallocate(xk);
        deallocate(triangle);
    }
}

int Kpoint::get_knum(const double xk[3],
                     const unsigned int nk[3]) const
{
    int i;
    double diff[3];
    double dnk[3];

    for (i = 0; i < 3; ++i) dnk[i] = static_cast<double>(nk[i]);
    for (i = 0; i < 3; ++i) diff[i] = static_cast<double>(nint(xk[i] * dnk[i])) - xk[i] * dnk[i];

    const auto norm = std::sqrt(diff[0] * diff[0]
                                + diff[1] * diff[1]
                                + diff[2] * diff[2]);

    if (norm >= eps12) return -1;

    const int iloc = nint(xk[0] * dnk[0] + 2.0 * dnk[0]) % nk[0];
    const int jloc = nint(xk[1] * dnk[1] + 2.0 * dnk[1]) % nk[1];
    const int kloc = nint(xk[2] * dnk[2] + 2.0 * dnk[2]) % nk[2];

    return kloc + nk[2] * jloc + nk[1] * nk[2] * iloc;
}

bool Kpoint::in_first_BZ(const double *xk_in) const
{
    int i;
    const auto nmax = 1;
    double tmp[3];

    for (i = 0; i < 3; ++i) tmp[i] = xk_in[i];

    rotvec(tmp, tmp, system->rlavec_p, 'T');

    const auto dist_min = std::sqrt(tmp[0] * tmp[0] + tmp[1] * tmp[1] + tmp[2] * tmp[2]);

    auto ncount = 0;

    auto iloc = ncount;

    for (i = -nmax; i <= nmax; ++i) {
        for (int j = -nmax; j <= nmax; ++j) {
            for (int k = -nmax; k <= nmax; ++k) {

                if (i == 0 && j == 0 && k == 0) continue;

                ++ncount;

                tmp[0] = xk_in[0] - static_cast<double>(i);
                tmp[1] = xk_in[1] - static_cast<double>(j);
                tmp[2] = xk_in[2] - static_cast<double>(k);

                rotvec(tmp, tmp, system->rlavec_p, 'T');
                const auto dist = std::sqrt(tmp[0] * tmp[0] + tmp[1] * tmp[1] + tmp[2] * tmp[2]);

                if (dist < dist_min) {
                    iloc = ncount;
                    break;
                }
            }
        }
    }

    return iloc == 0;
}

void Kpoint::get_symmetrization_matrix_at_k(const double *xk_in,
                                            std::vector<int> &sym_list,
                                            double S_avg[3][3]) const
{
    int i, j;
    double srot[3][3];
    double srot_inv[3][3], srot_inv_t[3][3];
    double xk_orig[3], xk_sym[3];
    double xk_diff[3];

    sym_list.clear();

    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            S_avg[i][j] = 0.0;
        }
    }

    for (i = 0; i < 3; ++i) xk_orig[i] = xk_in[i];

    for (auto isym = 0; isym < symmetry->nsym; ++isym) {

        for (i = 0; i < 3; ++i) {
            for (j = 0; j < 3; ++j) {
                srot[i][j] = static_cast<double>(symmetry->SymmList[isym].rot[i][j]);
            }
        }

        invmat3(srot_inv, srot);
        transpose3(srot_inv_t, srot_inv);
        rotvec(xk_sym, xk_orig, srot_inv_t);

        for (i = 0; i < 3; ++i) {
            xk_sym[i] = xk_sym[i] - nint(xk_sym[i]);
            xk_diff[i] = std::fmod(xk_sym[i] - xk_orig[i], 1.0);
        }

        if (std::sqrt(std::pow(xk_diff[0], 2)
                      + std::pow(xk_diff[1], 2)
                      + std::pow(xk_diff[2], 2)) < eps10) {
            sym_list.push_back(isym);

            for (i = 0; i < 3; ++i) {
                for (j = 0; j < 3; ++j) {
                    S_avg[i][j] += srot_inv_t[i][j];
                }
            }
        }
    }

    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            S_avg[i][j] /= static_cast<double>(sym_list.size());
        }
    }
}

void Kpoint::get_commensurate_kpoints(const double lavec_super[3][3],
                                      const double lavec_prim[3][3],
                                      std::vector<std::vector<double>> &klist) const
{
    int i, j;
    double inv_lavec_super[3][3];
    double convmat[3][3];

    invmat3(inv_lavec_super, lavec_super);
    matmul3(convmat, inv_lavec_super, lavec_prim);
    transpose3(convmat, convmat);

    const auto det = convmat[0][0] * (convmat[1][1] * convmat[2][2] - convmat[2][1] * convmat[1][2])
                     - convmat[1][0] * (convmat[0][1] * convmat[2][2] - convmat[2][1] * convmat[0][2])
                     + convmat[2][0] * (convmat[0][1] * convmat[1][2] - convmat[1][1] * convmat[0][2]);

    const auto nkmax = static_cast<int>(std::ceil(1.0 / det));

    const auto tol = 1.0e-6;
    const auto max_denom = 10000; // for safety
    int k;

    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {

            const auto frac = std::abs(convmat[i][j]);

            if (frac < tol) {
                convmat[i][j] = 0.0;
            } else {
                auto found_denom = false;
                for (k = 1; k < max_denom; ++k) {
                    if (std::abs(frac * static_cast<double>(k) - 1.0) < tol) {
                        found_denom = true;
                        break;
                    }
                }
                if (found_denom) {
                    const auto sign = (0.0 < convmat[i][j]) - (convmat[i][j] < 0.0);
                    convmat[i][j] = static_cast<double>(sign) / static_cast<double>(k);
                } else {
                    exit("get_commensurate_kpoints",
                         "The denominator of the conversion matrix > 10000");
                }
            }
        }
    }

    const auto nmax_cell = static_cast<int>(nkmax) / 2 + 1;
    std::vector<int> signs({-1, 1});
    std::vector<std::vector<int>> comb;
    comb.clear();
    for (i = 0; i < nmax_cell; ++i) {
        for (j = 0; j < nmax_cell; ++j) {
            for (k = 0; k < nmax_cell; ++k) {
                for (const auto &sx: signs) {
                    for (const auto &sy: signs) {
                        for (const auto &sz: signs) {
                            comb.push_back({i * sx, j * sy, k * sz});
                        }
                    }
                }
            }
        }
    }

    double qtmp[3], qdiff[3];

    for (const auto &p: comb) {
        for (i = 0; i < 3; ++i) qtmp[i] = p[i];
        rotvec(qtmp, qtmp, convmat);

        for (i = 0; i < 3; ++i) {
            qtmp[i] = std::fmod(qtmp[i], 1.0);
            if (qtmp[i] >= 0.5) {
                qtmp[i] -= 1.0;
            } else if (qtmp[i] < -0.5) {
                qtmp[i] += 1.0;
            }
        }

        auto new_entry = true;

        for (const auto &elem: klist) {

            for (i = 0; i < 3; ++i) {
                qdiff[i] = std::fmod(qtmp[i] - elem[i], 1.0);
                if (qdiff[i] >= 0.5) {
                    qdiff[i] -= 1.0;
                } else if (qdiff[i] < -0.5) {
                    qdiff[i] += 1.0;
                }
            }

            const auto norm = std::sqrt(qdiff[0] * qdiff[0]
                                        + qdiff[1] * qdiff[1]
                                        + qdiff[2] * qdiff[2]);

            if (norm < tol) {
                new_entry = false;
                break;
            }
        }
        if (new_entry) klist.push_back({qtmp[0], qtmp[1], qtmp[2]});

        if (klist.size() == nkmax) break;
    }
}

void KpointMeshUniform::get_unique_triplet_k(const int ik,
                                             const std::vector<SymmetryOperation> &symmlist,
                                             const bool use_triplet_symmetry,
                                             const bool use_permutation_symmetry,
                                             std::vector<KsListGroup> &triplet,
                                             const int sign) const
{
    // This function returns the irreducible set of (k2, k3) satisfying the momentum conservation.
    // When sign = -1 (default), pairs satisfying - k1 + k2 + k3 = G are returned.
    // When sign =  1, pairs satisfying k1 + k2 + k3 = G are returned.
    //

    int i;
    unsigned int num_group_k;
    int ks_in[2];
    const auto knum = kpoint_irred_all[ik][0].knum;
    bool *flag_found;
    std::vector<KsList> kslist;
    double xk0[3], xk1[3], xk2[3];

    allocate(flag_found, nk);

    if (use_triplet_symmetry) {
        num_group_k = small_group_of_k[ik].size();
    } else {
        num_group_k = 1;
    }

    for (i = 0; i < 3; ++i) xk0[i] = xk[knum][i];
    for (i = 0; i < nk; ++i) flag_found[i] = false;

    triplet.clear();

    for (auto ik1 = 0; ik1 < nk; ++ik1) {

        for (i = 0; i < 3; ++i) xk1[i] = xk[ik1][i];

        if (sign == -1) {
            for (i = 0; i < 3; ++i) xk2[i] = xk0[i] - xk1[i];
        } else if (sign == 1) {
            for (i = 0; i < 3; ++i) xk2[i] = -xk0[i] - xk1[i];
        } else {
            //exit("get_unituq_triplet_k", "Invalid sign");
        }

        const auto ik2 = get_knum(xk2);

        kslist.clear();

        if (ik1 > ik2 && use_permutation_symmetry) continue;

        // Add symmety-connected triplets to kslist
        for (auto isym = 0; isym < num_group_k; ++isym) {

            ks_in[0] = knum_sym(ik1, symmlist[small_group_of_k[ik][isym]].rot);
            ks_in[1] = knum_sym(ik2, symmlist[small_group_of_k[ik][isym]].rot);

            if (!flag_found[ks_in[0]]) {
                kslist.emplace_back(2, ks_in, small_group_of_k[ik][isym]);
                flag_found[ks_in[0]] = true;
            }

            if (ks_in[0] != ks_in[1] && use_permutation_symmetry && !flag_found[ks_in[1]]) {
                const auto tmp = ks_in[0];
                ks_in[0] = ks_in[1];
                ks_in[1] = tmp;

                kslist.emplace_back(2, ks_in, small_group_of_k[ik][isym]);
                flag_found[ks_in[0]] = true;
            }
        }
        if (!kslist.empty()) {
            triplet.emplace_back(kslist);
        }
    }

    deallocate(flag_found);
}

void KpointMeshUniform::get_unique_quartet_k(const int ik,
                                             const std::vector<SymmetryOperation> &symmlist,
                                             const bool use_quartet_symmetry,
                                             const bool use_permutation_symmetry,
                                             std::vector<KsListGroup> &quartet,
                                             const int sign) const
{
    // This function returns the irreducible set of (k2, k3, k4) satisfying the momentum conservation.
    // When sign = -1 (default), pairs satisfying - k1 + k2 + k3 + k4 = G are returned.
    // When sign =  1, pairs satisfying k1 + k2 + k3 + k4 = G are returned.
    //

    int i;
    unsigned int num_group_k;
    std::vector<int> ks_in(3);
    const auto knum = kpoint_irred_all[ik][0].knum;
    bool **flag_found;
    std::vector<KsList> kslist;
    double xk0[3], xk1[3], xk2[3], xk3[3];

    allocate(flag_found, nk, nk);

    if (use_quartet_symmetry) {
        num_group_k = small_group_of_k[ik].size();
    } else {
        num_group_k = 1;
    }

    for (i = 0; i < 3; ++i) xk0[i] = xk[knum][i];
    for (i = 0; i < nk; ++i) {
        for (int j = 0; j < nk; ++j) {
            flag_found[i][j] = false;
        }
    }

    quartet.clear();

    for (int ik1 = 0; ik1 < nk; ++ik1) {

        for (i = 0; i < 3; ++i) xk1[i] = xk[ik1][i];

        int ik2_start = 0;
        if (use_permutation_symmetry) ik2_start = ik1;

        for (int ik2 = ik2_start; ik2 < nk; ++ik2) {
            for (i = 0; i < 3; ++i) xk2[i] = this->xk[ik2][i];

            if (sign == -1) {
                for (i = 0; i < 3; ++i) xk3[i] = xk0[i] - xk1[i] - xk2[i];
            } else {
                for (i = 0; i < 3; ++i) xk3[i] = -xk0[i] - xk1[i] - xk2[i];
            }

            int ik3 = get_knum(xk3);

            kslist.clear();

            if (ik3 > ik2 && use_permutation_symmetry) continue;

            for (int isym = 0; isym < num_group_k; ++isym) {

                ks_in[0] = knum_sym(ik1, symmlist[small_group_of_k[ik][isym]].rot);
                ks_in[1] = knum_sym(ik2, symmlist[small_group_of_k[ik][isym]].rot);
                ks_in[2] = knum_sym(ik3, symmlist[small_group_of_k[ik][isym]].rot);

                if (!flag_found[ks_in[0]][ks_in[1]]) {
                    kslist.emplace_back(3, &ks_in[0], small_group_of_k[ik][isym]);
                    flag_found[ks_in[0]][ks_in[1]] = true;
                }

                if (use_permutation_symmetry) {
                    std::sort(ks_in.begin(), ks_in.end());
                    do {
                        if (!flag_found[ks_in[0]][ks_in[1]]) {
                            kslist.emplace_back(3, &ks_in[0], small_group_of_k[ik][isym]);
                            flag_found[ks_in[0]][ks_in[1]] = true;
                        }
                    } while (std::next_permutation(ks_in.begin(), ks_in.end()));
                }
            }

            if (!kslist.empty()) {
                quartet.emplace_back(kslist);
            }
        }
    }

    deallocate(flag_found);
}
