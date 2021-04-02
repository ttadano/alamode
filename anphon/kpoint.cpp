/*
kpoint.cpp

Copyright (c) 2014, 2015, 2016 Terumasa Tadano

This file is distributed under the terms of the MIT license.
Please see the file 'LICENCE.txt' in the root directory
or http://opensource.org/licenses/mit-license.php for information.
*/

#include "mpi_common.h"
#include "kpoint.h"
#include "constants.h"
#include "memory.h"
#include "error.h"
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
    xk = nullptr;
    kaxis = nullptr;
    kvec_na = nullptr;
    knum_minus = nullptr;
    nkx = 0;
    nky = 0;
    nkz = 0;
    nk = 0;
    kp_planes = nullptr;
    kp_planes_tri = nullptr;
    small_group_of_k = nullptr;
}


void Kpoint::deallocate_variables()
{
    if (xk) {
        memory->deallocate(xk);
    }
    if (kaxis) {
        memory->deallocate(kaxis);
    }
    if (kvec_na) {
        memory->deallocate(kvec_na);
    }
    if (knum_minus) {
        memory->deallocate(knum_minus);
    }
    if (kp_planes) {
        memory->deallocate(kp_planes);
    }
    if (kp_planes_tri) {
        memory->deallocate(kp_planes_tri);
    }
    if (small_group_of_k) {
        memory->deallocate(small_group_of_k);
    }
}


void Kpoint::kpoint_setups(const std::string mode)
{
    symmetry->symmetry_flag = true;

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

            setup_kpoint_given(kpInp, nk, xk, kvec_na);

            if (mympi->my_rank == 0) {
                std::cout << "  Number of k points : " << nk << std::endl << std::endl;
                std::cout << "  List of k points : " << std::endl;
                for (i = 0; i < nk; ++i) {
                    std::cout << std::setw(5) << i + 1 << ":";
                    for (j = 0; j < 3; ++j) {
                        std::cout << std::setw(15) << xk[i][j];
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

            setup_kpoint_band(kpInp, nk, xk, kvec_na, kaxis);

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
                std::cout << "  Number of k points : " << nk << std::endl << std::endl;

            }

            break;

        case 2:

            if (mympi->my_rank == 0) {
                std::cout << "  KPMODE = 2: Uniform grid" << std::endl;
            }

            setup_kpoint_mesh(kpInp, nk, nkx, nky, nkz,
                              xk,
                              kvec_na,
                              symmetry->symmetry_flag,
                              kpoint_irred_all);

            nk_irred = kpoint_irred_all.size();

            weight_k.clear();
            for (i = 0; i < kpoint_irred_all.size(); ++i) {
                weight_k.push_back(static_cast<double>(kpoint_irred_all[i].size()) / static_cast<double>(nk));
            }

            if (mympi->my_rank == 0) {
                std::cout << "  Gamma-centered uniform grid with the following mesh density: " << std::endl;
                std::cout << "  nk1:" << std::setw(4) << nkx << std::endl;
                std::cout << "  nk2:" << std::setw(4) << nky << std::endl;
                std::cout << "  nk3:" << std::setw(4) << nkz << std::endl;
                std::cout << std::endl;
                std::cout << "  Number of k points : " << nk << std::endl;
                std::cout << "  Number of irreducible k points : " << kpoint_irred_all.size() << std::endl << std::endl;
                std::cout << "  List of irreducible k points (reciprocal coordinate, weight) : " << std::endl;

                for (i = 0; i < kpoint_irred_all.size(); ++i) {
                    std::cout << "  " << std::setw(5) << i + 1 << ":";
                    for (j = 0; j < 3; ++j) {
                        std::cout << std::setprecision(5) << std::setw(14)
                                  << std::scientific << kpoint_irred_all[i][0].kval[j];
                    }
                    std::cout << std::setprecision(6) << std::setw(11)
                              << std::fixed << weight_k[i] << std::endl;
                }
                std::cout << std::endl;
            }

            memory->allocate(knum_minus, nk);
            gen_nkminus(nk, knum_minus, xk);

            for (i = 0; i < nk_irred; ++i) {
                for (j = 0; j < kpoint_irred_all[i].size(); ++j) {
                    kmap_to_irreducible.insert(std::map<int, int>::value_type(kpoint_irred_all[i][j].knum, i));
                }
            }
            // Compute small group of every irreducible k points for later use
            memory->allocate(small_group_of_k, nk_irred);
            calc_small_groups_k_irred(small_group_of_k);

            break;

        case 3:

            if (mympi->my_rank == 0) {
                std::cout << "  KPMODE = 3: Momentum-resolved final state amplitude" << std::endl;
            }

            setup_kpoint_plane(kpInp, nplanes, kp_planes);

            if (mympi->my_rank == 0) {
                std::cout << "  Number of planes : " << nplanes << std::endl;
                for (i = 0; i < nplanes; ++i) {
                    std::cout << "  The number of k points in plane " << std::setw(3) << i + 1 << " :";
                    std::cout << std::setw(8) << kp_planes[i].size() << std::endl;
                }
                std::cout << std::endl;
            }

            break;

        default:
            error->exit("setup_kpoints", "This cannot happen.");
    }
}

void Kpoint::setup_kpoint_given(const std::vector<KpointInp> &kpinfo,
                                unsigned int &n,
                                double **&k,
                                double **&kdirec) const
{
    int i;

    n = kpinfo.size();
    MPI_Bcast(&n, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

    memory->allocate(k, n, 3);
    memory->allocate(kdirec, n, 3);

    if (mympi->my_rank == 0) {
        int j = 0;
        for (const auto &it : kpinfo) {
            for (i = 0; i < 3; ++i) {
                k[j][i] = std::atof(it.kpelem[i].c_str());
            }

            rotvec(kdirec[j], k[j], system->rlavec_p, 'T');

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
}

void Kpoint::setup_kpoint_band(const std::vector<KpointInp> &kpinfo,
                               unsigned int &n,
                               double **&xk,
                               double **&kdirec,
                               double *&axis) const
{
    int j, k;


    if (mympi->my_rank == 0) {

        std::string **kp_symbol;
        unsigned int *nk_path;
        double **k_start, **k_end;

        const auto npath = kpinfo.size();

        memory->allocate(kp_symbol, npath, 2);
        memory->allocate(k_start, npath, 3);
        memory->allocate(k_end, npath, 3);
        memory->allocate(nk_path, npath);

        n = 0;
        int i = 0;

        for (const auto &it : kpinfo) {
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

        memory->allocate(xk, n, 3);
        memory->allocate(kdirec, n, 3);
        memory->allocate(axis, n);

        unsigned int ik = 0;
        double direc_tmp[3], tmp[3];

        for (i = 0; i < npath; ++i) {
            for (j = 0; j < 3; ++j) {
                direc_tmp[j] = k_end[i][j] - k_start[i][j];
            }

            rotvec(direc_tmp, direc_tmp, system->rlavec_p, 'T');
            auto norm = std::pow(direc_tmp[0], 2)
                        + std::pow(direc_tmp[1], 2)
                        + std::pow(direc_tmp[2], 2);
            norm = std::sqrt(norm);

            if (norm > eps) {
                for (j = 0; j < 3; ++j) direc_tmp[j] /= norm;
            }

            for (j = 0; j < nk_path[i]; ++j) {
                for (k = 0; k < 3; ++k) {
                    xk[ik][k] = k_start[i][k]
                                + (k_end[i][k] - k_start[i][k])
                                  * static_cast<double>(j) / static_cast<double>(nk_path[i] - 1);

                    kdirec[ik][k] = direc_tmp[k];
                }

                if (ik == 0) {
                    axis[ik] = 0.0;
                } else {
                    if (j == 0) {
                        axis[ik] = axis[ik - 1];
                    } else {
                        for (k = 0; k < 3; ++k) tmp[k] = xk[ik][k] - xk[ik - 1][k];
                        rotvec(tmp, tmp, system->rlavec_p, 'T');
                        axis[ik] = axis[ik - 1]
                                   + std::sqrt(tmp[0] * tmp[0]
                                               + tmp[1] * tmp[1]
                                               + tmp[2] * tmp[2]);
                    }
                }
                ++ik;
            }
        }
        memory->deallocate(nk_path);
        memory->deallocate(k_start);
        memory->deallocate(k_end);
        memory->deallocate(kp_symbol);
    }

    MPI_Bcast(&n, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

    if (mympi->my_rank > 0) {
        memory->allocate(xk, n, 3);
        memory->allocate(kdirec, n, 3);
    }

    MPI_Bcast(&xk[0][0], 3 * n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&kdirec[0][0], 3 * n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

void Kpoint::setup_kpoint_mesh(const std::vector<KpointInp> &kpinfo,
                               unsigned int &nk,
                               unsigned int &nkx,
                               unsigned int &nky,
                               unsigned int &nkz,
                               double **&xk,
                               double **&kdirec,
                               const bool usesym,
                               std::vector<std::vector<KpointList>> &kp_irreducible) const
{
    int j;
    unsigned int nk_tmp[3];

    if (mympi->my_rank == 0) {

        nkx = std::atoi(kpinfo[0].kpelem[0].c_str());
        nky = std::atoi(kpinfo[0].kpelem[1].c_str());
        nkz = std::atoi(kpinfo[0].kpelem[2].c_str());

        nk = nkx * nky * nkz;

        memory->allocate(xk, nk, 3);
        memory->allocate(kdirec, nk, 3);

        nk_tmp[0] = nkx;
        nk_tmp[1] = nky;
        nk_tmp[2] = nkz;

        gen_kmesh(usesym, nk_tmp, xk, kp_irreducible);

        for (int i = 0; i < nk; ++i) {
            for (j = 0; j < 3; ++j) kdirec[i][j] = xk[i][j];

            rotvec(kdirec[i], kdirec[i], system->rlavec_p, 'T');
            const auto norm = kdirec[i][0] * kdirec[i][0]
                              + kdirec[i][1] * kdirec[i][1]
                              + kdirec[i][2] * kdirec[i][2];

            if (norm > eps) {
                for (j = 0; j < 3; ++j) kdirec[i][j] /= std::sqrt(norm);
            }
        }
    }

    MPI_Bcast(&nk, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    MPI_Bcast(&nkx, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    MPI_Bcast(&nky, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    MPI_Bcast(&nkz, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

    if (mympi->my_rank > 0) {
        memory->allocate(xk, nk, 3);
        memory->allocate(kdirec, nk, 3);
    }

    MPI_Bcast(&xk[0][0], 3 * nk, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&kdirec[0][0], 3 * nk, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    mpi_broadcast_kpoint_vector(kp_irreducible);
}

void Kpoint::mpi_broadcast_kpoint_vector(std::vector<std::vector<KpointList>> &kp_irreducible) const
{
    int i, j, ik;
    double **xk_tmp;
    unsigned int *knum_tmp;
    unsigned int *kequiv_tmp;

    std::vector<KpointList> kp_group;
    std::vector<double> ktmp;

    // Broadcast kp_irredicible to all threads

    auto nk_irred = kp_irreducible.size();
    MPI_Bcast(&nk_irred, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

    memory->allocate(xk_tmp, nk, 3);
    memory->allocate(knum_tmp, nk);
    memory->allocate(kequiv_tmp, nk_irred);

    if (mympi->my_rank == 0) {
        ik = 0;

        for (i = 0; i < nk_irred; ++i) {
            kequiv_tmp[i] = kp_irreducible[i].size();

            for (j = 0; j < kequiv_tmp[i]; ++j) {
                knum_tmp[ik] = kp_irreducible[i][j].knum;
                for (int k = 0; k < 3; ++k) {
                    xk_tmp[ik][k] = kp_irreducible[i][j].kval[k];
                }
                ++ik;
            }
        }
    }

    MPI_Bcast(&knum_tmp[0], nk, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    MPI_Bcast(&kequiv_tmp[0], nk_irred, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    MPI_Bcast(&xk_tmp[0][0], 3 * nk, MPI_DOUBLE, 0, MPI_COMM_WORLD);


    if (mympi->my_rank > 0) {
        kp_irreducible.clear();
        kp_group.clear();

        ik = 0;

        for (i = 0; i < nk_irred; ++i) {
            kp_group.clear();
            for (j = 0; j < kequiv_tmp[i]; ++j) {
                ktmp.clear();
                ktmp.push_back(xk_tmp[ik][0]);
                ktmp.push_back(xk_tmp[ik][1]);
                ktmp.push_back(xk_tmp[ik][2]);
                kp_group.emplace_back(knum_tmp[ik], ktmp);
                ++ik;
            }
            kp_irreducible.push_back(kp_group);
        }
    }

    memory->deallocate(xk_tmp);
    memory->deallocate(knum_tmp);
    memory->deallocate(kequiv_tmp);
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

        memory->allocate(naxis, nkp, 2);
        memory->allocate(xk_plane, nkp, 3);

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
        memory->deallocate(naxis);
        memory->deallocate(xk_plane);
    }
}


void Kpoint::setup_kpoint_plane(const std::vector<KpointInp> &kpinfo,
                                unsigned int &nplane,
                                std::vector<KpointPlane> *&kp_plane)
{
    nplane = kpinfo.size();
    MPI_Bcast(&nplane, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

    memory->allocate(kp_plane, nplane);
    memory->allocate(kp_planes_tri, nplane);

    if (mympi->my_rank == 0) {
        gen_kpoints_plane(kpinfo, kp_plane, kp_planes_tri);
    }

    mpi_broadcast_kplane_vector(nplane, kp_plane);
}


void Kpoint::gen_kmesh(const bool usesym,
                       const unsigned int nk_in[3],
                       double **xk_out,
                       std::vector<std::vector<KpointList>> &kplist_out) const
{
    unsigned int ik;
    double **xkr;

    const auto nk_tot = nk_in[0] * nk_in[1] * nk_in[2];
    int nsym;


    memory->allocate(xkr, nk_tot, 3);

    for (unsigned int ix = 0; ix < nk_in[0]; ++ix) {
        for (unsigned int iy = 0; iy < nk_in[1]; ++iy) {
            for (unsigned int iz = 0; iz < nk_in[2]; ++iz) {
                ik = iz + iy * nk_in[2] + ix * nk_in[2] * nk_in[1];
                xkr[ik][0] = static_cast<double>(ix) / static_cast<double>(nk_in[0]);
                xkr[ik][1] = static_cast<double>(iy) / static_cast<double>(nk_in[1]);
                xkr[ik][2] = static_cast<double>(iz) / static_cast<double>(nk_in[2]);
            }
        }
    }

    if (usesym) {
        nsym = symmetry->SymmList.size();
    } else {
        nsym = 1;
    }
    reduce_kpoints(nsym, xkr, nk_in, kplist_out);

    for (ik = 0; ik < nk_tot; ++ik) {
        for (unsigned int i = 0; i < 3; ++i) {
            xk_out[ik][i] = xkr[ik][i] - static_cast<double>(nint(xkr[ik][i]));
        }
    }

    memory->deallocate(xkr);
}

void Kpoint::reduce_kpoints(const unsigned int nsym,
                            double **xkr,
                            const unsigned int nk_in[3],
                            std::vector<std::vector<KpointList>> &kplist_out) const
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

    memory->allocate(symop_k, nsym, 3, 3);

    for (isym = 0; isym < nsym; ++isym) {

        for (i = 0; i < 3; ++i) {
            for (j = 0; j < 3; ++j) {
                srot[i][j] = static_cast<double>(symmetry->SymmList[isym].rot[i][j]);
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

    kplist_out.clear();

    const auto nk_tot = nk_in[0] * nk_in[1] * nk_in[2];
    memory->allocate(k_found, nk_tot);

    for (ik = 0; ik < nk_tot; ++ik) k_found[ik] = false;

    for (ik = 0; ik < nk_tot; ++ik) {

        if (k_found[ik]) continue;

        k_group.clear();

        for (i = 0; i < 3; ++i) xk_orig[i] = xkr[ik][i];

        for (isym = 0; isym < nsym; ++isym) {

            rotvec(xk_sym, xk_orig, symop_k[isym]);

            for (i = 0; i < 3; ++i) xk_sym[i] = xk_sym[i] - nint(xk_sym[i]);

            int nloc = get_knum(xk_sym, nk_in);


            if (nloc == -1) {

                error->exit("reduce_kpoints", "Cannot find the kpoint");

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

            if (symmetry->time_reversal_sym) {

                for (i = 0; i < 3; ++i) xk_sym[i] *= -1.0;

                nloc = get_knum(xk_sym, nk_in);

                if (nloc == -1) {

                    error->exit("reduce_kpoints", "Cannot find the kpoint");

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
        kplist_out.push_back(k_group);
    }

    memory->deallocate(k_found);
    memory->deallocate(symop_k);
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
            error->exit("gen_kpoints_plane",
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

        memory->allocate(xk, N1 * N2, 3);

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

        memory->allocate(triangle, number_of_triangle_tiles, 3);

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

        memory->deallocate(xk);
        memory->deallocate(triangle);
    }
}

void Kpoint::gen_nkminus(const unsigned int nk,
                         unsigned int *minus_k,
                         double **xk_in) const
{
    for (unsigned int ik = 0; ik < nk; ++ik) {

        const auto ik_minus = get_knum(-xk_in[ik][0], -xk_in[ik][1], -xk_in[ik][2]);

        if (ik_minus == -1)
            error->exit("gen_nkminus",
                        "-xk doesn't exist on the mesh point.");
        if (ik_minus < ik) continue;

        minus_k[ik] = ik_minus;
        minus_k[ik_minus] = ik;
    }
}

int Kpoint::get_knum(const double kx,
                     const double ky,
                     const double kz) const
{
    double diff[3];
    const auto dkx = static_cast<double>(nkx);
    const auto dky = static_cast<double>(nky);
    const auto dkz = static_cast<double>(nkz);

    diff[0] = static_cast<double>(nint(kx * dkx)) - kx * dkx;
    diff[1] = static_cast<double>(nint(ky * dky)) - ky * dky;
    diff[2] = static_cast<double>(nint(kz * dkz)) - kz * dkz;

    const auto norm = std::sqrt(diff[0] * diff[0]
                                + diff[1] * diff[1]
                                + diff[2] * diff[2]);

    if (norm >= eps12) return -1;

    const int iloc = nint(kx * dkx + 2.0 * dkx) % nkx;
    const int jloc = nint(ky * dky + 2.0 * dky) % nky;
    const int kloc = nint(kz * dkz + 2.0 * dkz) % nkz;

    return kloc + nkz * jloc + nky * nkz * iloc;
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


void Kpoint::generate_irreducible_kmap(int *kequiv,
                                       unsigned int &nk_irreducible,
                                       std::vector<int> &k_irreducible,
                                       const unsigned int n1,
                                       const unsigned int n2,
                                       const unsigned int n3,
                                       double **xk_in,
                                       const int nsymop,
                                       int ***symrot) const
{
    int i, j;
    const int nktot = n1 * n2 * n3;

    double xk_orig[3], xk_sym[3];
    double srot[3][3], srot_inv[3][3], srot_inv_t[3][3];

    for (i = 0; i < nktot; ++i) kequiv[i] = i;

    for (i = 0; i < nktot; ++i) {

        if (kequiv[i] < i) continue;

        for (j = 0; j < 3; ++j) xk_orig[j] = xk_in[i][j];

        for (int isym = 0; isym < nsymop; ++isym) {

            for (j = 0; j < 3; ++j) {
                for (int k = 0; k < 3; ++k) {
                    srot[j][k] = symrot[isym][j][k];
                }
            }

            invmat3(srot_inv, srot);
            transpose3(srot_inv_t, srot_inv);

            rotvec(xk_sym, xk_orig, srot_inv_t);

            const auto ksym = get_knum(xk_sym[0], xk_sym[1], xk_sym[2]);

            if (ksym > kequiv[i]) {
                kequiv[ksym] = i;
            }
        }
    }

    nk_irreducible = 0;
    k_irreducible.clear();

    std::map<int, int> map_to_irreducible_index;

    for (i = 0; i < nktot; ++i) {
        if (kequiv[i] == i) {
            map_to_irreducible_index.insert(std::map<int, int>::value_type(i, nk_irreducible));
            ++nk_irreducible;
            k_irreducible.push_back(i);
        }
    }

    for (i = 0; i < nktot; ++i) {
        kequiv[i] = map_to_irreducible_index[kequiv[i]];
    }

    map_to_irreducible_index.clear();
}

void Kpoint::calc_small_groups_k_irred(std::vector<int> *small_group)
{
    for (int ik = 0; ik < nk_irred; ++ik) {
        small_group[ik] = get_small_group_of_k(kpoint_irred_all[ik][0].knum);
    }
}

std::vector<int> Kpoint::get_small_group_of_k(const int ik) const
{
    std::vector<int> small_group;
    small_group.clear();
    for (auto isym = 0; isym < symmetry->nsym; ++isym) {
        const auto ksym = knum_sym(ik, isym);
        if (ksym == ik) {
            small_group.push_back(isym);
        }
    }
    return small_group;
}

void Kpoint::get_small_group_k(const double *xk_in,
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

int Kpoint::knum_sym(const int ik_in,
                     const int symop_num) const
{
    // Returns kpoint index of S(symop_num)*xk[ik_in]
    // Works only for gamma-centered mesh calculations
    int i;

    double srot[3][3];
    double srot_inv[3][3], srot_inv_t[3][3];
    double xk_orig[3], xk_sym[3];

    if (symop_num < 0 || symop_num >= symmetry->nsym) {
        error->exit("knum_sym", "Invalid symop_num");
    }

    for (i = 0; i < 3; ++i) {
        for (auto j = 0; j < 3; ++j) {
            srot[i][j] = static_cast<double>(symmetry->SymmList[symop_num].rot[i][j]);
        }
    }

    invmat3(srot_inv, srot);
    transpose3(srot_inv_t, srot_inv);

    for (i = 0; i < 3; ++i) xk_orig[i] = xk[ik_in][i];

    rotvec(xk_sym, xk_orig, srot_inv_t);
    for (i = 0; i < 3; ++i) {
        xk_sym[i] = xk_sym[i] - nint(xk_sym[i]);
    }

    return get_knum(xk_sym[0], xk_sym[1], xk_sym[2]);
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
                    error->exit("get_commensurate_kpoints",
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
                for (const auto &sx : signs) {
                    for (const auto &sy : signs) {
                        for (const auto &sz : signs) {
                            comb.push_back({i * sx, j * sy, k * sz});
                        }
                    }
                }
            }
        }
    }

    double qtmp[3], qdiff[3];

    for (const auto &p : comb) {
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

        for (const auto &elem : klist) {

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


void Kpoint::get_unique_triplet_k(const int ik,
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
    int num_group_k;
    int ks_in[2];
    const int knum = kpoint_irred_all[ik][0].knum;
    bool *flag_found;
    std::vector<KsList> kslist;
    double xk[3], xk1[3], xk2[3];

    memory->allocate(flag_found, nk);

    if (use_triplet_symmetry) {
        num_group_k = small_group_of_k[ik].size();
    } else {
        num_group_k = 1;
    }

    for (i = 0; i < 3; ++i) xk[i] = this->xk[knum][i];
    for (i = 0; i < nk; ++i) flag_found[i] = false;

    triplet.clear();

    for (auto ik1 = 0; ik1 < nk; ++ik1) {

        for (i = 0; i < 3; ++i) xk1[i] = this->xk[ik1][i];

        if (sign == -1) {
            for (i = 0; i < 3; ++i) xk2[i] = xk[i] - xk1[i];
        } else if (sign == 1) {
            for (i = 0; i < 3; ++i) xk2[i] = -xk[i] - xk1[i];
        } else {
            error->exit("get_unituq_triplet_k", "Invalid sign");
        }

        const auto ik2 = get_knum(xk2[0], xk2[1], xk2[2]);

        kslist.clear();

        if (ik1 > ik2 && use_permutation_symmetry) continue;

        // Add symmety-connected triplets to kslist
        for (auto isym = 0; isym < num_group_k; ++isym) {

            ks_in[0] = knum_sym(ik1, small_group_of_k[ik][isym]);
            ks_in[1] = knum_sym(ik2, small_group_of_k[ik][isym]);

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

    memory->deallocate(flag_found);
}


void Kpoint::get_unique_quartet_k(const int ik,
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
    int num_group_k;
    std::vector<int> ks_in(3);
    int knum = kpoint_irred_all[ik][0].knum;
    bool **flag_found;
    std::vector<KsList> kslist;
    double xk[3], xk1[3], xk2[3], xk3[3];

    memory->allocate(flag_found, nk, nk);

    if (use_quartet_symmetry) {
        num_group_k = small_group_of_k[ik].size();
    } else {
        num_group_k = 1;
    }

    for (i = 0; i < 3; ++i) xk[i] = this->xk[knum][i];
    for (i = 0; i < nk; ++i) {
        for (int j = 0; j < nk; ++j) {
            flag_found[i][j] = false;
        }
    }

    quartet.clear();

    for (int ik1 = 0; ik1 < nk; ++ik1) {

        for (i = 0; i < 3; ++i) xk1[i] = this->xk[ik1][i];

        int ik2_start = 0;
        if (use_permutation_symmetry) ik2_start = ik1;

        for (int ik2 = ik2_start; ik2 < nk; ++ik2) {
            for (i = 0; i < 3; ++i) xk2[i] = this->xk[ik2][i];

            if (sign == -1) {
                for (i = 0; i < 3; ++i) xk3[i] = xk[i] - xk1[i] - xk2[i];
            } else {
                for (i = 0; i < 3; ++i) xk3[i] = -xk[i] - xk1[i] - xk2[i];
            }

            int ik3 = get_knum(xk3[0], xk3[1], xk3[2]);

            kslist.clear();

            if (ik3 > ik2 && use_permutation_symmetry) continue;

            for (int isym = 0; isym < num_group_k; ++isym) {

                ks_in[0] = knum_sym(ik1, small_group_of_k[ik][isym]);
                ks_in[1] = knum_sym(ik2, small_group_of_k[ik][isym]);
                ks_in[2] = knum_sym(ik3, small_group_of_k[ik][isym]);

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

    memory->deallocate(flag_found);
}
