/*
 ewald.cpp

 Copyright (c) 2015 Tatsuro Nishimoto
 Copyright (c) 2016 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#include "mpi_common.h"
#include "ewald.h"
#include "constants.h"
#include "dynamical.h"
#include "error.h"
#include "kpoint.h"
#include "mathfunctions.h"
#include "memory.h"
#include "parsephon.h"
#include "system.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <complex>
#include <boost/math/special_functions/erf.hpp>

using namespace PHON_NS;

Ewald::Ewald(PHON *phon) : Pointers(phon)
{
    set_default_variables();
}

Ewald::~Ewald()
{
    deallocate_variables();
}

void Ewald::set_default_variables()
{
    is_longrange = false;
    print_fc2_ewald = false;
    file_longrange = "";
    prec_ewald = 1.0e-12;
    rate_ab = 1.0;
    multiplicity = nullptr;
    Born_charge = nullptr;
    distall_ewald = nullptr;
    force_permutation_sym = true;
}

void Ewald::deallocate_variables()
{
    if (multiplicity) {
        deallocate(multiplicity);
    }
    if (Born_charge) {
        deallocate(Born_charge);
    }
    if (distall_ewald) {
        deallocate(distall_ewald);
    }
}

void Ewald::init()
{
    MPI_Bcast(&is_longrange, 1, MPI_CXX_BOOL, 0, MPI_COMM_WORLD);
    MPI_Bcast(&prec_ewald, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&rate_ab, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (is_longrange) {
        int nsize[3] = {1, 1, 1};

        allocate(multiplicity, system->nat, system->nat);
        allocate(Born_charge, system->natmin, 3, 3);
        allocate(distall_ewald, system->nat, system->nat);

        get_pairs_of_minimum_distance(system->nat, nsize, system->xr_s);

        for (int i = 0; i < system->natmin; ++i) {
            for (int j = 0; j < 3; ++j) {
                for (int k = 0; k < 3; ++k) {
                    Born_charge[i][j][k] = dynamical->borncharge[i][j][k];
                }
            }
        }

        prepare_Ewald(dynamical->dielec);
        prepare_G();
        compute_ewald_fcs();
    }
}

void Ewald::prepare_Ewald(const double dielectric[3][3])
{
    int icrd;

    double lavec_norm[3], rlavec_norm[3];
    double e_lavec[3], e_rlavec[3];
    double lavec_enorm, rlavec_enorm, lavec_min[2], rlavec_min[2];

    if (mympi->my_rank == 0) {
        std::cout << std::endl;
        std::cout << "  Preparing for the Ewald summation ..." << std::endl << std::endl;
    }

    double p = -std::log(prec_ewald);

    for (icrd = 0; icrd < 3; ++icrd) {
        for (int jcrd = 0; jcrd < 3; ++jcrd) {
            epsilon[icrd][jcrd] = dielectric[icrd][jcrd];
        }
    }

    // Calculating convergence parameters
    invmat3(epsilon_inv, epsilon);

    // For calculating Coulombic (dipole-dipole) FCs
    for (icrd = 0; icrd < 3; ++icrd) {
        lavec_norm[icrd] = std::sqrt(std::pow(system->lavec_s[icrd][0], 2.0)
                                     + std::pow(system->lavec_s[icrd][1], 2.0)
                                     + std::pow(system->lavec_s[icrd][2], 2.0));

        rlavec_norm[icrd] = std::sqrt(std::pow(system->rlavec_s[icrd][0], 2.0)
                                      + std::pow(system->rlavec_s[icrd][1], 2.0)
                                      + std::pow(system->rlavec_s[icrd][2], 2.0));

        rotvec(e_lavec, system->lavec_s[icrd], epsilon_inv);
        rotvec(e_rlavec, system->rlavec_s[icrd], epsilon);

        lavec_enorm = std::sqrt(system->lavec_s[icrd][0] * e_lavec[0]
                                + system->lavec_s[icrd][1] * e_lavec[1]
                                + system->lavec_s[icrd][2] * e_lavec[2]);

        rlavec_enorm = std::sqrt(system->rlavec_s[icrd][0] * e_rlavec[0]
                                 + system->rlavec_s[icrd][1] * e_rlavec[1]
                                 + system->rlavec_s[icrd][2] * e_rlavec[2]);

        if (lavec_enorm < lavec_min[0] || icrd == 0) {
            lavec_min[0] = lavec_enorm;
            lavec_min[1] = lavec_norm[icrd];
        }
        if (rlavec_enorm < rlavec_min[0] || icrd == 0) {
            rlavec_min[0] = rlavec_enorm;
            rlavec_min[1] = rlavec_norm[icrd];
        }
    }

    Lmax_sub = std::sqrt(2.0 * p / (lavec_min[0] * rlavec_min[0] * std::pow(rate_ab, 1.0 / 3.0))) * lavec_min[1];
    Gmax_sub = std::sqrt(2.0 * p / (lavec_min[0] * rlavec_min[0]) * std::pow(rate_ab, 1.0 / 3.0)) * rlavec_min[1];
    lambda_sub = std::sqrt(rlavec_min[0] / (2.0 * lavec_min[0]) * std::pow(rate_ab, 1.0 / 3.0));

    for (icrd = 0; icrd < 3; ++icrd) {
        nl_sub[icrd] = static_cast<int>(Lmax_sub / lavec_norm[icrd]) + 1;
        ng_sub[icrd] = static_cast<int>(Gmax_sub / rlavec_norm[icrd]) + 1;
    }
    num_l_sub = (2 * nl_sub[0] + 1) * (2 * nl_sub[1] + 1) * (2 * nl_sub[2] + 1);
    num_g_sub = (2 * ng_sub[0] + 1) * (2 * ng_sub[1] + 1) * (2 * ng_sub[2] + 1);


    // For calculating Coulombic (dipole-dipole) dynamical matrix
    for (icrd = 0; icrd < 3; ++icrd) {
        lavec_norm[icrd] = std::sqrt(std::pow(system->lavec_p[icrd][0], 2.0)
                                     + std::pow(system->lavec_p[icrd][1], 2.0)
                                     + std::pow(system->lavec_p[icrd][2], 2.0));

        rlavec_norm[icrd] = std::sqrt(std::pow(system->rlavec_p[icrd][0], 2.0)
                                      + std::pow(system->rlavec_p[icrd][1], 2.0)
                                      + std::pow(system->rlavec_p[icrd][2], 2.0));

        rotvec(e_lavec, system->lavec_p[icrd], epsilon_inv);
        lavec_enorm = std::sqrt(system->lavec_p[icrd][0] * e_lavec[0]
                                + system->lavec_p[icrd][1] * e_lavec[1]
                                + system->lavec_p[icrd][2] * e_lavec[2]);

        rotvec(e_rlavec, system->rlavec_p[icrd], epsilon);
        rlavec_enorm = std::sqrt(system->rlavec_p[icrd][0] * e_rlavec[0]
                                 + system->rlavec_p[icrd][1] * e_rlavec[1]
                                 + system->rlavec_p[icrd][2] * e_rlavec[2]);

        if (lavec_enorm < lavec_min[0] || icrd == 0) {
            lavec_min[0] = lavec_enorm;
            lavec_min[1] = lavec_norm[icrd];
        }
        if (rlavec_enorm < rlavec_min[0] || icrd == 0) {
            rlavec_min[0] = rlavec_enorm;
            rlavec_min[1] = rlavec_norm[icrd];
        }
    }

    Lmax = std::sqrt(2.0 * p / (lavec_min[0] * rlavec_min[0] * std::pow(rate_ab, 1.0 / 3.0))) * lavec_min[1];
    Gmax = std::sqrt(2.0 * p / (lavec_min[0] * rlavec_min[0]) * std::pow(rate_ab, 1.0 / 3.0)) * rlavec_min[1];
    lambda = std::sqrt(rlavec_min[0] / (2.0 * lavec_min[0]) * std::pow(rate_ab, 1.0 / 3.0));

    for (icrd = 0; icrd < 3; ++icrd) {
        nl[icrd] = static_cast<int>(Lmax / lavec_norm[icrd]) + 1;
        ng[icrd] = static_cast<int>(Gmax / rlavec_norm[icrd]) + 1;
    }
    num_l = (2 * nl[0] + 1) * (2 * nl[1] + 1) * (2 * nl[2] + 1);
    num_g = (2 * ng[0] + 1) * (2 * ng[1] + 1) * (2 * ng[2] + 1);

    det_epsilon = epsilon[0][0] * (epsilon[1][1] * epsilon[2][2] - epsilon[1][2] * epsilon[2][1])
                  - epsilon[0][1] * (epsilon[1][0] * epsilon[2][2] - epsilon[1][2] * epsilon[2][0])
                  + epsilon[0][2] * (epsilon[1][0] * epsilon[2][1] - epsilon[1][1] * epsilon[2][0]);

    if (mympi->my_rank == 0) {

        std::cout << "  Inverse dielectric tensor : " << std::endl;
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                std::cout << std::setw(15) << epsilon_inv[i][j];
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;

        std::cout << "  Determinant of epsilon: " << std::setw(15) << det_epsilon
                  << std::endl << std::endl;

        std::cout << "  Parameters for the Ewald summation :" << std::endl;
        std::cout << "  - Force constant" << std::endl;
        std::cout << "    Lambda : " << std::setw(15) << lambda_sub << std::endl;
        std::cout << "    Lmax   : " << std::setw(15) << Lmax_sub << std::endl;
        std::cout << "    Gmax   : " << std::setw(15) << Gmax_sub << std::endl;
        std::cout << "    Maximum number of real-space cells : "
                  << std::setw(3) << nl_sub[0] << "x" << std::setw(3) << nl_sub[1] << "x" << std::setw(3) << nl_sub[2]
                  << std
                  ::endl;
        std::cout << "    Maximum number of reciprocal cells : "
                  << std::setw(3) << ng_sub[0] << "x" << std::setw(3) << ng_sub[1] << "x" << std::setw(3) << ng_sub[2]
                  << std
                  ::endl;
        std::cout << std::endl;
        std::cout << "  - Dynamical matrix" << std::endl;
        std::cout << "    Lambda : " << std::setw(15) << lambda << std::endl;
        std::cout << "    Lmax   : " << std::setw(15) << Lmax << std::endl;
        std::cout << "    Gmax   : " << std::setw(15) << Gmax << std::endl;
        std::cout << "    Maximum number of real-space cells : "
                  << std::setw(3) << nl[0] << "x" << std::setw(3) << nl[1] << "x" << std::setw(3) << nl[2] << std::endl;
        std::cout << "    Maximum number of reciprocal cells : "
                  << std::setw(3) << ng[0] << "x" << std::setw(3) << ng[1] << "x" << std::setw(3) << ng[2] << std::endl;
        std::cout << std::endl;
    }
}

void Ewald::prepare_G()
{
    // Accumulate reciprocal lattice vectors

    int ix, iy, iz;
    double g_tmp[3], gnorm;

    G_vector_sub.clear();
    G_vector.clear();

    for (ix = -ng_sub[0]; ix <= ng_sub[0]; ++ix) {
        for (iy = -ng_sub[1]; iy <= ng_sub[1]; ++iy) {
            for (iz = -ng_sub[2]; iz <= ng_sub[2]; ++iz) {
                if (ix == 0 && iy == 0 && iz == 0) continue;
                g_tmp[0] = static_cast<double>(ix);
                g_tmp[1] = static_cast<double>(iy);
                g_tmp[2] = static_cast<double>(iz);
                rotvec(g_tmp, g_tmp, system->rlavec_s, 'T');
                gnorm = std::sqrt(g_tmp[0] * g_tmp[0] + g_tmp[1] * g_tmp[1] + g_tmp[2] * g_tmp[2]);
                if (gnorm <= Gmax_sub) {
                    G_vector_sub.emplace_back(g_tmp);
                }
            }
        }
    }

    for (ix = -ng[0]; ix <= ng[0]; ++ix) {
        for (iy = -ng[1]; iy <= ng[1]; ++iy) {
            for (iz = -ng[2]; iz <= ng[2]; ++iz) {
                if (ix == 0 && iy == 0 && iz == 0) continue;
                g_tmp[0] = static_cast<double>(ix);
                g_tmp[1] = static_cast<double>(iy);
                g_tmp[2] = static_cast<double>(iz);
                rotvec(g_tmp, g_tmp, system->rlavec_p, 'T');
                gnorm = std::sqrt(g_tmp[0] * g_tmp[0] + g_tmp[1] * g_tmp[1] + g_tmp[2] * g_tmp[2]);
                if (gnorm <= Gmax) {
                    G_vector.emplace_back(g_tmp);
                }
            }
        }
    }
}

void Ewald::get_pairs_of_minimum_distance(const int nat,
                                          const int nsize[3],
                                          double **xf) const
{
    // Get pairs and multiplicities

    int icell = 0;
    int iat, jat;
    double dist_tmp;
    double ***xcrd; // fractional coordinate

    int ncell = (2 * nsize[0] + 1) * (2 * nsize[1] + 1) * (2 * nsize[2] + 1);

    allocate(xcrd, ncell, nat, 3);

    for (iat = 0; iat < nat; ++iat) {
        for (int icrd = 0; icrd < 3; ++icrd) {
            xcrd[0][iat][icrd] = xf[iat][icrd];
        }
        rotvec(xcrd[0][iat], xcrd[0][iat], system->lavec_s);
    }

    for (int isize = -nsize[0]; isize <= nsize[0]; ++isize) {
        for (int jsize = -nsize[1]; jsize <= nsize[1]; ++jsize) {
            for (int ksize = -nsize[2]; ksize <= nsize[2]; ++ksize) {
                if (isize == 0 && jsize == 0 && ksize == 0) continue;
                ++icell;
                for (iat = 0; iat < nat; ++iat) {
                    xcrd[icell][iat][0] = xf[iat][0] + static_cast<double>(isize);
                    xcrd[icell][iat][1] = xf[iat][1] + static_cast<double>(jsize);
                    xcrd[icell][iat][2] = xf[iat][2] + static_cast<double>(ksize);
                    rotvec(xcrd[icell][iat], xcrd[icell][iat], system->lavec_s);
                }

            }
        }
    }

    for (iat = 0; iat < nat; ++iat) {
        for (jat = 0; jat < nat; ++jat) {
            for (icell = 0; icell < ncell; ++icell) {
                dist_tmp = std::sqrt(std::pow(xcrd[0][iat][0] - xcrd[icell][jat][0], 2.0)
                                     + std::pow(xcrd[0][iat][1] - xcrd[icell][jat][1], 2.0)
                                     + std::pow(xcrd[0][iat][2] - xcrd[icell][jat][2], 2.0));

                distall_ewald[iat][jat].emplace_back(icell, dist_tmp);
            }
            std::sort(distall_ewald[iat][jat].begin(), distall_ewald[iat][jat].end());
        }
    }
    double dist_hold = -1.0;
    for (iat = 0; iat < nat; ++iat) {
        for (jat = 0; jat < nat; ++jat) {
            multiplicity[iat][jat] = 0;

            for (auto it = distall_ewald[iat][jat].begin();
                 it != distall_ewald[iat][jat].end(); ++it) {
                if (it == distall_ewald[iat][jat].begin()) dist_hold = (*it).dist;
                dist_tmp = (*it).dist;
                if (std::abs(dist_tmp - dist_hold) < 1.0e-3) multiplicity[iat][jat] += 1;
            }

        }
    }
    deallocate(xcrd);
}

void Ewald::compute_ewald_fcs()
{
    int i, j;
    int iat, jat;
    int icrd, jcrd;
    int atm_s;
    int nat = system->nat;
    int natmin = system->natmin;
    double **fcs_ewald;
    double **fc_ewald_short, **fc_ewald_long;
    double **fcs_total, **fcs_other;
    std::string file_fcs_ewald = input->job_title + ".fc2_ewald";

    if (mympi->my_rank == 0) {
        std::cout << " Calculating long-range (dipole-dipole) FCs in the supercell ...";
    }

    allocate(fcs_ewald, 3 * natmin, 3 * nat);
    allocate(fc_ewald_short, 3, 3);
    allocate(fc_ewald_long, 3, 3);

    for (iat = 0; iat < natmin; ++iat) {
        atm_s = system->map_p2s[iat][0];
        for (jat = 0; jat < nat; ++jat) {
            calc_short_term_ewald_fcs(atm_s, jat, fc_ewald_short);
            calc_long_term_ewald_fcs(atm_s, jat, fc_ewald_long);

            for (icrd = 0; icrd < 3; ++icrd) {
                for (jcrd = 0; jcrd < 3; ++jcrd) {
                    fcs_ewald[3 * iat + icrd][3 * jat + jcrd]
                            = fc_ewald_short[icrd][jcrd] + fc_ewald_long[icrd][jcrd];
                }
            }
        }
    }

    deallocate(fc_ewald_short);
    deallocate(fc_ewald_long);

    allocate(fcs_total, 3 * natmin, 3 * nat);
    allocate(fcs_other, 3 * natmin, 3 * nat);

    for (i = 0; i < 3 * natmin; ++i) {
        for (j = 0; j < 3 * nat; ++j) {
            fcs_total[i][j] = 0.0;
        }
    }

    for (const auto &it: fcs_phonon->fc2_ext) {
        fcs_total[3 * it.atm1 + it.xyz1][3 * it.atm2 + it.xyz2] += it.fcs_val;
    }

    for (i = 0; i < 3 * natmin; ++i) {
        for (j = 0; j < 3 * nat; ++j) {
            fcs_other[i][j] = fcs_total[i][j] - fcs_ewald[i][j];
        }
    }
    FcsClassExtent fcext_tmp;

    for (iat = 0; iat < natmin; ++iat) {
        atm_s = system->map_p2s[iat][0];

        for (icrd = 0; icrd < 3; ++icrd) {
            for (jat = 0; jat < nat; ++jat) {
                for (jcrd = 0; jcrd < 3; ++jcrd) {

                    fcext_tmp.atm1 = iat;
                    fcext_tmp.xyz1 = icrd;
                    fcext_tmp.atm2 = jat;
                    fcext_tmp.xyz2 = jcrd;

                    int nmulti = multiplicity[atm_s][jat];
                    fcext_tmp.fcs_val = fcs_other[3 * iat + icrd][3 * jat + jcrd] / static_cast<double>(nmulti);

                    if (std::abs(fcext_tmp.fcs_val) > eps15) {

                        for (int icell = 0; icell < nmulti; ++icell) {
                            fcext_tmp.cell_s = distall_ewald[atm_s][jat][icell].cell;
                            fc2_without_dipole.push_back(fcext_tmp);
                        }
                    }
                }
            }
        }
    }

    if (mympi->my_rank == 0) {
        if (print_fc2_ewald) {

            std::ofstream ofs_fcs_ewald;

            ofs_fcs_ewald.open(file_fcs_ewald.c_str(), std::ios::out);
            if (!ofs_fcs_ewald) exit("compute_ewald_fcs", "cannot open file PREFIX.fcs_ewald");

            ofs_fcs_ewald << "# Harmonic force constants" << std::endl;
            ofs_fcs_ewald << "# atom1, xyz1, atom2, xyz2, fc2 original, fc2 dipole-dipole, fc2_orig - fc2_dipole" << std
            ::endl;

            for (iat = 0; iat < natmin; ++iat) {
                for (icrd = 0; icrd < 3; ++icrd) {
                    for (jat = 0; jat < nat; ++jat) {
                        for (jcrd = 0; jcrd < 3; ++jcrd) {
                            ofs_fcs_ewald << std::setw(5) << iat + 1;
                            ofs_fcs_ewald << std::setw(5) << icrd + 1;
                            ofs_fcs_ewald << std::setw(5) << jat + 1;
                            ofs_fcs_ewald << std::setw(5) << jcrd + 1;
                            ofs_fcs_ewald << std::setw(15) << fcs_total[3 * iat + icrd][3 * jat + jcrd];
                            ofs_fcs_ewald << std::setw(15) << fcs_ewald[3 * iat + icrd][3 * jat + jcrd];
                            ofs_fcs_ewald << std::setw(15) << fcs_other[3 * iat + icrd][3 * jat + jcrd];
                            ofs_fcs_ewald << std::endl;
                        }
                    }
                }
            }
            ofs_fcs_ewald.close();
        }
    }

    deallocate(fcs_ewald);
    deallocate(fcs_total);
    deallocate(fcs_other);

    if (mympi->my_rank == 0) {
        std::cout << " done." << std::endl;
        if (print_fc2_ewald) {
            std::cout << std::endl;
            std::cout << " FC2_EWALD = 1: Dipole-dipole and short-ranged components of harmonic " << std::endl;
            std::cout << "                force constants are printed in " << file_fcs_ewald << std::endl;
        }
    }
}

void Ewald::compute_ewald_fcs2()
{
    int i, j;
    int iat, jat;
    int icrd, jcrd;
    int atm_s;
    int nat = system->nat;
    int natmin = system->natmin;
    double **fcs_ewald;
    double **fc_ewald_short, **fc_ewald_long;
    double **fcs_total, **fcs_other;
    std::string file_fcs_ewald = input->job_title + ".fc2_ewald";

    if (mympi->my_rank == 0) {
        std::cout << " Calculating long-range (dipole-dipole) FCs in the supercell ...";
    }

    std::vector<std::vector<double>> k_commensurate;
    kpoint->get_commensurate_kpoints(system->lavec_s,
                                     system->lavec_p,
                                     k_commensurate);

    allocate(fcs_ewald, 3 * natmin, 3 * nat);
    allocate(fc_ewald_short, 3, 3);
    allocate(fc_ewald_long, 3, 3);

    for (iat = 0; iat < natmin; ++iat) {
        atm_s = system->map_p2s[iat][0];
        for (jat = 0; jat < nat; ++jat) {
            calc_short_term_ewald_fcs(atm_s, jat, fc_ewald_short);
            calc_long_term_ewald_fcs(atm_s, jat, fc_ewald_long);

            for (icrd = 0; icrd < 3; ++icrd) {
                for (jcrd = 0; jcrd < 3; ++jcrd) {
                    fcs_ewald[3 * iat + icrd][3 * jat + jcrd]
                            = fc_ewald_short[icrd][jcrd] + fc_ewald_long[icrd][jcrd];
                }
            }
        }
    }

    deallocate(fc_ewald_short);
    deallocate(fc_ewald_long);

    allocate(fcs_total, 3 * natmin, 3 * nat);
    allocate(fcs_other, 3 * natmin, 3 * nat);

    for (i = 0; i < 3 * natmin; ++i) {
        for (j = 0; j < 3 * nat; ++j) {
            fcs_total[i][j] = 0.0;
        }
    }

    for (const auto &it: fcs_phonon->fc2_ext) {
        fcs_total[3 * it.atm1 + it.xyz1][3 * it.atm2 + it.xyz2] += it.fcs_val;
    }

    for (i = 0; i < 3 * natmin; ++i) {
        for (j = 0; j < 3 * nat; ++j) {
            fcs_other[i][j] = fcs_total[i][j] - fcs_ewald[i][j];
        }
    }
    FcsClassExtent fcext_tmp;

    for (iat = 0; iat < natmin; ++iat) {
        atm_s = system->map_p2s[iat][0];

        for (icrd = 0; icrd < 3; ++icrd) {
            for (jat = 0; jat < nat; ++jat) {
                for (jcrd = 0; jcrd < 3; ++jcrd) {

                    fcext_tmp.atm1 = iat;
                    fcext_tmp.xyz1 = icrd;
                    fcext_tmp.atm2 = jat;
                    fcext_tmp.xyz2 = jcrd;

                    int nmulti = multiplicity[atm_s][jat];
                    fcext_tmp.fcs_val = fcs_other[3 * iat + icrd][3 * jat + jcrd] / static_cast<double>(nmulti);

                    if (std::abs(fcext_tmp.fcs_val) > eps15) {

                        for (int icell = 0; icell < nmulti; ++icell) {
                            fcext_tmp.cell_s = distall_ewald[atm_s][jat][icell].cell;
                            fc2_without_dipole.push_back(fcext_tmp);
                        }
                    }
                }
            }
        }
    }

    if (mympi->my_rank == 0) {
        if (print_fc2_ewald) {

            std::ofstream ofs_fcs_ewald;

            ofs_fcs_ewald.open(file_fcs_ewald.c_str(), std::ios::out);
            if (!ofs_fcs_ewald) exit("compute_ewald_fcs", "cannot open file PREFIX.fcs_ewald");

            ofs_fcs_ewald << "# Harmonic force constants" << std::endl;
            ofs_fcs_ewald << "# atom1, xyz1, atom2, xyz2, fc2 original, fc2 dipole-dipole, fc2_orig - fc2_dipole" << std
            ::endl;

            for (iat = 0; iat < natmin; ++iat) {
                for (icrd = 0; icrd < 3; ++icrd) {
                    for (jat = 0; jat < nat; ++jat) {
                        for (jcrd = 0; jcrd < 3; ++jcrd) {
                            ofs_fcs_ewald << std::setw(5) << iat + 1;
                            ofs_fcs_ewald << std::setw(5) << icrd + 1;
                            ofs_fcs_ewald << std::setw(5) << jat + 1;
                            ofs_fcs_ewald << std::setw(5) << jcrd + 1;
                            ofs_fcs_ewald << std::setw(15) << fcs_total[3 * iat + icrd][3 * jat + jcrd];
                            ofs_fcs_ewald << std::setw(15) << fcs_ewald[3 * iat + icrd][3 * jat + jcrd];
                            ofs_fcs_ewald << std::setw(15) << fcs_other[3 * iat + icrd][3 * jat + jcrd];
                            ofs_fcs_ewald << std::endl;
                        }
                    }
                }
            }
            ofs_fcs_ewald.close();
        }
    }

    deallocate(fcs_ewald);
    deallocate(fcs_total);
    deallocate(fcs_other);

    if (mympi->my_rank == 0) {
        std::cout << " done." << std::endl;
        if (print_fc2_ewald) {
            std::cout << std::endl;
            std::cout << " FC2_EWALD = 1: Dipole-dipole and short-ranged components of harmonic " << std::endl;
            std::cout << "                force constants are printed in " << file_fcs_ewald << std::endl;
        }
    }
}

void Ewald::calc_short_term_ewald_fcs(const int iat,
                                      const int jat,
                                      double **fc_l_out)
{
    // Real lattice sum part for FCs
    // iat : atom index in the supercell (should be in the center primitive cell)
    // jat : atom index in the supercell

    int i;
    int icrd, jcrd;
    int icell, jcell, kcell;
    int kat, kkd;
    double xnorm;
    double x_tmp[3], trans[3];
    std::vector<std::vector<double>> func_L(3, std::vector<double>(3, 0.0));

    for (icrd = 0; icrd < 3; ++icrd) {
        for (jcrd = 0; jcrd < 3; ++jcrd) {
            fc_l_out[icrd][jcrd] = 0.0;
        }
    }

    if (iat == jat) {
        // (k,l) = (k',l')

        for (icell = -nl_sub[0]; icell <= nl_sub[0]; ++icell) {
            for (jcell = -nl_sub[1]; jcell <= nl_sub[1]; ++jcell) {
                for (kcell = -nl_sub[2]; kcell <= nl_sub[2]; ++kcell) {

                    if (icell == 0 && jcell == 0 && kcell == 0) {
                        // l'' = l = 0

                        for (kat = 0; kat < system->nat; ++kat) {

                            if (kat == iat) continue; // k'' = k

                            kkd = system->map_s2p[kat].atom_num;
                            for (i = 0; i < 3; ++i) {
                                x_tmp[i] = system->xr_s[iat][i] - system->xr_s[kat][i];
                            }
                            rotvec(x_tmp, x_tmp, system->lavec_s);
                            xnorm = std::sqrt(x_tmp[0] * x_tmp[0]
                                              + x_tmp[1] * x_tmp[1]
                                              + x_tmp[2] * x_tmp[2]);

                            if (xnorm < Lmax_sub) {

                                calc_realspace_sum(iat, kat, x_tmp, lambda_sub, func_L);

                                if (force_permutation_sym) {
                                    for (icrd = 0; icrd < 3; ++icrd) {
                                        for (jcrd = 0; jcrd < 3; ++jcrd) {
                                            fc_l_out[icrd][jcrd] += 0.5 * (func_L[icrd][jcrd]
                                                                           + func_L[jcrd][icrd]);
                                        }
                                    }
                                } else {
                                    for (icrd = 0; icrd < 3; ++icrd) {
                                        for (jcrd = 0; jcrd < 3; ++jcrd) {
                                            fc_l_out[icrd][jcrd] += func_L[icrd][jcrd];
                                        }
                                    }
                                }
                            }
                        }

                    } else {
                        // l'' != 0

                        trans[0] = static_cast<double>(icell);
                        trans[1] = static_cast<double>(jcell);
                        trans[2] = static_cast<double>(kcell);
                        rotvec(trans, trans, system->lavec_s);

                        for (kat = 0; kat < system->nat; ++kat) {
                            kkd = system->map_s2p[kat].atom_num;
                            for (i = 0; i < 3; ++i) {
                                x_tmp[i] = system->xr_s[iat][i] - system->xr_s[kat][i];
                            }
                            rotvec(x_tmp, x_tmp, system->lavec_s);
                            for (i = 0; i < 3; ++i) {
                                x_tmp[i] -= trans[i];
                            }
                            xnorm = std::sqrt(x_tmp[0] * x_tmp[0]
                                              + x_tmp[1] * x_tmp[1]
                                              + x_tmp[2] * x_tmp[2]);

                            if (xnorm < Lmax_sub) {
                                calc_realspace_sum(iat, kat, x_tmp, lambda_sub, func_L);

                                if (force_permutation_sym) {
                                    for (icrd = 0; icrd < 3; ++icrd) {
                                        for (jcrd = 0; jcrd < 3; ++jcrd) {
                                            fc_l_out[icrd][jcrd] += 0.5 * (func_L[icrd][jcrd]
                                                                           + func_L[jcrd][icrd]);
                                        }
                                    }
                                } else {
                                    for (icrd = 0; icrd < 3; ++icrd) {
                                        for (jcrd = 0; jcrd < 3; ++jcrd) {
                                            fc_l_out[icrd][jcrd] += func_L[icrd][jcrd];
                                        }
                                    }
                                }
                            }
                        }

                        for (i = 0; i < 3; ++i) {
                            x_tmp[i] = system->xr_s[iat][i] - system->xr_s[jat][i];
                        }
                        rotvec(x_tmp, x_tmp, system->lavec_s);
                        for (i = 0; i < 3; ++i) {
                            x_tmp[i] -= trans[i];
                        }
                        xnorm = std::sqrt(x_tmp[0] * x_tmp[0]
                                          + x_tmp[1] * x_tmp[1]
                                          + x_tmp[2] * x_tmp[2]);

                        if (xnorm < Lmax_sub) {
                            calc_realspace_sum(iat, jat, x_tmp, lambda_sub, func_L);

                            for (icrd = 0; icrd < 3; ++icrd) {
                                for (jcrd = 0; jcrd < 3; ++jcrd) {
                                    fc_l_out[icrd][jcrd] -= func_L[icrd][jcrd];
                                }
                            }
                        }
                    }
                }
            }
        }

    } else {

        // case of i != j
        // (k,l) != (k',l')

        for (icell = -nl_sub[0]; icell <= nl_sub[0]; ++icell) {
            for (jcell = -nl_sub[1]; jcell <= nl_sub[1]; ++jcell) {
                for (kcell = -nl_sub[2]; kcell <= nl_sub[2]; ++kcell) {

                    trans[0] = static_cast<double>(icell);
                    trans[1] = static_cast<double>(jcell);
                    trans[2] = static_cast<double>(kcell);

                    for (i = 0; i < 3; ++i) {
                        x_tmp[i] = system->xr_s[iat][i] - system->xr_s[jat][i] - trans[i];
                    }
                    rotvec(x_tmp, x_tmp, system->lavec_s);
                    xnorm = std::sqrt(x_tmp[0] * x_tmp[0]
                                      + x_tmp[1] * x_tmp[1]
                                      + x_tmp[2] * x_tmp[2]);

                    if (xnorm < Lmax_sub) {
                        calc_realspace_sum(iat, jat, x_tmp, lambda_sub, func_L);

                        for (icrd = 0; icrd < 3; ++icrd) {
                            for (jcrd = 0; jcrd < 3; ++jcrd) {
                                fc_l_out[icrd][jcrd] -= func_L[icrd][jcrd];
                            }
                        }
                    }

                }
            }
        }

    }
}

void Ewald::calc_long_term_ewald_fcs(const int iat,
                                     const int jat,
                                     double **fc_g_out)
{
    // Reciprocal lattice sum part for FCs

    int i;
    int icrd, jcrd;
    int acrd, bcrd;
    double gnorm2;
    double x_tmp[3], g_tmp[3], epsilon_gvector[3];
    double common_tmp;

    for (icrd = 0; icrd < 3; ++icrd) {
        for (jcrd = 0; jcrd < 3; ++jcrd) {
            fc_g_out[icrd][jcrd] = 0.0;
        }
    }
    double volume = system->volume(system->lavec_s[0],
                                   system->lavec_s[1],
                                   system->lavec_s[2]);

    int ikd = system->map_s2p[iat].atom_num;
    int jkd = system->map_s2p[jat].atom_num;

    double factor = 4.0 * pi / volume;

    if (iat == jat) {

        for (const auto &it: G_vector_sub) {
            for (i = 0; i < 3; ++i) g_tmp[i] = it.vec[i];
            rotvec(epsilon_gvector, g_tmp, epsilon);
            gnorm2 = g_tmp[0] * epsilon_gvector[0]
                     + g_tmp[1] * epsilon_gvector[1]
                     + g_tmp[2] * epsilon_gvector[2];

            for (int kat = 0; kat < system->nat; ++kat) {
                int kkd = system->map_s2p[kat].atom_num;

                for (i = 0; i < 3; ++i) {
                    x_tmp[i] = system->xr_s[iat][i] - system->xr_s[kat][i];
                }
                rotvec(x_tmp, x_tmp, system->lavec_s);

                common_tmp = factor * std::exp(-0.25 * gnorm2 / std::pow(lambda_sub, 2.0)) / gnorm2
                             * std::cos(g_tmp[0] * x_tmp[0]
                                        + g_tmp[1] * x_tmp[1]
                                        + g_tmp[2] * x_tmp[2]);

                for (icrd = 0; icrd < 3; ++icrd) {
                    for (jcrd = 0; jcrd < 3; ++jcrd) {
                        for (acrd = 0; acrd < 3; ++acrd) {
                            for (bcrd = 0; bcrd < 3; ++bcrd) {
                                if (force_permutation_sym) {
                                    fc_g_out[icrd][jcrd] -= g_tmp[acrd] * g_tmp[bcrd] * common_tmp
                                                            *
                                                            (Born_charge[ikd][acrd][icrd] * Born_charge[kkd][bcrd][jcrd]
                                                             +
                                                             Born_charge[jkd][acrd][jcrd] *
                                                             Born_charge[kkd][bcrd][icrd]);
                                } else {
                                    fc_g_out[icrd][jcrd] -= g_tmp[acrd] * g_tmp[bcrd] * common_tmp * 2.0
                                                            * (Born_charge[ikd][acrd][icrd] *
                                                               Born_charge[kkd][bcrd][jcrd]);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    for (const auto &it: G_vector_sub) {
        for (i = 0; i < 3; ++i) g_tmp[i] = it.vec[i];
        rotvec(epsilon_gvector, g_tmp, epsilon);
        gnorm2 = g_tmp[0] * epsilon_gvector[0]
                 + g_tmp[1] * epsilon_gvector[1]
                 + g_tmp[2] * epsilon_gvector[2];

        for (i = 0; i < 3; ++i) {
            x_tmp[i] = system->xr_s[iat][i] - system->xr_s[jat][i];
        }
        rotvec(x_tmp, x_tmp, system->lavec_s);

        common_tmp = 2.0 * factor * std::exp(-0.25 * gnorm2 / std::pow(lambda_sub, 2.0)) / gnorm2
                     * std::cos(g_tmp[0] * x_tmp[0] + g_tmp[1] * x_tmp[1] + g_tmp[2] * x_tmp[2]);

        for (icrd = 0; icrd < 3; ++icrd) {
            for (jcrd = 0; jcrd < 3; ++jcrd) {
                for (acrd = 0; acrd < 3; ++acrd) {
                    for (bcrd = 0; bcrd < 3; ++bcrd) {
                        fc_g_out[icrd][jcrd] += g_tmp[acrd] * g_tmp[bcrd] * common_tmp
                                                * Born_charge[ikd][acrd][icrd] * Born_charge[jkd][bcrd][jcrd];
                    }
                }
            }
        }
    }
}

void Ewald::add_longrange_matrix(const double *xk_in,
                                 const double *kvec_in,
                                 std::complex<double> **dymat_k_out)
{
    int icrd, jcrd, iat, jat;
    int natmin = system->natmin;
    int neval = 3 * system->natmin;
    double xk[3];
    std::complex<double> **dymat_tmp_l, **dymat_tmp_g;

    allocate(dymat_tmp_l, 3, 3);
    allocate(dymat_tmp_g, 3, 3);

    rotvec(xk, xk_in, system->rlavec_p, 'T');

    for (int i = 0; i < neval; ++i) {
        for (int j = 0; j < neval; ++j) {
            dymat_k_out[i][j] = std::complex<double>(0.0, 0.0);
        }
    }

    for (iat = 0; iat < natmin; ++iat) {
        for (jat = 0; jat < natmin; ++jat) {
            calc_short_term_dynamical_matrix(iat, jat, xk, dymat_tmp_l);
            calc_long_term_dynamical_matrix(iat, jat, xk, kvec_in, dymat_tmp_g);
            for (icrd = 0; icrd < 3; ++icrd) {
                for (jcrd = 0; jcrd < 3; ++jcrd) {
                    dymat_k_out[3 * iat + icrd][3 * jat + jcrd] = dymat_tmp_l[icrd][jcrd]
                                                                  + dymat_tmp_g[icrd][jcrd];
                }
            }

        }
    }
    deallocate(dymat_tmp_l);
    deallocate(dymat_tmp_g);


//    for (iat = 0; iat < natmin; ++iat) {
//        for (icrd = 0; icrd < 3; ++icrd) {
//            for (jat = 0; jat < natmin; ++jat) {
//                for (jcrd = 0; jcrd < 3; ++jcrd) {
//
//                    // Hermiticity
//                    std::complex<double> check = dymat_k_out[3 * iat + icrd][3 * jat + jcrd]
//                                                 - std::conj(dymat_k_out[3 * jat + jcrd][3 * iat + icrd]);
//                    if (std::abs(check) > eps10) {
//                        std::cout << std::endl;
//                        exit("add_longrange_matrix",
//                                    "Hermiticity of Dynamical matrix is broken.");
//                    }
//                }
//            }
//        }
//    }
}

void Ewald::calc_short_term_dynamical_matrix(const int iat,
                                             const int jat,
                                             double *xk_in,
                                             std::complex<double> **mat_out)
{
    // Real lattice sum part for a dynamical matrix
    // iat : atom index in the primitive cell
    // jat : atom index in the primitive cell

    int i;
    int icrd, jcrd, kat;
    int atm_s3;
    int icell, jcell, kcell;
    double xnorm, phase;
    double x_tmp[3], trans[3];
    std::complex<double> im(0.0, 1.0);
    std::vector<std::vector<double>> func_L(3, std::vector<double>(3, 0.0));

    // Substitute quantities into variables
    int atm_s1 = system->map_p2s[iat][0];
    int atm_s2 = system->map_p2s[jat][0];

    for (i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            mat_out[i][j] = std::complex<double>(0.0, 0.0);
        }
    }

    if (iat == jat) {
        // k = k'

        for (icell = -nl[0]; icell <= nl[0]; ++icell) {
            for (jcell = -nl[1]; jcell <= nl[1]; ++jcell) {
                for (kcell = -nl[2]; kcell <= nl[2]; ++kcell) {

                    trans[0] = static_cast<double>(icell);
                    trans[1] = static_cast<double>(jcell);
                    trans[2] = static_cast<double>(kcell);
                    rotvec(trans, trans, system->lavec_p);

                    // Lattice vector = 0
                    if (icell == 0 && jcell == 0 && kcell == 0) {

                        for (kat = 0; kat < system->natmin; ++kat) {
                            if (kat == iat) continue;

                            atm_s3 = system->map_p2s[kat][0];
                            for (i = 0; i < 3; ++i) {
                                x_tmp[i] = system->xr_s[atm_s1][i] - system->xr_s[atm_s3][i];
                            }
                            rotvec(x_tmp, x_tmp, system->lavec_s);
                            xnorm = std::sqrt(x_tmp[0] * x_tmp[0]
                                              + x_tmp[1] * x_tmp[1]
                                              + x_tmp[2] * x_tmp[2]);

                            if (xnorm < Lmax) {
                                calc_realspace_sum(atm_s1, atm_s3, x_tmp, lambda, func_L);

                                if (force_permutation_sym) {
                                    for (icrd = 0; icrd < 3; ++icrd) {
                                        for (jcrd = 0; jcrd < 3; ++jcrd) {
                                            mat_out[icrd][jcrd] += 0.5 * (func_L[icrd][jcrd]
                                                                          + func_L[jcrd][icrd]);
                                        }
                                    }
                                } else {
                                    for (icrd = 0; icrd < 3; ++icrd) {
                                        for (jcrd = 0; jcrd < 3; ++jcrd) {
                                            mat_out[icrd][jcrd] += func_L[icrd][jcrd];
                                        }
                                    }
                                }
                            }
                        }

                    } else {

                        for (kat = 0; kat < system->natmin; ++kat) {
                            atm_s3 = system->map_p2s[kat][0];
                            for (i = 0; i < 3; ++i) {
                                x_tmp[i] = system->xr_s[atm_s1][i] - system->xr_s[atm_s3][i];
                            }
                            rotvec(x_tmp, x_tmp, system->lavec_s);
                            for (i = 0; i < 3; ++i) {
                                x_tmp[i] -= trans[i];
                            }
                            xnorm = std::sqrt(x_tmp[0] * x_tmp[0]
                                              + x_tmp[1] * x_tmp[1]
                                              + x_tmp[2] * x_tmp[2]);

                            if (xnorm < Lmax) {

                                calc_realspace_sum(atm_s1, atm_s3, x_tmp, lambda, func_L);

                                if (force_permutation_sym) {
                                    for (icrd = 0; icrd < 3; ++icrd) {
                                        for (jcrd = 0; jcrd < 3; ++jcrd) {
                                            mat_out[icrd][jcrd] += 0.5 * (func_L[icrd][jcrd]
                                                                          + func_L[jcrd][icrd]);
                                        }
                                    }
                                } else {
                                    for (icrd = 0; icrd < 3; ++icrd) {
                                        for (jcrd = 0; jcrd < 3; ++jcrd) {
                                            mat_out[icrd][jcrd] += func_L[icrd][jcrd];
                                        }
                                    }
                                }
                            }
                        }

                        for (i = 0; i < 3; ++i) {
                            x_tmp[i] = system->xr_s[atm_s1][i] - system->xr_s[atm_s2][i];
                        }
                        rotvec(x_tmp, x_tmp, system->lavec_s);
                        for (i = 0; i < 3; ++i) {
                            x_tmp[i] -= trans[i];
                        }
                        xnorm = std::sqrt(x_tmp[0] * x_tmp[0]
                                          + x_tmp[1] * x_tmp[1]
                                          + x_tmp[2] * x_tmp[2]);

                        if (xnorm < Lmax) {
                            calc_realspace_sum(atm_s1, atm_s2, x_tmp, lambda, func_L);
                            phase = xk_in[0] * trans[0] + xk_in[1] * trans[1] + xk_in[2] * trans[2];

                            for (icrd = 0; icrd < 3; ++icrd) {
                                for (jcrd = 0; jcrd < 3; ++jcrd) {
                                    mat_out[icrd][jcrd] -= func_L[icrd][jcrd] * std::exp(im * phase);
                                }
                            }
                        }
                    }
                }
            }
        }

    } else {

        for (icell = -nl[0]; icell <= nl[0]; ++icell) {
            for (jcell = -nl[1]; jcell <= nl[1]; ++jcell) {
                for (kcell = -nl[2]; kcell <= nl[2]; ++kcell) {

                    trans[0] = static_cast<double>(icell);
                    trans[1] = static_cast<double>(jcell);
                    trans[2] = static_cast<double>(kcell);
                    rotvec(trans, trans, system->lavec_p);

                    for (i = 0; i < 3; ++i) {
                        x_tmp[i] = system->xr_s[atm_s1][i] - system->xr_s[atm_s2][i];
                    }
                    rotvec(x_tmp, x_tmp, system->lavec_s);
                    for (i = 0; i < 3; ++i) {
                        x_tmp[i] -= trans[i];
                    }
                    xnorm = std::sqrt(x_tmp[0] * x_tmp[0]
                                      + x_tmp[1] * x_tmp[1]
                                      + x_tmp[2] * x_tmp[2]);

                    if (xnorm < Lmax) {
                        calc_realspace_sum(atm_s1, atm_s2, x_tmp, lambda, func_L);
                        phase = xk_in[0] * trans[0] + xk_in[1] * trans[1] + xk_in[2] * trans[2];

                        for (icrd = 0; icrd < 3; ++icrd) {
                            for (jcrd = 0; jcrd < 3; ++jcrd) {
                                mat_out[icrd][jcrd] -= func_L[icrd][jcrd] * std::exp(im * phase);
                            }
                        }
                    }
                }
            }
        }

    }

    const auto mi = system->mass[atm_s1];
    const auto mj = system->mass[atm_s2];
    for (icrd = 0; icrd < 3; ++icrd) {
        for (jcrd = 0; jcrd < 3; ++jcrd) {
            mat_out[icrd][jcrd] /= std::sqrt(mi * mj);
        }
    }
}

void Ewald::calc_long_term_dynamical_matrix(const int iat,
                                            const int jat,
                                            const double *xk_in,
                                            const double *kvec_in,
                                            std::complex<double> **mat_out)
{
    // Real lattice sum part for a dynamical matrix

    int i, j;
    int icrd, jcrd, acrd, bcrd;
    double vec[3], e_kvec[3];
    std::complex<double> im(0.0, 1.0);
    double tmp;

    int atm_s1 = system->map_p2s[iat][0];
    int atm_s2 = system->map_p2s[jat][0];
    double mi = system->mass[atm_s1];
    double mj = system->mass[atm_s2];
    double vol_p = system->volume_p;
    for (i = 0; i < 3; ++i) {
        vec[i] = system->xr_s[atm_s1][i] - system->xr_s[atm_s2][i];
    }
    rotvec(vec, vec, system->lavec_s);
    double phase = xk_in[0] * vec[0] + xk_in[1] * vec[1] + xk_in[2] * vec[2];
    rotvec(e_kvec, xk_in, epsilon);

    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            mat_out[i][j] = std::complex<double>(0.0, 0.0);
        }
    }
    double kd = xk_in[0] * e_kvec[0] + xk_in[1] * e_kvec[1] + xk_in[2] * e_kvec[2];

    std::complex<double> exp_phase = std::exp(im * phase);

    if (std::sqrt(kd) > eps10) {

        // constant analytic term

        for (icrd = 0; icrd < 3; ++icrd) {
            for (jcrd = 0; jcrd < 3; ++jcrd) {
                tmp = 0.0;

                for (acrd = 0; acrd < 3; ++acrd) {
                    for (bcrd = 0; bcrd < 3; ++bcrd) {
                        tmp += xk_in[acrd] * xk_in[bcrd]
                               * Born_charge[iat][acrd][icrd] * Born_charge[jat][bcrd][jcrd];
                    }
                }
                mat_out[icrd][jcrd] += 2.0 * tmp / kd * exp_phase
                                       * std::exp(-0.25 * kd / std::pow(lambda, 2.0));
            }
        }

    } else {

        // Treat non-analytic term

        double kdirec[3], e_kdirec[3];
        for (i = 0; i < 3; ++i) kdirec[i] = kvec_in[i];
        rotvec(e_kdirec, kdirec, epsilon);
        double norm = kdirec[0] * e_kdirec[0] + kdirec[1] * e_kdirec[1] + kdirec[2] * e_kdirec[2];

        if (norm > eps) {
            for (icrd = 0; icrd < 3; ++icrd) {
                for (jcrd = 0; jcrd < 3; ++jcrd) {
                    tmp = 0.0;

                    for (acrd = 0; acrd < 3; ++acrd) {
                        for (bcrd = 0; bcrd < 3; ++bcrd) {
                            tmp += kdirec[acrd] * kdirec[bcrd]
                                   * Born_charge[iat][acrd][icrd] * Born_charge[jat][bcrd][jcrd];
                        }
                    }
                    mat_out[icrd][jcrd] += 2.0 * tmp / norm * exp_phase;
                }
            }
        }
    }

    double g[3], gk[3], vecl[3], g_tmp[3], gk_tmp[3];
    double common;

    for (auto &it: G_vector) {
        for (int l = 0; l < 3; ++l) {
            g[l] = it.vec[l];
            gk[l] = g[l] + xk_in[l];
        }

        if (iat == jat) {

            rotvec(g_tmp, g, epsilon);
            double gd = g[0] * g_tmp[0] + g[1] * g_tmp[1] + g[2] * g_tmp[2];
            common = std::exp(-0.25 * gd / std::pow(lambda, 2.0)) / gd;

            for (int kat = 0; kat < system->natmin; ++kat) {
                int atm_s3 = system->map_p2s[kat][0];

                for (i = 0; i < 3; ++i) {
                    vecl[i] = system->xr_s[atm_s1][i] - system->xr_s[atm_s3][i];
                }
                rotvec(vecl, vecl, system->lavec_s);
                double phase_g1 = g[0] * vecl[0] + g[1] * vecl[1] + g[2] * vecl[2];
                exp_phase = std::exp(im * phase_g1);

                for (icrd = 0; icrd < 3; ++icrd) {
                    for (jcrd = 0; jcrd < 3; ++jcrd) {
                        tmp = 0.0;

                        for (acrd = 0; acrd < 3; ++acrd) {
                            for (bcrd = 0; bcrd < 3; ++bcrd) {
                                if (force_permutation_sym) {
                                    tmp += g[acrd] * g[bcrd]
                                           * (Born_charge[iat][acrd][icrd] * Born_charge[kat][bcrd][jcrd]
                                              + Born_charge[jat][acrd][jcrd] * Born_charge[kat][bcrd][icrd]);
                                } else {
                                    tmp += g[acrd] * g[bcrd] * 2.0
                                           * (Born_charge[iat][acrd][icrd] * Born_charge[kat][bcrd][jcrd]);
                                }
                            }
                        }
                        mat_out[icrd][jcrd] -= tmp * common * exp_phase;
                    }
                }
            }
        }

        rotvec(gk_tmp, gk, epsilon);
        double gkd = gk[0] * gk_tmp[0] + gk[1] * gk_tmp[1] + gk[2] * gk_tmp[2];
        double phase_g2 = gk[0] * vec[0] + gk[1] * vec[1] + gk[2] * vec[2];

        common = 2.0 * std::exp(-0.25 * gkd / std::pow(lambda, 2.0)) / gkd;
        exp_phase = std::exp(im * phase_g2);

        for (icrd = 0; icrd < 3; ++icrd) {
            for (jcrd = 0; jcrd < 3; ++jcrd) {
                tmp = 0.0;

                for (acrd = 0; acrd < 3; ++acrd) {
                    for (bcrd = 0; bcrd < 3; ++bcrd) {
                        tmp += gk[acrd] * gk[bcrd]
                               * Born_charge[iat][acrd][icrd] * Born_charge[jat][bcrd][jcrd];
                    }
                }
                mat_out[icrd][jcrd] += tmp * common * exp_phase;
            }
        }
    }

    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            mat_out[i][j] *= 4.0 * pi / (vol_p * std::sqrt(mi * mj));
        }
    }
}

void Ewald::calc_realspace_sum(const int iat,
                               const int jat,
                               const double xdist[3],
                               const double lambda_in,
                               std::vector<std::vector<double>> &ret)
{
    // iat : atom index in the supercell
    // jat : atom index in the supercell
    // xdist: distance between the two atomic sites
    unsigned int icrd, jcrd, acrd, bcrd;
    double tmp;
    double **hmat_tmp;
    const double lambda3 = std::pow(lambda_in, 3.0);
    allocate(hmat_tmp, 3, 3);

    calc_anisotropic_hmat(lambda_in, xdist, hmat_tmp);

    const auto ikd = system->map_s2p[iat].atom_num;
    const auto jkd = system->map_s2p[jat].atom_num;

    for (icrd = 0; icrd < 3; ++icrd) {
        for (jcrd = 0; jcrd < 3; ++jcrd) {
            ret[icrd][jcrd] = 0.0;
        }
    }

    for (icrd = 0; icrd < 3; ++icrd) {
        for (jcrd = 0; jcrd < 3; ++jcrd) {
            tmp = 0.0;
            for (acrd = 0; acrd < 3; ++acrd) {
                for (bcrd = 0; bcrd < 3; ++bcrd) {
                    tmp += hmat_tmp[acrd][bcrd]
                           * Born_charge[ikd][acrd][icrd]
                           * Born_charge[jkd][bcrd][jcrd];
                }
            }
            ret[icrd][jcrd] = 2.0 * tmp * lambda3;
        }
    }
    deallocate(hmat_tmp);
}

void Ewald::calc_anisotropic_hmat(const double lambda_in,
                                  const double *x,
                                  double **hmat_out)
{
    // Compute H_ab(0\kappa;\ell'\kappa')
    int icrd, jcrd;
    double common_tmp[2];
    double x_tmp[3], y_tmp[3];

    for (icrd = 0; icrd < 3; ++icrd) {
        for (jcrd = 0; jcrd < 3; ++jcrd) {
            hmat_out[icrd][jcrd] = 0.0;
        }
    }

    for (int i = 0; i < 3; ++i) {
        y_tmp[i] = x[i] * lambda_in;
    }

    rotvec(x_tmp, y_tmp, epsilon_inv);

    double yd = std::sqrt(x_tmp[0] * y_tmp[0] + x_tmp[1] * y_tmp[1] + x_tmp[2] * y_tmp[2]);
    if (yd == 0.0) {
        exit("ewald->calc_anisotropic_hmat", "components of hmat diverge.");
    }
    double yd_inv = 1.0 / yd;
    double yd2 = std::pow(yd, 2.0);
    double yd2_inv = yd_inv * yd_inv;
    double erfc_y = boost::math::erfc(yd);
    double exp_y2 = std::exp(-yd2);
    double two_over_sqrtpi = 2.0 / std::sqrt(pi);

    common_tmp[0] = (3.0 * yd_inv * yd2_inv * erfc_y + two_over_sqrtpi * (3.0 * yd2_inv + 2.0) * exp_y2)
                    * yd2_inv / std::sqrt(det_epsilon);
    common_tmp[1] = (yd_inv * yd2_inv * erfc_y + two_over_sqrtpi * yd2_inv * exp_y2) / std::sqrt(det_epsilon);

    for (icrd = 0; icrd < 3; ++icrd) {
        for (jcrd = 0; jcrd < 3; ++jcrd) {
            hmat_out[icrd][jcrd] = x_tmp[icrd] * x_tmp[jcrd] * common_tmp[0]
                                   - epsilon_inv[icrd][jcrd] * common_tmp[1];
        }
    }
}
