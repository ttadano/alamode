/*
fcs_phonon.cpp

Copyright (c) 2014, 2015, 2016 Terumasa Tadano

This file is distributed under the terms of the MIT license.
Please see the file 'LICENCE.txt' in the root directory 
or http://opensource.org/licenses/mit-license.php for information.
*/

#include "mpi_common.h"
#include "fcs_phonon.h"
#include "constants.h"
#include "dynamical.h"
#include "error.h"
#include "gruneisen.h"
#include "memory.h"
#include "phonons.h"
#include "anharmonic_core.h"
#include "system.h"
#include "thermodynamics.h"
#include "mathfunctions.h"
#include <string>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include <Eigen/LU>

using namespace PHON_NS;

Fcs_phonon::Fcs_phonon(PHON *phon) : Pointers(phon)
{
    set_default_variables();
}

Fcs_phonon::~Fcs_phonon()
{
    deallocate_variables();
}

void Fcs_phonon::set_default_variables()
{
    maxorder = 0;
    file_fcs = "";
    file_fc2 = "";
    file_fc3 = "";
    file_fc4 = "";

    update_fc2 = false;
    force_constant_with_cell = nullptr;
}

void Fcs_phonon::deallocate_variables()
{
    if (force_constant_with_cell) {
        deallocate(force_constant_with_cell);
    }
}

void Fcs_phonon::setup(std::string mode)
{
    if (mympi->my_rank == 0) {
        std::cout << " Force constant" << std::endl;
        std::cout << " ==============" << std::endl << std::endl;
    }

    MPI_Bcast(&anharmonic_core->quartic_mode, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&gruneisen->print_gruneisen, 1, MPI_CXX_BOOL, 0, MPI_COMM_WORLD);
    MPI_Bcast(&thermodynamics->calc_FE_bubble, 1, MPI_CXX_BOOL, 0, MPI_COMM_WORLD);

    if (mode == "PHONONS") {
        require_cubic = false;
        require_quartic = false;
        maxorder = 1;

        if (gruneisen->print_gruneisen || thermodynamics->calc_FE_bubble) {
            require_cubic = true;
            maxorder = 2;
        }
        if (gruneisen->print_newfcs) {
            require_cubic = true;
            maxorder = 2;

            if (anharmonic_core->quartic_mode > 0) {
                require_quartic = true;
                maxorder = 3;
            }
        }

    } else if (mode == "RTA") {
        require_cubic = true;

        if (anharmonic_core->quartic_mode > 0) {
            maxorder = 3;
            require_quartic = true;
        } else {
            maxorder = 2;
            require_quartic = false;
        }
    } else if (mode == "SCPH") {
        require_cubic = true;
        require_quartic = true;
        maxorder = 3;
        anharmonic_core->quartic_mode = 1;
    }

    allocate(force_constant_with_cell, maxorder);

    if (mympi->my_rank == 0) {
        load_fc2_xml();
        load_fcs_xml();

        for (auto i = 0; i < maxorder; ++i) {
            std::cout << "  Number of non-zero IFCs for " << i + 2 << " order: ";
            std::cout << force_constant_with_cell[i].size() << std::endl;
        }
        std::cout << std::endl;

        std::cout << "  Maximum deviation from the translational invariance: " << std::endl;
        for (auto i = 0; i < maxorder; ++i) {
            const auto maxdev = examine_translational_invariance(i,
                                                                 system->get_supercell(i).number_of_atoms,
                                                                 system->get_primcell().number_of_atoms,
                                                                 system->get_mapping_super_alm(i).from_true_primitive,
                                                                 force_constant_with_cell[i]);
            std::cout << "   Order " << i + 2 << " : " << std::setw(12)
                      << std::scientific << maxdev << std::endl;
        }
        std::cout << std::endl;
    }

    MPI_Bcast_fc2_ext();
    MPI_Bcast_fcs_array(maxorder);
    replicate_force_constants(maxorder);

}

void Fcs_phonon::replicate_force_constants(const int maxorder_in)
{
    std::vector<FcsArrayWithCell> force_constant_replicate;
    std::vector<Eigen::Vector3d> relvecs;
    Eigen::Vector3d relvec_tmp;
    std::vector<std::vector<unsigned int>> map_trans;
    Eigen::Vector3d xshift, x0, x_shifted;
    Eigen::Vector3d xdiff;
    std::vector<unsigned int> map_now;

    for (auto order = 0; order < maxorder_in; ++order) {

        force_constant_replicate.clear();
        map_trans.clear();

        const auto map_tmp = system->get_mapping_super_alm(order);
        const auto cell_tmp = system->get_supercell(order);
        const auto ntran_tmp = map_tmp.from_true_primitive[0].size();

        for (auto itran = 0; itran < ntran_tmp; ++itran) {
            xshift = cell_tmp.x_fractional.row(map_tmp.from_true_primitive[0][itran])
                     - cell_tmp.x_fractional.row(map_tmp.from_true_primitive[0][0]);
            map_now.clear();

            for (auto iat = 0; iat < cell_tmp.number_of_atoms; ++iat) {
                auto kat = -1;
                for (auto jat = 0; jat < cell_tmp.number_of_atoms; ++jat) {
                    xdiff = cell_tmp.x_fractional.row(jat) - cell_tmp.x_fractional.row(iat);
                    xdiff = (xdiff + xshift).unaryExpr([](const double x) { return x - static_cast<double>(nint(x)); });
                    if (xdiff.norm() < eps6) {
                        kat = jat;
                        break;
                    }
                }
                if (kat == -1) {
                    exit("replicate_force_constants", "Equivalent atom cound not be found.");
                } else {
                    map_now.emplace_back(kat);
                }
            }
            map_trans.emplace_back(map_now);
        }

        std::vector<AtomCellSuper> pairs_tmp(order + 2);
        std::vector<unsigned int> atom_super(order + 2), atom_super_tran(order + 2);
        std::vector<unsigned int> atom_new_prim(order + 2), atom_new_super(order + 2);

        for (const auto &it: force_constant_with_cell[order]) {
            for (auto i = 0; i < order + 2; ++i) {
                atom_super[i] = map_tmp.from_true_primitive[it.pairs[i].index / 3][it.pairs[i].tran];
            }

            for (const auto &it_trans: map_trans) {
                for (auto i = 0; i < order + 2; ++i) {
                    atom_super_tran[i] = it_trans[atom_super[i]];
                }

                if (system->get_map_s2p(order)[atom_super_tran[0]].tran_num != 0) {
                    continue;
                }

                for (auto i = 0; i < order + 2; ++i) {
                    atom_new_prim[i] = system->get_map_s2p(order)[atom_super_tran[i]].atom_num;
                    pairs_tmp[i].index = 3 * atom_new_prim[i] + it.pairs[i].index % 3;
                }

                relvecs.clear();
                for (auto i = 0; i < order + 1; ++i) {
                    for (auto j = 0; j < 3; ++j) {
                        relvec_tmp[j] = it.relvecs[i][j] + cell_tmp.x_cartesian(atom_super_tran[0], j)
                                        - cell_tmp.x_cartesian(system->get_map_p2s(order)[atom_new_prim[i + 1]][0], j);
                    }
                    relvec_tmp = system->get_primcell().lattice_vector.inverse() * relvec_tmp;
                    relvecs.emplace_back(relvec_tmp);
                }
                force_constant_replicate.emplace_back(it.fcs_val, pairs_tmp, relvecs);
            }
        }

        force_constant_with_cell[order].clear();
        std::copy(force_constant_replicate.begin(),
                  force_constant_replicate.end(),
                  std::back_inserter(force_constant_with_cell[order]));
    }
}

void Fcs_phonon::load_fc2_xml()
{
    using namespace boost::property_tree;

    unsigned int atm1, atm2, xyz1, xyz2, cell_s;
    ptree pt;
    std::stringstream ss1, ss2;
    FcsClassExtent fcext_tmp;

    if (update_fc2) {
        try {
            read_xml(file_fc2, pt);
        }
        catch (std::exception &e) {
            auto str_error = "Cannot open file FC2XML ( " + file_fc2 + " )";
            exit("load_fc2_xml", str_error.c_str());
        }
    } else {
        try {
            read_xml(file_fcs, pt);
        }
        catch (std::exception &e) {
            auto str_error = "Cannot open file FCSXML ( " + file_fcs + " )";
            exit("load_fc2_xml", str_error.c_str());
        }
    }

    fc2_ext.clear();

    BOOST_FOREACH (const ptree::value_type &child_, pt.get_child("Data.ForceConstants.HARMONIC")) {
                    const auto &child = child_.second;
                    const auto str_p1 = child.get<std::string>("<xmlattr>.pair1");
                    const auto str_p2 = child.get<std::string>("<xmlattr>.pair2");

                    ss1.str("");
                    ss2.str("");
                    ss1.clear();
                    ss2.clear();

                    ss1 << str_p1;
                    ss2 << str_p2;

                    ss1 >> atm1 >> xyz1;
                    ss2 >> atm2 >> xyz2 >> cell_s;

                    fcext_tmp.atm1 = atm1 - 1;
                    fcext_tmp.xyz1 = xyz1 - 1;
                    fcext_tmp.atm2 = atm2 - 1;
                    fcext_tmp.xyz2 = xyz2 - 1;
                    fcext_tmp.cell_s = cell_s - 1;
                    fcext_tmp.fcs_val = boost::lexical_cast<double>(child.data());

                    fc2_ext.push_back(fcext_tmp);
                }
    pt.clear();
}

void Fcs_phonon::load_fcs_xml() const
{
    using namespace boost::property_tree;
    ptree pt;
    std::string str_tag;
    unsigned int atmn, xyz, cell_s;

    std::vector<Triplet> tri_vec;

    std::stringstream ss;

    AtomCellSuper ivec_tmp{};
    std::vector<AtomCellSuper> ivec_with_cell, ivec_copy;
    std::vector<Eigen::Vector3d> relvecs;

    Eigen::Vector3d relvec_tmp;
    std::vector<int> atoms_prim_tmp, coords_tmp;

    const auto xf_image = dynamical->get_xrs_image();

    std::cout << "  Reading force constants from the XML file ... ";

    try {
        read_xml(file_fcs, pt);
    }
    catch (std::exception &e) {
        auto str_error = "Cannot open file FCSXML ( " + fcs_phonon->file_fcs + " )";
        exit("load_fcs_xml", str_error.c_str());
    }

    for (unsigned int order = 0; order < maxorder; ++order) {

        if (order == 0) {
            str_tag = "Data.ForceConstants.HARMONIC";
        } else {
            str_tag = "Data.ForceConstants.ANHARM" + std::to_string(order + 2);
        }

        const auto map_tmp = system->get_mapping_super_alm(order);
        const auto xf_tmp = system->get_supercell(order).x_fractional;

        auto child_ = pt.get_child_optional(str_tag);

        if (!child_) {
            auto str_tmp = str_tag + " flag not found in the XML file";
            exit("load_fcs_xml", str_tmp.c_str());
        }

        BOOST_FOREACH (const ptree::value_type &child_, pt.get_child(str_tag)) {
                        const auto &child = child_.second;

                        auto fcs_val = boost::lexical_cast<double>(child.data());

                        ivec_with_cell.clear();

                        for (auto i = 0; i < order + 2; ++i) {
                            auto str_attr = "<xmlattr>.pair" + std::to_string(i + 1);
                            auto str_pairs = child.get<std::string>(str_attr);

                            ss.str("");
                            ss.clear();
                            ss << str_pairs;

                            if (i == 0) {
                                ss >> atmn >> xyz;
                                ivec_tmp.index = 3 * map_tmp.from_true_primitive[atmn - 1][0] + xyz - 1;
                                ivec_tmp.cell_s = 0;
                                ivec_tmp.tran = 0; // dummy
                                ivec_with_cell.push_back(ivec_tmp);
                            } else {
                                ss >> atmn >> xyz >> cell_s;
                                ivec_tmp.index = 3 * (atmn - 1) + xyz - 1;
                                ivec_tmp.cell_s = cell_s - 1;
                                ivec_tmp.tran = 0; // dummy
                                ivec_with_cell.push_back(ivec_tmp);
                            }
                        }

                        if (std::abs(fcs_val) > eps) {
                            do {
                                ivec_copy.clear();
                                for (auto i = 0; i < ivec_with_cell.size(); ++i) {
                                    atmn = ivec_with_cell[i].index / 3;
                                    xyz = ivec_with_cell[i].index % 3;
                                    ivec_tmp.index = 3 * map_tmp.to_true_primitive[atmn].atom_num + xyz;
                                    ivec_tmp.cell_s = ivec_with_cell[i].cell_s;
                                    ivec_tmp.tran = map_tmp.to_true_primitive[atmn].tran_num;
                                    ivec_copy.push_back(ivec_tmp);
                                }
                                force_constant_with_cell[order].emplace_back(fcs_val, ivec_copy);
                            } while (std::next_permutation(ivec_with_cell.begin() + 1, ivec_with_cell.end()));
                        }
                    }

        // Register relative vector information for later use
        for (auto &it: force_constant_with_cell[order]) {
            relvecs.clear();
            for (auto i = 1; i < order + 2; ++i) {
                const auto atom1_s = map_tmp.from_true_primitive[it.pairs[i].index / 3][it.pairs[i].tran];
                const auto atom2_s = map_tmp.from_true_primitive[it.pairs[0].index / 3][0];
                for (auto j = 0; j < 3; ++j) {
                    relvec_tmp[j] = xf_tmp(atom1_s, j) + xf_image[it.pairs[i].cell_s][j]
                                    - xf_tmp(atom2_s, j);
                }
                relvec_tmp = system->get_supercell(order).lattice_vector * relvec_tmp;
                relvecs.emplace_back(relvec_tmp);
            }
            it.relvecs = relvecs;
        }
    }

    std::cout << "done !" << std::endl;
}

void Fcs_phonon::MPI_Bcast_fc2_ext()
{
    unsigned int i;
    double *fcs_tmp;
    unsigned int **ind;
    FcsClassExtent fcext_tmp;

    auto nfcs = fc2_ext.size();
    MPI_Bcast(&nfcs, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

    allocate(fcs_tmp, nfcs);
    allocate(ind, nfcs, 5);

    if (mympi->my_rank == 0) {
        for (i = 0; i < nfcs; ++i) {
            fcs_tmp[i] = fc2_ext[i].fcs_val;
            ind[i][0] = fc2_ext[i].atm1;
            ind[i][1] = fc2_ext[i].xyz1;
            ind[i][2] = fc2_ext[i].atm2;
            ind[i][3] = fc2_ext[i].xyz2;
            ind[i][4] = fc2_ext[i].cell_s;
        }
    }
    MPI_Bcast(&fcs_tmp[0], nfcs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&ind[0][0], nfcs * 5, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

    if (mympi->my_rank != 0) {
        for (i = 0; i < nfcs; ++i) {
            fcext_tmp.atm1 = ind[i][0];
            fcext_tmp.xyz1 = ind[i][1];
            fcext_tmp.atm2 = ind[i][2];
            fcext_tmp.xyz2 = ind[i][3];
            fcext_tmp.cell_s = ind[i][4];
            fcext_tmp.fcs_val = fcs_tmp[i];
            fc2_ext.push_back(fcext_tmp);
        }
    }
    deallocate(fcs_tmp);
    deallocate(ind);
}

double Fcs_phonon::examine_translational_invariance(const int order,
                                                    const unsigned int nat,
                                                    const unsigned int natmin,
                                                    const std::vector<std::vector<unsigned int>> &map_p2s_in,
                                                    const std::vector<FcsArrayWithCell> &fc_in) const
{
    int j, k, l, m;

    double dev;

    double ret = 0.0;

    if (order == 0) {
        double **sum2;
        allocate(sum2, 3 * natmin, 3);

        for (j = 0; j < 3 * natmin; ++j) {
            for (k = 0; k < 3; ++k) {
                sum2[j][k] = 0.0;
            }
        }

        for (const auto &it: fc_in) {
            j = it.pairs[0].index;
            k = it.pairs[1].index % 3;
            sum2[j][k] += it.fcs_val;
        }

        for (j = 0; j < 3 * natmin; ++j) {
            for (k = 0; k < 3; ++k) {
                dev = std::abs(sum2[j][k]);
                if (ret < dev) ret = dev;
            }
        }
        deallocate(sum2);

    } else if (order == 1) {

        double ***sum3;
        allocate(sum3, 3 * natmin, 3 * nat, 3);

        for (j = 0; j < 3 * natmin; ++j) {
            for (k = 0; k < 3 * nat; ++k) {
                for (l = 0; l < 3; ++l) {
                    sum3[j][k][l] = 0.0;
                }
            }
        }

        for (const auto &it: fc_in) {
            j = it.pairs[0].index;
            k = 3 * map_p2s_in[it.pairs[1].index / 3][it.pairs[1].tran] + it.pairs[1].index % 3;
            l = it.pairs[2].index % 3;
            sum3[j][k][l] += it.fcs_val;
        }
        for (j = 0; j < 3 * natmin; ++j) {
            for (k = 0; k < 3 * nat; ++k) {
                for (l = 0; l < 3; ++l) {
                    dev = std::abs(sum3[j][k][l]);
                    if (ret < dev) ret = dev;
                }
            }
        }
        deallocate(sum3);

    } else if (order == 2) {

        double ****sum4;
        allocate(sum4, 3 * natmin, 3 * nat, 3 * nat, 3);

        for (j = 0; j < 3 * natmin; ++j) {
            for (k = 0; k < 3 * nat; ++k) {
                for (l = 0; l < 3 * nat; ++l) {
                    for (m = 0; m < 3; ++m) {
                        sum4[j][k][l][m] = 0.0;
                    }
                }
            }
        }

        for (const auto &it: fc_in) {
            j = it.pairs[0].index;
            k = 3 * map_p2s_in[it.pairs[1].index / 3][it.pairs[1].tran] + it.pairs[1].index % 3;
            l = 3 * map_p2s_in[it.pairs[2].index / 3][it.pairs[2].tran] + it.pairs[2].index % 3;
            m = it.pairs[3].index % 3;
            sum4[j][k][l][m] += it.fcs_val;
        }

        for (j = 0; j < 3 * natmin; ++j) {
            for (k = 0; k < 3 * nat; ++k) {
                for (l = 0; l < 3 * nat; ++l) {
                    for (m = 0; m < 3; ++m) {
                        dev = std::abs(sum4[j][k][l][m]);
                        if (ret < dev) ret = dev;

                    }
                }
            }
        }
        deallocate(sum4);

    }

    return ret;
}

void Fcs_phonon::MPI_Bcast_fcs_array(const unsigned int N) const
{
    int j, k;
    double *fcs_tmp;
    unsigned int ***ind;
    double ***relative_vector_tmp;

    AtomCellSuper ivec_tmp;
    std::vector<AtomCellSuper> ivec_array;

    std::vector<Eigen::Vector3d> relvecs;
    Eigen::Vector3d relvec_tmp;

    for (unsigned int i = 0; i < N; ++i) {

        int len = force_constant_with_cell[i].size();
        int nelem = i + 2;

        MPI_Bcast(&len, 1, MPI_INT, 0, MPI_COMM_WORLD);

        allocate(fcs_tmp, len);
        allocate(ind, len, nelem, 3);
        allocate(relative_vector_tmp, len, nelem - 1, 3);

        if (mympi->my_rank == 0) {
            for (j = 0; j < len; ++j) {
                fcs_tmp[j] = force_constant_with_cell[i][j].fcs_val;
                for (k = 0; k < nelem; ++k) {
                    ind[j][k][0] = force_constant_with_cell[i][j].pairs[k].index;
                    ind[j][k][1] = force_constant_with_cell[i][j].pairs[k].tran;
                    ind[j][k][2] = force_constant_with_cell[i][j].pairs[k].cell_s;
                }
                for (k = 0; k < nelem - 1; ++k) {
                    for (auto l = 0; l < 3; ++l) {
                        relative_vector_tmp[j][k][l] = force_constant_with_cell[i][j].relvecs[k][l];
                    }
                }
            }
        }

        MPI_Bcast(&fcs_tmp[0], len, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&ind[0][0][0], 3 * nelem * len, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
        MPI_Bcast(&relative_vector_tmp[0][0][0], 3 * len * (nelem - 1), MPI_DOUBLE, 0, MPI_COMM_WORLD);

        if (mympi->my_rank > 0) {
            force_constant_with_cell[i].clear();

            for (j = 0; j < len; ++j) {

                ivec_array.clear();

                for (k = 0; k < nelem; ++k) {
                    ivec_tmp.index = ind[j][k][0];
                    ivec_tmp.tran = ind[j][k][1];
                    ivec_tmp.cell_s = ind[j][k][2];

                    ivec_array.push_back(ivec_tmp);
                }

                relvecs.clear();

                for (k = 0; k < nelem - 1; ++k) {
                    for (auto l = 0; l < 3; ++l) {
                        relvec_tmp[l] = relative_vector_tmp[j][k][l];
                    }
                    relvecs.emplace_back(relvec_tmp);
                }

                force_constant_with_cell[i].emplace_back(fcs_tmp[j],
                                                         ivec_array,
                                                         relvecs);
            }
        }

        deallocate(fcs_tmp);
        deallocate(ind);
        deallocate(relative_vector_tmp);
    }
}
