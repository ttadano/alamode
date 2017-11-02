/*
fcs_phonon.cpp

Copyright (c) 2014, 2015, 2016 Terumasa Tadano

This file is distributed under the terms of the MIT license.
Please see the file 'LICENCE.txt' in the root directory 
or http://opensource.org/licenses/mit-license.php for information.
*/

#include "mpi_common.h"
#include "dynamical.h"
#include "kpoint.h"
#include "fcs_phonon.h"
#include "system.h"
#include "memory.h"
#include "error.h"
#include "phonons.h"
#include "relaxation.h"
#include "constants.h"
#include "gruneisen.h"
#include "xml_parser.h"
#include <string>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/foreach.hpp>
#include <boost/optional.hpp>
#include <boost/lexical_cast.hpp>

using namespace PHON_NS;

Fcs_phonon::Fcs_phonon(PHON *phon): Pointers(phon)
{
}

Fcs_phonon::~Fcs_phonon()
{
    memory->deallocate(force_constant_with_cell);
}

void Fcs_phonon::setup(std::string mode)
{
    unsigned int i;

    if (mympi->my_rank == 0) {
        std::cout << " Force constant" << std::endl;
        std::cout << " ==============" << std::endl << std::endl;
    }

    MPI_Bcast(&relaxation->quartic_mode, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&gruneisen->print_gruneisen, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD);

    if (mode == "PHONONS") {
        require_cubic = false;
        require_quartic = false;
        maxorder = 1;

        if (gruneisen->print_gruneisen) {
            require_cubic = true;
            maxorder = 2;
        }
        if (gruneisen->print_newfcs) {
            require_cubic = true;
            maxorder = 2;

            if (relaxation->quartic_mode > 0) {
                require_quartic = true;
                maxorder = 3;
            }
        }

    } else if (mode == "RTA") {
        require_cubic = true;

        if (relaxation->quartic_mode > 0) {
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
        relaxation->quartic_mode = 1;
    }

    memory->allocate(force_constant_with_cell, maxorder);

    if (mympi->my_rank == 0) {
        double *maxdev;

        load_fc2_xml();
        load_fcs_xml();

        for (i = 0; i < maxorder; ++i) {
            std::cout << "  Number of non-zero IFCs for " << i + 2 << " order: ";
            if (i == 0) {
                std::cout << fc2_ext.size() << std::endl;
            } else {
                std::cout << force_constant_with_cell[i].size() << std::endl;
            }
        }
        std::cout << std::endl;

        memory->allocate(maxdev, maxorder);
        examine_translational_invariance(maxorder, system->nat_anharm, system->natmin,
                                         maxdev, fc2_ext, force_constant_with_cell);

        std::cout << "  Maximum deviation from the translational invariance: " << std::endl;
        for (i = 0; i < maxorder; ++i) {
            std::cout << "   Order " << i + 2 << " : " << std::setw(12)
                << std::scientific << maxdev[i] << std::endl;
        }
        std::cout << std::endl;
        memory->deallocate(maxdev);
    }

    MPI_Bcast_fc2_ext();
    MPI_Bcast_fcs_array(maxorder);
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
            std::string str_error = "Cannot open file FC2XML ( " + file_fc2 + " )";
            error->exit("load_fc2_xml", str_error.c_str());
        }
    } else {
        try {
            read_xml(file_fcs, pt);
        }
        catch (std::exception &e) {
            std::string str_error = "Cannot open file FCSXML ( " + file_fcs + " )";
            error->exit("load_fc2_xml", str_error.c_str());
        }
    }

    fc2_ext.clear();

    BOOST_FOREACH (const ptree::value_type& child_, pt.get_child("Data.ForceConstants.HARMONIC")) {
            const ptree &child = child_.second;
            const std::string str_p1 = child.get<std::string>("<xmlattr>.pair1");
            const std::string str_p2 = child.get<std::string>("<xmlattr>.pair2");

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
}

void Fcs_phonon::load_fcs_xml()
{
    using namespace boost::property_tree;
    ptree pt;
    unsigned int order;
    std::string str_tag;
    unsigned int i;
    unsigned int atmn, xyz, cell_s;

    double fcs_val;

    std::vector<Triplet> tri_vec;

    std::stringstream ss;
    std::string str_pairs;
    std::string str_attr;

    AtomCellSuper ivec_tmp;
    std::vector<AtomCellSuper> ivec_with_cell, ivec_copy;


    std::cout << "  Reading force constants from the XML file ... ";

    try {
        read_xml(file_fcs, pt);
    }
    catch (std::exception &e) {
        std::string str_error = "Cannot open file FCSXML ( " + fcs_phonon->file_fcs + " )";
        error->exit("load_fcs_xml", str_error.c_str());
    }

    for (order = 0; order < maxorder; ++order) {

        if (order == 0) {
            str_tag = "Data.ForceConstants.HARMONIC";
        } else {
            str_tag = "Data.ForceConstants.ANHARM" + boost::lexical_cast<std::string>(order + 2);
        }

        boost::optional<ptree&> child_ = pt.get_child_optional(str_tag);

        if (!child_) {
            std::string str_tmp = str_tag + " flag not found in the XML file";
            error->exit("load_fcs_xml", str_tmp.c_str());
        }

        BOOST_FOREACH (const ptree::value_type& child_, pt.get_child(str_tag)) {
                const ptree &child = child_.second;

                fcs_val = boost::lexical_cast<double>(child.data());

                ivec_with_cell.clear();

                for (i = 0; i < order + 2; ++i) {
                    str_attr = "<xmlattr>.pair" + boost::lexical_cast<std::string>(i + 1);
                    str_pairs = child.get<std::string>(str_attr);

                    ss.str("");
                    ss.clear();
                    ss << str_pairs;

                    if (i == 0) {

                        ss >> atmn >> xyz;
                        if (update_fc2) {
                            ivec_tmp.index = 3 * system->map_p2s_anharm_orig[atmn - 1][0] + xyz - 1;
                        } else {
                            ivec_tmp.index = 3 * system->map_p2s_anharm[atmn - 1][0] + xyz - 1;
                        }
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

                        for (i = 0; i < ivec_with_cell.size(); ++i) {
                            atmn = ivec_with_cell[i].index / 3;
                            xyz = ivec_with_cell[i].index % 3;
                            ivec_tmp.index = 3 * system->map_s2p_anharm[atmn].atom_num + xyz;
                            ivec_tmp.cell_s = ivec_with_cell[i].cell_s;
                            ivec_tmp.tran = system->map_s2p_anharm[atmn].tran_num;
                            ivec_copy.push_back(ivec_tmp);
                        }

                        force_constant_with_cell[order].push_back(FcsArrayWithCell(fcs_val, ivec_copy));

                    } while (std::next_permutation(ivec_with_cell.begin() + 1, ivec_with_cell.end()));
                }
            }
    }

    std::cout << "done !" << std::endl;
}

void Fcs_phonon::MPI_Bcast_fc_class(const unsigned int N)
{
    unsigned int i;
    int j, k;
    int len;
    int nelem;
    double *fcs_tmp;
    unsigned int ***ind;

    Triplet tri_tmp;
    std::vector<Triplet> tri_vec;

    for (i = 0; i < N; ++i) {

        len = force_constant[i].size();
        nelem = i + 2;

        MPI_Bcast(&len, 1, MPI_INT, 0, MPI_COMM_WORLD);

        memory->allocate(fcs_tmp, len);
        memory->allocate(ind, len, nelem, 3);

        if (mympi->my_rank == 0) {
            for (j = 0; j < len; ++j) {
                fcs_tmp[j] = force_constant[i][j].fcs_val;
                for (k = 0; k < nelem; ++k) {
                    ind[j][k][0] = force_constant[i][j].elems[k].atom;
                    ind[j][k][1] = force_constant[i][j].elems[k].cell;
                    ind[j][k][2] = force_constant[i][j].elems[k].xyz;
                }
            }
        }

        MPI_Bcast(&fcs_tmp[0], len, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&ind[0][0][0], 3 * nelem * len, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

        if (mympi->my_rank > 0) {
            force_constant[i].clear();

            for (j = 0; j < len; ++j) {

                tri_vec.clear();

                for (k = 0; k < nelem; ++k) {
                    tri_tmp.atom = ind[j][k][0];
                    tri_tmp.cell = ind[j][k][1];
                    tri_tmp.xyz = ind[j][k][2];

                    tri_vec.push_back(tri_tmp);
                }
                force_constant[i].push_back(FcsClass(fcs_tmp[j], tri_vec));
            }
        }

        memory->deallocate(fcs_tmp);
        memory->deallocate(ind);
    }
}

void Fcs_phonon::MPI_Bcast_fc2_ext()
{
    unsigned int i;
    double *fcs_tmp;
    unsigned int **ind;
    unsigned int nfcs;
    FcsClassExtent fcext_tmp;

    nfcs = fc2_ext.size();
    MPI_Bcast(&nfcs, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

    memory->allocate(fcs_tmp, nfcs);
    memory->allocate(ind, nfcs, 5);

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
    memory->deallocate(fcs_tmp);
    memory->deallocate(ind);
}


void Fcs_phonon::examine_translational_invariance(const int n, const unsigned int nat,
                                                  const unsigned int natmin, double *ret,
                                                  std::vector<FcsClassExtent> &fc2,
                                                  std::vector<FcsArrayWithCell> *fcs)
{
    int i, j, k, l, m;

    double dev;
    double **sum2;
    double ***sum3;
    double ****sum4;

    bool force_asr = false;
    FcsClassExtent fc2_tmp;
    std::vector<FcsClassExtent>::iterator it_target;


    for (i = 0; i < n; ++i) ret[i] = 0.0;

    for (i = 0; i < n; ++i) {

        if (i == 0) {
            memory->allocate(sum2, 3 * natmin, 3);

            for (j = 0; j < 3 * natmin; ++j) {
                for (k = 0; k < 3; ++k) {
                    sum2[j][k] = 0.0;
                }
            }
            for (auto it = fc2.cbegin(); it != fc2.cend(); ++it) {
                sum2[3 * (*it).atm1 + (*it).xyz1][(*it).xyz2] += (*it).fcs_val;
            }

            if (force_asr) {
                std::cout << "  force_asr = true: Modify harmonic force constans so that the ASR is satisfied." << std::endl;
                for (j = 0; j < natmin; ++j) {
                    for (k = 0; k < 3; ++k) {
                        for (m = 0; m < 3; ++m) {
                            fc2_tmp.atm1 = j;
                            fc2_tmp.xyz1 = k;
                            fc2_tmp.atm2 = system->map_p2s[j][0];
                            fc2_tmp.xyz2 = m;
                            fc2_tmp.cell_s = 0;
                            fc2_tmp.fcs_val = sum2[3 * j + k][m];
                            it_target = std::find(fc2.begin(), fc2.end(), fc2_tmp);
                            if (std::abs(fc2_tmp.fcs_val) > eps12) {
                                if (it_target != fc2.end()) {
                                    fc2[it_target - fc2.begin()].fcs_val -= fc2_tmp.fcs_val;
                                } else {
                                    error->exit("examine_translational_invariance",
                                                "Corresponding IFC not found.");
                                }
                            }
                        }
                    }
                }
                for (j = 0; j < 3 * natmin; ++j) {
                    for (k = 0; k < 3; ++k) {
                        sum2[j][k] = 0.0;
                    }
                }

                for (auto it = fc2.cbegin(); it != fc2.cend(); ++it) {
                    sum2[3 * (*it).atm1 + (*it).xyz1][(*it).xyz2] += (*it).fcs_val;
                }
            }

            for (j = 0; j < 3 * natmin; ++j) {
                for (k = 0; k < 3; ++k) {
                    dev = std::abs(sum2[j][k]);
                    if (ret[i] < dev) ret[i] = dev;
                }
            }
            memory->deallocate(sum2);

        } else if (i == 1) {

            memory->allocate(sum3, 3 * natmin, 3 * nat, 3);

            for (j = 0; j < 3 * natmin; ++j) {
                for (k = 0; k < 3 * nat; ++k) {
                    for (l = 0; l < 3; ++l) {
                        sum3[j][k][l] = 0.0;
                    }
                }
            }

            for (auto it = fcs[i].cbegin(); it != fcs[i].cend(); ++it) {
                j = (*it).pairs[0].index;
                k = 3 * (natmin * (*it).pairs[1].tran + (*it).pairs[1].index / 3) + (*it).pairs[1].index % 3;
                l = (*it).pairs[2].index % 3;
                sum3[j][k][l] += (*it).fcs_val;
            }
            for (j = 0; j < 3 * natmin; ++j) {
                for (k = 0; k < 3 * nat; ++k) {
                    for (l = 0; l < 3; ++l) {
                        dev = std::abs(sum3[j][k][l]);
                        if (ret[i] < dev) ret[i] = dev;
                    }
                }
            }

            memory->deallocate(sum3);

        } else if (i == 2) {

            memory->allocate(sum4, 3 * natmin, 3 * nat, 3 * nat, 3);

            for (j = 0; j < 3 * natmin; ++j) {
                for (k = 0; k < 3 * nat; ++k) {
                    for (l = 0; l < 3 * nat; ++l) {
                        for (m = 0; m < 3; ++m) {
                            sum4[j][k][l][m] = 0.0;
                        }
                    }
                }
            }

            for (auto it = fcs[i].cbegin(); it != fcs[i].cend(); ++it) {
                j = (*it).pairs[0].index;
                k = 3 * system->map_p2s_anharm[(*it).pairs[1].index / 3][(*it).pairs[1].tran]
                    + (*it).pairs[1].index % 3;
                l = 3 * system->map_p2s_anharm[(*it).pairs[2].index / 3][(*it).pairs[2].tran]
                    + (*it).pairs[2].index % 3;
                m = (*it).pairs[3].index % 3;

                sum4[j][k][l][m] += (*it).fcs_val;
            }

            for (j = 0; j < 3 * natmin; ++j) {
                for (k = 0; k < 3 * nat; ++k) {
                    for (l = 0; l < 3 * nat; ++l) {
                        for (m = 0; m < 3; ++m) {
                            dev = std::abs(sum4[j][k][l][m]);
                            if (ret[i] < dev) ret[i] = dev;

                        }
                    }
                }
            }

            memory->deallocate(sum4);

        }

    }
}


void Fcs_phonon::MPI_Bcast_fcs_array(const unsigned int N)
{
    unsigned int i;
    int j, k;
    int len;
    int nelem;
    double *fcs_tmp;
    unsigned int ***ind;

    AtomCellSuper ivec_tmp;
    std::vector<AtomCellSuper> ivec_array;

    for (i = 0; i < N; ++i) {

        len = force_constant_with_cell[i].size();
        nelem = i + 2;

        MPI_Bcast(&len, 1, MPI_INT, 0, MPI_COMM_WORLD);

        memory->allocate(fcs_tmp, len);
        memory->allocate(ind, len, nelem, 3);

        if (mympi->my_rank == 0) {
            for (j = 0; j < len; ++j) {
                fcs_tmp[j] = force_constant_with_cell[i][j].fcs_val;
                for (k = 0; k < nelem; ++k) {
                    ind[j][k][0] = force_constant_with_cell[i][j].pairs[k].index;
                    ind[j][k][1] = force_constant_with_cell[i][j].pairs[k].tran;
                    ind[j][k][2] = force_constant_with_cell[i][j].pairs[k].cell_s;
                }
            }
        }

        MPI_Bcast(&fcs_tmp[0], len, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&ind[0][0][0], 3 * nelem * len, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

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
                force_constant_with_cell[i].push_back(FcsArrayWithCell(fcs_tmp[j],
                                                                       ivec_array));
            }
        }

        memory->deallocate(fcs_tmp);
        memory->deallocate(ind);
    }
}
