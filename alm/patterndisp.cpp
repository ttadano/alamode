/*
patterndisp.cpp

Copyright (c) 2014, 2015, 2016 Terumasa Tadano

This file is distributed under the terms of the MIT license.
Please see the file 'LICENCE.txt' in the root directory 
or http://opensource.org/licenses/mit-license.php for information.
*/

#include <iostream>
#include <iomanip>
#include "patterndisp.h"
#include "memory.h"
#include "error.h"
#include "system.h"
#include "fcs.h"
#include "interaction.h"
#include "symmetry.h"
#include "constants.h"
#include "mathfunctions.h"
#include "constraint.h"
#include <map>
#include <boost/bimap.hpp>

using namespace ALM_NS;

Displace::Displace(ALM *alm) : Pointers(alm)
{
}

Displace::~Displace()
{
    if (alm->mode == "suggest") {
        memory->deallocate(pattern_all);
    }
}

void Displace::gen_displacement_pattern()
{
    int i, j, m, order;
    int maxorder = interaction->maxorder;
    std::string preferred_basis;
    std::vector<int> group_tmp;
    std::vector<ConstraintClass> *constsym;

    std::vector<int> pairs;
    std::set<int> *include_set;
    std::set<DispAtomSet> *dispset;

    std::vector<int> *nequiv;
    std::vector<FcProperty> *fc_table, *fc_zeros;

    std::vector<ConstraintTypeFix> *const_fix_tmp;
    std::vector<ConstraintTypeRelate> *const_relate_tmp;
    boost::bimap<int, int> *index_bimap_tmp;

    std::cout << " DISPLACEMENT PATTERN" << std::endl;
    std::cout << " ====================" << std::endl << std::endl;

    // Decide preferred basis (Cartesian or Lattice)
    int ncompat_cart = 0;
    int ncompat_latt = 0;
    for (auto it = symmetry->SymmData.begin(); it != symmetry->SymmData.end(); ++it) {
        if ((*it).compatible_with_cartesian) ++ncompat_cart;
        if ((*it).compatible_with_lattice) ++ncompat_latt;
    }
    if (ncompat_cart >= ncompat_latt) {
        preferred_basis = "Cartesian";
    } else {
        preferred_basis = "Lattice";
    }

    memory->allocate(fc_table, maxorder);
    memory->allocate(fc_zeros, maxorder);
    memory->allocate(constsym, maxorder);
    memory->allocate(nequiv, maxorder);
    memory->allocate(const_fix_tmp, maxorder);
    memory->allocate(const_relate_tmp, maxorder);
    memory->allocate(index_bimap_tmp, maxorder);

    for (order = 0; order < maxorder; ++order) {
        fcs->generate_force_constant_table(order, interaction->pairs[order],
                                           symmetry->SymmData, preferred_basis,
                                           fc_table[order], nequiv[order],
                                           fc_zeros[order], false);

        constraint->get_symmetry_constraint(order, interaction->pairs[order],
                                            symmetry->SymmData, preferred_basis,
                                            fc_table[order], nequiv[order],
                                            constsym[order]);
    }

    constraint->get_mapping_constraint(maxorder, nequiv, constsym, const_fix_tmp,
                                       const_relate_tmp, index_bimap_tmp, true);

    for (order = 0; order < maxorder; ++order) {
        std::cout << "  Number of free" << std::setw(9)
            << interaction->str_order[order] << " FCs : "
            << index_bimap_tmp[order].size() << std::endl;
    }
    std::cout << std::endl;

    memory->deallocate(constsym);
    memory->deallocate(const_fix_tmp);
    memory->deallocate(const_relate_tmp);
    memory->deallocate(fc_zeros);

    memory->allocate(include_set, maxorder);

    for (order = 0; order < maxorder; ++order) {
        include_set[order].clear();

        for (boost::bimap<int, int>::const_iterator it = index_bimap_tmp[order].begin();
             it != index_bimap_tmp[order].end(); ++it) {
            include_set[order].insert((*it).right);
        }
    }

    memory->deallocate(index_bimap_tmp);


    std::cout << "  Generating displacement patterns in ";
    //    if (disp_basis[0] == 'C') {
    std::cout << "Cartesian coordinate... ";
    //    } else {
    //        std::cout << "fractional coordinate...";
    //    }

    memory->allocate(dispset, maxorder);

    for (order = 0; order < maxorder; ++order) {

        m = 0;

        for (i = 0; i < nequiv[order].size(); ++i) {

            if (include_set[order].find(i) != include_set[order].end()) {

                group_tmp.clear();

                // Store first order + 1 indexes as a necessary displacement pattern.
                // Here, duplicate entries will be removed. 
                // For example, (iij) will be reduced to (ij).
                for (j = 0; j < order + 1; ++j) {
                    group_tmp.push_back(fc_table[order][m].elems[j]);
                }
                group_tmp.erase(std::unique(group_tmp.begin(), group_tmp.end()),
                                group_tmp.end());

                // Avoid equivalent entries using set.
                dispset[order].insert(DispAtomSet(group_tmp));

            }

            m += nequiv[order][i];
        }
    }
    memory->deallocate(include_set);
    memory->deallocate(nequiv);
    memory->deallocate(fc_table);

    memory->allocate(pattern_all, maxorder);
    generate_pattern_all(maxorder, pattern_all,
                         dispset, preferred_basis);

    memory->deallocate(dispset);


    std::cout << " done!" << std::endl;
    std::cout << std::endl;

    for (order = 0; order < maxorder; ++order) {
        std::cout << "  Number of disp. patterns for " << std::setw(9)
            << interaction->str_order[order] << " : "
            << pattern_all[order].size() << std::endl;
    }
    std::cout << std::endl;
}

void Displace::generate_pattern_all(const int N,
                                    std::vector<AtomWithDirection> *pattern,
                                    std::set<DispAtomSet> *dispset_in,
                                    const std::string preferred_basis)
{
    int i, j;
    int order;
    int atom_tmp, natom_disp;
    double sign_double;
    double disp_tmp[3];

    std::vector<int> atoms, vec_tmp, nums;
    std::vector<double> directions, directions_copy;
    std::vector<std::vector<int>> *sign_prod, sign_reduced;


    memory->allocate(sign_prod, N);

    for (order = 0; order < N; ++order) {
        vec_tmp.clear();
        generate_signvecs(order + 1, sign_prod[order], vec_tmp);
    }

    for (order = 0; order < N; ++order) {

        pattern[order].clear();

        for (auto it = dispset_in[order].begin(); it != dispset_in[order].end(); ++it) {

            atoms.clear();
            directions.clear();
            nums.clear();

            for (i = 0; i < (*it).atomset.size(); ++i) {

                atom_tmp = (*it).atomset[i] / 3;

                nums.push_back((*it).atomset[i]);
                atoms.push_back(atom_tmp);

                for (j = 0; j < 3; ++j) {
                    disp_tmp[j] = 0.0;
                }
                disp_tmp[(*it).atomset[i] % 3] = 1.0;

                for (j = 0; j < 3; ++j) directions.push_back(disp_tmp[j]);
            }

            natom_disp = atoms.size();

            if (trim_dispsign_for_evenfunc) {
                find_unique_sign_pairs(natom_disp, sign_prod[natom_disp - 1],
                                       nums, sign_reduced, preferred_basis);
            } else {
                sign_reduced.clear();
                std::copy(sign_prod[natom_disp - 1].begin(),
                          sign_prod[natom_disp - 1].end(),
                          std::back_inserter(sign_reduced));
            }

            directions_copy.clear();
            std::copy(directions.begin(), directions.end(),
                      std::back_inserter(directions_copy));

            for (auto it2 = sign_reduced.cbegin(); it2 != sign_reduced.cend(); ++it2) {
                directions.clear();

                for (i = 0; i < (*it2).size(); ++i) {
                    sign_double = static_cast<double>((*it2)[i]);

                    for (j = 0; j < 3; ++j) {
                        disp_tmp[j] = directions_copy[3 * i + j] * sign_double;
                    }

                    //                    if (disp_basis[0] == 'F') {
                    //                        rotvec(disp_tmp, disp_tmp, system->rlavec);
                    //                        for (j = 0; j < 3; ++j) {
                    //                            disp_tmp[j] /= 2.0 * pi;
                    //                        }
                    //                    }

                    if (preferred_basis == "Lattice") {
                        rotvec(disp_tmp, disp_tmp, system->lavec);
                        double norm = 0.0;
                        for (j = 0; j < 3; ++j) norm += disp_tmp[j] * disp_tmp[j];
                        norm = std::sqrt(norm);
                        for (j = 0; j < 3; ++j) disp_tmp[j] /= norm;
                    }

                    for (j = 0; j < 3; ++j) {
                        directions.push_back(disp_tmp[j]);
                    }
                }
                pattern[order].push_back(AtomWithDirection(atoms, directions));
            }
        }
    }

    memory->deallocate(sign_prod);
}

void Displace::generate_signvecs(const int N,
                                 std::vector<std::vector<int>> &sign,
                                 std::vector<int> vec)
{
    // returns the product of signs ('+','-')

    if (N == 0) {
        sign.push_back(vec);
        vec.clear();
    } else {
        std::vector<int> vec_tmp;

        vec_tmp.clear();
        std::copy(vec.begin(), vec.end(), std::back_inserter(vec_tmp));

        vec_tmp.push_back(1);
        generate_signvecs(N - 1, sign, vec_tmp);

        vec_tmp.clear();
        std::copy(vec.begin(), vec.end(), std::back_inserter(vec_tmp));
        vec_tmp.push_back(-1);
        generate_signvecs(N - 1, sign, vec_tmp);
    }
}

void Displace::find_unique_sign_pairs(const int N,
                                      std::vector<std::vector<int>> sign_in,
                                      std::vector<int> pair_in,
                                      std::vector<std::vector<int>> &sign_out,
                                      const std::string preferred_basis)
{
    int isym, i, j, k;
    int mapped_atom;
    int mapped_index;
    int nat = system->nat;

    bool flag_avail;

    double disp_tmp;
    double **disp, **disp_sym;

    std::vector<int> symnum_vec;
    std::vector<int>::iterator loc;
    std::vector<int> atom_tmp, pair_tmp;
    std::vector<std::vector<int>> sign_found;
    std::vector<int> sign_tmp;
    std::vector<int> list_disp_atom;
    std::vector<IndexWithSign> index_for_sort;

    memory->allocate(disp, nat, 3);
    memory->allocate(disp_sym, nat, 3);

    sign_out.clear();
    symnum_vec.clear();
    list_disp_atom.clear();

    for (i = 0; i < pair_in.size(); ++i) {
        list_disp_atom.push_back(pair_in[i] / 3);
    }
    list_disp_atom.erase(std::unique(list_disp_atom.begin(), list_disp_atom.end()),
                         list_disp_atom.end());

    for (i = 0; i < nat; ++i) {
        for (j = 0; j < 3; ++j) {
            disp[i][j] = 0.0;
        }
    }

    for (i = 0; i < pair_in.size(); ++i) {
        disp[pair_in[i] / 3][pair_in[i] % 3] = 1.0;
    }

    // Find symmetry operations which can be used to
    // reduce the number of sign patterns (+, -) of displacements

    for (isym = 0; isym < symmetry->nsym; ++isym) {

        flag_avail = true;
        pair_tmp.clear();

        for (i = 0; i < N; ++i) {
            if (!flag_avail) break;

            mapped_atom = symmetry->map_sym[pair_in[i] / 3][isym];
            mapped_index = 3 * mapped_atom + pair_in[i] % 3;

            loc = std::find(pair_in.begin(), pair_in.end(), mapped_index);
            if (loc == pair_in.end()) {
                flag_avail = false;
            }

            pair_tmp.push_back(mapped_index);
        }

        if (!flag_avail) continue;

        // The symmetry may be useful only when the pair doesn't change
        // by the symmetry operation.
        std::sort(pair_tmp.begin(), pair_tmp.end());

        if (pair_tmp == pair_in) {

            pair_tmp.clear();

            for (i = 0; i < list_disp_atom.size(); ++i) {
                mapped_atom = symmetry->map_sym[list_disp_atom[i]][isym];

                for (j = 0; j < 3; ++j) {
                    disp_sym[mapped_atom][j] = 0.0;

                    if (preferred_basis == "Cartesian") {
                        for (k = 0; k < 3; ++k) {
                            disp_sym[mapped_atom][j]
                                += symmetry->SymmData[isym].rotation_cart[j][k] * disp[list_disp_atom[i]][k];
                        }
                    } else if (preferred_basis == "Lattice") {
                        for (k = 0; k < 3; ++k) {
                            disp_sym[mapped_atom][j]
                                += static_cast<double>(symmetry->SymmData[isym].rotation[j][k])
                                * disp[list_disp_atom[i]][k];
                        }
                    } else {
                        error->exit("find_unique_sign_pairs",
                                    "Invalid basis. This cannot happen.");
                    }

                    disp_tmp = disp_sym[mapped_atom][j];
                    if (std::abs(disp_tmp) > eps) {
                        pair_tmp.push_back(3 * mapped_atom + j);
                    }
                }
            }

            std::sort(pair_tmp.begin(), pair_tmp.end());

            if (pair_tmp == pair_in) {
                symnum_vec.push_back(isym);
            }
        }
    }

    // Now find unique pairs of displacement directions 

    sign_found.clear();

    for (std::vector<std::vector<int>>::const_iterator it = sign_in.begin();
         it != sign_in.end(); ++it) {

        // if the sign has already been found before, cycle the loop.
        // else, add the current sign pairs to the return variable.
        if (std::find(sign_found.begin(), sign_found.end(), (*it)) != sign_found.end()) {
            continue;
        } else {
            sign_out.push_back(*it);
        }

        for (i = 0; i < system->nat; ++i) {
            for (j = 0; j < 3; ++j) {
                disp[i][j] = 0.0;
            }
        }

        for (i = 0; i < N; ++i) {
            disp[pair_in[i] / 3][pair_in[i] % 3] = static_cast<double>((*it)[i]);
        }

        for (isym = 0; isym < symnum_vec.size(); ++isym) {

            index_for_sort.clear();

            for (i = 0; i < list_disp_atom.size(); ++i) {
                mapped_atom = symmetry->map_sym[list_disp_atom[i]][symnum_vec[isym]];

                for (j = 0; j < 3; ++j) {
                    disp_sym[mapped_atom][j] = 0.0;
                    if (preferred_basis == "Cartesian") {
                        for (k = 0; k < 3; ++k) {
                            disp_sym[mapped_atom][j]
                                += symmetry->SymmData[symnum_vec[isym]].rotation_cart[j][k]
                                * disp[list_disp_atom[i]][k];
                        }
                    } else if (preferred_basis == "Lattice") {
                        for (k = 0; k < 3; ++k) {
                            disp_sym[mapped_atom][j]
                                += static_cast<double>(symmetry->SymmData[symnum_vec[isym]].rotation[j][k])
                                * disp[list_disp_atom[i]][k];
                        }
                    } else {
                        error->exit("find_unique_sign_pairs",
                                    "Invalid basis. This cannot happen.");
                    }

                    disp_tmp = disp_sym[mapped_atom][j];

                    if (std::abs(disp_tmp) > eps) {

                        if (disp_tmp < 0.0) {
                            index_for_sort.push_back(IndexWithSign(3 * mapped_atom + j, -1));
                        } else {
                            index_for_sort.push_back(IndexWithSign(3 * mapped_atom + j, 1));
                        }
                    }
                }

            }
            std::sort(index_for_sort.begin(), index_for_sort.end());

            sign_tmp.clear();

            for (i = 0; i < index_for_sort.size(); ++i) {
                sign_tmp.push_back(index_for_sort[i].sign);
            }

            if ((sign_tmp.size() == N) &&
                (std::find(sign_found.begin(), sign_found.end(), sign_tmp) == sign_found.end())) {
                sign_found.push_back(sign_tmp);
                std::sort(sign_found.begin(), sign_found.end());
            }
        }
    }

    memory->deallocate(disp);
    memory->deallocate(disp_sym);
}
