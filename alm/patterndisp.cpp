/*
patterndisp.cpp

Copyright (c) 2014 Terumasa Tadano

This file is distributed under the terms of the MIT license.
Please see the file 'LICENCE.txt' in the root directory 
or http://opensource.org/licenses/mit-license.php for information.
*/

#include <iostream>
#include "patterndisp.h"
#include "memory.h"
#include "error.h"
#include "system.h"
#include "fcs.h"
#include "interaction.h"
#include "symmetry.h"
#include "constants.h"
#include "mathfunctions.h"
#include <map>

using namespace ALM_NS;

Displace::Displace(ALM *alm) : Pointers(alm) {}
Displace::~Displace() {
    memory->deallocate(pattern_all);
}

void Displace::gen_displacement_pattern()
{
    int i, j, m, order;
    int maxorder = interaction->maxorder;
    std::vector<int> group_tmp;

    std::cout << " Generating displacement patterns in ";
    if (disp_basis[0] == 'C') {
        std::cout << "Cartesian coordinate...";
    } else {
        std::cout << "fractional coordinate...";
    }

    memory->allocate(dispset, maxorder);

    for (order = 0; order < maxorder; ++order) {

        m = 0;

        for (i = 0; i < fcs->ndup[order].size(); ++i) {

            group_tmp.clear();

            // Store first order + 1 indexes as a necessarily displacement pattern.
            // Here, duplicate entries will be removed. 
            // For example, (iij) will be reduced to (ij).
            for (j = 0; j < order + 1; ++j) {
                group_tmp.push_back(fcs->fc_set[order][m].elems[j]);
            }
            group_tmp.erase(std::unique(group_tmp.begin(), group_tmp.end()), group_tmp.end());

            // Avoid equivalent entries using set.
            dispset[order].insert(DispAtomSet(group_tmp));
            m += fcs->ndup[order][i];
        }
    }

    memory->allocate(pattern_all, maxorder);

    generate_pattern_all(maxorder, pattern_all);

    memory->deallocate(dispset);

    std::cout << " done!" << std::endl;
}

void Displace::generate_pattern_all(const int N, std::vector<AtomWithDirection> *pattern)
{
    int i, j;
    int order;
    int atom_tmp;
    double disp_tmp[3];
    double norm;
    double *direc_tmp;

    std::vector<int> atoms;
    std::vector<double> directions;
    std::vector<std::vector<int> > *sign_prod, sign_reduced;
    std::vector<int> vec_tmp;

    int natom_disp;
    double sign_double;
    std::vector<int> nums;
    std::vector<double> directions_copy;

    memory->allocate(sign_prod, N);

    for (order = 0; order < N; ++order) {
        vec_tmp.clear();
        generate_signvecs(order + 1, sign_prod[order], vec_tmp);
    }

    for (order = 0; order < N; ++order) {

        pattern[order].clear();

        for (std::set<DispAtomSet>::iterator it = dispset[order].begin(); it != dispset[order].end(); ++it) {

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
                find_unique_sign_pairs(natom_disp, sign_prod[natom_disp - 1], nums, sign_reduced);
            } else {
                sign_reduced.clear();
                std::copy(sign_prod[natom_disp - 1].begin(), sign_prod[natom_disp - 1].end(), 
                    std::back_inserter(sign_reduced));
            }

            directions_copy.clear();
            std::copy(directions.begin(), directions.end(), std::back_inserter(directions_copy));

            for (std::vector<std::vector<int> >::const_iterator it = sign_reduced.begin(); 
                it != sign_reduced.end(); ++it) {
                    directions.clear();

                    for (i = 0; i < (*it).size(); ++i) {
                        sign_double = static_cast<double>((*it)[i]);

                        for (j = 0; j < 3; ++j) {
                            disp_tmp[j] = directions_copy[3 * i + j] * sign_double;
                        }

                        if (disp_basis[0] == 'F') { 
                            rotvec(disp_tmp, disp_tmp, system->rlavec);
                            for (j = 0; j < 3; ++j) {
                                disp_tmp[j] /= 2.0 * pi;
                            }
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

void Displace::generate_signvecs(const int N, std::vector<std::vector<int> > &sign, std::vector<int> vec)
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

void Displace::find_unique_sign_pairs(const int N, std::vector<std::vector<int> > sign_in, 
                                      std::vector<int> pair_in, std::vector<std::vector<int> > &sign_out) 
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
    std::vector<std::vector<int> > sign_found;
    std::vector<int> sign_tmp;
    std::vector<int> list_disp_atom;
    std::vector<IndexWithSign> index_for_sort;

    memory->allocate(disp, nat, 3);
    memory->allocate(disp_sym, nat, 3);

    sign_out.clear();
    symnum_vec.clear();
    list_disp_atom.clear();

    for (i = 0; i < pair_in.size(); ++i) {
        list_disp_atom.push_back(pair_in[i]/ 3);
    }
    list_disp_atom.erase(std::unique(list_disp_atom.begin(), list_disp_atom.end() ), list_disp_atom.end());

    for (i = 0; i < nat; ++i) {
        for (j = 0; j < 3; ++j) {
            disp[i][j] = 0.0;
        }
    }

    for (i = 0; i < pair_in.size(); ++i) {
        disp[pair_in[i]/3][pair_in[i]%3] = 1.0;
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
                    for (k = 0; k < 3; ++k) {
                        disp_sym[mapped_atom][j] += symmetry->symrel[isym][j][k] * disp[list_disp_atom[i]][k];
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

    for (std::vector<std::vector<int> >::const_iterator it = sign_in.begin(); it != sign_in.end(); ++it) {

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
            disp[pair_in[i]/3][pair_in[i]%3] = static_cast<double>((*it)[i]);
        }

        for (isym = 0; isym < symnum_vec.size(); ++isym) {

            index_for_sort.clear();

            for (i = 0; i < list_disp_atom.size(); ++i) {
                mapped_atom = symmetry->map_sym[list_disp_atom[i]][symnum_vec[isym]];

                for (j = 0; j < 3; ++j) {
                    disp_sym[mapped_atom][j] = 0.0;
                    for (k = 0; k < 3; ++k) {
                        disp_sym[mapped_atom][j] += symmetry->symrel[symnum_vec[isym]][j][k] * disp[list_disp_atom[i]][k];
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

            if (sign_tmp.size() == N && std::find(sign_found.begin(), sign_found.end(), sign_tmp) == sign_found.end()) {
                sign_found.push_back(sign_tmp);
                std::sort(sign_found.begin(), sign_found.end());
            }
        }
    }

    memory->deallocate(disp);
    memory->deallocate(disp_sym);
}

// void Displace::estimate_best_direction_harmonic(std::vector<DispDirectionHarmonic> &disp_harm)
// {
//     int i, j, k, l;
//     int nsym_disp, nsym_max;
//     int ii, jj;
//     int direc;
//     int nat = system->nat;
//     int **flag_disp;
// 
//     double factor = 1.0e-5;
//     double dprod, norm1, norm2;
//     double direc_tmp[3];
//     double disp_tmp[3], disp_tmp_frac[3], disp_best[3];
//     double **pos_disp;
// 
//     std::vector<DirectionVec> directionlist_tmp, direction_new;
// 
//     memory->allocate(flag_disp, nat, 3);
// 
//     for (i = 0; i < nat; ++i) {
//         for (j = 0; j < 3; ++j) {
//             flag_disp[i][j] = 0;
//         }
//     }
// 
//     for (std::set<DispAtomSet>::iterator it = dispset[0].begin(); it != dispset[0].end(); ++it) {
//         int index_tmp = (*it).atomset[0];
//         flag_disp[index_tmp / 3][index_tmp % 3] = 1;
//     }
// 
//     for (i = 0; i < nat; ++i) {
// 
//         directionlist_tmp.clear();
// 
//         for (j = 0; j < 3; ++j) {
//             if (flag_disp[i][j] == 1) {
//                 for (k = 0; k < 3; ++k) {
//                     direc_tmp[k] = 0.0;
//                 }
//                 direc_tmp[j] = 1.0;
//                 directionlist_tmp.push_back(direc_tmp);
//             }
//         }
//         if (directionlist_tmp.size() > 0) {
//             disp_harm.push_back(DispDirectionHarmonic(i, directionlist_tmp));
//         }
//     }
//     memory->deallocate(flag_disp);
// 
//     memory->allocate(pos_disp, nat, 3);
// 
//     for (i = 0; i < disp_harm.size(); ++i) {
// 
//         direction_new.clear();
// 
//         for (j = 0; j < disp_harm[i].directionlist.size(); ++j) {
// 
//             for (k = 0; k < 3; ++k) {
//                 disp_tmp[k] = disp_harm[i].directionlist[j].direction[k];
//                 if (std::abs(disp_tmp[k]) > eps) direc = k;
//             }
//             rotvec(disp_tmp_frac, disp_tmp, system->rlavec);
// 
// 
//             for (k = 0; k < nat; ++k) {
//                 for (l = 0; l < 3; ++l) {
//                     pos_disp[k][l] = system->xcoord[k][l];
//                 }
//             }
// 
//             for (l = 0; l < 3; ++l) pos_disp[disp_harm[i].atom][l] += factor * disp_tmp_frac[l];
// 
//             nsym_max = symmetry->numsymop(nat, pos_disp, symmetry->tolerance);
// 
//             for (k = 0; k < 3; ++k) disp_best[k] = disp_tmp[k];
// 
//             for (ii = 1; ii >= -1; --ii) {
//                 for (jj = 1; jj >= -1; --jj) {
//                     disp_tmp[direc] = 1.0;
//                     disp_tmp[(direc + 1) % 3] = static_cast<double>(ii);
//                     disp_tmp[(direc + 2) % 3] = static_cast<double>(jj);
// 
//                     rotvec(disp_tmp_frac, disp_tmp, system->rlavec);
// 
//                     for (k = 0; k < nat; ++k) {
//                         for (l = 0; l < 3; ++l) {
//                             pos_disp[k][l] = system->xcoord[k][l];
//                         }
//                     }
// 
//                     for (l = 0; l < 3; ++l) pos_disp[disp_harm[i].atom][l] += factor * disp_tmp_frac[l];
//                     nsym_disp = symmetry->numsymop(nat, pos_disp, symmetry->tolerance);
// 
//                     if (nsym_disp > nsym_max) {
//                         bool isok = true;
// 
//                         for (k = 0; k < direction_new.size(); ++k) {
//                             dprod = 0.0;
//                             norm1 = 0.0;
//                             norm2 = 0.0;
// 
//                             for (l = 0; l < 3; ++l) {
//                                 dprod = direction_new[k].direction[l] * disp_tmp[l];
//                                 norm1 = std::pow(direction_new[k].direction[l], 2);
//                                 norm2 = std::pow(disp_tmp[l], 2);
//                             }
//                             if (std::abs(dprod - std::sqrt(norm1*norm2)) < eps10) {
//                                 isok = false;
//                             }
//                         }
// 
//                         if (isok) {
//                             nsym_max = nsym_disp;
//                             for (l = 0; l < 3; ++l) disp_best[l] = disp_tmp[l];
//                         }
//                     }
//                 }
//             }
//             direction_new.push_back(disp_best);
//         }
// 
// 
//         disp_harm_best.push_back(DispDirectionHarmonic(disp_harm[i].atom, direction_new));
//     }
// 
//     memory->deallocate(pos_disp);
// 
//     std::copy(disp_harm_best.begin(), disp_harm_best.end(), disp_harm.begin());
// 
//     disp_harm_best.clear();
// }
