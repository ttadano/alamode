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

using namespace ALM_NS;

Displace::Displace(ALM *alm) : Pointers(alm) {}
Displace::~Displace() {}

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

            for (j = 0; j < order + 1; ++j) {
                group_tmp.push_back(fcs->fc_set[order][m].elems[j]);
            }
            dispset[order].insert(DispAtomSet(group_tmp));
            m += fcs->ndup[order][i];
        }
    }

    estimate_best_direction_harmonic(disp_harm);

    memory->allocate(pattern_all, maxorder);
    generate_pattern_all(maxorder, pattern_all);

    std::cout << " done!" << std::endl;
}

void Displace::estimate_best_direction_harmonic(std::vector<DispDirectionHarmonic> &disp_harm)
{
    int i, j, k, l;
    int nsym_disp, nsym_max;
    int ii, jj;
    int direc;
    int nat = system->nat;
    int **flag_disp;

    double factor = 1.0e-5;
    double dprod, norm1, norm2;
    double direc_tmp[3];
    double disp_tmp[3], disp_tmp_frac[3], disp_best[3];
    double **pos_disp;

    std::vector<DirectionVec> directionlist_tmp, direction_new;

    memory->allocate(flag_disp, nat, 3);

    for (i = 0; i < nat; ++i) {
        for (j = 0; j < 3; ++j) {
            flag_disp[i][j] = 0;
        }
    }

    for (std::set<DispAtomSet>::iterator it = dispset[0].begin(); it != dispset[0].end(); ++it) {
        int index_tmp = (*it).atomset[0];
        flag_disp[index_tmp / 3][index_tmp % 3] = 1;
    }

    for (i = 0; i < nat; ++i) {

        directionlist_tmp.clear();

        for (j = 0; j < 3; ++j) {
            if (flag_disp[i][j] == 1) {
                for (k = 0; k < 3; ++k) {
                    direc_tmp[k] = 0.0;
                }
                direc_tmp[j] = 1.0;
                directionlist_tmp.push_back(direc_tmp);
            }
        }
        if (directionlist_tmp.size() > 0) {
            disp_harm.push_back(DispDirectionHarmonic(i, directionlist_tmp));
        }
    }
    memory->deallocate(flag_disp);

    memory->allocate(pos_disp, nat, 3);

    for (i = 0; i < disp_harm.size(); ++i) {

        direction_new.clear();

        for (j = 0; j < disp_harm[i].directionlist.size(); ++j) {

            for (k = 0; k < 3; ++k) {
                disp_tmp[k] = disp_harm[i].directionlist[j].direction[k];
                if (std::abs(disp_tmp[k]) > eps) direc = k;
            }
            rotvec(disp_tmp_frac, disp_tmp, system->rlavec);


            for (k = 0; k < nat; ++k) {
                for (l = 0; l < 3; ++l) {
                    pos_disp[k][l] = system->xcoord[k][l];
                }
            }

            for (l = 0; l < 3; ++l) pos_disp[disp_harm[i].atom][l] += factor * disp_tmp_frac[l];

            nsym_max = symmetry->numsymop(nat, pos_disp, symmetry->tolerance);

            for (k = 0; k < 3; ++k) disp_best[k] = disp_tmp[k];

            for (ii = 1; ii >= -1; --ii) {
                for (jj = 1; jj >= -1; --jj) {
                    disp_tmp[direc] = 1.0;
                    disp_tmp[(direc + 1) % 3] = static_cast<double>(ii);
                    disp_tmp[(direc + 2) % 3] = static_cast<double>(jj);

                    rotvec(disp_tmp_frac, disp_tmp, system->rlavec);

                    for (k = 0; k < nat; ++k) {
                        for (l = 0; l < 3; ++l) {
                            pos_disp[k][l] = system->xcoord[k][l];
                        }
                    }

                    for (l = 0; l < 3; ++l) pos_disp[disp_harm[i].atom][l] += factor * disp_tmp_frac[l];
                    nsym_disp = symmetry->numsymop(nat, pos_disp, symmetry->tolerance);

                    if (nsym_disp > nsym_max) {
                        bool isok = true;

                        for (k = 0; k < direction_new.size(); ++k) {
                            dprod = 0.0;
                            norm1 = 0.0;
                            norm2 = 0.0;

                            for (l = 0; l < 3; ++l) {
                                dprod = direction_new[k].direction[l] * disp_tmp[l];
                                norm1 = std::pow(direction_new[k].direction[l], 2);
                                norm2 = std::pow(disp_tmp[l], 2);
                            }
                            if (std::abs(dprod - std::sqrt(norm1*norm2)) < eps10) {
                                isok = false;
                            }
                        }

                        if (isok) {
                            nsym_max = nsym_disp;
                            for (l = 0; l < 3; ++l) disp_best[l] = disp_tmp[l];
                        }
                    }
                }
            }
            direction_new.push_back(disp_best);
        }


        disp_harm_best.push_back(DispDirectionHarmonic(disp_harm[i].atom, direction_new));
    }

    memory->deallocate(pos_disp);

    std::copy(disp_harm_best.begin(), disp_harm_best.end(), disp_harm.begin());

    disp_harm_best.clear();
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

    std::vector<std::vector<int> > *sign_prod;
    std::vector<int> vec_tmp;

    memory->allocate(sign_prod, N);

    for (order = 0; order < N; ++order) {
        vec_tmp.clear();
        generate_signvecs(order + 1, sign_prod[order], vec_tmp);
    }

    for (order = 0; order < N; ++order) {

        pattern[order].clear();

        if (order == 0) {

            // Special treatment for harmonic terms

            for (std::vector<DispDirectionHarmonic>::iterator it  = disp_harm.begin(); 
                it != disp_harm.end(); ++it) {
                    DispDirectionHarmonic disp_now = *it;


                    atom_tmp = disp_now.atom;

                    for (std::vector<DirectionVec>::iterator it2  = disp_now.directionlist.begin(); 
                        it2 != disp_now.directionlist.end(); ++it2) {

                            atoms.clear();
                            directions.clear();

                            atoms.push_back(disp_now.atom);

                            for (i = 0; i < 3; ++i) {
                                disp_tmp[i] = (*it2).direction[i];
                            }
                            norm = disp_tmp[0] * disp_tmp[0] + disp_tmp[1] * disp_tmp[1] + disp_tmp[2] * disp_tmp[2];
                            for (i = 0; i < 3; ++i) disp_tmp[i] /= std::sqrt(norm);

                            if (disp_basis[0] == 'F') {
                                rotvec(disp_tmp, disp_tmp, system->rlavec);
                                for (i = 0; i < 3; ++i) disp_tmp[i] /= 2.0 * pi;
                            }

                            for (i = 0; i < sign_prod[0].size(); ++i) {
                                directions.clear();

                                for (j = 0; j < 3; ++j) {
                                    directions.push_back(disp_tmp[j] * static_cast<double>(sign_prod[order][i][0]));
                                }
                                pattern[order].push_back(AtomWithDirection(atoms, directions));		
                            }
                    }			
            }

        } else {

            // Anharmonic terms

            int natom_disp;
            std::vector<int>::iterator loc;
            std::vector<int> nums;
            std::vector<double> directions_copy;
            double sign_double;


            for (std::set<DispAtomSet>::iterator it = dispset[order].begin(); it != dispset[order].end(); ++it) {

                atoms.clear();
                directions.clear();

                nums.clear();

                for (i = 0; i < (*it).atomset.size(); ++i) {

                    atom_tmp = (*it).atomset[i] / 3;

                    loc = std::find(nums.begin(), nums.end(), (*it).atomset[i]);

                    if (loc != nums.end()) continue;

                    nums.push_back((*it).atomset[i]);
                    atoms.push_back(atom_tmp);

                    for (j = 0; j < 3; ++j) {
                        disp_tmp[j] = 0.0;
                    }
                    disp_tmp[(*it).atomset[i] % 3] = 1.0;

                    if (disp_basis[0] == 'F') {
                        rotvec(disp_tmp, disp_tmp, system->rlavec);
                        for (j = 0; j < 3; ++j) disp_tmp[j] /= 2.0 * pi;
                    }

                    for (j = 0; j < 3; ++j) directions.push_back(disp_tmp[j]);
                }

                natom_disp = atoms.size();

                directions_copy.clear();
                std::copy(directions.begin(), directions.end(), std::back_inserter(directions_copy));
                
                for (std::vector<std::vector<int> >::const_iterator it = sign_prod[natom_disp - 1].begin(); 
                    it != sign_prod[natom_disp - 1].end(); ++it) {

                        directions.clear();

                        for (i = 0; i < (*it).size(); ++i) {
                            sign_double = static_cast<double>((*it)[i]);
                            for (j = 0; j < 3; ++j) {
                                directions.push_back(directions_copy[3 * i + j] * sign_double);
                            }
                        }
                        pattern[order].push_back(AtomWithDirection(atoms, directions));
                }               
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