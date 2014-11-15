/*
interaction.cpp

Copyright (c) 2014 Terumasa Tadano

This file is distributed under the terms of the MIT license.
Please see the file 'LICENCE.txt' in the root directory 
or http://opensource.org/licenses/mit-license.php for information.
*/

#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <set>
#include <boost/lexical_cast.hpp>
#include "interaction.h"
#include "memory.h"
#include "system.h"
#include "error.h"
#include "symmetry.h"
#include "combination.h"
#include "constants.h"
#include "listcomparison.h"
#include "files.h"
#include "timer.h"
#include "fcs.h"
#include <cmath>

using namespace ALM_NS;

Interaction::Interaction(ALM *alm) : Pointers(alm) {
    nsize[0] = nsize[1] = nsize[2] = 1;
}

Interaction::~Interaction() 
{
    memory->deallocate(xcrd);
    memory->deallocate(str_order);
    memory->deallocate(nbody_include);
    memory->deallocate(pairs);
    memory->deallocate(mindist_pairs);
    memory->deallocate(interaction_pair);
    memory->deallocate(mindist_cluster);
}

void Interaction::init()
{
    int i, j, k;
    int nat = system->nat;
    int nkd = system->nkd;

    std::cout << " INTERACTION" << std::endl;
    std::cout << " ===========" << std::endl << std::endl;

    memory->allocate(str_order, maxorder);
    set_ordername();

    std::cout << "  +++ Cutoff Radii Matrix in Bohr Unit (NKD x NKD matrix) +++" << std::endl;

    for (i = 0; i < maxorder; ++i) {
        std::cout << "  " <<  std::setw(9) << str_order[i] << std::endl; 
        for (j = 0; j < nkd; ++j) {
            for (k = 0; k < nkd; ++k) {
                if (rcs[i][j][k] < 0) {
                    std::cout << std::setw(9) << "None";
                } else {
                    std::cout << std::setw(9) << rcs[i][j][k];
                }
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }

    for (i = 0; i < 3; ++i) {
        if(!is_periodic[i]) nsize[i] = 0;
    }

    nneib = (2 * nsize[0] + 1) * (2 * nsize[1] + 1) * (2 * nsize[2] + 1);
    memory->allocate(mindist_pairs, nat, nat);
    memory->allocate(xcrd, nneib, nat, 3);
    memory->allocate(interaction_pair, maxorder, symmetry->natmin);
    memory->allocate(pairs, maxorder);
    memory->allocate(mindist_cluster, maxorder, symmetry->natmin);

    get_pairs_of_minimum_distance(nat, system->xcoord, mindist_pairs);
    print_neighborlist(mindist_pairs);
    search_interactions(interaction_pair, pairs);
    calc_mindist_clusters(interaction_pair, mindist_pairs, mindist_cluster);

    timer->print_elapsed();
    std::cout << " --------------------------------------------------------------" << std::endl;
    std::cout << std::endl;
}

double Interaction::distance(double *x1, double *x2)
{
    double dist;    
    dist = std::pow(x1[0] - x2[0], 2) + std::pow(x1[1] - x2[1], 2) + std::pow(x1[2] - x2[2], 2);
    dist = std::sqrt(dist);

    return dist;
}

void Interaction::get_pairs_of_minimum_distance(int nat, double **xf, std::vector<DistInfo> **mindist_pairs)
{
    int icell = 0;
    int i, j, k;
    int isize, jsize, ksize;
    double dist_tmp;
    double vec[3];
    std::vector<DistInfo> **distall;

    memory->allocate(distall, nat, nat);

    for (i = 0; i < nat; ++i) {
        for (j = 0; j < 3; ++j) {
            xcrd[0][i][j] = xf[i][j];
        }
    }

    for (isize = -nsize[0]; isize <= nsize[0] ; ++isize) {
        for (jsize = -nsize[1]; jsize <= nsize[1] ; ++jsize) {
            for (ksize = -nsize[2]; ksize <= nsize[2] ; ++ksize) {

                if (isize == 0 && jsize == 0 && ksize == 0) continue;

                ++icell;
                for (i = 0; i < nat; ++i) {
                    xcrd[icell][i][0] = xf[i][0] + static_cast<double>(isize);
                    xcrd[icell][i][1] = xf[i][1] + static_cast<double>(jsize);
                    xcrd[icell][i][2] = xf[i][2] + static_cast<double>(ksize);
                }
            }
        }
    }

    for (icell = 0; icell < nneib; ++icell) system->frac2cart(xcrd[icell]);

    for (i = 0; i < nat; ++i) {
        for (j = 0; j < nat; ++j) {
            for (icell = 0; icell < nneib; ++icell) {

                dist_tmp = distance(xcrd[0][i], xcrd[icell][j]);

                for (k = 0; k < 3; ++k) vec[k] = xcrd[icell][j][k] - xcrd[0][i][k];

                distall[i][j].push_back(DistInfo(icell, dist_tmp, vec));
            }
            std::sort(distall[i][j].begin(), distall[i][j].end());
        }
    }

    // Construct pairs of minimum distance.

    double dist_min;
    for (i = 0; i < nat; ++i) {
        for (j = 0; j < nat; ++j) {
            mindist_pairs[i][j].clear();

            dist_min = distall[i][j][0].dist;
            for (std::vector<DistInfo>::const_iterator it = distall[i][j].begin(); it != distall[i][j].end(); ++it) {
                if (std::abs((*it).dist - dist_min) < eps8) {
                    mindist_pairs[i][j].push_back(DistInfo(*it));
                }
            }
        }
    }

    memory->deallocate(distall);
}

void Interaction::print_neighborlist(std::vector<DistInfo> **mindist)
{
    int i, j, k;
    int iat;
    int nat = system->nat;
    int icount;

    double dist_tmp;
    std::vector<DistList> *neighborlist;

    memory->allocate(neighborlist, symmetry->natmin);

    for (i = 0; i < symmetry->natmin; ++i) {
        neighborlist[i].clear();

        iat = symmetry->map_p2s[i][0];

        for (j = 0; j < nat; ++j) {
            neighborlist[i].push_back(DistList(j, mindist[iat][j][0].dist));
        }
        std::sort(neighborlist[i].begin(), neighborlist[i].end());
    }

    std::cout << std::endl;
    std::cout << "  List of neighboring atoms below." << std::endl;
    std::cout << "  Format [N th-nearest shell, distance in Bohr (Number of atoms on the shell)]" << std::endl << std::endl;

    int nthnearest;
    std::vector<int> atomlist;

    for (i = 0; i < symmetry->natmin; ++i) {

        nthnearest = 0;
        atomlist.clear();

        iat = symmetry->map_p2s[i][0];
        std::cout << std::setw(5) << iat + 1 << " (" << std::setw(3) << system->kdname[system->kd[iat] - 1] << "): ";

        dist_tmp = 0.0;

        for (j = 0; j < nat; ++j) {

            if (neighborlist[i][j].dist < eps8) continue; // distance is zero

            if (std::abs(neighborlist[i][j].dist - dist_tmp) > eps6) {

                if (atomlist.size() > 0) {
                    nthnearest += 1;

                    if (nthnearest > 1) std::cout << std::setw(13) << " ";

                    std::cout << std::setw(3) << nthnearest << std::setw(10) << dist_tmp 
                        << " (" << std::setw(3) << atomlist.size() << ") -";

                    icount = 0;
                    for (k = 0; k < atomlist.size(); ++k) {

                        if (icount % 4 == 0 && icount > 0) {
                            std::cout << std::endl;
                            std::cout << std::setw(34) << " ";
                        }
                        ++icount;

                        std::cout << std::setw(4) << atomlist[k] + 1; 
                        std::cout <<  "(" << std::setw(3) << system->kdname[system->kd[atomlist[k]] - 1] << ")";

                    }
                    std::cout << std::endl;
                }


                dist_tmp = neighborlist[i][j].dist;
                atomlist.clear();
                atomlist.push_back(neighborlist[i][j].atom);
            } else {
                atomlist.push_back(neighborlist[i][j].atom);
            }
        }
        std::cout << std::endl;
    }
    atomlist.clear();
    memory->deallocate(neighborlist);
}

void Interaction::search_interactions(std::vector<int> **interaction_list_out, std::set<IntList> *pair_out)
{
    int i, j;
    int order;
    int natmin = symmetry->natmin;
    int iat, jat;
    int nat = system->nat;
    int ikd, jkd;

    double cutoff_tmp;
    std::set<IntList> *interacting_atom_pairs;
    std::vector<int> intlist;

    memory->allocate(interacting_atom_pairs, maxorder);

    for (order = 0; order < maxorder; ++order) {
        for (i = 0; i < natmin; ++i) {
            interaction_list_out[order][i].clear();

            iat = symmetry->map_p2s[i][0];
            ikd = system->kd[iat] - 1;

            for (jat = 0; jat < nat; ++jat) {

                jkd = system->kd[jat] - 1;

                cutoff_tmp = rcs[order][ikd][jkd];

                if (cutoff_tmp < 0.0) {

                    interaction_list_out[order][i].push_back(jat);

                } else {

                    if (mindist_pairs[iat][jat][0].dist <= cutoff_tmp) {
                        interaction_list_out[order][i].push_back(jat);
                    }
                }
            }            
        }
    }

    intlist.clear();
    std::cout << std::endl;
    std::cout << "  List of interacting atom pairs considered for each order:" << std::endl;
    for (order = 0; order < maxorder; ++order) {

        interacting_atom_pairs[order].clear();

        std::cout << std::endl << "   ***" << str_order[order] << "***" << std::endl;

        for (i = 0; i < natmin; ++i) {

            if (interaction_list_out[order][i].size() == 0) {
                std::cout << "   No interacting atoms! Skipped." << std::endl;
                continue; // no interaction
            }

            iat = symmetry->map_p2s[i][0];


            intlist.clear();
            for (std::vector<int>::const_iterator it = interaction_list_out[order][i].begin(); 
                it != interaction_list_out[order][i].end(); ++it) {
                    intlist.push_back((*it));
            }
            std::sort(intlist.begin(), intlist.end()); // Necessarily to sort here

            // write atoms inside the cutoff radius
            int id = 0;
            std::cout << "    Atom " << std::setw(5) << iat + 1
                << "(" << std::setw(3) << system->kdname[system->kd[iat]-1] << ")" << " interacts with atoms ... " << std::endl;

            for (int id = 0; id < intlist.size(); ++id) {
                if (id%6 == 0) {
                    if (id == 0) {
                        std::cout << "   ";
                    } else {
                        std::cout << std::endl;
                        std::cout << "   ";
                    }
                }
                std::cout << std::setw(5) << intlist[id] + 1 << "(" 
                    << std::setw(3) << system->kdname[system->kd[intlist[id]]-1] << ")";
            }

            std::cout << std::endl << std::endl;
            std::cout << "    Number of total interaction pairs = " 
                << interaction_list_out[order][i].size() << std::endl << std::endl;

            int *intarr;        
            memory->allocate(intarr, order + 2);

            if (intlist.size() > 0) {
                if (order == 0) {
                    for (unsigned int ielem = 0; ielem < intlist.size(); ++ielem) {
                        intarr[0] = iat;
                        intarr[1] = intlist[ielem];
                        insort(2, intarr);

                        interacting_atom_pairs[0].insert(IntList(2, intarr));
                    }
                } else if (order > 0) {
                    CombinationWithRepetition<int> g(intlist.begin(), intlist.end(), order + 1);
                    do {
                        std::vector<int> data = g.now();
                        intarr[0] = iat;
                        intarr[1] = data[0];
                        for (unsigned int isize = 1; isize < data.size() ; ++isize) {
                            intarr[isize + 1] = data[isize];
                        }

                        if (!is_incutoff(order+2, intarr)) continue;
                        insort(order+2, intarr);

                        interacting_atom_pairs[order].insert(IntList(order + 2, intarr));

                    } while(g.next());
                }
            }
            intlist.clear();
            memory->deallocate(intarr);
        }
    }

    std::cout << std::endl;
    int *pair_tmp;

    for (order = 0; order < maxorder; ++order) {

        pair_out[order].clear();

        if (order + 2 > nbody_include[order]) {
            std::cout << "  For " << std::setw(8) << interaction->str_order[order] << ", ";
            std::cout << "interactions related to more than" << std::setw(2) << nbody_include[order];
            std::cout << " atoms will be neglected." << std::endl;
        }

        memory->allocate(pair_tmp, order + 2);

        for (std::set<IntList>::const_iterator it = interacting_atom_pairs[order].begin(); 
            it != interacting_atom_pairs[order].end(); ++it) {
                for (j = 0; j < order + 2; ++j) {
                    pair_tmp[j] = (*it).iarray[j];
                }

                // Ignore many-body case 
                if (nbody(order + 2, pair_tmp) > nbody_include[order]) continue;

                pair_out[order].insert(IntList(order + 2, pair_tmp));
        }
        memory->deallocate(pair_tmp);
    }
    memory->deallocate(interacting_atom_pairs);
}

bool Interaction::is_incutoff(const int n, int *atomnumlist) 
{
    int i, j;
    int iat, jat, kat;
    int jkd, kkd;
    int ncheck = n - 1;
    int order = n - 2;   
    bool in_cutoff_tmp;
    double cutoff_tmp;
    double dist_tmp;
    std::vector<DistInfo>::const_iterator it, it2;

    iat = atomnumlist[0];

    for (i = 0; i < ncheck; ++i) {

        jat = atomnumlist[i + 1];
        jkd = system->kd[jat] - 1;

        for (j = i + 1; j < ncheck; ++j) {

            kat = atomnumlist[j + 1];
            kkd = system->kd[kat] - 1;

            cutoff_tmp = rcs[order][jkd][kkd];

            if (cutoff_tmp >= 0.0) {

                in_cutoff_tmp = false;

                for (it = mindist_pairs[iat][jat].begin(); it != mindist_pairs[iat][jat].end(); ++it) {
                    for (it2 = mindist_pairs[iat][kat].begin(); it2 != mindist_pairs[iat][kat].end(); ++it2) {
                        dist_tmp = distance(xcrd[(*it).cell][jat], xcrd[(*it2).cell][kat]);

                        if (dist_tmp <= cutoff_tmp) {
                            in_cutoff_tmp = true;
                        }
                    }
                }
                if (!in_cutoff_tmp) return false;
            }
        }
    }

    return true;
}

void Interaction::set_ordername()
{
    std::string strnum;

    str_order[0] = "HARMONIC";

    for (int i = 1;  i < maxorder; ++i) {
        strnum = boost::lexical_cast<std::string>(i+2);
        str_order[i] = "ANHARM" + strnum;
    }
}

int Interaction::nbody(const int n, const int *arr)
{
    std::vector<int> v;
    v.clear();
    int ret;

    for (unsigned int i = 0; i < n; ++i) {
        v.push_back(arr[i]);
    }
    std::sort(v.begin(), v.end());
    v.erase(std::unique(v.begin(), v.end()), v.end());

    ret = v.size();
    v.clear();

    return ret;
}

void Interaction::calc_mindist_clusters(std::vector<int> **interaction_pair_in, std::vector<DistInfo> **mindist_pair_in, 
                                        std::set<MinimumDistanceCluster> **mindist_cluster_out)
{
    std::vector<MinDistList> distance_list;

    int natmin = symmetry->natmin;
    int i, j, k;
    int iat, jat;
    int order;
    int ii;
    std::vector<int> intlist;

    double dist_tmp;
    double sum_dist_min, sum_dist;

    std::vector<int> cell_vector;
    std::vector<double> dist_vector;
    std::vector<std::vector<int> > pairs_icell, comb_cell, comb_cell_min;
    std::vector<int> accum_tmp;
    std::vector<int> atom_tmp, cell_tmp;


    for (order = 0; order < maxorder; ++order) {
        for (i = 0; i < symmetry->natmin; ++i) {

            mindist_cluster_out[order][i].clear();

            iat = symmetry->map_p2s[i][0];

            intlist.clear();
            for (std::vector<int>::const_iterator it  = interaction_pair_in[order][i].begin(); 
                                                  it != interaction_pair_in[order][i].end(); ++it) {
                    intlist.push_back((*it));
            }
            std::sort(intlist.begin(), intlist.end()); // Necessarily to sort here

            if (intlist.size() > 0) {
                if (order == 0){
                    for (unsigned int ielem = 0; ielem < intlist.size(); ++ielem) {

                        comb_cell_min.clear();
                        atom_tmp.clear();

                        jat = intlist[ielem];
                        atom_tmp.push_back(jat);

                        for (j = 0; j < mindist_pair_in[iat][jat].size(); ++j) {
                            cell_tmp.clear();
                            cell_tmp.push_back(mindist_pair_in[iat][jat][j].cell);
                            comb_cell_min.push_back(cell_tmp);
                        }
                        mindist_cluster_out[order][i].insert(MinimumDistanceCluster(atom_tmp, comb_cell_min));
                    }
                } else if (order > 0) {
                    CombinationWithRepetition<int> g(intlist.begin(), intlist.end(), order + 1);
                    do {
                        std::vector<int> data = g.now();

                        pairs_icell.clear();
                        for (j = 0; j < order + 1; ++j) {
                            jat = data[j];
                            cell_vector.clear();

                            for (ii = 0; ii < mindist_pair_in[iat][jat].size(); ++ii) {
                                cell_vector.push_back(mindist_pair_in[iat][jat][ii].cell);
                            }
                            pairs_icell.push_back(cell_vector);
                        }

                        accum_tmp.clear();
                        comb_cell.clear();
                        cell_combination(pairs_icell, 0, accum_tmp, comb_cell);

                        distance_list.clear();
                        for (j = 0; j < comb_cell.size(); ++j) {
                            dist_vector.clear();

                            for (k = 0; k < comb_cell[j].size(); ++k) {
                                for (ii = k + 1; ii < comb_cell[j].size(); ++ii) {
                                    dist_tmp = distance(xcrd[comb_cell[j][k]][data[k]], xcrd[comb_cell[j][ii]][data[ii]]);
                                    dist_vector.push_back(dist_tmp);
                                }
                            }
                            distance_list.push_back(MinDistList(comb_cell[j], dist_vector));
                        }

                        std::sort(distance_list.begin(), distance_list.end(), MinDistList::compare_sum_distance);

                        comb_cell_min.clear();

                        sum_dist_min = 0.0;
                        for (j = 0; j < distance_list[0].dist.size(); ++j) {
                            sum_dist_min += distance_list[0].dist[j];
                        }

                        for (j = 0; j < distance_list.size(); ++j) {
                            sum_dist = 0.0;

                            for (k = 0; k < distance_list[j].dist.size(); ++k) {
                                sum_dist += distance_list[j].dist[k];
                            }

                            // In the following, only pairs having minimum sum of distances
                            // are stored. However, we found that this treatment didn't
                            // return a reliable value of phonon linewidth.
                            
                            // if (std::abs(sum_dist - sum_dist_min) < eps8) {
                            //    comb_cell_min.push_back(distance_list[j].cell);
                            // } else {
                            //    break;
                            // }

                            // Therefore, we consider all duplicate pairs
                             
                            comb_cell_min.push_back(distance_list[j].cell);

                        }

                        mindist_cluster_out[order][i].insert(MinimumDistanceCluster(data, comb_cell_min));

                    } while(g.next());
                }
            }
            intlist.clear();
        }
    }
}


void Interaction::cell_combination(std::vector<std::vector<int> > array, int i, std::vector<int> accum, std::vector<std::vector<int> > &comb)
{
    if (i == array.size())  {
        comb.push_back(accum); 
    } else  {
        std::vector<int> row = array[i];
        for (int j = 0; j < row.size(); ++j) {
            std::vector<int> tmp(accum);
            tmp.push_back(row[j]);
            cell_combination(array,i+1,tmp, comb);
        }
    }
}
