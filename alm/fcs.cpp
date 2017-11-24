/*
 fcs.cpp

 Copyright (c) 2014, 2015, 2016 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory
 or http://opensource.org/licenses/mit-license.php for information.
*/

#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>
#include "../external/combination.hpp"
#include <boost/lexical_cast.hpp>
#include "files.h"
#include "interaction.h"
#include "error.h"
#include "memory.h"
#include "fcs.h"
#include "symmetry.h"
#include "system.h"
#include "timer.h"
#include "constants.h"

using namespace ALM_NS;

Fcs::Fcs(ALM *alm) : Pointers(alm)
{
};

Fcs::~Fcs()
{
    if (alm->mode == "fitting") {
        memory->deallocate(fc_table);
        memory->deallocate(nequiv);
        memory->deallocate(fc_zeros);
    }
};

void Fcs::init()
{
    int i;
    int maxorder = interaction->maxorder;

    std::cout << " FORCE CONSTANT" << std::endl;
    std::cout << " ==============" << std::endl << std::endl;

    memory->allocate(fc_table, maxorder);
    memory->allocate(nequiv, maxorder);
    memory->allocate(fc_zeros, maxorder);

    for (i = 0; i < maxorder; ++i) {
        generate_force_constant_table(i, interaction->pairs[i],
                                      symmetry->SymmData, "Cartesian",
                                      fc_table[i], nequiv[i], fc_zeros[i], true);
    }

    std::cout << std::endl;
    for (i = 0; i < maxorder; ++i) {
        std::cout << "  Number of " << std::setw(9)
            << interaction->str_order[i]
            << " FCs : " << nequiv[i].size();
        std::cout << std::endl;
    }
    std::cout << std::endl;


    //
    //    std::vector<FcProperty> fc_test, fc_zeros_tmp;
    //    std::vector<int> nmulti;
    //    generate_force_constant_table(1, interaction->pairs[1],
    //                                  symmetry->SymmData, "Lattice",
    //                                  fc_test, nmulti, fc_zeros_tmp, true);
    //    std::cout << "Nonzero independent IFCs:" << std::setw(5) << nmulti.size() << std::endl;
    //    std::cout << "Zero IFCs:" << std::endl;
    //    for (auto it = fc_zeros_tmp.begin(); it != fc_zeros_tmp.end(); ++it) {
    //        for (i = 0; i < (*it).elems.size(); ++i) {
    //            std::cout << std::setw(4) << (*it).elems[i];
    //        }
    //        std::cout << std::setw(5) << (*it).mother;
    //        std::cout << std::endl;
    //    }
    //    std::cout << std::endl;


    timer->print_elapsed();
    std::cout << " -------------------------------------------------------------------" << std::endl;
    std::cout << std::endl;
}

void Fcs::generate_force_constant_table(const int order,
                                        const std::set<IntList> pairs,
                                        const std::vector<SymmetryOperation> symmop,
                                        std::string basis,
                                        std::vector<FcProperty> &fc_vec,
                                        std::vector<int> &ndup,
                                        std::vector<FcProperty> &fc_zeros,
                                        const bool store_zeros)
{
    int i, j;
    int i1, i2;
    int i_prim;
    int *atmn, *atmn_mapped;
    int *ind, *ind_mapped;
    int *ind_tmp, *ind_mapped_tmp;
    int nxyz;
    unsigned int isym;

    double c_tmp;

    int **xyzcomponent;

    int nmother;
    int nat = system->nat;
    int nsym = symmop.size();
    int nsym_in_use;

    bool is_zero;
    bool *is_searched;
    int counter;
    int **map_sym;
    double ***rotation;

    if (order < 0) return;

    memory->allocate(rotation, nsym, 3, 3);
    memory->allocate(map_sym, nat, nsym);
    nsym_in_use = 0;
    counter = 0;
    if (basis == "Cartesian") {

        for (auto it = symmop.begin(); it != symmop.end(); ++it) {
            if ((*it).compatible_with_cartesian) {
                for (i = 0; i < 3; ++i) {
                    for (j = 0; j < 3; ++j) {
                        rotation[nsym_in_use][i][j] = (*it).rotation_cart[i][j];
                    }
                }
                for (i = 0; i < nat; ++i) {
                    map_sym[i][nsym_in_use] = symmetry->map_sym[i][counter];
                }
                ++nsym_in_use;
            }
            ++counter;
        }

    } else if (basis == "Lattice") {

        for (auto it = symmop.begin(); it != symmop.end(); ++it) {
            if ((*it).compatible_with_lattice) {
                for (i = 0; i < 3; ++i) {
                    for (j = 0; j < 3; ++j) {
                        rotation[nsym_in_use][i][j]
                            = static_cast<double>((*it).rotation[i][j]);
                    }
                }
                for (i = 0; i < nat; ++i) {
                    map_sym[i][nsym_in_use] = symmetry->map_sym[i][counter];
                }
                ++nsym_in_use;
            }
            ++counter;
        }


    } else {
        memory->deallocate(rotation);
        memory->deallocate(map_sym);
        error->exit("generate_force_constant_table", "Invalid basis inpout");
    }

    memory->allocate(atmn, order + 2);
    memory->allocate(atmn_mapped, order + 2);
    memory->allocate(ind, order + 2);
    memory->allocate(ind_mapped, order + 2);
    memory->allocate(ind_tmp, order);
    memory->allocate(ind_mapped_tmp, order + 2);
    memory->allocate(is_searched, 3 * nat);

    fc_vec.clear();
    ndup.clear();
    fc_zeros.clear();
    nmother = 0;

    nxyz = static_cast<int>(std::pow(3.0, order + 2));

    memory->allocate(xyzcomponent, nxyz, order + 2);
    get_xyzcomponent(order + 2, xyzcomponent);

    std::set<IntList> list_found;

    for (auto iter = pairs.begin(); iter != pairs.end(); ++iter) {

        for (i = 0; i < order + 2; ++i) atmn[i] = (*iter).iarray[i];

        for (i1 = 0; i1 < nxyz; ++i1) {
            for (i = 0; i < order + 2; ++i) ind[i] = 3 * atmn[i] + xyzcomponent[i1][i];

            if (!is_ascending(order + 2, ind)) continue;

            i_prim = min_inprim(order + 2, ind);
            std::swap(ind[0], ind[i_prim]);
            sort_tail(order + 2, ind);

            is_zero = false;

            if (list_found.find(IntList(order + 2, ind)) != list_found.end()) continue; // Already exits!

            // Search symmetrically-dependent parameter set

            int ndeps = 0;

            for (isym = 0; isym < nsym_in_use; ++isym) {

                for (i = 0; i < order + 2; ++i)
                    atmn_mapped[i] = map_sym[atmn[i]][isym];

                if (!is_inprim(order + 2, atmn_mapped)) continue;

                for (i2 = 0; i2 < nxyz; ++i2) {
                    c_tmp = coef_sym(order + 2, rotation[isym], xyzcomponent[i1], xyzcomponent[i2]);
                    if (std::abs(c_tmp) > eps12) {
                        for (i = 0; i < order + 2; ++i)
                            ind_mapped[i] = 3 * atmn_mapped[i] + xyzcomponent[i2][i];

                        i_prim = min_inprim(order + 2, ind_mapped);
                        std::swap(ind_mapped[0], ind_mapped[i_prim]);
                        sort_tail(order + 2, ind_mapped);

                        if (!is_zero) {
                            bool zeroflag = true;
                            for (i = 0; i < order + 2; ++i) {
                                zeroflag = zeroflag & (ind[i] == ind_mapped[i]);
                            }
                            zeroflag = zeroflag & (std::abs(c_tmp + 1.0) < eps8);
                            is_zero = zeroflag;
                        }

                        // Add to found list (set) and fcset (vector) if the created is new one.

                        if (list_found.find(IntList(order + 2, ind_mapped)) == list_found.end()) {
                            list_found.insert(IntList(order + 2, ind_mapped));

                            fc_vec.push_back(FcProperty(order + 2, c_tmp,
                                                        ind_mapped, nmother));
                            ++ndeps;

                            // Add equivalent interaction list (permutation) if there are two or more indices
                            // which belong to the primitive cell.
                            // This procedure is necessary for fitting.

                            for (i = 0; i < 3 * nat; ++i) is_searched[i] = false;
                            is_searched[ind_mapped[0]] = true;
                            for (i = 1; i < order + 2; ++i) {
                                if ((!is_searched[ind_mapped[i]]) && is_inprim(ind_mapped[i])) {

                                    for (j = 0; j < order + 2; ++j) ind_mapped_tmp[j] = ind_mapped[j];
                                    std::swap(ind_mapped_tmp[0], ind_mapped_tmp[i]);
                                    sort_tail(order + 2, ind_mapped_tmp);
                                    fc_vec.push_back(FcProperty(order + 2, c_tmp,
                                                                ind_mapped_tmp, nmother));

                                    ++ndeps;

                                    is_searched[ind_mapped[i]] = true;
                                }
                            }


                        }
                    }
                }
            } // close symmetry loop

            if (is_zero) {
                if (store_zeros) {
                    for (auto it = fc_vec.rbegin(); it != fc_vec.rbegin() + ndeps; ++it) {
                        (*it).mother = -1;
                        fc_zeros.push_back(*it);
                    }
                }
                for (i = 0; i < ndeps; ++i) fc_vec.pop_back();
            } else {
                ndup.push_back(ndeps);
                ++nmother;
            }

        } // close xyz component loop
    } // close atom number loop (iterator)

    memory->deallocate(xyzcomponent);
    list_found.clear();
    memory->deallocate(atmn);
    memory->deallocate(atmn_mapped);
    memory->deallocate(ind);
    memory->deallocate(ind_mapped);
    memory->deallocate(ind_tmp);
    memory->deallocate(ind_mapped_tmp);
    memory->deallocate(is_searched);
    memory->deallocate(rotation);
    memory->deallocate(map_sym);

    // sort fc_vec

    if (ndup.size() > 0) {
        std::sort(fc_vec.begin(), fc_vec.begin() + ndup[0]);
        int nbegin = ndup[0];
        int nend;
        for (int mm = 1; mm < ndup.size(); ++mm) {
            nend = nbegin + ndup[mm];
            std::sort(fc_vec.begin() + nbegin, fc_vec.begin() + nend);
            nbegin += ndup[mm];
        }
    }
}


double Fcs::coef_sym(const int n,
                     const int symnum,
                     const int *arr1,
                     const int *arr2)
{
    double tmp = 1.0;
    int i;

    for (i = 0; i < n; ++i) {
        tmp *= symmetry->SymmData[symnum].rotation_cart[arr2[i]][arr1[i]];
    }
    return tmp;
}

double Fcs::coef_sym(const int n,
                     double **rot,
                     const int *arr1,
                     const int *arr2)
{
    double tmp = 1.0;
    int i;

    for (i = 0; i < n; ++i) {
        tmp *= rot[arr2[i]][arr1[i]];
    }
    return tmp;
}

bool Fcs::is_ascending(const int n, const int *arr)
{
    int i;
    for (i = 0; i < n - 1; ++i) {
        if (arr[i] > arr[i + 1]) return false;
    }
    return true;
}

int Fcs::min_inprim(const int n, const int *arr)
{
    int i, j, atmnum;
    int natmin = symmetry->nat_prim;
    int minloc;
    int *ind;

    memory->allocate(ind, n);

    for (i = 0; i < n; ++i) {

        ind[i] = 3 * system->nat;
        atmnum = arr[i] / 3;

        for (j = 0; j < natmin; ++j) {
            if (symmetry->map_p2s[j][0] == atmnum) {
                ind[i] = arr[i];
                continue;
            }
        }
    }

    int minval = ind[0];
    minloc = 0;

    for (i = 0; i < n; ++i) {
        if (ind[i] < minval) {
            minval = ind[i];
            minloc = i;
        }
    }

    memory->deallocate(ind);
    return minloc;
}

bool Fcs::is_inprim(const int n, const int *arr)
{
    int i, j;
    int natmin = symmetry->nat_prim;

    for (i = 0; i < n; ++i) {
        for (j = 0; j < natmin; ++j) {
            if (symmetry->map_p2s[j][0] == arr[i]) return true;
        }
    }
    return false;
}

bool Fcs::is_inprim(const int n)
{
    int i, atmn;
    int natmin = symmetry->nat_prim;

    atmn = n / 3;

    for (i = 0; i < natmin; ++i) {
        if (symmetry->map_p2s[i][0] == atmn) return true;
    }

    return false;
}

void Fcs::get_xyzcomponent(int n, int **xyz)
{
    // Return xyz component for the given order using boost algorithm library

    std::vector<int> v;
    int i;

    for (i = 0; i < n; ++i) {
        v.push_back(0);
        v.push_back(1);
        v.push_back(2);
    }

    std::sort(v.begin(), v.end());

    int m = 0;

    do {
        xyz[m][0] = v[0];
        for (i = 1; i < n; ++i) xyz[m][i] = v[i];
        ++m;
    } while (boost::next_partial_permutation(v.begin(), v.begin() + n, v.end()));
}

void Fcs::sort_tail(const int n, int *arr)
{
    int i, m;

    m = n - 1;
    int *ind_tmp;

    memory->allocate(ind_tmp, m);

    for (i = 0; i < m; ++i) {
        ind_tmp[i] = arr[i + 1];
    }

    interaction->insort(m, ind_tmp);

    for (i = 0; i < m; ++i) {
        arr[i + 1] = ind_tmp[i];
    }

    memory->deallocate(ind_tmp);
}

std::string Fcs::easyvizint(const int n)
{
    int atmn;
    int crdn;
    atmn = n / 3 + 1;
    crdn = n % 3;
    std::string str_crd[3] = {"x", "y", "z"};
    std::string str_tmp;

    str_tmp = boost::lexical_cast<std::string>(atmn);
    str_tmp += str_crd[crdn];

    return str_tmp;
}
