#include <iostream>
#include <string>
#include <cmath>
#include <boost/algorithm/combination.hpp>
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

Fcs::Fcs(ALM *alm) : Pointers(alm){};
Fcs::~Fcs() {};

void Fcs::init(){

    int i;
    int maxorder = interaction->maxorder;
    memory->allocate(nints, maxorder);
    memory->allocate(pairs, maxorder);
    memory->allocate(nzero, maxorder);


    for(i = 0; i < maxorder; ++i) nzero[i] = 0;

    read_pairs(maxorder);

    fc_set = new std::vector<FcProperty> [maxorder];
    ndup = new std::vector<int> [maxorder];
    generate_fclists(maxorder);

    std::cout << std::endl;
    for(i = 0; i < maxorder; ++i){
        std::cout << "Number of " << std::setw(9) << interaction->str_order[i] << " FCs (nzero): " << ndup[i].size();
        std::cout << " ( " << nzero[i] << " ) " << std::endl;
    }
    std::cout << std::endl;

    // sort fc_set

    std::cout << "Sorting interaction arrays ...";
    for(int order = 0; order < maxorder; ++order){
        if(ndup[order].size() > 0) {
            std::sort(fc_set[order].begin(), fc_set[order].begin() + ndup[order][0]);
            int nbegin = ndup[order][0];
            int nend;
            for(int mm = 1; mm < ndup[order].size(); ++mm){
                nend  = nbegin + ndup[order][mm];
                std::sort(fc_set[order].begin() + nbegin, fc_set[order].begin() + nend);
                nbegin += ndup[order][mm];
            }
        }
    }
    std::cout << " done." << std::endl;

    memory->deallocate(nints);
    memory->deallocate(pairs);
    memory->deallocate(nzero);
    timer->print_elapsed();
}

void Fcs::read_pairs(int maxorder)
{

    int i, j, order;
    int *pair_tmp;

    files->ifs_int.open(files->file_int.c_str(), std::ios::out);
    if(!files->ifs_int) error->exit("read_pairs", "cannot open file_int");

    for(order = 0; order < maxorder; ++order){
        files->ifs_int >> nints[order];
        if(nints[order] == 0) continue;

        memory->allocate(pair_tmp, order + 2);

        for(i = 0; i < nints[order]; ++i){
            for(j = 0; j < order + 2; ++j){
                files->ifs_int >> pair_tmp[j];
            }
            pairs[order].insert(IntList(order + 2, pair_tmp));
        }
        memory->deallocate(pair_tmp);
    }

    files->ifs_int.close();
}

void Fcs::generate_fclists(int maxorder)
{
    int i, j;
    int i1, i2;
    int order;
    int i_prim;
    int *atmn, *atmn_mapped;
    int *ind, *ind_mapped;
    int *ind_tmp, *ind_mapped_tmp;
    int nxyz;
    IntList list_tmp;

    double c_tmp;

    int **xyzcomponent;

    int nmother;
    int nat = system->nat;

    bool is_zero;
    bool *is_searched;

    std::cout << "Generating Symmetrically-Independent Parameters ..." << std::endl;

    memory->allocate(atmn, maxorder + 1);
    memory->allocate(atmn_mapped, maxorder + 1);
    memory->allocate(ind, maxorder + 1);
    memory->allocate(ind_mapped, maxorder + 1);
    memory->allocate(ind_tmp, maxorder - 1);
    memory->allocate(ind_mapped_tmp, maxorder + 1);
    memory->allocate(is_searched, 3 * nat);

    for(order = 0; order < maxorder; ++order){

        std::cout << std::setw(8) << interaction->str_order[order] << " ...";

        fc_set[order].clear();
        ndup[order].clear();
        nmother = 0;

        nxyz = static_cast<int>(pow(static_cast<long double>(3), order + 2));

        memory->allocate(xyzcomponent, nxyz, order + 2);
        get_xyzcomponent(order + 2, xyzcomponent);

        std::set<IntList> list_found;

        for (std::set<IntList>::iterator iter = pairs[order].begin(); iter != pairs[order].end(); ++iter){

            IntList list_tmp = *iter;
            for (i = 0; i < order + 2; ++i) atmn[i] = list_tmp.iarray[i];

            for (i1 = 0; i1 < nxyz; ++i1){
                for (i = 0; i < order + 2; ++i) ind[i] = 3 * atmn[i] + xyzcomponent[i1][i];

                if (!is_ascending(order + 2, ind)) continue;

                i_prim = min_inprim(order + 2, ind);
                std::swap(ind[0], ind[i_prim]);
                sort_tail(order + 2, ind);

                is_zero = false;

                if(list_found.find(IntList(order + 2, ind)) != list_found.end()) continue; // already exits!

                // search symmetrically-dependent parameter set

                int ndeps = 0;

                for (int isym = 0; isym < symmetry->nsym; ++isym){

                    for (i = 0; i < order + 2; ++i) atmn_mapped[i] = symmetry->map_sym[atmn[i]][isym];
                    if (!is_inprim(order + 2, atmn_mapped)) continue;

                    for (i2 = 0; i2 < nxyz; ++i2){
                        c_tmp = coef_sym(order + 2, isym, xyzcomponent[i1], xyzcomponent[i2]);
                        if (std::abs(c_tmp) > eps12) {
                            for (i = 0; i < order + 2; ++i) ind_mapped[i] = 3 * atmn_mapped[i] + xyzcomponent[i2][i];

                            i_prim = min_inprim(order + 2, ind_mapped);
                            std::swap(ind_mapped[0], ind_mapped[i_prim]);
                            sort_tail(order + 2, ind_mapped);

                            if (!is_zero){
                                bool zeroflag = true;
                                for (i = 0; i < order + 2; ++i){
                                    zeroflag = zeroflag & (ind[i] == ind_mapped[i]);
                                }
                                zeroflag = zeroflag & (std::abs(c_tmp + 1.0) < eps8);
                                is_zero = zeroflag;
                            }

                            // add to found list (set) and fcset (vector) if the created is new one.

                            if (list_found.find(IntList(order + 2, ind_mapped)) == list_found.end()){
                                list_found.insert(IntList(order + 2, ind_mapped));

                                fc_set[order].push_back(FcProperty(order + 2, c_tmp, ind_mapped, nmother));
                                ++ndeps;

                                // Add equivalent interaction list (permutation) if there are two or more indices
                                // which belong to the primitive cell.
                                // This procedure is necessary for fitting.

                                for (i = 0; i < 3 * nat; ++i) is_searched[i] = false;
                                is_searched[ind_mapped[0]] = true;
                                for (i = 1; i < order + 2; ++i){
                                    if((!is_searched[ind_mapped[i]]) && is_inprim(ind_mapped[i])){

                                        for (j = 0; j < order + 2; ++j) ind_mapped_tmp[j] = ind_mapped[j];
                                        std::swap(ind_mapped_tmp[0], ind_mapped_tmp[i]);
                                        sort_tail(order + 2, ind_mapped_tmp);
                                        fc_set[order].push_back(FcProperty(order + 2, c_tmp, ind_mapped_tmp, nmother));
                                        ++ndeps;

                                        is_searched[ind_mapped[i]] = true;
                                    }
                                }


                            }
                        }
                    }
                } // close symmetry loop

                if(is_zero){
                    for (i = 0; i < ndeps; ++i) fc_set[order].pop_back();
                    ++nzero[order];
                } else {             
                    ndup[order].push_back(ndeps);
                    ++nmother;
                }

            } // close xyz component loop
        } // close atom number loop (iterator)

        /*std::cout << "ORDER: " << order << " Size: " << ndup[order].size() << std::endl;
        for(unsigned int m = 0; m < ndup[order].size(); ++m){
        std::cout << "ORDER: " << order << std::setw(5) << ndup[order][m] << std::endl;
        }
        std::cout << "ORDER: " << order << " nzero " << nzero[order]<< std::endl;

        std::cout << "**ORDER = " << order << " **" << std::endl;
        std::cout << "Number of Parameters = " << ndup[order].size() << std::endl;

        int mmm = 0;
        for (unsigned int m = 0; m < ndup[order].size(); ++m){
        std::cout << "#" << m << " ndup = " << ndup[order][m] << std::endl;
        for (unsigned int mm = 0; mm < ndup[order][m]; ++mm){
        for (i = 0; i < order + 2; ++i){
        std::cout << std::setw(6) << easyvizint(fc_set[order][mmm].elems[i]);    
        }
        std::cout << " " << fc_set[order][mmm].mother << std::endl;
        ++mmm;
        } 
        std::cout << std::endl;
        }
        */
        memory->deallocate(xyzcomponent);
        list_found.clear();
        std::cout << ".. done. " << std::endl;
    } //close order loop

    memory->deallocate(atmn);
    memory->deallocate(atmn_mapped);
    memory->deallocate(ind);
    memory->deallocate(ind_mapped);
    memory->deallocate(ind_tmp);
    memory->deallocate(ind_mapped_tmp);
    memory->deallocate(is_searched);

    std::cout << "Finished !" << std::endl;
}

double Fcs::coef_sym(const int n, const int symnum, const int *arr1, const int *arr2)
{
    double tmp = 1.0;
    int i;

    for (i = 0; i < n; ++i){
        tmp *= symmetry->symrel[symnum][arr2[i]][arr1[i]];
    }
    return tmp;
}

bool Fcs::is_ascending(const int n, const int *arr)
{
    int i;
    for (i = 0; i < n - 1; ++i){
        if(arr[i] > arr[i+1]) return false;
    }
    return true;
}

int Fcs::min_inprim(const int n, const int *arr)
{
    int i, j, atmnum;
    int natmin = symmetry->natmin;
    int minloc;
    int *ind;

    memory->allocate(ind, n);

    for (i = 0; i < n; ++i){

        ind[i] = system->nat;
        atmnum = arr[i] / 3;

        for (j = 0; j < natmin; ++j){
            if (symmetry->map_p2s[j][0] == atmnum) {
                ind[i] = arr[i];
                continue;
            }
        }
    }

    int minval = ind[0];
    minloc = 0;

    for (i = 0; i < n; ++i){
        if(ind[i] < minval){
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
    int natmin = symmetry->natmin;

    for (i = 0; i < n; ++i){
        for (j = 0; j < natmin; ++j){
            if(symmetry->map_p2s[j][0] == arr[i]) return true;
        }
    }
    return false;
}

bool Fcs::is_inprim(const int n){

    int i, atmn;
    int natmin = symmetry->natmin;

    atmn = n / 3;

    for (i = 0; i < natmin; ++i){
        if(symmetry->map_p2s[i][0] == atmn) return true;
    }

    return false;
}

void Fcs::get_xyzcomponent(int n, int **xyz)
{
    // return xyz component for the given order using boost algorithm library

    std::vector<int> v;
    int i;

    for(i = 0; i < n; ++i){
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
    } while(boost::next_partial_permutation(v.begin(), v.begin() + n, v.end()));
}

void Fcs::sort_tail(const int n, int *arr){

    int i, m;

    m = n - 1;
    int *ind_tmp;

    memory->allocate(ind_tmp, m);

    for (i = 0; i < m; ++i){
        ind_tmp[i] = arr[i + 1];
    }

    interaction->insort(m, ind_tmp);

    for (i = 0; i < m; ++i){
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

    //   str_tmp = std::to_string(static_cast<long double>(atmn));
    str_tmp = boost::lexical_cast<std::string>(atmn);
    str_tmp += str_crd[crdn];

    return  str_tmp;
}
