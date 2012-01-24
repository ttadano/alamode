#include "files.h"
#include "interaction.h"
#include "error.h"
#include "memory.h"
#include <iostream>
#include "fcs.h"
#include "symmetry.h"
#include <cmath>
#include <boost/algorithm/combination.hpp>

using namespace ALM_NS;

Fcs::Fcs(ALM *alm) : Pointers(alm){};
Fcs::~Fcs() {};

void Fcs::init(){

    int i, j;
    int maxorder = interaction->maxorder;

    nints = new int [maxorder];
    //    pairs = new std::vector<Pairs>[maxorder];
    pairs = new std::set<IntList>[maxorder];
    read_pairs(maxorder);
    generate_fclists(maxorder);
}

void Fcs::read_pairs(int maxorder)
{

    int i, j, order;
    int *pair_tmp;

    files->ifs_int.open(files->file_int, std::ios::out);
    if(!files->ifs_int) error->exit("read_pairs", "cannot open file_int");

    for(order = 0; order < maxorder; ++order){
        files->ifs_int >> nints[order];
        if(nints[order] == 0) continue;

        pair_tmp = new int [order + 2];
        for(i = 0; i < nints[order]; ++i){
            for(j = 0; j < order + 2; ++j){
                files->ifs_int >> pair_tmp[j];
            }
            //       pairs[order].push_back(Pairs(order + 2, pair_tmp));
            pairs[order].insert(IntList(order + 2, pair_tmp));
        }
        delete pair_tmp;
        std::cout << pairs[order].size() << std::endl;
    }
}

void Fcs::generate_fclists(int maxorder)
{
    int i, j;
    int order;
    int i_prim;
    int *atmn, *atmn_mapped;
    int *ind, *ind_mapped;
    int nxyz;
    IntList list_tmp;

    double c_tmp;

    int **xyzcomponent;

    atmn = new int [maxorder];
    atmn_mapped = new int [maxorder];
    ind = new int  [maxorder];
    ind_mapped = new int [maxorder];

    for(order = 0; order < maxorder; ++order){

        nxyz = static_cast<int>(pow(static_cast<long double>(3), order + 2));

        memory->allocate(xyzcomponent, nxyz, order + 2);
        for (i = 0; i < nxyz; ++i){
            std::cout << std::setw(10) << &xyzcomponent[i][0] << std::setw(10) << &xyzcomponent[i][1] << std::endl;
        }
        get_xyzcomponent(order + 2, xyzcomponent);
        memory->deallocate(xyzcomponent);
        exit(1);
        std::cout << "OK" << std::endl;

        for(std::set<IntList>::iterator iter = pairs[order].begin(); iter != pairs[order].end(); ++iter){

            list_tmp = *iter;
            for (i = 0; i < order + 2; ++i){
                atmn[i] = list_tmp.iarray[i];
            }

            for (int i1 = 0; i1 < nxyz; ++i1){
                for (i = 0; i < order + 2; ++i){
                    ind[i] = 3 * atmn[i] + xyzcomponent[i1][i];
                }

                if(!is_ascending(order + 2, ind)) continue;
                i_prim = min_inprim(order + 2, ind);
                std::swap(ind[0], ind[i_prim]);

                // search symmetrycally-dependent parameter set

                int ndeps = 0;

                for(int isym = 0; isym < symmetry->nsym; ++isym){
                    for(i = 0; i < order + 2; ++i){
                        atmn_mapped[i] = symmetry->map_sym[atmn[i]][isym];
                    }

                    if(!is_inprim(order + 2, atmn_mapped)) continue;

                    for (int i2 = 0; i2 < nxyz; ++i2){
                        c_tmp = coef_sym(order + 2, isym, xyzcomponent[i1], xyzcomponent[i2]);
                        if(abs(c_tmp) > 1.0e-8) {
                            for (i = 0; i < order + 2; ++i){
                                ind_mapped[i] = 3 * atmn_mapped[i] + xyzcomponent[i2][i];
                            }

                            i_prim = min_inprim(order + 2, ind_mapped);
                            std::swap(ind_mapped[0], ind_mapped[i_prim]);
                            ++ndeps;
                        }
                    }
                }
                //           std::cout << "ndeps" << ndeps << std::endl;
            }
        }
        memory->deallocate(xyzcomponent);
        std::cout << "OK2" << std::endl;
    }
}

double Fcs::coef_sym(const int n, const int symnum, const int *arr1, const int *arr2)
{
    double tmp = 1.0;
    int i;

    for (i = 0; i < n; ++i){
        tmp *= symmetry->symrel[symnum][arr2[i]][arr1[i]];
    }
    return tmp;

    //   if(symnum == 0) std::cout << symmetry->symrel[0][0][0] * symmetry->symrel[0][0][0] << std::endl;
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

    for (i = 0; i < n; ++i){
        atmnum = arr[i] / 3;
        for (j = 0; j < natmin; ++j){
            if(symmetry->map_p2s[j][0] == atmnum) return i; 
        }
    }
    // this cannot happen
    error->exit("min_inprim", "no indecis in the primitive cell");
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
        for (i = 0; i < n; ++i) {
            std::cout << std::setw(10) << &xyz[m][i];
        }
        std::cout << std::endl;
        ++m;
    } while(boost::next_partial_permutation(v.begin(), v.begin() + n, v.end()));
}
