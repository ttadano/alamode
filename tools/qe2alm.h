/*
 qe2alm.h

 Copyright (c) 2015 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#pragma once

#include <string>
#include <fstream>
#include <stdlib.h>
#include <new>
#include <iostream>
#include <vector>

std::string file_fc2qe, file_xmlout;
std::string flag_asr;
std::ifstream ifs_fc2qe;

int i, j, k;
int natmin, nkd, ibrav;
int nsize, nat;
double celldm[6];
double lavec[3][3], rlavec[3][3];

int *kd, **kd_super;
int **map_p2s;
double **xcrd, **xcrd_super;
std::string *kd_symbol;

bool include_na;
bool consider_asr;
    ;
int nq[3];
double dielec[3][3];
double ***born;
double ***fc2;

class DistInfo {
public:
    int cell;
    double dist;
    double relvec[3];
    
    DistInfo();
    DistInfo(const int n, const double d, const double x[3]) {
        cell = n;
        dist = d;
        for (int i = 0; i < 3; ++i) relvec[i] = x[i];
    }

    DistInfo(const DistInfo &obj) {
        cell = obj.cell;
        dist = obj.dist;
        for (int i = 0; i < 3; ++i) relvec[i] = obj.relvec[i];
    }
};

inline bool operator<(const DistInfo a, const DistInfo b) {
    return a.dist < b.dist;
}

std::vector<DistInfo> **mindist;

void calc_lattice_vector(const int, double [6], double [3][3]);
void recips(double [3][3], double[3][3]);
std::string double2string(const double );
void get_pairs_of_minimum_distance(const int, const int, int **, double **, std::vector<DistInfo> **);
double distance(double *, double*);

// memsize calculator
    
unsigned long memsize_in_MB(const int size_of_one, const int n1){
    unsigned long n = n1 * size_of_one;
    return n / 1000000;
}
unsigned long memsize_in_MB(const int size_of_one, const int n1, const int n2){
    unsigned long n = n1 * n2 * size_of_one;
    return n / 1000000;
}
unsigned long memsize_in_MB(const int size_of_one, const int n1, const int n2, const int n3){
    unsigned long n = n1 * n2 * n3 * size_of_one;
    return n / 1000000;
}
unsigned long memsize_in_MB(const int size_of_one, const int n1, const int n2, const int n3, const int n4){
    unsigned long n = n1 * n2 * n3 * n4 * size_of_one;
    return n / 1000000;
}

template <typename T>
T *allocate(T *&arr, const unsigned int n1){
    try{
        arr = new T [n1];
    }
    catch (std::bad_alloc &ba)
    {
        std::cout << " Caught an exception when trying to allocate 1-dimensional array" << std::endl;
        std::cout << " " << ba.what() << " : Array length = " << n1 << std::endl;
        std::cout << " " << ba.what() << " : Array size (MB) = " << memsize_in_MB(sizeof(T), n1) << std::endl;
        exit(EXIT_FAILURE);
    }
    return arr;
}

template <typename T>
T **allocate(T **&arr, const unsigned int n1, const unsigned int n2){
    try{
        arr = new T *[n1];
        arr[0] = new T [n1 * n2];
        for (unsigned int i = 1; i < n1; ++i){
            arr[i] = arr[0] + i * n2;
        }
    }
    catch (std::bad_alloc &ba)
    {
        std::cout << " Caught an exception when trying to allocate 2-dimensional array" << std::endl;
        std::cout << " " << ba.what() << " : Array length = " << n1 << "x" << n2 << std::endl;
        std::cout << " " << ba.what() << " : Array size (MB) = " << memsize_in_MB(sizeof(T), n1, n2) << std::endl;
        exit(EXIT_FAILURE);
    }
    return arr;
}

template <typename T>
T ***allocate(T ***&arr, const unsigned int n1, const unsigned int n2, const unsigned int n3){
    try{
        arr = new T **[n1];
        arr[0] = new T *[n1 * n2];
        arr[0][0] = new T [n1 * n2 * n3];
        for (unsigned int i = 0; i < n1; ++i){
            arr[i] = arr[0] + i * n2;
            for (unsigned int j = 0; j < n2; ++j){
                arr[i][j] = arr[0][0] + i * n2 * n3 + j * n3;
            }
        }
    }
    catch(std::bad_alloc &ba)
    {
        std::cout << " Caught an exception when trying to allocate 3-dimensional array" << std::endl;
        std::cout << " " << ba.what() << " : Array length = " << n1 << "x" << n2 << "x" << n3 << std::endl;
        std::cout << " " << ba.what() << " : Array size (MB) = " << memsize_in_MB(sizeof(T), n1, n2, n3) << std::endl;
        exit(EXIT_FAILURE);
    }
    return arr;
}

template <typename T>
T ****allocate(T ****&arr, const unsigned int n1, const unsigned int n2, const unsigned int n3, const unsigned int n4){
    try{
        arr = new T ***[n1];
        arr[0] = new T **[n1 * n2];
        arr[0][0] = new T *[n1 * n2 * n3];
        arr[0][0][0] = new T [n1 * n2 * n3 * n4];
        
        for (unsigned int i = 0; i < n1; ++i){
            arr[i] = arr[0] + i * n2;
            for (unsigned int j = 0; j < n2; ++j){
                arr[i][j] = arr[0][0] + i * n2 * n3 + j * n3;
                for (unsigned int k = 0; k < n3; ++k){
                    arr[i][j][k] = arr[0][0][0] + i * n2 * n3 * n4 + j * n3 * n4 + k * n4;
                }
            }
        }
    }
    catch(std::bad_alloc &ba)
    {
        std::cout << " Caught an exception when trying to allocate 4-dimensional array" << std::endl;
        std::cout << " " << ba.what() << " : Array length = " << n1 << "x" << n2 << "x" << n3 << "x" << n4 << std::endl;
        std::cout << " " << ba.what() << " : Array size (MB) = " << memsize_in_MB(sizeof(T), n1, n2, n3, n4) << std::endl;
        exit(EXIT_FAILURE);
    }
    return arr;
}

// deallocator

template <typename T>
void deallocate(T *&arr){
    delete [] arr;
}

template <typename T>
void deallocate(T **&arr){
    delete [] arr[0];
    delete [] arr;
}

template <typename T>
void deallocate(T ***&arr){
    delete [] arr[0][0];
    delete [] arr[0];
    delete [] arr;
}

template <typename T>
void deallocate(T ****&arr){
    delete [] arr[0][0][0];
    delete [] arr[0][0];
    delete [] arr[0];
    delete [] arr;
}


