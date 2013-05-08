#pragma once

#include <string>
#include <cctype>
#include <algorithm>
#include <functional>
#include <fstream>
#include "../alm_c++/constants.h"


int locate_tag(std::string);
std::ifstream ifs;

int i, j, k, l;

int nat, nkd, ns;
int nkx, nky, nkz;
int nk;
int nt;
double *temp, tmin, tmax, dt;
double volume;

double **omega, ***tau;
double ****vel;
int *n_weight;

std::string calc;

int beg_k, end_k;
int beg_s, end_s;

void calc_tau(int);
void calc_tau_temp(int, int);
void calc_kappa();
void calc_kappa_size(double, double, int, int [3]);
double Cv(double, double);

static const double Ryd_to_kayser = Hz_to_kayser / time_ry;
static const double kayser_to_Ryd = 1.0 / Ryd_to_kayser;
static const double T_to_Ryd = k_Boltzmann / Ryd;


// trim from start
static inline std::string &ltrim(std::string &s) {
	s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
	return s;
}

// trim from end
static inline std::string &rtrim(std::string &s) {
	s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
	return s;
}

// trim from both ends
static inline std::string &trim(std::string &s) {
	return ltrim(rtrim(s));
}



// memsize calculator

int memsize_in_MB(int size_of_one, int n1){
	unsigned long n = n1 * size_of_one;
	return n / 1000000;
}
int memsize_in_MB(int size_of_one, int n1, int n2){
	unsigned long n = n1 * n2 * size_of_one;
	return n / 1000000;
}
int memsize_in_MB(int size_of_one, int n1, int n2, int n3){
	unsigned long n = n1 * n2 * n3 * size_of_one;
	return n / 1000000;
}
int memsize_in_MB(int size_of_one, int n1, int n2, int n3, int n4){
	unsigned long n = n1 * n2 * n3 * n4 * size_of_one;
	return n / 1000000;
}


template <typename T>
T *allocate(T *&arr, int n1){

	try{
		arr = new T [n1];
	}
	catch (std::bad_alloc &ba)
	{
		std::cout << "Caught an exception when trying to allocate 1-dimensional array" << std::endl;
		std::cout << ba.what() << " : Array size (MB) = " << memsize_in_MB(sizeof(T), n1) << std::endl;
		exit(EXIT_FAILURE);
	}
	return arr;
}

template <typename T>
T **allocate(T **&arr, int n1, int n2){
	try{
		arr = new T *[n1];
		arr[0] = new T [n1 * n2];
		for (int i = 1; i < n1; ++i){
			arr[i] = arr[0] + i * n2;
		}
	}
	catch (std::bad_alloc &ba)
	{
		std::cout << "Caught an exception when trying to allocate 2-dimensional array" << std::endl;
		std::cout << ba.what() << " : Array size (MB) = " << memsize_in_MB(sizeof(T), n1, n2) << std::endl;
		exit(EXIT_FAILURE);
	}
	return arr;
}

template <typename T>
T ***allocate(T ***&arr, int n1, int n2, int n3){
	try{
		arr = new T **[n1];
		arr[0] = new T *[n1 * n2];
		arr[0][0] = new T [n1 * n2 * n3];
		for (int i = 0; i < n1; ++i){
			arr[i] = arr[0] + i * n2;
			for (int j = 0; j < n2; ++j){
				arr[i][j] = arr[0][0] + i * n2 * n3 + j * n3;
			}
		}
	}
	catch(std::bad_alloc &ba)
	{
		std::cout << "Caught an exception when trying to allocate 3-dimensional array" << std::endl;
		std::cout << ba.what() << " : Array size (MB) = " << memsize_in_MB(sizeof(T), n1, n2, n3) << std::endl;
		exit(EXIT_FAILURE);
	}
	return arr;
}

template <typename T>
T ****allocate(T ****&arr, int n1, int n2, int n3, int n4){

	try{
		arr = new T ***[n1];
		arr[0] = new T **[n1 * n2];
		arr[0][0] = new T *[n1 * n2 * n3];
		arr[0][0][0] = new T [n1 * n2 * n3 * n4];

		for (int i = 0; i < n1; ++i){
			arr[i] = arr[0] + i * n2;
			for (int j = 0; j < n2; ++j){
				arr[i][j] = arr[0][0] + i * n2 * n3 + j * n3;
				for (int k = 0; k < n3; ++k){
					arr[i][j][k] = arr[0][0][0] + i * n2 * n3 * n4 + j * n3 * n4 + k * n4;
				}
			}
		}
	}
	catch(std::bad_alloc &ba)
	{
		std::cout << "Caught an exception when trying to allocate 3-dimensional array" << std::endl;
		std::cout << ba.what() << " : Array size (MB) = " << memsize_in_MB(sizeof(T), n1, n2, n3, n4) << std::endl;
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
