/*
 mathfunctions.h

 Copyright (c) 2014 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#pragma once

#include <iostream>
#include <cstdlib>

template <typename T>
inline void matmul3(T ret[3][3], const T amat[3][3], const T bmat[3][3]) {
	int i, j, k;

	T ret_tmp[3][3];

	for (i = 0; i < 3; ++i) {
		for (j = 0; j < 3; ++j) {
			ret_tmp[i][j] = 0.0;
			for (k = 0; k < 3; ++k) ret_tmp[i][j] += amat[i][k] * bmat[k][j]; 	        
		}
	}

	for (i = 0; i < 3; ++i) {
		for (j = 0; j < 3; ++j) {
			ret[i][j] = ret_tmp[i][j];
		}
	}
}

inline void transpose3(double ret[3][3], const double mat[3][3]) 
{
	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 3; ++j) {
			ret[i][j] = mat[j][i];
		}
	}
}

inline void rotvec(double vec_out[3], double vec_in[3], double mat[3][3], char mode = 'N')
{
	// Perform matrix x vector multiplication. 
	//
	// vec_out = mat      * vec_in   (mode = 'N')
	//          (mat)^{t} * vec_in   (mode = 'T')
	//

	unsigned int i;
	double vec_tmp[3];

	for (i = 0; i < 3; ++i){
		vec_tmp[i] = vec_in[i];
	}

	if (mode == 'N') {
		for (i = 0; i < 3; ++i){
			vec_out[i] = mat[i][0] * vec_tmp[0] + mat[i][1] * vec_tmp[1] + mat[i][2] * vec_tmp[2];
		}
	} else if (mode == 'T'){
		for (i = 0; i < 3; ++i){
			vec_out[i] = mat[0][i] * vec_tmp[0] + mat[1][i] * vec_tmp[1] + mat[2][i] * vec_tmp[2];
		}
	} else {
		std::cout << "Invalid mode " << mode << std::endl;
		exit(1);
	}
}

inline void rotvec(double vec_out[3], double vec_in[3], double **mat, char mode = 'N')
{
	// Perform matrix x vector multiplication. 
	//
	// vec_out = mat      * vec_in   (mode = 'N')
	//          (mat)^{t} * vec_in   (mode = 'T')
	//

	unsigned int i;
	double vec_tmp[3];

	for (i = 0; i < 3; ++i){
		vec_tmp[i] = vec_in[i];
	}

	if (mode == 'N') {
		for (i = 0; i < 3; ++i){
			vec_out[i] = mat[i][0] * vec_tmp[0] + mat[i][1] * vec_tmp[1] + mat[i][2] * vec_tmp[2];
		}
	} else if (mode == 'T'){
		for (i = 0; i < 3; ++i){
			vec_out[i] = mat[0][i] * vec_tmp[0] + mat[1][i] * vec_tmp[1] + mat[2][i] * vec_tmp[2];
		}
	} else {
		std::cout << "Invalid mode " << mode << std::endl;
		exit(1);
	}
}

inline void invmat3(double invmat[3][3], double mat[3][3])
{
	unsigned int i, j;
	double det;
	double mat_tmp[3][3];

	for (i = 0; i < 3; ++i) {
		for (j = 0; j < 3; ++j) {
			mat_tmp[i][j] = mat[i][j];
		}
	}

	det = mat_tmp[0][0] * mat_tmp[1][1] * mat_tmp[2][2] 
	+ mat_tmp[1][0] * mat_tmp[2][1] * mat_tmp[0][2] 
	+ mat_tmp[2][0] * mat_tmp[0][1] * mat_tmp[1][2]
	- mat_tmp[0][0] * mat_tmp[2][1] * mat_tmp[1][2] 
	- mat_tmp[2][0] * mat_tmp[1][1] * mat_tmp[0][2]
	- mat_tmp[1][0] * mat_tmp[0][1] * mat_tmp[2][2];

	if(std::abs(det) < 1.0e-12) {
		std::cout << "invmat3: Given matrix is singular" << std::endl;
		exit(1);
	}

	double factor = 1.0 / det;

	invmat[0][0] = (mat_tmp[1][1] * mat_tmp[2][2] - mat_tmp[1][2] * mat_tmp[2][1]) * factor;
	invmat[0][1] = (mat_tmp[0][2] * mat_tmp[2][1] - mat_tmp[0][1] * mat_tmp[2][2]) * factor;
	invmat[0][2] = (mat_tmp[0][1] * mat_tmp[1][2] - mat_tmp[0][2] * mat_tmp[1][1]) * factor;

	invmat[1][0] = (mat_tmp[1][2] * mat_tmp[2][0] - mat_tmp[1][0] * mat_tmp[2][2]) * factor;
	invmat[1][1] = (mat_tmp[0][0] * mat_tmp[2][2] - mat_tmp[0][2] * mat_tmp[2][0]) * factor;
	invmat[1][2] = (mat_tmp[0][2] * mat_tmp[1][0] - mat_tmp[0][0] * mat_tmp[1][2]) * factor;

	invmat[2][0] = (mat_tmp[1][0] * mat_tmp[2][1] - mat_tmp[1][1] * mat_tmp[2][0]) * factor;
	invmat[2][1] = (mat_tmp[0][1] * mat_tmp[2][0] - mat_tmp[0][0] * mat_tmp[2][1]) * factor;
	invmat[2][2] = (mat_tmp[0][0] * mat_tmp[1][1] - mat_tmp[0][1] * mat_tmp[1][0]) * factor;
}

inline void invmat3_i(int invmat[3][3], int mat[3][3])
{
	int det;

	det = mat[0][0] * mat[1][1] * mat[2][2] 
	+ mat[1][0] * mat[2][1] * mat[0][2] 
	+ mat[2][0] * mat[0][1] * mat[1][2]
	- mat[0][0] * mat[2][1] * mat[1][2] 
	- mat[2][0] * mat[1][1] * mat[0][2]
	- mat[1][0] * mat[0][1] * mat[2][2];

	if(std::abs(det) == 0) {
		std::cout << "invmat3_i: Given matrix is singular" << std::endl;
		exit(1);
	}

	invmat[0][0] = (mat[1][1] * mat[2][2] - mat[1][2] * mat[2][1]) / det;
	invmat[0][1] = (mat[0][2] * mat[2][1] - mat[0][1] * mat[2][2]) / det;
	invmat[0][2] = (mat[0][1] * mat[1][2] - mat[0][2] * mat[1][1]) / det;

	invmat[1][0] = (mat[1][2] * mat[2][0] - mat[1][0] * mat[2][2]) / det;
	invmat[1][1] = (mat[0][0] * mat[2][2] - mat[0][2] * mat[2][0]) / det;
	invmat[1][2] = (mat[0][2] * mat[1][0] - mat[0][0] * mat[1][2]) / det;

	invmat[2][0] = (mat[1][0] * mat[2][1] - mat[1][1] * mat[2][0]) / det;
	invmat[2][1] = (mat[0][1] * mat[2][0] - mat[0][0] * mat[2][1]) / det;
	invmat[2][2] = (mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0]) / det;

}

inline int nint(double x)
{
	return int(x + 0.5 - (x < 0.0));
}