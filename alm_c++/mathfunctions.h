#pragma once

#include <iostream>
#include <cstdlib>

template <typename T>
inline void matmul3(T ret[3][3], T amat[3][3], T bmat[3][3]) {
	int i, j, k;

	for (i = 0; i < 3; ++i) {
		for (j = 0; j < 3; ++j) {
			ret[i][j] = 0.0;
			for (k = 0; k < 3; ++k) ret[i][j] += amat[i][k] * bmat[k][j]; 	        
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

inline int nint(double x)
{
	return int(x + 0.5 - (x < 0.0));
}