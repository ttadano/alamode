#include "memory.h"

using namespace ALM_NS;

Memory::Memory() {};

double **Memory::create_d2_array(int n1, int n2, double **arr){
	arr = new double*[n1];
	for (int i = 0; i < n1; i++){
		arr[i] = new double[n2];
	}
	return arr;
}

Memory::~Memory() {};