// main.cpp 

#include "alamode.h"
#include "input.h"
#include "files.h"
#include "memory.h"
#include <iostream>
#include <string>
#include <vector>

using namespace ALM_NS;
using namespace std;

int main(int argc, char **argv)
{
	string str2, disp_file, force_file;
	int nat, nkd, nsym, nnp, ndata;
	int *kd;
	bool multiply_data, constraint;
	double eps;
	double lavec[3][3];
	double **rcs, **xeq;
	string *kdname;
	double *masskd, *mass;
	// Read Job prefix
	cin >> job_prefix;
	// Read nat and nkd
	cin >> nat >> nkd;
	cin >> nsym >> nnp;
	// Read lattice vector in Bohr unit
	cin >> lavec[0][0] >> lavec[0][1] >> lavec[0][2];
	cin >> lavec[1][0] >> lavec[1][1] >> lavec[1][2];
	cin >> lavec[2][0] >> lavec[2][1] >> lavec[2][2];
	// Read Cutoff Radius for each species
	//	vector<vector<double>> rcs(nkd, vector<int>(3));
//	rcs = new double*[nkd];
//	for (int i = 0; i < nkd; i++){
//		rcs[i] = new double[3];
//	}

	Memory memory;
	//rcs = mem1.create_d2_array(nkd, 3, rcs);
	rcs = memory.create(nkd, 3);
	for (int i = 0; i < nkd; i++){
		cin >> rcs[i][0] >> rcs[i][1] >> rcs[i][2];
	}
	cin >> ndata;
	cin >> eps;
	cin >> disp_file;
	cin >> force_file;
	cin >> multiply_data >> constraint;
	// Read species mass
	kdname = new string[nkd];
	masskd = new double[nkd];
	for (int i = 0; i < nkd; i++){
		cin >> kdname[i] >> masskd[i];
	}
	// Read atomic coordinates
	kd = new int[nat];
	xeq = new double*[3];
	for (int i = 0; i < 3; i++){
		xeq[i] = new double[nat];
	}
	for (int i = 0; i < nat; i++){
		cin >> kd[i] >> xeq[0][i] >> xeq[1][i] >> xeq[2][i];
	}
	cin >> str2;
	cout << str2 << endl;
	cout << rcs[0][0] << endl;
	cin.get();
	cin.get();
	cin.get();
	return 0;
}