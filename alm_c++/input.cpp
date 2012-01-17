#include "input.h"
#include <iostream>
#include <string>
#include "memory.h"
#include "files.h"
#include "interaction.h"
#include "system.h"
#include "symmetry.h"

using namespace ALM_NS;

Input::Input(ALM *alm, int narg, char **arg): Pointers(alm) {}

Input::~Input() {}

void Input::sparce_input()
{
    using namespace std;
    string job_title;
    string disp_file, force_file;
	int nat, nkd, nsym, nnp, ndata;
	int *kd;
    bool is_periodic[3];
	bool multiply_data, constraint;
	double eps;
	double lavec[3][3];
	double **rcs, **xeq;
	string *kdname;
	double *masskd, *mass;
	// Read Job prefix
	cin >> job_title;
	// Read nat and nkd
	cin >> nat >> nkd;
	cin >> nsym >> nnp;
	// Read lattice vector in Bohr unit
	cin >> lavec[0][0] >> lavec[0][1] >> lavec[0][2];
	cin >> lavec[1][0] >> lavec[1][1] >> lavec[1][2];
	cin >> lavec[2][0] >> lavec[2][1] >> lavec[2][2];
	// Read Cutoff Radius for each species
    memory->allocate(rcs,nkd,3);
    memory->allocate(interaction->rcs, nkd, 3);
	for (int i = 0; i < nkd; i++){
		cin >> rcs[i][0] >> rcs[i][1] >> rcs[i][2];
	}
	cin >> ndata;
	cin >> eps;
	cin >> disp_file;
	cin >> force_file;
	cin >> multiply_data >> constraint;
    cin >> is_periodic[0] >> is_periodic[1] >> is_periodic[2];
	// Read species mass
	kdname = new string[nkd];
	masskd = new double[nkd];
	for (int i = 0; i < nkd; i++){
		cin >> kdname[i] >> masskd[i];
	}
	// Read atomic coordinates
	kd = new int[nat];
    memory->allocate(xeq, nat, 3);
	for (int i = 0; i < nat; i++){
		cin >> kd[i] >> xeq[i][0] >> xeq[i][1] >> xeq[i][2];
	}

    files->job_title = job_title;
    system->nat = nat;
    system->nkd = nkd;
    symmetry->nsym = nsym;
    symmetry->nnp = nnp;
    files->file_disp = disp_file;
    files->file_force = force_file;
    system->ndata = ndata;
    for (int i = 0; i < nkd; i++){
        for (int j = 0; j < 3; j++){
        interaction->rcs[i][j] = rcs[i][j];
        }
    }
    symmetry->multiply_data = multiply_data;

    for (int i = 0; i < 3; i++) interaction->is_periodic[i] = is_periodic[i];

    for (int i = 0; i < 3; i++){
        for (int j = 0; j < 3; j++){
        system->lavec[i][j] = lavec[i][j];
        }
    }
    system->kdname = new string[nkd];
    system->mass_kd = new double[nkd];
    for (int i = 0; i < nkd; i++){
        system->kdname[i] = kdname[i];
        system->mass_kd[i] = masskd[i];
    }

    system->kd = new int[nat];
    memory->allocate(system->xcoord, nat, 3);
    for (int i = 0; i < nat; i++){
        system->kd[i] = kd[i];
        for (int j = 0; j < 3; j++){
            system->xcoord[i][j] = xeq[i][j];
        }
    }
     delete masskd;
     delete kd;
     memory->deallocate(xeq);
}