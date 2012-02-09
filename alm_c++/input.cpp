#include "input.h"
#include <iostream>
#include <string>
#include "memory.h"
#include "files.h"
#include "interaction.h"
#include "system.h"
#include "symmetry.h"
#include "error.h"
#include "fitting.h"

using namespace ALM_NS;

Input::Input(ALM *alm, int narg, char **arg): Pointers(alm) {}

Input::~Input() {}

void Input::sparce_input()
{
    using namespace std;
    string job_title;
    string disp_file, force_file, fc2_file;
	int nat, nkd, nsym, nnp, ndata;
	int *kd;
    bool is_periodic[3];
	bool multiply_data;
    int constraint_flag;
	double lavec[3][3];
	double **rcs, **xeq;
	string *kdname;
	double *masskd;
    int maxorder;
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
    cin >> maxorder;
    if(maxorder < 1) error->exit("sparce_input", "maxorder has to be a positive integer");
    interaction->maxorder = maxorder;
    memory->allocate(rcs,nkd,maxorder);
    memory->allocate(interaction->rcs, nkd, maxorder);
	for (int i = 0; i < nkd; ++i){
        for (int j = 0; j < maxorder; ++j){
        cin >> rcs[i][j];
        }
	}
	cin >> ndata;
	cin >> disp_file;
	cin >> force_file;
	cin >> multiply_data >> constraint_flag;

    if(constraint_flag == 2) {
        cin >> fc2_file;
        fitting->fc2_file = fc2_file;
    }

    cin >> is_periodic[0] >> is_periodic[1] >> is_periodic[2];
	// Read species mass
    memory->allocate(kdname, nkd);
    memory->allocate(masskd, nkd);
	for (int i = 0; i < nkd; i++){
		cin >> kdname[i] >> masskd[i];
	}
	// Read atomic coordinates
    memory->allocate(kd, nat);
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
        for (int j = 0; j < maxorder; j++){
        interaction->rcs[i][j] = rcs[i][j];
        }
    }
    symmetry->multiply_data = multiply_data;
    fitting->constraint = constraint_flag;

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
     memory->deallocate(masskd);
     memory->deallocate(kd);
     memory->deallocate(xeq);
}