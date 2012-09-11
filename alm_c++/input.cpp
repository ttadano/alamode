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
#include "constraint.h"
#include "fcs.h"

using namespace ALM_NS;

Input::Input(ALM *alm, int narg, char **arg): Pointers(alm) {}

Input::~Input() {}

void Input::parce_input()
{
    int i, j, k;

    using namespace std;
    string job_title;
    string disp_file, force_file, fc2_file, refsys_file;
    string rotation_axis;
    int nat, nkd, nsym, nnp;
    int ndata, nstart, nend, nskip;
    int *kd;
    bool is_periodic[3];
    int multiply_data;
    int constraint_flag;
    double lavec[3][3];
    double ***rcs, **xeq;
    string *kdname;
    double *masskd;
    int maxorder;
    unsigned int nboot;

    int *nbody_include;

    // Read Job prefix
    cin >> job_title;

    // Read nat and nkd
    cin >> nat >> nkd;
    cin >> nsym >> nnp;

    // Read lattice vector in Bohr unit
    /*
    In this program, lattice vectors a1, a2, a3 and 
    matrix lavec are related to each other as follows:

    lavec = (a1, a2, a3)  

      (a11 a12 a13)
    = (a21 a22 a23)
      (a31 a32 a33) , that is

    a1 = t(a11, a21, a31) etc.
    */
    cin >> lavec[0][0] >> lavec[1][0] >> lavec[2][0]; // a1
    cin >> lavec[0][1] >> lavec[1][1] >> lavec[2][1]; // a2
    cin >> lavec[0][2] >> lavec[1][2] >> lavec[2][2]; // a3

    // Read Cutoff Radius for each species
    cin >> maxorder;
    if(maxorder < 1) error->exit("parce_input", "maxorder has to be a positive integer");
    interaction->maxorder = maxorder;
    memory->allocate(rcs, maxorder, nkd, nkd);
    memory->allocate(interaction->rcs, maxorder, nkd, nkd);

    for (i = 0; i < maxorder; ++i){
        for (j = 0; j < nkd; ++j){
            for (k = 0; k < nkd; ++k){
                cin >> rcs[i][j][k];
            }
        }
        for (j = 0; j < nkd; ++j) {
            for (k = j + 1; k < nkd; ++k){
                if (rcs[i][j][k] != rcs[i][k][j]) error->exit("input", "Inconsistent cutoff radius rcs for order =", k + 1);
            }
        }
    }

    memory->allocate(nbody_include, maxorder);
    memory->allocate(fcs->nbody_include, maxorder);

    for (i = 0; i < maxorder; ++i) {
       cin >> nbody_include[i];
    }
    if (nbody_include[0] != 2) {
        error->warn("parce_input", "Harmonic interaction is always 2 body (except on-site 1 body)");
    }

    for (i = 0; i < maxorder; ++i){
       fcs->nbody_include[i] = nbody_include[i];    
    }

    memory->deallocate(nbody_include);

    // Read data info used for fitting
    cin >> ndata >> nstart >> nend >> nskip;
    if(ndata <= 0 || nstart <= 0 || nend <= 0 
        || nstart > ndata || nend > ndata || nstart > nend) {
            error->exit("parce_input", "ndata, nstart, nend are not consistent with each other");
    }
    if(nskip < -1) error->exit("parce_input", "nskip has to be larger than -2.");
    if(nskip == -1) {
        cin >> nboot;
        if(nboot <= 0) error->exit("parce_input", "nboot has to be a positive integer");
        fitting->nboot = nboot;
    }
    cin >> disp_file;
    cin >> force_file;

    cin >> multiply_data;
    if (multiply_data == 2) cin >> refsys_file;
    cin >> constraint_flag;
    if(constraint_flag == 2 || constraint_flag == 4 || constraint_flag == 6) {
        cin >> fc2_file;
        constraint->fc2_file = fc2_file;
    }

    if(constraint_flag >= 3) {
        cin >> rotation_axis;
        constraint->rotation_axis = rotation_axis;
    }

    cin >> is_periodic[0] >> is_periodic[1] >> is_periodic[2];

    // Read species mass
    memory->allocate(kdname, nkd);
    memory->allocate(masskd, nkd);
    for (i = 0; i < nkd; ++i){
        cin >> kdname[i] >> masskd[i];
    }

    // Read atomic coordinates
    memory->allocate(kd, nat);
    memory->allocate(xeq, nat, 3);
    for (i = 0; i < nat; ++i){
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
    system->nstart = nstart;
    system->nend = nend;
    system->nskip = nskip;

    for (i = 0; i < maxorder; ++i){
        for (j = 0; j < nkd; ++j) {
            for (k = 0; k < nkd; ++k){
                interaction->rcs[i][j][k] = rcs[i][j][k];
            }
        }
    }
    symmetry->multiply_data = multiply_data;
    if(multiply_data == 2) symmetry->refsys_file = refsys_file;
    constraint->constraint_mode = constraint_flag;

    for (i = 0; i < 3; ++i) interaction->is_periodic[i] = is_periodic[i];

    for (i = 0; i < 3; ++i){
        for (j = 0; j < 3; ++j){
            system->lavec[i][j] = lavec[i][j];
        }
    }
    system->kdname = new string[nkd];
    system->mass_kd = new double[nkd];
    for (int i = 0; i < nkd; ++i){
        system->kdname[i] = kdname[i];
        system->mass_kd[i] = masskd[i];
    }

    system->kd = new int[nat];
    memory->allocate(system->xcoord, nat, 3);
    for (i = 0; i < nat; ++i){
        system->kd[i] = kd[i];
        for (j = 0; j < 3; ++j){
            system->xcoord[i][j] = xeq[i][j];
        }
    }
    memory->deallocate(masskd);
    memory->deallocate(kd);
    memory->deallocate(xeq);
}
