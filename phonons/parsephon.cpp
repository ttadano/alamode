#include "parsephon.h"
#include "error.h"
#include "system.h"
#include "kpoint.h"
#include "fcs_phonon.h"
#include "dynamical.h"
#include "write_phonons.h"
#include <iostream>
#include <string>

using namespace PHON_NS;

Input::Input(PHON *phon, int narg, char **arg): Pointers(phon) {}

Input::~Input() {}

void Input::parce_input()
{
    using namespace std;
    string mode;

    cin >> job_title;
    cin >> mode;
    if(mode == "phonons") {
        cout << "Calculation of PHONONS" << endl << endl;
        read_input_phonons();
    } else if (mode == "boltzmann"){
        cout << "Calculation of Thermal Conductivity" << endl << endl;
        error->exit("parse_input", "Sorry :( Boltzmann is not supported yet");
    } else {
        error->exit("parse_input", "invalid mode");
    }
}

void Input::read_input_phonons()
{
    using namespace std;

    string file_fcs;
    double lavec[3][3];
    int kpoint_mode;
    unsigned int nkmax;
    bool eigenvectors, writeanime, nonanalytic;

    string file_born;
    double na_sigma;

    unsigned int i, j;

    for (i = 0; i < 3; ++i){
        cin >> lavec[i][0] >> lavec[i][1] >> lavec[i][2];
    }
    cin >> file_fcs;
    cin >> kpoint_mode >> nkmax;

    cin >> eigenvectors >> writeanime >> nonanalytic;
    if(nonanalytic) cin >> file_born >> na_sigma;

    // distribute input parameters to each class

    for (i = 0; i < 3; ++i){
        for(j = 0; j < 3; ++j){
            system->lavec_p[i][j] = lavec[i][j];
        }
    }

    dynamical->eigenvectors = eigenvectors;
    writes->writeanime = writeanime;
    dynamical->nonanalytic = nonanalytic;

    if(nonanalytic) {
        dynamical->file_born = file_born;
        dynamical->na_sigma = na_sigma;
    }

    fcs_phonon->file_fcs = file_fcs;
    kpoint->kpoint_mode = kpoint_mode;
    kpoint->nkmax = nkmax;

}