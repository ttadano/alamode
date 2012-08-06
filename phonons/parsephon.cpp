#include "parsephon.h"
#include "error.h"
#include "system.h"
#include "kpoint.h"
#include "fcs_phonon.h"
#include "dynamical.h"
#include "write_phonons.h"
#include <iostream>
#include <string>
#include <algorithm>
#include <boost/algorithm/string.hpp>

using namespace PHON_NS;

Input::Input(PHON *phon, int narg, char **arg): Pointers(phon) {}

Input::~Input() {}

void Input::parce_input()
{
    using namespace std;
    string mode;

    cin >> job_title;
    cin >> mode;
    boost::to_lower(mode);
//    transform(mode.begin(), mode.end(), mode.begin(), tolower);

    phon->mode = mode;

    if(mode == "phonons") {
        cout << "Calculation of PHONONS" << endl << endl;
        read_input_phonons();
    } else if (mode == "boltzmann"){
        cout << "Calculation of Thermal Conductivity" << endl << endl;
        error->warn("parse_input", "Sorry :( Boltzmann is still under implementation");
        read_input_boltzmann();
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
    bool eigenvectors, writeanime, nonanalytic;
    int nbands;

    string file_born;
    double na_sigma;

    unsigned int i, j;

    for (i = 0; i < 3; ++i){
        cin >> lavec[0][i] >> lavec[1][i] >> lavec[2][i];
    }
    cin >> file_fcs;

    cin >> eigenvectors >> writeanime >> nonanalytic;
    cin >> nbands;
    cin >> kpoint_mode;
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
    writes->nbands = nbands;

    if(nonanalytic) {
        dynamical->file_born = file_born;
        dynamical->na_sigma = na_sigma;
    }

    fcs_phonon->file_fcs = file_fcs;
    kpoint->kpoint_mode = kpoint_mode;

}

void Input::read_input_boltzmann()
{
    using namespace std;

    string file_fcs;
    double lavec[3][3];

    unsigned int i, j;

    for (i = 0; i < 3; ++i) {
        cin >> lavec[0][i] >> lavec[1][i] >> lavec[2][i];
    }
    cin >> file_fcs;

    for (i = 0; i < 3; ++i){
        for (j = 0; j < 3; ++j){
        system->lavec_p[i][j] = lavec[i][j];
        }
    }
    fcs_phonon->file_fcs = file_fcs;
}
