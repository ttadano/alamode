#include "mpi_common.h"
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
#include "relaxation.h"
#include "conductivity.h"
#include "symmetry_core.h"

using namespace PHON_NS;

Input::Input(PHON *phon, int narg, char **arg): Pointers(phon) {}

Input::~Input() {}

void Input::parce_input()
{
    using namespace std;
    string mode;

    if (mympi->my_rank == 0) {
        cin >> job_title;
        cin >> mode;
        boost::to_lower(mode);
    }
/*
    int len = job_title.length();
    MPI_Bcast(&len, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(const_cast<char *>(job_title.data()), len, MPI_CHAR, 0, MPI_COMM_WORLD);
    std::cout << job_title.data() << std::endl;
*/

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
    double Tmin, Tmax, dT;

    unsigned int i, j;
    unsigned int nsym, nnp;

    cin >> nsym >> nnp;
    for (i = 0; i < 3; ++i){
        cin >> lavec[0][i] >> lavec[1][i] >> lavec[2][i];
    }
    cin >> file_fcs;

    cin >> eigenvectors >> writeanime >> nonanalytic;
    cin >> nbands;
    cin >> Tmin >> Tmax >> dT;
    cin >> kpoint_mode;
    if(nonanalytic) cin >> file_born >> na_sigma;

    // distribute input parameters to each class

    symmetry->nsym = nsym;
    symmetry->nnp = nnp;

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

    if (Tmin > Tmax) error->exit("read_input_phonon", "Tmin is larger than Tmax");

    system->Tmin = Tmin;
    system->Tmax = Tmax;
    system->dT = dT;

}

void Input::read_input_boltzmann()
{
    using namespace std;

    string file_fcs;
    double lavec[3][3];
    double epsilon;
    double Tmin, Tmax, dT;

    unsigned int i, j;
    unsigned int nsym, nnp;

    cin >> nsym >> nnp;
    for (i = 0; i < 3; ++i) {
        cin >> lavec[0][i] >> lavec[1][i] >> lavec[2][i];
    }
    cin >> file_fcs;
    cin >> epsilon;
    cin >> Tmin >> Tmax >> dT;

    if (Tmin > Tmax) {
        error->exit("read_input_boltzmann", "Tmin is bigger than Tmax");
    }

    symmetry->nsym = nsym;
    symmetry->nnp = nnp;

    for (i = 0; i < 3; ++i){
        for (j = 0; j < 3; ++j){
            system->lavec_p[i][j] = lavec[i][j];
        }
    }
    fcs_phonon->file_fcs = file_fcs;
    relaxation->epsilon = epsilon;

    system->Tmin = Tmin;
    system->Tmax = Tmax;
    system->dT = dT;
}
