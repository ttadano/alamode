#include "write_phonons.h"
#include "system.h"
#include "dynamical.h"
#include "kpoint.h"
#include "parsephon.h"
#include "error.h"
#include "phonon_dos.h"
#include "../alm_c++/constants.h"
#include <iomanip>
#include <fstream>

using namespace PHON_NS;

Writes::Writes(PHON *phon): Pointers(phon){
    Ry_to_kayser = std::pow(Hz_to_kayser, 2) / (amu_ry * std::pow(time_ry, 2));
};

Writes::~Writes(){};

void Writes::write_phonon_info()
{
    if(kpoint->kpoint_mode == 1){
        write_phonon_bands();
    }

    if(dos->flag_dos) {
        write_phonon_dos();
    }

    if(writeanime) {
        write_mode_anime();
    }
}

void Writes::write_phonon_bands()
{
    std::ofstream ofs_bands;

    file_bands = input->job_title + ".bands";
    ofs_bands.open(file_bands.c_str(), std::ios::out);
    if(!ofs_bands) error->exit("write_phonon_bands", "cannot open file_bands");

    unsigned int i, j;

    unsigned int nk = kpoint->nk;

    double *kaxis = kpoint->kaxis;
    double **eval = dynamical->eval_phonon;

    unsigned int neval = 3 * system->natmin;

    ofs_bands << "# k-axis, Eigenvalues [cm^-1]" << std::endl;
    ofs_bands.setf(std::ios::fixed);

    for (i = 0; i < nk; ++i){
        ofs_bands << std::setw(8) << kaxis[i];
        for (j = 0; j < neval; ++j){
            ofs_bands << std::setw(12) << in_kayser(eval[i][j]);
        }
        ofs_bands << std::endl;
    }

    ofs_bands.close();
}

void Writes::write_phonon_dos()
{
    std::ofstream ofs_dos;

    file_bands = input->job_title + ".dos";
    ofs_dos.open(file_bands.c_str(), std::ios::out);
    if(!ofs_dos) error->exit("write_phonon_dos", "cannot open file_dos");
    ofs_dos.close();

}

void Writes::write_mode_anime()
{
    std::ofstream ofs_anime;

    file_anime = input->job_title + ".axsf";
    ofs_anime.open(file_anime.c_str(), std::ios::out);
    if(!ofs_anime) error->exit("write_mode_anime", "cannot open file_anime");
    ofs_anime.close();

}

double Writes::in_kayser(const double x)
{
    double val = x;
    val *= Ry_to_kayser;

    if(val < 0.0){
        return -std::sqrt(-val);
    } else {
        return std::sqrt(val);
    }
}