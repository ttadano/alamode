#pragma once

#include "pointers.h"

namespace PHON_NS
{
    class Dos: protected Pointers{
    public:
        Dos(class PHON *);
        ~Dos();

        void setup();
        bool flag_dos;
        void calc_dos();

        double emin, emax, delta_e;
        int n_energy;
        double *energy_dos;
        double *dos_phonon;

    private:
        void prepare_tetrahedron(const int, const int, const int);
        double dos_integration(const unsigned int, const unsigned int, double **,const double);

        unsigned int ntetra;
        int **tetras;
    };

}