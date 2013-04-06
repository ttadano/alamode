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
        void calc_tdos();

        double emin, emax, delta_e;
        int n_energy;
        double *energy_dos;
        double *dos_phonon;
        double **pdos_phonon;
    };
}
