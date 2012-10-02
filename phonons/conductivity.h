#pragma once

#include "pointers.h"

namespace PHON_NS {

    class Conductivity: protected Pointers {
    public:
        Conductivity(class PHON *);
        ~Conductivity();
        void setup_kl();
        void finish_kl();
        void calc_kl();

        double **tau;
        double kl[3][3];

    private:
        double ***vel;
        void calc_kl_at_T(const double);
        void calc_kl_mpi(const unsigned int, unsigned int *, double *, unsigned int *, unsigned int **, const unsigned int, double *, double ***); 
        unsigned int nk, ns;
    };
}
