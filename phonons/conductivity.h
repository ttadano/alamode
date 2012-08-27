#pragma once

#include "pointers.h"

namespace PHON_NS {

    class Conductivity: protected Pointers {
    public:
        Conductivity(class PHON *);
        ~Conductivity();
        void setup_kl();
        void calc_kl();

        double **tau;
        double kl[3][3];
        double Tmin, Tmax, dT;

    private:
        double ***vel;
        void gen_tau(const double);
        void calc_kl_at_T(const double);
        unsigned int nk, ns;
        double *func;
    };
}