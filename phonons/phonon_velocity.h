#pragma once

#include "pointers.h"

namespace PHON_NS {
    class Phonon_velocity: protected Pointers {
    public:
        Phonon_velocity(class PHON *);
        ~Phonon_velocity();

       void calc_phonon_vel_band();
       void phonon_vel_k(double *, double **);

       double **phvel;
   
    private:
        double diff(double *, const unsigned int, double);
        double freq(const double);
    };
}
