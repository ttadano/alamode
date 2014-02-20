#pragma once

#include "pointers.h"

namespace PHON_NS {
    class Phonon_velocity: protected Pointers {
    public:
        Phonon_velocity(class PHON *);
        ~Phonon_velocity();

        void calc_group_velocity(const int);
        void calc_phonon_vel_band();
        void phonon_vel_k(double *, double **);

        bool printvel;
        double **phvel;

    private:
        double diff(double *, const unsigned int, double);
    };
}
