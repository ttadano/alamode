#pragma once

#include "pointers.h"
#include <complex>

namespace PHON_NS {




    class Relaxation: protected Pointers {
    public:
        Relaxation(class PHON *);
        ~Relaxation();

        void setup_relaxation();
        void calc_V3();

    private:
        std::complex<double> ***V3;
    };
}