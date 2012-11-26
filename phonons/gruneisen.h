#pragma once

#include "pointers.h"
#include <string>
#include <complex>

namespace PHON_NS {
    class Gruneisen: protected Pointers {
    public:
        Gruneisen(class PHON *);
        ~Gruneisen();
        
        double delta_a;
        void setup();
        std::complex<double> **gruneisen;
        void calc_gruneisen();
        void calc_gruneisen2();
        void finish_gruneisen();

    private:
         double ****fc2_plus, ****fc2_minus;
         double ****dfc2;
         void prepare_delta_fc2();
         void prepare_newfc2();
         void calc_dfc2_reciprocal(std::complex<double> **, double *);
        
    };
}