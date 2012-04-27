#pragma once

#include "pointers.h"
#include <string>

namespace ALM_NS {
    class System: protected Pointers {
    public:
        System(class ALM *);
        ~System();
        void init();
        void recips(double [3][3], double [3][3]);
        void frac2cart(double **);
     
        
        int nat, nkd;
        int ndata, nstart, nend, nskip;
        int *kd;
        double lavec[3][3], rlavec[3][3];
        double *mass_kd, *mass;
        double **xcoord; // fractional coordinate
        double **x_cartesian;
        std::string *kdname;

    };
}