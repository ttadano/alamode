#pragma once

#include "pointers.h"
#include <string>

namespace ALM_NS {
    class Ewald: protected Pointers{
    public:
        Ewald(class ALM *);
        ~Ewald();

        bool is_longrange;
        std::string file_longrange;

        void init();
        double *Q;
        double alpha;

        void ewald_force(const int, double **, double **);
    
    private:
        void load_charge();
        void ewald_force_short(const int, double **, double **);
        void ewald_force_long(const int, double **, double **);
    };
}