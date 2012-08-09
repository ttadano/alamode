#pragma once

#include "pointers.h"

namespace PHON_NS {

    class Integration: protected Pointers {
    public:
        Integration(class PHON *);
        ~Integration();

        void setup_integration();
        void finish_integration();
        double do_tetrahedron(double *, double *, const double);

    private:
        unsigned int ntetra;
        int **tetras;
        void prepare_tetrahedron(const int, const int, const int);
        inline double fij(const double, const double, const double);
    };

}