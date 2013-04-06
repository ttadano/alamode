#pragma once

#include "pointers.h"
#include <vector>

namespace PHON_NS {

    struct tetra_pair {
        double e;
        double f;
    };

    inline bool operator<(const tetra_pair &a, const tetra_pair &b){
        return a.e < b.e;
    }

    class Integration: protected Pointers {
    public:
        Integration(class PHON *);
        ~Integration();

        void setup_integration();
        void finish_integration();
        double do_tetrahedron(double *, double *, const double);
        double dos_integration(double *, const double);

    private:
        unsigned int ntetra;
        int **tetras;
        void prepare_tetrahedron(const int, const int, const int);
        inline double fij(const double, const double, const double);
        inline double volume(int *);
        std::vector<tetra_pair> tetra_data;
        inline double refold(double);
    };
}
