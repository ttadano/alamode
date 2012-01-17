#ifndef ALM_INTERACTION_HEADER
#define ALM_INTERACTION_HEADER

#include "pointers.h"
#include <Eigen/core>

namespace ALM_NS {
    class Interaction: protected Pointers {
    public:
        Interaction(class ALM *);
        ~Interaction();
        void init();

        double **rcs;
        double **distlist;

        double distance(double *, double *);

        bool is_periodic[3];

        int nneib;
        double ***xcrd;

    private:
        int nsize[3];
        void calc_distlist(int, double **);
    };
}

#endif