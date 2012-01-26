#ifndef ALM_FITTING_HEADER
#define ALM_FITTING_HEADER

#include "pointers.h"

namespace ALM_NS {
    class Fitting: protected Pointers {
    public:
        Fitting(class ALM *);
        ~Fitting();

        void fitmain();

        double **amat;
        double *fsum;
    private:
        int inprim_index(const int);
    };

}
#endif