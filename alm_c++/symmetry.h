#ifndef ALM_SYMMETRY_HEADER
#define ALM_SYMMETRY_HEADER

#include "pointers.h"

namespace ALM_NS {
    class Symmetry: protected Pointers {
    public:
        Symmetry(class ALM *);
        ~Symmetry();

        int nsym, nnp;
        
        int ***symrel_int;
        double ***symrel;
        double **tnons;


    };
}
#endif