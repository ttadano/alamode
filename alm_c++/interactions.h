#ifndef ALM_INTERACTION_HEADER
#define ALM_INTERACTION_HEADER

#include "pointers.h"

namespace ALM_NS {
    class Interaction: protected Pointers {
    public:
        Interaction(class ALM *);
        ~Interaction();
    
        double **rcs;
        double **distance;
    };
}

#endif