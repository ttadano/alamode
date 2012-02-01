#ifndef ALM_SYSTEM_HEADER
#define ALM_SYSTEM_HEADER

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
        int ndata;
        int *kd;
        double lavec[3][3], rlavec[3][3];
        double *mass_kd, *mass;
        double **xcoord;
        std::string *kdname;

    };
}
#endif