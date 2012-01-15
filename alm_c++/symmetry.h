#ifndef ALM_SYMMETRY_HEADER
#define ALM_SYMMETRY_HEADER

#include "pointers.h"

namespace ALM_NS {
    class Symmetry: protected Pointers {
    public:
        Symmetry(class ALM *);
        ~Symmetry();
        
        void init();
        void gensym(int, int, int, double[3][3], double[3][3], 
            double **, int *, int ***, double **);
        void findsym(int, int, int *, double[3][3], double[3][3],
            double **, int, int ***, int **);
        bool is_ortho(int[3][3], double[3][3], double[3][3]);
        bool is_invariant(int[3][3], int, int*, double **, int[3], int);

        int nsym, nnp;
        int ntran;
        int maxsym;
        
        int ***symrel_int;
        int *symnum_tran;
        double ***symrel;
        double **tnons;


    };
}
#endif