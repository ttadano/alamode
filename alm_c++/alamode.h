/* Declaration of pointers used in the whole program. */

#ifndef ALM_ALAMODE_HEADER
#define ALM_ALAMODE_HEADER
#include "mpi.h"

namespace ALM_NS {

    class ALM {
    public:
        class Memory *memory;
        ALM(int, char **, MPI_Comm);
        ~ALM();

    };
}

#endif