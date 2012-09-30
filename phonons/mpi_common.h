#pragma once

#ifdef _WIN32
#include <mpi.h>
#else
#include "mpi.h"
#endif

#include "pointers.h"

namespace PHON_NS {
     
    class MyMPI: protected Pointers {
    public:
        MyMPI(class PHON *, MPI_Comm);
        ~MyMPI();

        int my_rank;
        int nprocs;
    };
    
}