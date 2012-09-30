#include "mpi_common.h"
#include <iostream>

using namespace PHON_NS;

MyMPI::MyMPI(PHON *phon, MPI_Comm comm): Pointers(phon) {

    MPI_Comm_rank(comm, &my_rank);
    MPI_Comm_size(comm, &nprocs);

    std::cout << "my_rank = " << my_rank << std::endl;
}

MyMPI::~MyMPI() {}