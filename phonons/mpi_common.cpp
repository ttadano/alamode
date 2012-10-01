#include "mpi_common.h"
#include <iostream>
#include <string>
#include <cstring>

using namespace PHON_NS;

MyMPI::MyMPI(PHON *phon, MPI_Comm comm): Pointers(phon) {

    MPI_Comm_rank(comm, &my_rank);
    MPI_Comm_size(comm, &nprocs);
}

MyMPI::~MyMPI() {}

void MyMPI::MPI_Bcast_string(std::string &str, int root, MPI_Comm comm)
{
    int len;
    len = str.length();
    MPI_Bcast(&len, 1, MPI_INT, 0, comm);
    
    char ctmp[len + 1];
    std::strcpy(ctmp, str.c_str());
    MPI_Bcast(&ctmp, len + 1, MPI_CHAR, 0, comm);
    str = std::string(ctmp);
}
