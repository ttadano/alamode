#include "phonons.h"
#include "mpi_common.h"
#include <iostream>

using namespace PHON_NS;

int main(int argc, char **argv)
{
    std::cout << "Phonons program C++ version 0.1" << std::endl;

    MPI_Init(&argc, &argv);
 
    /*
    MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
    MPI_Comm_size(MPI_COMM_WORLD,&mpi_size);

    std::cout << "my_rank = " << my_rank << std::endl;

    MPI_Finalize();
    exit(1);
    */

    PHON *phon = new PHON(argc, argv, MPI_COMM_WORLD);

    MPI_Finalize();
    
    std::cout << "Bye! :)" << std::endl;
    return 0;
}