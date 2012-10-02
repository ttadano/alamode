#include "mpi_common.h"
#include "phonons.h"
#include <iostream>
//#include <boost/mpi.hpp>

using namespace PHON_NS;

//namespace mpi = boost::mpi;

int main(int argc, char **argv)
{
    std::cout << "Phonons program version 1.0 (MPI)" << std::endl;

    MPI_Init(&argc, &argv);

    PHON *phon = new PHON(argc, argv, MPI_COMM_WORLD);

    MPI_Finalize();

    std::cout << "Bye! :)" << std::endl;
    return 0;
}
