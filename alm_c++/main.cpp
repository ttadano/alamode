// main.cpp 

#include "alamode.h"
#include "mpi.h"

using namespace ALM_NS;

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    ALM *alm = new ALM(argc, argv, MPI_COMM_WORLD);

    MPI_Finalize();
    return 0;
}