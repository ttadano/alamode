// main.cpp 

#include "mpi.h"
#include "alamode.h"
#include "input.h"

using namespace ALM_NS;

int main(int argc, char **argv)
{
	MPI_Init(&argc, &argv);
	ALM *alm = new ALM(argc, argv, MPI_COMM_WORLD);
	alm->memory;
	MPI_Finalize();
	return 0;
}