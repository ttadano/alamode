#include "mpi_common.h"
#include "phonons.h"
#include <iostream>

using namespace PHON_NS;

int main(int argc, char **argv)
{
	MPI_Init(&argc, &argv);

	PHON *phon = new PHON(argc, argv, MPI_COMM_WORLD);

	delete phon;

	MPI_Finalize();

	return 0;
}
