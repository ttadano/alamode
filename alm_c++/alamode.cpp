#include "memory.h"
#include "alamode.h"
#include "mpi.h"
#include <iostream>

using namespace ALM_NS;

ALM::ALM(int narg, char **arg, MPI_Comm communicator)
{
   	memory = new Memory(this);
    std::cout << "OK" << std::endl;
}