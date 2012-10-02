#include <iostream>
#include <cstdlib>
#include <string>
#include "error.h"
#include "mpi_common.h"

using namespace PHON_NS;

Error::Error(PHON *phon): Pointers(phon) {}

Error::~Error() {}

void Error::warn(const char *file, const char *message)
{
    std::cout << "WARNING in " << file << "  MESSAGE: " << message << std::endl;
}

void Error::exit(const char *file, const char *message)
{
    std::cout << "ERROR in " << file << "  MESSAGE: " << message << std::endl;
    std::exit(EXIT_FAILURE);
}

void Error::exit(const char *file, const char *message, int info)
{
    MPI_Finalize();
    std::cout << "ERROR in " << file << "  MESSAGE: " << message << info << std::endl;
    std::exit(EXIT_FAILURE);
}

void Error::exitall(const char *file, const char *message)
{
    MPI_Finalize();
    std::cout << "ERROR in " << file << "  MESSAGE: " << message << std::endl;
    std::exit(EXIT_FAILURE);
}
