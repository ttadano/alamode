/*
 error.cpp

 Copyright (c) 2014 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#include "mpi_common.h"
#include <iostream>
#include <cstdlib>
#include <string>
#include "error.h"

using namespace PHON_NS;

Error::Error(PHON *phon) : Pointers(phon) {}

Error::~Error() {}

void Error::warn(const char *file,
                 const char *message) const
{
    std::cout << "WARNING in " << file << "  MESSAGE: " << message << std::endl;
}

void Error::exit(const char *file,
                 const char *message) const
{
    std::cout << "ERROR in " << file << "  MESSAGE: " << message << std::endl;
    std::exit(EXIT_FAILURE);
}

void Error::exit(const char *file,
                 const char *message,
                 int info) const
{
    MPI_Finalize();
    std::cout << "ERROR in " << file << "  MESSAGE: " << message << info << std::endl;
    std::exit(EXIT_FAILURE);
}

void Error::exitall(const char *file,
                    const char *message) const
{
    MPI_Finalize();
    std::cout << "ERROR in " << file << "  MESSAGE: " << message << std::endl;
    std::exit(EXIT_FAILURE);
}

void Error::exit(const char *file,
                 const char *message,
                 const char *info) const
{
    std::cout << "ERROR in " << file << "  MESSAGE: " << message << info << std::endl;
    std::exit(EXIT_FAILURE);
}
