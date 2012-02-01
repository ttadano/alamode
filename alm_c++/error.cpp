#include <iostream>
#include <cstdlib>
#include <string>
#include "error.h"

using namespace ALM_NS;

Error::Error(ALM *alm): Pointers(alm) {}

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
    std::cout << "ERROR in " << file << "  MESSAGE: " << message << info << std::endl;
    std::exit(EXIT_FAILURE);
}