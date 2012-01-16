#include "error.h"
#include <cstdlib>
#include <iostream>

using namespace ALM_NS;

Error::Error(ALM *alm): Pointers(alm) {}

Error::~Error() {}

void Error::warn(const char *file, const char *message)
{
    std::cout << "WARNING in " << file << "    Message:" << message << std::endl;
}


void Error::exit(const char *file, const char *message)
{
    std::cout << "ERROR in " << file << "    Message:" << message << std::endl;
    std::exit(1);
}

void Error::exit(const char *file, const char *message, int info)
{
    std::cout << "ERROR in " << file << "    Message:" << message << info << std::endl;
    std::exit(1);
}