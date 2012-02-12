#include "phonons.h"
#include <iostream>

using namespace PHON_NS;

int main(int argc, char **argv)
{
    std::cout << "Phonons program C++ version 0.1" << std::endl;

    PHON *phon = new PHON(argc, argv);
    
    std::cout << "Bye! :)" << std::endl;
    return 0;
}