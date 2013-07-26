// main.cpp 

#include "alamode.h"
#include <iostream>

using namespace ALM_NS;

int main(int argc, char **argv)
{
    std::cout << "Alamode C++ version 0.2: svn rev. 319" << std::endl << std::endl;

    ALM *alm = new ALM(argc, argv);
    
    std::cout << "Bye! :)" << std::endl;
    return 0;
}
