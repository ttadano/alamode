// main.cpp 

#include "alamode.h"
#include <iostream>

using namespace ALM_NS;

int main(int argc, char **argv)
{
    std::cout << "Program ALM (version 0.4)" << std::endl << std::endl;

    ALM *alm = new ALM(argc, argv);

    std::cout << "Bye! :)" << std::endl;
    return 0;
}
