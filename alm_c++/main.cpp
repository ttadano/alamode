// main.cpp 

#include "alamode.h"
#include <iostream>

using namespace ALM_NS;

int main(int argc, char **argv)
{
    std::cout << "Alamode C++ version 0.1" << std::endl;

    ALM *alm = new ALM(argc, argv);
    
    std::cout << "Bye! :)";
    return 0;
}