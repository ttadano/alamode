// main.cpp 

#define _CRTDBG_MAP_ALLOC

#include <stdlib.h>

#include <crtdbg.h>

#include "alamode.h"
#include <iostream>

using namespace ALM_NS;

int main(int argc, char **argv)
{
//    _CrtSetDbgFlag(_CrtSetDbgFlag(0) | _CRTDBG_LEAK_CHECK_DF);
    _CrtSetDbgFlag ( _CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF );

    ALM *alm = new ALM(argc, argv);

    _CrtDumpMemoryLeaks();
    return 0;
}
    