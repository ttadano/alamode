// main.cpp 

#include <stdlib.h>
#include <iostream>
#include "alamode.h"

using namespace ALM_NS;

int main(int argc, char **argv)
{
//    _CrtSetDbgFlag ( _CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF );

    ALM *alm = new ALM(argc, argv);

//    _CrtDumpMemoryLeaks();
    return 0;
}
    