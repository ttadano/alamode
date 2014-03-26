/*
 main.cpp

 Copyright (c) 2014 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

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
    