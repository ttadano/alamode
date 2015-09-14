/*
 memory.cpp

 Copyright (c) 2014 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#include "mpi_common.h"
#include "memory.h"

using namespace PHON_NS;

Memory::Memory(PHON *phon) : Pointers(phon) {}

Memory::~Memory() {};
