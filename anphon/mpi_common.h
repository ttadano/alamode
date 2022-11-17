/*
 mpi_common.h

 Copyright (c) 2014 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#pragma once

#ifdef _WIN32
#include <mpi.h>
#else

#include "mpi.h"

#endif

#include <string>
#include "pointers.h"

namespace PHON_NS {
class MyMPI : protected Pointers {
public:
    MyMPI(class PHON *,
          MPI_Comm);

    ~MyMPI();

    void MPI_Bcast_string(std::string &,
                          int,
                          MPI_Comm) const;

    int my_rank;
    int nprocs;
};
}
