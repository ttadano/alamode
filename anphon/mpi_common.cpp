/*
 mpi_common.cpp

 Copyright (c) 2014 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#include "mpi_common.h"
#include <string>
#include <cstring>

using namespace PHON_NS;

MyMPI::MyMPI(PHON *phon,
             MPI_Comm comm) : Pointers(phon)
{
    MPI_Comm_rank(comm, &my_rank);
    MPI_Comm_size(comm, &nprocs);
}

MyMPI::~MyMPI() {}

void MyMPI::MPI_Bcast_string(std::string &str,
                             int root,
                             MPI_Comm comm) const
{
    int len = str.length();
    MPI_Bcast(&len, 1, MPI_INT, 0, comm);

    // limited to 512 characters
    char ctmp[512];
    std::strcpy(ctmp, str.c_str());
    MPI_Bcast(&ctmp, len + 1, MPI_CHAR, 0, comm);
    str = std::string(ctmp);
}
