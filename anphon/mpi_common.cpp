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
#include <Eigen/Core>
#include <iostream>

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
    MPI_Bcast(&len, 1, MPI_INT, root, comm);

    // limited to 512 characters
    char ctmp[512];
    std::strcpy(ctmp, str.c_str());
    MPI_Bcast(&ctmp, len + 1, MPI_CHAR, root, comm);
    str = std::string(ctmp);
}


void MyMPI::MPI_Bcast_CellClass(Cell &cell,
                                int root,
                                MPI_Comm comm) const
{
    MPI_Bcast(cell.lattice_vector.data(), 9, MPI_DOUBLE, root, comm);

    int natoms;
    if (my_rank == root) {
        natoms = cell.x_fractional.rows();
    }
    MPI_Bcast(&natoms, 1, MPI_INT, root, comm);

    if (my_rank != root) {
        cell.x_fractional.resize(natoms, 3);
        cell.kind.resize(natoms);
    }
    MPI_Bcast(cell.x_fractional.data(), 3 * natoms, MPI_DOUBLE, root, comm);
    MPI_Bcast(cell.kind.data(), natoms, MPI_INT, root, comm);
}

void MyMPI::MPI_Bcast_SpinClass(Spin &spin,
                                int root,
                                MPI_Comm comm) const
{
    MPI_Bcast(&spin.lspin, 1, MPI_INT, root, comm);
    MPI_Bcast(&spin.noncollinear, 1, MPI_INT, root, comm);
    MPI_Bcast(&spin.time_reversal_symm, 1, MPI_INT, root, comm);

    int natoms;
    if (my_rank == root) {
        natoms = spin.magmom.size();
    }
    MPI_Bcast(&natoms, 1, MPI_INT, root, comm);
    if (my_rank != root) {
        spin.magmom.resize(natoms, std::vector<double>(3));
    }
    MPI_Bcast(spin.magmom.data(), 3 * natoms, MPI_DOUBLE, root, comm);
}

void MyMPI::MPI_Bcast_MappingTable(MappingTable &mapping,
                                   int root,
                                   MPI_Comm comm) const
{
    int natmin_tmp, ntran_tmp;
    if (my_rank == root) {
        natmin_tmp = mapping.from_true_primitive.size();
        ntran_tmp = mapping.from_true_primitive[0].size();
    }
    MPI_Bcast(&natmin_tmp, 1, MPI_INT, 0, comm);
    MPI_Bcast(&ntran_tmp, 1, MPI_INT, 0, comm);

    if (mympi->my_rank != root) {
        mapping.from_true_primitive.resize(natmin_tmp, std::vector<unsigned int>(ntran_tmp));
        mapping.to_true_primitive.resize(natmin_tmp * ntran_tmp);
    }
    MPI_Bcast(&mapping.from_true_primitive[0][0],
              natmin_tmp * ntran_tmp, MPI_UNSIGNED, root, comm);
    MPI_Bcast(&mapping.to_true_primitive[0],
              natmin_tmp * ntran_tmp * sizeof mapping.to_true_primitive[0],
              MPI_BYTE, root, comm);
}

