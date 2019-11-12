/*
 dielec.cpp

 Copyright (c) 2019 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#include "mpi_common.h"
#include "dynamical.h"
#include "error.h"
#include "mathfunctions.h"
#include "memory.h"
#include "system.h"
#include "write_phonons.h"
#include "phonon_dos.h"
#include "fcs_phonon.h"
#include "dielec.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <complex>

using namespace PHON_NS;

Dielec::Dielec(PHON *phon): Pointers(phon)
{
    set_default_variables();
}

Dielec::~Dielec()
{
    deallocate_variables();
}

void Dielec::set_default_variables()
{
    calc_dielectric_constant = 0;
    dielec = nullptr;
}

void Dielec::deallocate_variables()
{
    if (dielec) {
        memory->deallocate(dielec);
    }
}

void Dielec::init()
{
    MPI_Bcast(&calc_dielectric_constant, 1, MPI_INT, 0, MPI_COMM_WORLD);
}

void Dielec::compute_dielectric_constant()
{
    unsigned int ns;
    std::cout << "ok\n";
    ns = dynamical->neval;

        std::cout << "ok\n";


    double *xk, *kvec;
    double *eval;
    std::complex<double> **evec;



    memory->allocate(xk, 3);
    memory->allocate(kvec, 3);
    memory->allocate(eval, ns);
    memory->allocate(evec, ns, ns);

    std::cout << "ok\n";

    for (auto i = 0; i < 3; ++i) {
        xk[i] = 0.0;
        kvec[i] = 0.0;
    }

    // dynamical->eval_k(xk, kvec, fcs_phonon->fc2_ext, eval, evec, true);

    // for (auto i = 0; i < ns; ++i) {
    //     std::cout << "eval = " << eval[i] << std::endl;
    //     for (auto j = 0; j < ns; ++j) {
    //         std::cout << std::setw(15) << evec[i][j].real();
    //         std::cout << std::setw(15) << evec[i][j].imag();
    //         std::cout << std::endl;
    //     }
    //     std::cout << std::endl;
    // }

}