/*
 dielec.cpp

 Copyright (c) 2019 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#include "mpi_common.h"
#include "dielec.h"
#include "constants.h"
#include "dynamical.h"
#include "error.h"
#include "mathfunctions.h"
#include "memory.h"
#include "system.h"
#include "write_phonons.h"
#include "parsephon.h"
#include "phonon_dos.h"
#include "fcs_phonon.h"
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
    std::cout << "DIELEC = " << calc_dielectric_constant << std::endl;
}

void Dielec::compute_dielectric_constant()
{
    const auto ns = dynamical->neval;
    const auto nomega = dos->n_energy;
    const auto freq_dyn = dos->energy_dos;
    const auto zstar = dynamical->borncharge;

    double *xk, *kdirec;
    double *eval;
    std::complex<double> **evec;

    memory->allocate(xk, 3);
    memory->allocate(kdirec, 3);
    memory->allocate(eval, ns);
    memory->allocate(evec, ns, ns);

    for (auto i = 0; i < 3; ++i) xk[i] = 0.0;

    rotvec(kdirec, xk, system->rlavec_p, 'T');
    auto norm = kdirec[0] * kdirec[0] + kdirec[1] * kdirec[1] + kdirec[2] * kdirec[2];
    if (norm > eps) {
        for (auto i = 0; i < 3; ++i) kdirec[i] /= std::sqrt(norm);
    }

//    kdirec[0] = 1.0;
    dynamical->eval_k(xk, kdirec, fcs_phonon->fc2_ext, eval, evec, true);

    for (auto i = 0; i < ns; ++i) {
        std::cout << "eval = " << eval[i] << std::endl;
        for (auto j = 0; j < ns; ++j) {
            std::cout << std::setw(15) << evec[i][j].real();
            std::cout << std::setw(15) << evec[i][j].imag();
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }

    for (auto i = 0; i < ns; ++i) {
        for (auto j = 0; j < ns; ++j) {
            evec[i][j] /= std::sqrt(system->mass[system->map_p2s[j / 3][0]]);
        }
    }

    for (auto i = 0; i < ns; ++i) {
        std::cout << "U = " << eval[i] << std::endl;
        for (auto j = 0; j < ns; ++j) {
            std::cout << std::setw(15) << real(evec[i][j]);
            std::cout << std::setw(15) << evec[i][j].imag();
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }

    double ***s_born;
    double **zstar_u;

    memory->allocate(zstar_u, 3, ns);
    std::cout << "Zstar_u:\n";

    for (auto i = 0; i < 3; ++i) {
        for (auto is = 0; is < ns; ++is) {
            zstar_u[i][is] = 0.0;

            for (auto j = 0; j < ns; ++j) {
                zstar_u[i][is] += zstar[j / 3][i][j % 3] * evec[is][j].real();
            }
        }
    }

    for (auto is = 0; is < ns; ++is) {
        for (auto i = 0; i < 3; ++i) {
            std::cout << std::setw(15) << zstar_u[i][is];
        }
        std::cout << '\n';
    }
    std::cout << std::endl;

    memory->allocate(s_born, 3, 3, ns);

    std::cout << "S_born:\n";

    for (auto i = 0; i < 3; ++i) {
        for (auto j = 0; j < 3; ++j) {
            for (auto is = 0; is < ns; ++is) {
                s_born[i][j][is] = zstar_u[i][is] * zstar_u[j][is];
            }
        }
    }

    for (auto is = 0; is < ns; ++is) {
        for (auto i = 0; i < 3; ++i) {
            for (auto j = 0; j < 3; ++j) {
                std::cout << std::setw(15) << s_born[i][j][is];
            }
            std::cout << '\n';
        }
        std::cout << '\n';
    }

    memory->allocate(dielec, nomega, 3, 3);

    auto freq_conv_factor = time_ry * time_ry / (Hz_to_kayser * Hz_to_kayser);
    auto factor = 8.0 * pi / system->volume_p;
    double w2_tmp;
    for (auto iomega = 0; iomega < nomega; ++iomega) {
        w2_tmp = freq_dyn[iomega]*freq_dyn[iomega] * freq_conv_factor;

        for (auto i = 0; i < 3; ++i) { 
            for (auto j = 0; j < 3; ++j) {
                dielec[iomega][i][j] = 0.0;

                for (auto is = 3; is < ns; ++is) {
                    dielec[iomega][i][j] += s_born[i][j][is] / (eval[is] - w2_tmp);
                }
                dielec[iomega][i][j] *= factor;
            }
        }
    }

    std::ofstream ofs_dielec;
    auto file_dielec = input->job_title + ".dielec";

    ofs_dielec.open(file_dielec.c_str(), std::ios::out);
    if (!ofs_dielec) error->exit("write_phonon_vel", "cannot open file_vel");

    for (auto iomega = 0; iomega < nomega; ++iomega) {
        ofs_dielec << std::setw(10) << freq_dyn[iomega];
        for (auto i = 0; i < 3; ++i) {
                ofs_dielec << std::setw(15) << dielec[iomega][i][i];
        }
        ofs_dielec << '\n';
    }
    ofs_dielec << std::endl;
    ofs_dielec.close();

    memory->deallocate(xk);
    memory->deallocate(kdirec);
    memory->deallocate(eval);
    memory->deallocate(evec);

    memory->deallocate(zstar_u);
    memory->deallocate(s_born);



}