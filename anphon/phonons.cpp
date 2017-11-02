/*
 phonons.cpp

 Copyright (c) 2014, 2015, 2016 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#include "mpi_common.h"
#include <iostream>
#include <iomanip>
#include "phonons.h"
#include "timer.h"
#include "parsephon.h"
#include "memory.h"
#include "error.h"
#include "gruneisen.h"
#include "system.h"
#include "symmetry_core.h"
#include "kpoint.h"
#include "fcs_phonon.h"
#include "dynamical.h"
#include "phonon_velocity.h"
#include "thermodynamics.h"
#include "write_phonons.h"
#include "phonon_dos.h"
#include "integration.h"
#include "relaxation.h"
#include "conductivity.h"
#include "isotope.h"
#include "selfenergy.h"
#include "version.h"
#include "scph.h"
#include "ewald.h"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace PHON_NS;

PHON::PHON(int narg, char **arg, MPI_Comm comm)
{
    mympi = new MyMPI(this, comm);
    input = new Input(this);

    create_pointers();

    if (mympi->my_rank == 0) {
        std::cout << " +-----------------------------------------------------------------+" << std::endl;
        std::cout << " +                         Program ANPHON                          +" << std::endl;
        std::cout << " +                             Ver.";
        std::cout << std::setw(7) << ALAMODE_VERSION;
        std::cout << "                         +" << std::endl;
        std::cout << " +-----------------------------------------------------------------+" << std::endl;

        std::cout << std::endl;
        std::cout << " Job started at " << timer->DateAndTime() << std::endl;
        std::cout << " The number of MPI processes: " << mympi->nprocs << std::endl;
#ifdef _OPENMP
        std::cout << " The number of OpenMP threads: " 
            << omp_get_max_threads() << std::endl;
#endif
        std::cout << std::endl;

        input->parce_input(narg, arg);
        writes->write_input_vars();
    }

    mympi->MPI_Bcast_string(input->job_title, 0, MPI_COMM_WORLD);
    mympi->MPI_Bcast_string(mode, 0, MPI_COMM_WORLD);


    if (mode == "PHONONS") {

        execute_phonons();

    } else if (mode == "RTA") {

        execute_RTA();

    } else if (mode == "SCPH") {

        execute_self_consistent_phonon();

    } else {
        error->exit("phonons", "invalid mode: ", mode.c_str());
    }

    if (mympi->my_rank == 0) {
        std::cout << std::endl << " Job finished at "
            << timer->DateAndTime() << std::endl;
    }
    destroy_pointers();
}

PHON::~PHON()
{
    delete input;
    delete mympi;
}

void PHON::create_pointers()
{
    memory = new Memory(this);
    timer = new Timer(this);
    error = new Error(this);
    system = new System(this);
    symmetry = new Symmetry(this);
    kpoint = new Kpoint(this);
    fcs_phonon = new Fcs_phonon(this);
    dynamical = new Dynamical(this);
    integration = new Integration(this);
    phonon_velocity = new Phonon_velocity(this);
    thermodynamics = new Thermodynamics(this);
    relaxation = new Relaxation(this);
    selfenergy = new Selfenergy(this);
    conductivity = new Conductivity(this);
    writes = new Writes(this);
    dos = new Dos(this);
    gruneisen = new Gruneisen(this);
    isotope = new Isotope(this);
    scph = new Scph(this);
    ewald = new Ewald(this);
}

void PHON::destroy_pointers()
{
    delete memory;
    delete timer;
    delete error;
    delete system;
    delete symmetry;
    delete kpoint;
    delete fcs_phonon;
    delete dynamical;
    delete integration;
    delete phonon_velocity;
    delete thermodynamics;
    delete relaxation;
    delete selfenergy;
    delete conductivity;
    delete writes;
    delete dos;
    delete gruneisen;
    delete isotope;
    delete scph;
    delete ewald;
}

void PHON::setup_base()
{
    system->setup();
    symmetry->setup_symmetry();
    kpoint->kpoint_setups(mode);
    fcs_phonon->setup(mode);
    dynamical->setup_dynamical(mode);
    dos->setup();
    ewald->init();
    if (mympi->my_rank == 0) {
        std::cout << " Now, move on to phonon calculations." << std::endl;
    }
}

void PHON::execute_phonons()
{
    if (mympi->my_rank == 0) {
        std::cout << "                      MODE = phonons                         " << std::endl;
        std::cout << "                                                             " << std::endl;
        std::cout << "      Phonon calculation within harmonic approximation       " << std::endl;
        std::cout << "      Harmonic force constants will be used.                 " << std::endl;

        if (gruneisen->print_gruneisen) {
            std::cout << std::endl;
            std::cout << "      GRUNEISEN = 1 : Cubic force constants are necessary." << std::endl;
        }
        std::cout << std::endl;
    }

    setup_base();

    dynamical->diagonalize_dynamical_all();
    phonon_velocity->calc_group_velocity(kpoint->kpoint_mode);

    if (dos->flag_dos) {
        integration->setup_integration();
        dos->calc_dos_all();
    }

    gruneisen->setup();
    if (gruneisen->print_gruneisen) {
        gruneisen->calc_gruneisen();
    }

    if (mympi->my_rank == 0) {
        writes->print_phonon_energy();
        writes->write_phonon_info();

        if (gruneisen->print_newfcs) {
            gruneisen->write_new_fcsxml_all();
        }
    }

    dynamical->finish_dynamical();
    gruneisen->finish_gruneisen();

    if (dos->flag_dos) {
        integration->finish_integration();
    }
}

void PHON::execute_RTA()
{
    if (mympi->my_rank == 0) {
        std::cout << "                        MODE = RTA                           " << std::endl;
        std::cout << "                                                             " << std::endl;
        std::cout << "      Calculation of phonon line width (lifetime) and        " << std::endl;
        std::cout << "      lattice thermal conductivity within the RTA            " << std::endl;
        std::cout << "      (relaxation time approximation).                       " << std::endl;
        std::cout << "      Harmonic and anharmonic force constants will be used.  " << std::endl;
        std::cout << std::endl;

        if (restart_flag) {
            std::cout << std::endl;
            std::cout << "      Restart mode is switched on!                                  " << std::endl;
            std::cout << "      The calculation will be restart from the existing result file." << std::endl;
            std::cout << "      If you want to start a calculation from scratch,              " << std::endl;
            std::cout << "      please set RESTART = 0 in the input file                      " << std::endl;
            std::cout << std::endl;
        }
    }

    setup_base();

    if (kpoint->kpoint_mode < 3) {
        dynamical->diagonalize_dynamical_all();
    }
    if (kpoint->kpoint_mode == 2) {
        integration->setup_integration();
    }
    isotope->setup_isotope_scattering();
    isotope->calc_isotope_selfenergy_all();

    //   relaxation->setup_mode_analysis();
    relaxation->setup_relaxation();
    selfenergy->setup_selfenergy();

    if (relaxation->ks_analyze_mode) {
        relaxation->perform_mode_analysis();
    } else {
        writes->setup_result_io();
        conductivity->setup_kappa();
        conductivity->prepare_restart();
        conductivity->calc_anharmonic_imagself();
        conductivity->compute_kappa();
        writes->write_kappa();
        writes->write_selfenergy_isotope();
    }

    if (kpoint->kpoint_mode == 2) {
        integration->finish_integration();
    }

    dynamical->finish_dynamical();
    relaxation->finish_relaxation();

    if (!relaxation->ks_analyze_mode) {
        conductivity->finish_kappa();
    }
}

void PHON::execute_self_consistent_phonon()
{
    if (mympi->my_rank == 0) {
        std::cout << "                        MODE = SCPH                           " << std::endl;
        std::cout << "                                                             " << std::endl;
        std::cout << "      Self-consistent phonon calculation to estimate         " << std::endl;
        std::cout << "      anharmonic phonon frequencies.                         " << std::endl;
        std::cout << "      Harmonic and quartic force constants will be used.  " << std::endl;
        std::cout << std::endl;
    }

    setup_base();
    dos->setup();

    dynamical->diagonalize_dynamical_all();

    if (kpoint->kpoint_mode == 2) {
        integration->setup_integration();
    }

    scph->setup_scph();
    scph->exec_scph();

    if (kpoint->kpoint_mode == 2) {
        integration->finish_integration();
    }

    dynamical->finish_dynamical();

    scph->finish_scph();
}
