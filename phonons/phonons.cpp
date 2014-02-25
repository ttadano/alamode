#include "mpi_common.h"
#include <iostream>
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
#include "phonon_thermodynamics.h"
#include "write_phonons.h"
#include "phonon_dos.h"
#include "integration.h"
#include "interpolation.h"
#include "relaxation.h"
#include "conductivity.h"
#include <omp.h>
#include "isotope.h"
#include "selfenergy.h"

using namespace PHON_NS;

PHON::PHON(int narg, char **arg, MPI_Comm comm)
{
    mympi = new MyMPI(this, comm);
    input = new Input(this);

    create_pointers();

    if (mympi->my_rank == 0) {
        std::cout << " +------------------------------------------------------------+" << std::endl;
        std::cout << " +                      Program PHONONS                       +" << std::endl;
        std::cout << " +                           Ver. 1.1.0                       +" << std::endl;
        std::cout << " +------------------------------------------------------------+" << std::endl;

        std::cout << std::endl;
        std::cout << " Job started at " << timer->DateAndTime() <<  std::endl;
        std::cout << " The number of MPI threads: " << mympi->nprocs << std::endl;
#ifdef _OPENMP
        std::cout << " The number of OpenMP threads: " << omp_get_max_threads() << std::endl;
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

    } else if (mode == "INTERPOLATION") {

        execute_interpolation();

    } else if (mode == "GRUNEISEN") {

        execute_gruneisen();

    } else {

        error->exit("phonons", "invalid mode");

    }

    if (mympi->my_rank == 0) {
        std::cout << std::endl << " Job finished at " << timer->DateAndTime() << std::endl;
    }
    destroy_pointers();
}

PHON::~PHON(){
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
    phonon_thermodynamics = new Phonon_thermodynamics(this);
    relaxation = new Relaxation(this);
    selfenergy = new Selfenergy(this);
    conductivity = new Conductivity(this);
    interpolation = new Interpolation(this);
    writes = new Writes(this);
    dos = new Dos(this);
    gruneisen = new Gruneisen(this);
    isotope = new Isotope(this);
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
    delete phonon_thermodynamics;
    //	delete relaxation;
    //	delete selfenergy;
    delete interpolation;
    delete conductivity;
    delete writes;
    delete dos;
    delete gruneisen;
    delete isotope;
}


void PHON::setup_base()
{
    system->setup();
    symmetry->setup_symmetry();
    kpoint->kpoint_setups(mode);
    fcs_phonon->setup(mode);
    dynamical->setup_dynamical(mode);

}

void PHON::execute_phonons()
{
    if (mympi->my_rank == 0) {
        std::cout << "                      MODE = phonons                         " << std::endl;
        std::cout << "                                                             " << std::endl;
        std::cout << "      Phonon calculation within harmonic approximation       " << std::endl;
        std::cout << "      Harmonic force constants will be used.                 " << std::endl;
        std::cout << std::endl;
    }

    setup_base();

    dos->setup();
    dynamical->diagonalize_dynamical_all();

    if (kpoint->kpoint_mode == 1)  {
        phonon_velocity->calc_phonon_vel_band();
    }

    if (dos->flag_dos) {
        integration->setup_integration();
        dos->calc_dos2();
    }

    if (mympi->my_rank == 0) {
        writes->write_phonon_info();
    }

    dynamical->finish_dynamical();

    if (kpoint->kpoint_mode == 1) {
        memory->deallocate(phonon_velocity->phvel);
    }

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

    dos->setup();

    if (kpoint->kpoint_mode < 3) {
        dynamical->diagonalize_dynamical_all();
    }
    relaxation->setup_mode_analysis();

    if (!relaxation->ks_analyze_mode) {
        writes->setup_result_io();
    }

    if (kpoint->kpoint_mode == 2) {
        integration->setup_integration();
    }

    relaxation->setup_relaxation();
    selfenergy->setup_selfenergy();
    isotope->setup_isotope_scattering();
    isotope->calc_isotope_selfenergy_all();

    if (relaxation->ks_analyze_mode) {
        relaxation->compute_mode_tau();
    } else {
        conductivity->setup_kappa();
        conductivity->prepare_restart();
        conductivity->calc_anharmonic_tau();
        conductivity->compute_kappa();
        writes->write_kappa();
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

void PHON::execute_interpolation()
{
    setup_base();

    dos->setup();
    dynamical->diagonalize_dynamical_all();

    interpolation->prepare_interpolation();
    interpolation->exec_interpolation();
    interpolation->finish_interpolation();

    // 
    // 		integration->setup_integration();
    // 		relaxation->setup_relaxation();
}

void PHON::execute_gruneisen()
{

    if (mympi->my_rank == 0) {
        std::cout << "                      MODE = Gruneisen                       " << std::endl;
        std::cout << "                                                             " << std::endl;
        std::cout << "             Calculation of Gruneisen parameters.            " << std::endl;
        std::cout << "      Harmonic and anharmonic force constants will be used.  " << std::endl;
        std::cout << std::endl;
    }

    setup_base();

    dos->setup();
    dynamical->diagonalize_dynamical_all();

    if (mympi->my_rank == 0) {
        gruneisen->setup();
        gruneisen->calc_gruneisen();
        // gruneisen->calc_gruneisen2();
        writes->write_gruneisen();
        gruneisen->finish_gruneisen();
    }
}
