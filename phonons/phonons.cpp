#include "phonons.h"
#include <iostream>
#include "timer.h"
#include "parsephon.h"
#include "memory.h"
#include "error.h"
#include "system.h"
#include "kpoint.h"
#include "fcs_phonon.h"
#include "dynamical.h"
#include "phonon_velocity.h"
#include "phonon_thermodynamics.h"
#include "write_phonons.h"
#include "phonon_dos.h"
#include "integration.h"
#include "relaxation.h"

using namespace PHON_NS;

PHON::PHON(int narg, char **arg)
{
    std::cout << std::endl << "Job started at " << timer->DataAndTime() <<  std::endl;
    input = new Input(this, narg, arg);
    create_pointers();
    input->parce_input();

    if (mode == "phonons") {

        system->setup();
        kpoint->kpoint_setups();
        fcs_phonon->setup(mode);
        dos->setup();

        dynamical->setup_dynamical(mode);
        dynamical->diagonalize_dynamical_all();

        // Calculate the group velocity of phonons along given direction in
        // the reciprocal space.
        if (kpoint->kpoint_mode == 1) phonon_velocity->calc_phonon_vel_band();

        if (dos->flag_dos) {
            integration->setup_integration();
            dos->calc_dos();
        }

        writes->write_phonon_info();

        memory->deallocate(dynamical->evec_phonon);
        memory->deallocate(dynamical->eval_phonon);
        if (kpoint->kpoint_mode == 1) {
            memory->deallocate(phonon_velocity->phvel);
        }

        integration->finish_integration();

    } else if (mode == "boltzmann") {

        system->setup();
        kpoint->kpoint_setups();
        fcs_phonon->setup(mode);
        dynamical->setup_dynamical(mode);

        integration->setup_integration();
        relaxation->setup_relaxation();

        dynamical->diagonalize_dynamical_all();
        relaxation->calc_ReciprocalV();
     //   relaxation->calc_selfenergy(100.0);
        phonon_thermodynamics->test_fB(100.0);

    } else {
        error->exit("phonons", "invalid mode");
    }

    destroy_pointers();

    std::cout << std::endl << "Job finished at " << timer->DataAndTime() << std::endl;
}

PHON::~PHON(){
    delete input;
}

void PHON::create_pointers()
{
    memory = new Memory(this);
    error = new Error(this);
    system = new System(this);
    kpoint = new Kpoint(this);
    fcs_phonon = new Fcs_phonon(this);
    dynamical = new Dynamical(this);
    integration = new Integration(this);
    phonon_velocity = new Phonon_velocity(this);
    phonon_thermodynamics = new Phonon_thermodynamics(this);
    relaxation = new Relaxation(this);
    writes = new Writes(this);
    dos = new Dos(this);
}

void PHON::destroy_pointers()
{
    delete memory;
    delete error;
    delete system;
    delete kpoint;
    delete fcs_phonon;
    delete dynamical;
    delete integration;
    delete phonon_velocity;
    delete phonon_thermodynamics;
    delete relaxation;
    delete writes;
    delete dos;
}
