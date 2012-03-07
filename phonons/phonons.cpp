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
#include "write_phonons.h"

using namespace PHON_NS;

PHON::PHON(int narg, char **arg)
{
    std::cout << std::endl << "Job started at " << timer->DataAndTime() <<  std::endl;
    input = new Input(this, narg, arg);
    create_pointers();
    input->parce_input();

    system->setup();
    kpoint->kpoint_setups();
    fcs_phonon->setup();

    dynamical->calc_dynamical_matrix();
    dynamical->diagonalize_dynamical();

    writes->write_phonon_info();

    destroy_pointers();

    std::cout << std::endl << "Job finished at " << timer->DataAndTime() << std::endl;
}
PHON::~PHON(){}

void PHON::create_pointers()
{
  memory = new Memory(this);
  error = new Error(this);
  system = new System(this);
  kpoint = new Kpoint(this);
  fcs_phonon = new Fcs_phonon(this);
  dynamical = new Dynamical(this);
  writes = new Writes(this);
}

void PHON::destroy_pointers()
{
    delete memory;
    delete error;
    delete system;
    delete kpoint;
    delete fcs_phonon;
    delete dynamical;
    delete writes;
}