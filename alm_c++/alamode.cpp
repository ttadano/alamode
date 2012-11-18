#include <iostream>
#include "interaction.h"
#include "symmetry.h"
#include "input.h"
#include "system.h"
#include "files.h"
#include "memory.h"
#include "timer.h"
#include "fcs.h"
#include "fitting.h"
#include "constraint.h"
#include "timer.h"
#include "writes.h"
#include "ewald.h"

using namespace ALM_NS;

ALM::ALM(int narg, char **arg)
{
    timer = new Timer(this);

    std::cout << std::endl << "Job started at " << timer->DataAndTime() <<  std::endl;
    input = new Input(this, narg, arg);
    create();
    input->parce_input();
    initialize();
    constraint->setup();
    fitting->fitmain();
    writes->writeall();
    finalize();
    std::cout << std::endl << "Job finished at " << timer->DataAndTime() << std::endl;
}

void ALM::create()
{
    memory = new Memory(this);
    files = new Files(this);
    system = new System(this);
    interaction = new Interaction(this);
    fcs = new Fcs(this);
    symmetry = new Symmetry(this);
    fitting = new Fitting(this);
    constraint = new Constraint(this);
    ewald = new Ewald(this);
    writes = new Writes(this);
}

void ALM::initialize()
{
    system->init();
    files->init();
    symmetry->init();
    interaction->init();
    fcs->init();
    ewald->init();
}
ALM::~ALM()
{
    finalize();
    delete input;
    delete timer;
}

void ALM::finalize()
{
    delete memory;
    delete files;
    delete interaction;
    delete fcs;
    delete symmetry;
    delete system;
    delete fitting;
    delete constraint;
    delete ewald;
    delete writes;
}