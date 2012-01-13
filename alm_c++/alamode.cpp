#include "memory.h"
#include "alamode.h"
#include "interactions.h"
#include "symmetry.h"
#include "input.h"
#include "system.h"
#include "files.h"
#include <iostream>

using namespace ALM_NS;

ALM::ALM(int narg, char **arg)
{
    input = new Input(this, narg, arg);
    create();
    input->sparce_input();
    initialize();
    finalize();
}

void ALM::create()
{
    memory = new Memory(this);
    files = new Files(this);
    system = new System(this);
    interaction = new Interaction(this);
    symmetry = new Symmetry(this);
}

void ALM::initialize()
{
    system->init();
}
ALM::~ALM()
{
    finalize();
    delete input;
}

void ALM::finalize()
{
    delete memory;
    delete files;
    delete interaction;
    delete symmetry;
    delete system;
}
