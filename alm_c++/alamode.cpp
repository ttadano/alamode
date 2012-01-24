#include "interaction.h"
#include "symmetry.h"
#include "input.h"
#include "system.h"
#include "files.h"
#include "memory.h"
#include "timer.h"
#include "fcs.h"

using namespace ALM_NS;

ALM::ALM(int narg, char **arg)
{
    timer = new Timer(this);
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
    fcs = new Fcs(this);
    symmetry = new Symmetry(this);
}

void ALM::initialize()
{
    system->init();
    files->init();
    symmetry->init();
    interaction->init();
    fcs->init();
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
}
