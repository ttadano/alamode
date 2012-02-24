#include "phonons.h"
#include <iostream>
#include "timer.h"
#include "parsephon.h"
#include "memory.h"
#include "error.h"

using namespace PHON_NS;

PHON::PHON(int narg, char **arg)
{
    std::cout << std::endl << "Job started at " << timer->DataAndTime() <<  std::endl;
    input = new Input(this, narg, arg);
    create_pointers();
    input->parce_input();
    destroy_pointers();

    std::cout << std::endl << "Job finished at " << timer->DataAndTime() << std::endl;
}
PHON::~PHON(){}

void PHON::create_pointers()
{
  memory = new Memory(this);
  error = new Error(this);
}

void PHON::destroy_pointers()
{
    delete memory;
    delete error;
}