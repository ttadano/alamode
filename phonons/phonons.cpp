#include "phonons.h"
#include <iostream>
#include "timer.h"

using namespace PHON_NS;

PHON::PHON(int narg, char **arg)
{
    std::cout << std::endl << "Job started at " << timer->DataAndTime() <<  std::endl;
    std::cout << std::endl << "Job finished at " << timer->DataAndTime() << std::endl;
}
PHON::~PHON(){}