#include "relaxation.h"

using namespace PHON_NS;

Relaxation::Relaxation(PHON *phon): Pointers(phon) {}
Relaxation::~Relaxation(){};


void Relaxation::setup_relaxation()
{
    parse_anharmonic_fcs();
}

void Relaxation::parse_anharmonic_fcs()
{

}