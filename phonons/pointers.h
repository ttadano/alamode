#pragma once

#include "phonons.h"

namespace PHON_NS {

    class Pointers {
    public:
        Pointers(PHON *ptr) :
          phon(ptr),
              memory(ptr->memory),
          error(ptr->error),
          input(ptr->input),
          system(ptr->system),
          kpoint(ptr->kpoint),
          fcs_phonon(ptr->fcs_phonon),
          dynamical(ptr->dynamical),
          timer(ptr->timer) {}
          virtual ~Pointers(){}
    protected:
        PHON *phon;
        Memory *&memory;
        Error *&error;
        Input *&input;
        System *&system;
        Kpoint *&kpoint;
        Fcs_phonon *&fcs_phonon;
        Dynamical *&dynamical;
        Timer *&timer;
    };
}