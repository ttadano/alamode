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
          phonon_velocity(ptr->phonon_velocity),
          writes(ptr->writes),
          dos(ptr->dos),
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
        Phonon_velocity *&phonon_velocity;
        Writes *&writes;
        Dos *&dos;
        Timer *&timer;
    };
}