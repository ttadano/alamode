#pragma once

#include "alamode.h"

namespace ALM_NS {

    class Pointers {
    public:
        Pointers(ALM *ptr) :
          alm(ptr),
              memory(ptr->memory),
              input(ptr->input),
              system(ptr->system),
          interaction(ptr->interaction),
          symmetry(ptr->symmetry),
          fitting(ptr->fitting),
          constraint(ptr->constraint),
          files(ptr->files),
          error(ptr->error),
          ewald(ptr->ewald),
		  displace(ptr->displace),
          fcs(ptr->fcs),
          writes(ptr->writes),
          timer(ptr->timer) {}
          virtual ~Pointers(){}
    protected:
        ALM *alm;
        Memory *&memory;
        Input *&input;
        System *&system;
        Interaction *&interaction;
        Fcs *&fcs;
        Symmetry *&symmetry;
        Fitting *&fitting;
        Constraint *&constraint;
        Files *&files;
        Ewald *&ewald;
		Displace *&displace;
        Writes *&writes;
        Error *&error;
        Timer *&timer;
    };
}
