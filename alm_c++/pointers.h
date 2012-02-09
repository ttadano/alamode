#ifndef ALM_POINTERS_HEADER
#define ALM_POINTERS_HEADER

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
          files(ptr->files),
          error(ptr->error),
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
        Files *&files;
        Writes *&writes;
        Error *&error;
        Timer *&timer;
    };
}
#endif