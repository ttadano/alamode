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
          files(ptr->files),
          error(ptr->error),
          timer(ptr->timer) {}
          virtual ~Pointers(){}
    protected:
        ALM *alm;
        Memory *&memory;
        Input *&input;
        System *&system;
        Interaction *&interaction;
        Symmetry *&symmetry;
        Files *&files;
        Error *&error;
        Timer *&timer;
    };
}
#endif