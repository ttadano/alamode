#ifndef PHON_POINTERS_HEADER
#define PHON_POINTERS_HEADER

#include "phonons.h"

namespace PHON_NS {

    class Pointers {
    public:
        Pointers(PHON *ptr) :
          phon(ptr),
              memory(ptr->memory),
          error(ptr->error),
          timer(ptr->timer) {}
          virtual ~Pointers(){}
    protected:
        PHON *phon;
        Memory *&memory;
        Error *&error;
        Timer *&timer;
    };
}
#endif