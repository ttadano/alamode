#ifndef ALM_POINTERS_HEADER
#define ALM_POINTERS_HEADER

#include "alamode.h"

namespace ALM_NS {

    class Pointers {
    public:
        Pointers(ALM *ptr) :
          alm(ptr),
              memory(ptr->memory) {}
          virtual ~Pointers(){}
    protected:
        ALM *alm;
        Memory *&memory;
    };
}
#endif