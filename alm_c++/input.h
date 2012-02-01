#ifndef ALM_INPUT_HEADER
#define ALM_INPUT_HEADER

#include "pointers.h"

namespace ALM_NS {
    class Input: protected Pointers {
    public:
        Input(class ALM *, int, char **);
        ~Input();
        void sparce_input();
    };
}

#endif