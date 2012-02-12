#ifndef PHON_INPUT_HEADER
#define PHON_INPUT_HEADER

#include "pointers.h"

namespace PHON_NS {
    class Input: protected Pointers {
    public:
        Input(class PHON *, int, char **);
        ~Input();
        void sparce_input();
    };
}

#endif