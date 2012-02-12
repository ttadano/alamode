#pragma once

#include "pointers.h"

namespace PHON_NS {
    class Input: protected Pointers {
    public:
        Input(class PHON *, int, char **);
        ~Input();
        void sparce_input();
    };
}