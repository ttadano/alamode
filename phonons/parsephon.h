#pragma once

#include "pointers.h"

namespace PHON_NS {
    class Input: protected Pointers {
    public:
        Input(class PHON *, int, char **);
        ~Input();
        void parce_input();
    private:
        void read_input_dispersion();
    };
}