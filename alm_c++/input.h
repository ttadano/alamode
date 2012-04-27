#pragma once

#include "pointers.h"

namespace ALM_NS {
    class Input: protected Pointers {
    public:
        Input(class ALM *, int, char **);
        ~Input();
        void parce_input();
    };
}