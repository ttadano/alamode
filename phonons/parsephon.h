#pragma once

#include "pointers.h"
#include <string>

namespace PHON_NS {
    class Input: protected Pointers {
    public:
        Input(class PHON *, int, char **);
        ~Input();
        void parce_input();

        std::string job_title;

    private:
        void read_input_phonons();
    };
}