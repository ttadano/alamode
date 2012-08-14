#pragma once

#include "pointers.h"

namespace PHON_NS {
    class Relaxation: protected Pointers {
    public:
        Relaxation(class PHON *);
        ~Relaxation();

        void setup_relaxation();

    private:
        void parse_anharmonic_fcs();
    };
}