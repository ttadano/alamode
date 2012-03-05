#pragma once

#include "pointers.h"
#include <string>

namespace PHON_NS {
    class Fcs_phonon: protected Pointers {
    public:
        Fcs_phonon(class PHON *);
        ~Fcs_phonon();

        void setup();

        std::string file_fcs;
        double ****fc2;

    private:
        void load_fc2();
        unsigned int coordinate_index(const char);
    };
}