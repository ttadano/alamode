#pragma once
#include <string>

namespace PHON_NS
{
    class PHON {
    public:
        class Memory *memory;
        class Error *error;
        class Timer *timer;
        class Input *input;
        class System *system;
        class Kpoint *kpoint;
        class Integration *integration;
        class Fcs_phonon *fcs_phonon;
        class Dynamical *dynamical;
        class Phonon_velocity *phonon_velocity;
        class Phonon_thermodynamics *phonon_thermodynamics;
        class Writes *writes;
        class Dos *dos;

        PHON(int, char**);
        ~PHON();

        void create_pointers();
        void destroy_pointers();

        std::string mode;
    };
}