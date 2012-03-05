#pragma once

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
        class Fcs_phonon *fcs_phonon;

        PHON(int, char**);
        ~PHON();

        void create_pointers();
        void destroy_pointers();
    };
}