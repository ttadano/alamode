/* Declaration of pointers used in the whole program. */

#pragma once
#include <string>

namespace ALM_NS {

    class ALM {
    public:
        class Memory *memory;
        class Input *input;
        class System *system;
        class Interaction *interaction;
        class Fcs *fcs;
        class Symmetry *symmetry;
        class Fitting *fitting;
        class Constraint *constraint;
        class Files *files;
        class Ewald *ewald;
        class Displace *displace;
        class Writes *writes;
        class Error *error;
        class Timer *timer;
        ALM(int, char **);
        ~ALM();
        void create();
        void initialize();
        void finalize();

        std::string mode;
    };
}
