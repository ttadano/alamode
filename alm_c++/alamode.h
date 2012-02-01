/* Declaration of pointers used in the whole program. */

#ifndef ALM_ALAMODE_HEADER
#define ALM_ALAMODE_HEADER

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
        class Files *files;
        class Error *error;
        class Timer *timer;
        ALM(int, char **);
        ~ALM();
        void create();
        void initialize();
        void finalize();

    };
}

#endif