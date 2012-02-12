#ifndef PHON_PHONONS_HEADER
#define PHON_PHONONS_HEADER

namespace PHON_NS
{
    class PHON {
    public:
        class Memory *memory;
        class Error *error;
        class Timer *timer;
        class Input *input;

        PHON(int, char**);
        ~PHON();

        void create_pointers();
        void destroy_pointers();
    };
}
#endif