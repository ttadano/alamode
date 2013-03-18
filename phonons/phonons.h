#pragma once

#ifdef _WIN32
#include <mpi.h>
#else
#include "mpi.h"
#endif

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
        class Symmetry *symmetry;
        class Kpoint *kpoint;
        class Integration *integration;
        class Fcs_phonon *fcs_phonon;
        class Dynamical *dynamical;
        class Phonon_velocity *phonon_velocity;
        class Phonon_thermodynamics *phonon_thermodynamics;
        class Relaxation *relaxation;
        class Conductivity *conductivity;
        class Writes *writes;
        class Dos *dos;
        class Gruneisen *gruneisen;
        class MyMPI *mympi;

        PHON(int, char**, MPI_Comm);
        ~PHON();

        void create_pointers();
        void destroy_pointers();

        std::string mode;
		bool restart_flag;
    };
}
