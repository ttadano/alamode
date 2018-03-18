/*
 phonons.h

 Copyright (c) 2014 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#pragma once

#ifdef _WIN32
#include <mpi.h>
#else
#include "mpi.h"
#endif

#include <string>

namespace PHON_NS
{
    class PHON
    {
    public:
        PHON(int,
             char **,
             MPI_Comm);

        virtual ~PHON();

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
        class Thermodynamics *thermodynamics;
        class AnharmonicCore *anharmonic_core;
        class ModeAnalysis *mode_analysis;
        class Selfenergy *selfenergy;
        class Conductivity *conductivity;
        class Writes *writes;
        class Dos *dos;
        class Gruneisen *gruneisen;
        class MyMPI *mympi;
        class Isotope *isotope;
        class Scph *scph;
        class Ewald *ewald;

        void create_pointers();
        void destroy_pointers();

        std::string mode;
        bool restart_flag;

        void execute_phonons();
        void execute_RTA();
        void execute_self_consistent_phonon();
        void setup_base();
    };
}
