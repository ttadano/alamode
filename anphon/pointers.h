/*
 pointers.h

 Copyright (c) 2014 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#pragma once

#include "phonons.h"

namespace PHON_NS
{
    class Pointers
    {
    public:
        Pointers(PHON *ptr) :
            phon(ptr),
            memory(ptr->memory),
            error(ptr->error),
            input(ptr->input),
            system(ptr->system),
            symmetry(ptr->symmetry),
            kpoint(ptr->kpoint),
            integration(ptr->integration),
            fcs_phonon(ptr->fcs_phonon),
            dynamical(ptr->dynamical),
            phonon_velocity(ptr->phonon_velocity),
            thermodynamics(ptr->thermodynamics),
            relaxation(ptr->relaxation),
            selfenergy(ptr->selfenergy),
            conductivity(ptr->conductivity),
            writes(ptr->writes),
            dos(ptr->dos),
            gruneisen(ptr->gruneisen),
            mympi(ptr->mympi),
            isotope(ptr->isotope),
            scph(ptr->scph),
            ewald(ptr->ewald),
            timer(ptr->timer) {}

        virtual ~Pointers() {}

    protected:
        PHON *phon;
        Memory *&memory;
        Error *&error;
        Input *&input;
        System *&system;
        Symmetry *&symmetry;
        Kpoint *&kpoint;
        Integration *&integration;
        Fcs_phonon *&fcs_phonon;
        Dynamical *&dynamical;
        Phonon_velocity *&phonon_velocity;
        Thermodynamics *&thermodynamics;
        Relaxation *&relaxation;
        Selfenergy *&selfenergy;
        Conductivity *&conductivity;
        Writes *&writes;
        Dos *&dos;
        Gruneisen *&gruneisen;
        MyMPI *&mympi;
        Isotope *&isotope;
        Scph *&scph;
        Ewald *&ewald;
        Timer *&timer;
    };
}

