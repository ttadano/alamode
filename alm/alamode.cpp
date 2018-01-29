/*
 alamode.cpp

 Copyright (c) 2014, 2015, 2016 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#include <iostream>
#include <iomanip>
#include "interaction.h"
#include "symmetry.h"
#include "input.h"
#include "system.h"
#include "files.h"
#include "memory.h"
#include "timer.h"
#include "fcs.h"
#include "fitting.h"
#include "constraint.h"
#include "timer.h"
#include "writes.h"
#include "patterndisp.h"
#include "version.h"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace ALM_NS;

ALM::ALM(int narg, char **arg)
{
    std::cout << " +-----------------------------------------------------------------+" << std::endl;
    std::cout << " +                         Program ALM                             +" << std::endl;
    std::cout << " +                             Ver.";
    std::cout << std::setw(7) << ALAMODE_VERSION;
    std::cout << "                         +" << std::endl;
    std::cout << " +-----------------------------------------------------------------+" << std::endl;
    std::cout << std::endl;

    timer = new Timer(this);

#ifdef _OPENMP
    std::cout << " Number of OpenMP threads = " 
        << omp_get_max_threads() << std::endl << std::endl;
#endif
    std::cout << " Job started at " << timer->DateAndTime() << std::endl;

    input = new Input(this, narg, arg);
    create();
    input->parse_input(narg, arg);
    writes->write_input_vars();
    initialize();

    if (mode == "fitting") {

        fcs->init();
        constraint->setup();
        fitting->fitmain();
        writes->writeall();

    } else if (mode == "suggest") {

        displace->gen_displacement_pattern();
        writes->write_displacement_pattern();

    }

    finalize();

    std::cout << std::endl << " Job finished at "
        << timer->DateAndTime() << std::endl;
}

void ALM::create()
{
    memory = new Memory(this);
    files = new Files(this);
    system = new System(this);
    interaction = new Interaction(this);
    fcs = new Fcs(this);
    symmetry = new Symmetry(this);
    fitting = new Fitting(this);
    constraint = new Constraint(this);
    displace = new Displace(this);
    writes = new Writes(this);
}

void ALM::initialize()
{
    system->init();
    files->init();
    symmetry->init();
    interaction->init();
}

ALM::~ALM()
{
    delete input;
    delete timer;
}

void ALM::finalize()
{
    delete files;
    delete interaction;
    delete fcs;
    delete symmetry;
    delete system;
    delete fitting;
    delete constraint;
    delete displace;
    delete writes;
    delete memory;
}
