#include <iostream>
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
#include "ewald.h"
#include "patterndisp.h"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace ALM_NS;

ALM::ALM(int narg, char **arg)
{
	timer = new Timer(this);

#ifdef _OPENMP
	std::cout << "Number of OpenMP threads = " << omp_get_max_threads() << std::endl << std::endl;
#endif
	std::cout << "Job started at " << timer->DateAndTime() << std::endl;

	input = new Input(this, narg, arg);
	create();
	input->parce_input();
	writes->write_input_vars();
	initialize();

	if (mode == "fitting") {

		constraint->setup();
		fitting->fitmain();
		writes->writeall();

	} else if (mode == "suggest") {

		displace->gen_displacement_pattern();
		writes->write_displacement_pattern();

	}

	finalize();

	std::cout << std::endl << "Job finished at " << timer->DateAndTime() << std::endl;
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
	ewald = new Ewald(this);
	displace = new Displace(this);
	writes = new Writes(this);
}

void ALM::initialize()
{
	system->init();
	files->init();
	symmetry->init();
	interaction->init();
	fcs->init();
	ewald->init();
}
ALM::~ALM()
{
	finalize();
	delete input;
	delete timer;
}

void ALM::finalize()
{
	delete memory;
	delete files;
	delete interaction;
	delete fcs;
	delete symmetry;
	delete system;
	delete fitting;
	delete constraint;
	delete ewald;
	delete displace;
	delete writes;
}
