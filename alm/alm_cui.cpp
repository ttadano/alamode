/*
 alm_cui.cpp

 Copyright (c) 2014, 2015, 2016 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory
 or http://opensource.org/licenses/mit-license.php for information.
*/

#include "alm_cui.h"
#include "alm.h"
#include "input_parser.h"
#include "timer.h"
#include "version.h"
#include "writer.h"
#include <iostream>
#include <iomanip>

#ifdef _OPENMP

#include <omp.h>

#endif

using namespace ALM_NS;

ALMCUI::ALMCUI() {}

ALMCUI::~ALMCUI() {}

void ALMCUI::run(const int narg,
                 char **arg) const
{
    auto alm = new ALM();

    // alm->mode is set herein.
    auto input_parser = new InputParser();
    input_parser->run(alm, narg, arg);
    auto run_mode = input_parser->get_run_mode();
    delete input_parser;

    if (alm->get_verbosity() > 0) {
        std::cout << " +-----------------------------------------------------------------+" << std::endl;
        std::cout << " +                         Program ALM                             +" << std::endl;
        std::cout << " +                             Ver.";
        std::cout << std::setw(7) << ALAMODE_VERSION;
        std::cout << "                         +" << std::endl;
        std::cout << " +-----------------------------------------------------------------+" << std::endl;
        std::cout << std::endl;
#ifdef _OPENMP
        std::cout << " Number of OpenMP threads = "
                  << omp_get_max_threads() << std::endl << std::endl;
#endif

        std::cout << " Job started at " << alm->timer->DateAndTime() << std::endl;
    }

    auto writer = new Writer();

    if (alm->get_verbosity() > 0) {
        writer->write_input_vars(alm, run_mode);
    }

    alm->init_fc_table();

    if (run_mode == "optimize") {
        alm->run_optimize();
        if (alm->get_optimizer_control().linear_model == 1 ||
            (alm->get_optimizer_control().linear_model >= 2
             && alm->get_optimizer_control().cross_validation == 0)) {
            writer->writeall(alm);
        }
    } else if (run_mode == "suggest") {
        alm->run_suggest();
        writer->write_displacement_pattern(alm);
    }
    delete writer;

    if (alm->get_verbosity() > 0) {
        std::cout << std::endl << " Job finished at "
                  << alm->timer->DateAndTime() << std::endl;
    }

    delete alm;
}
