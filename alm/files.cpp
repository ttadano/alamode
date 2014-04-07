/*
 files.cpp

 Copyright (c) 2014 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#include <iostream>
#include "files.h"
#include "error.h"
#include "memory.h"
#include "interaction.h"
#include <boost/lexical_cast.hpp>

using namespace ALM_NS;

Files::Files(ALM *alm): Pointers(alm) {}
Files::~Files() {
    if (alm->mode == "fitting") {
        closefiles();
    }
}

void Files::setfilenames()
{
    int i;

    file_fcs = job_title + ".fcs";
    file_info = job_title + ".info";

    if (alm->mode == "suggest") {

        memory->allocate(file_disp_pattern, interaction->maxorder);

        for (i = 0; i < interaction->maxorder; ++i) {
            if (i == 0) {
                file_disp_pattern[i] = job_title + ".pattern_HARMONIC";
            } else {
                file_disp_pattern[i] = job_title + ".pattern_ANHARM" + boost::lexical_cast<std::string>(i+2);
            }
        }
    }
}

void Files::openfiles()
{
    ifs_disp.open(file_disp.c_str(), std::ios::in);
    if(!ifs_disp) error->exit("openfiles", "cannot open disp file");
    ifs_force.open(file_force.c_str(), std::ios::in);
    if(!ifs_force) error->exit("openfiles", "cannot open force file");
}

void Files::closefiles()
{
    ifs_disp.close();
    ifs_force.close();
}

void Files::init()
{
    setfilenames();

    if (alm->mode == "fitting") openfiles();
}
