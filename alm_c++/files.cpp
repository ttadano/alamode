#include "files.h"
#include <iostream>

using namespace ALM_NS;

Files::Files(ALM *alm): Pointers(alm) {}
Files::~Files() {}

void Files::setfilenames()
{
    file_log = job_title + ".alm";
    file_fcs = job_title + ".fcs";
    file_premd = job_title + ".premd";
}

void Files::openfiles()
{
    ofs_log.open(file_log, std::ios::out);
    ofs_fcs.open(file_fcs, std::ios::out);
    ofs_premd.open(file_premd, std::ios::out);
   
    ifs_disp.open(file_disp);
    ifs_force.open(file_force);
}

void Files::closefiles()
{
    ofs_log.close();
    ofs_fcs.close();
    ofs_premd.close();

    ifs_disp.close();
    ifs_force.close();
}

void Files::init()
{
    setfilenames();
    openfiles();

    std::cout << "Input Filenames" << std::endl;
    std::cout << "  Displacement: " << file_disp << std::endl;
    std::cout << "  Force       : " << file_force << std::endl;
}