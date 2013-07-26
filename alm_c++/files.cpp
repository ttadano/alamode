#include <iostream>
#include "files.h"
#include "error.h"


using namespace ALM_NS;

Files::Files(ALM *alm): Pointers(alm) {}
Files::~Files() {}

void Files::setfilenames()
{
    file_fcs = job_title + ".fcs";
    file_info = job_title + ".info";
    file_disp_sym = file_disp + ".SYM";
    file_force_sym = file_force + ".SYM";
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
    openfiles();
}
