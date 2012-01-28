#include "files.h"
#include "error.h"
#include <iostream>

using namespace ALM_NS;

Files::Files(ALM *alm): Pointers(alm) {}
Files::~Files() {}

void Files::setfilenames()
{
    file_fcs = job_title + ".fcs";
    file_premd = job_title + ".premd";
    file_disp_sym = file_disp + ".SYM";
    file_force_sym = file_force + ".SYM";
    file_int = job_title + ".pairs";
}

void Files::openfiles()
{
    ofs_fcs.open(file_fcs, std::ios::out);
    if(!ofs_fcs) error->exit("openfiles", "cannot open fcs file");
    ofs_premd.open(file_premd, std::ios::out);
    if(!ofs_premd) error->exit("openfiles", "cannot open premd file");
   
    ifs_disp.open(file_disp, std::ios::in);
    if(!ifs_disp) error->exit("openfiles", "cannot open disp file");
    ifs_force.open(file_force, std::ios::in);
    if(!ifs_force) error->exit("openfiles", "cannot open force file");
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