/*
 files.h

 Copyright (c) 2014 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#pragma once

#include <string>
#include <fstream>
#include "pointers.h"

namespace ALM_NS
{
    class Files : protected Pointers
    {
    public:
        Files(class ALM *);
        ~Files();

        void init();

        std::string job_title;
        std::string file_fcs, file_info;

        std::string file_disp, file_force;

        std::string file_int;
        std::string *file_disp_pattern;

        std::ofstream ofs_int;
        std::ifstream ifs_int;

        std::ifstream ifs_disp, ifs_force;

        std::ofstream *ofs_disp_pattern;

    private:
        void openfiles();
        void closefiles();
        void setfilenames();
    };
}

