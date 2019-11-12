/*
 dielec.h

 Copyright (c) 2019 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#pragma once

#include "pointers.h"

namespace PHON_NS
{
    class Dielec : protected Pointers
    {
    public:
        Dielec(class PHON *);
        ~Dielec();

        void init();
        void compute_dielectric_constant();
        int calc_dielectric_constant;

    private:

        void set_default_variables();
        void deallocate_variables();
        double ***dielec;
    };
}
