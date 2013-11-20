#pragma once

#include "pointers.h"
#include <string>
#include <fstream>

namespace ALM_NS{
    class Writes: protected Pointers{
    public:
        Writes(class ALM *);
        ~Writes();
        
        void writeall();
		void write_input_vars();
		void write_displacement_pattern();

    private:
        void wrtfcs();
        void wrtmisc();

      std::ofstream ofs_info;

    };
}
