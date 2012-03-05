#pragma once

#include "pointers.h"
#include <string>

namespace PHON_NS {
    class Kpoint: protected Pointers {
    public:
        Kpoint(class PHON *);
        ~Kpoint();

      void read_kpoints();

      int kpoint_mode;
      unsigned int nkmax;
      unsigned int nkx, nky, nkz;
      
      unsigned int npath, nk;

      double **xk;

    private:
        std::string **kp_symbol;
        double ***kp_bound;
        unsigned int *nkp;
    };
}