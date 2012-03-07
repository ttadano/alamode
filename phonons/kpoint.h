#pragma once

#include "pointers.h"
#include <string>

namespace PHON_NS {
    class Kpoint: protected Pointers {
    public:
        Kpoint(class PHON *);
        ~Kpoint();

      void kpoint_setups();

      int kpoint_mode;
      unsigned int nkmax;
      unsigned int nkx, nky, nkz;
      
      unsigned int npath, nk;

      double **xk;

    private:
        void gen_kpoints_band();
        void gen_kmesh();
        std::string **kp_symbol;
        double ***kp_bound;
        unsigned int *nkp;
    };
}