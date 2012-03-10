#pragma once

#include "pointers.h"
#include <string>

namespace PHON_NS
{
    class Writes: protected Pointers {
    public:

        Writes(class PHON *);
        ~Writes();
        void write_phonon_info();
        bool writeanime;

        double in_kayser(const double);

    private:

        void write_phonon_bands();
        void write_phonon_dos();
        void write_mode_anime();
	void write_eigenvectors();

        double Ry_to_kayser;

        std::string file_bands, file_dos;
        std::string file_anime, file_evec;
    
    };
}
