#pragma once

#include "mpi_common.h"
#include "pointers.h"
#include <string>
#include <fstream>

namespace PHON_NS
{
    class Writes: protected Pointers {
    public:

        Writes(class PHON *);
        ~Writes();
        void write_phonon_info();
        void write_selfenergy();
        void write_gruneisen();
		void setup_result_io();
		void finalize_result_io();

        bool writeanime;

	//	MPI_File fp_result;
        double in_kayser(const double);
        int nbands;

	//	char *file_result;
		std::string file_result;
		std::fstream fs_result;

    private:

        void write_phonon_bands();
        void write_phonon_vel();
        void write_phonon_vel_all();
        void write_phonon_dos();
        void write_mode_anime();
        void write_eigenvectors();
        void write_thermodynamics();

        double Ry_to_kayser;

        std::string file_bands, file_dos;
        std::string file_anime, file_evec;
        std::string file_vel;
    
    };
}
