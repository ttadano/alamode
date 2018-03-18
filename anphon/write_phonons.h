/*
 write_phonons.h

 Copyright (c) 2014 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#pragma once

#include "mpi_common.h"
#include "pointers.h"
#include <string>
#include <fstream>

namespace PHON_NS
{
    class Writes : protected Pointers
    {
    public:

        Writes(class PHON *);

        ~Writes();

        void write_phonon_info();
        void print_phonon_energy();
        void write_gruneisen();
        void setup_result_io();
        void write_input_vars();
        void write_kappa();
        void write_selfenergy_isotope();

        bool print_xsf;
        bool print_msd;
        bool print_anime;

        unsigned int anime_cellsize[3];
        double anime_kpoint[3];

        double in_kayser(const double);

        int nbands;

        std::string file_result;
        std::string anime_format;
        std::fstream fs_result;

    private:

        void write_phonon_bands();
        void write_phonon_vel();
        void write_phonon_vel_all();
        void write_phonon_dos();
        void write_two_phonon_dos();
        void write_scattering_phase_space();
        void write_scattering_amplitude();
        void write_normal_mode_direction();

        void write_normal_mode_animation(const double [3],
                                         const unsigned int [3]);

        void write_eigenvectors();
        void write_thermodynamics();
        void write_msd();
        void write_participation_ratio();

        double Ry_to_kayser;
    };
}
