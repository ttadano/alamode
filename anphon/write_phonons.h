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

namespace PHON_NS {
    class Writes : protected Pointers {
    public:

        Writes(class PHON *);

        ~Writes();

        void write_phonon_info();

        void print_phonon_energy() const;

        void write_gruneisen();

        void setup_result_io();

        void write_input_vars();

        void write_kappa() const;

        void write_selfenergy_isotope() const;

        double in_kayser(const double) const;

        void setWriteOptions(const bool print_msd_,
                              const bool print_xsf_,
                              const bool print_anime_,
                              const std::string &anime_format_,
                              const int anime_steps_,
                              const unsigned int anime_cellsize_[3],
                              const double anime_kpoint_[3],
                              const bool print_ucorr_,
                              const int shift_ucorr_[3]);

        void write_scph_energy(double ***,
                               const int bubble = 0) const;

        void write_scph_bands(double ***,
                              const int bubble = 0) const;

        void write_scph_dos(double **,
                            const int bubble = 0) const;

        void write_scph_thermodynamics(double *heat_capacity,
                                       double *FE_QHA,
                                       double *dFE_scph) const;

        void write_scph_msd(double **, const int bubble = 0) const;

        void write_scph_ucorr(double ***ucorr_scph, const int bubble = 0) const;

        void write_scph_dielec(double ****dielec_scph) const;

        unsigned int getVerbosity() const;
        void setVerbosity(unsigned int verbosity_in);

        bool getPrintMSD() const;
        bool getPrintUcorr() const;
        std::array<int, 3> getShiftUcorr() const;

        std::fstream fs_result;
        std::string file_result;
        int nbands;

    private:

        void write_phonon_bands() const;

        void write_phonon_vel() const;

        void write_phonon_vel_all() const;

        void write_phonon_dos() const;

        void write_two_phonon_dos() const;

        void write_scattering_phase_space() const;

        void write_scattering_amplitude() const;

        void write_normal_mode_direction() const;

        void write_normal_mode_animation(const double [3],
                                         const unsigned int [3]) const;

        void write_eigenvectors() const;

#ifdef _HDF5
        void write_eigenvectors_HDF5() const;
#endif

        void write_thermodynamics() const;

        void write_msd() const;

        void write_disp_correlation() const;

        void write_participation_ratio() const;

        void write_dielectric_function() const;

        double Ry_to_kayser;
        unsigned int verbosity;

        int anime_frames;

        bool print_xsf;
        bool print_msd;
        bool print_ucorr;
        bool print_anime;

        unsigned int anime_cellsize[3];
        double anime_kpoint[3];
        int shift_ucorr[3];

        std::string anime_format;


    public:

    };
}
