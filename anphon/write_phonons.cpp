/*
write_phonons.cpp

Copyright (c) 2014, 2015, 2016 Terumasa Tadano

This file is distributed under the terms of the MIT license.
Please see the file 'LICENCE.txt' in the root directory 
or http://opensource.org/licenses/mit-license.php for information.
*/

#include "mpi_common.h"
#include <iomanip>
#include <sys/stat.h>
#include "constants.h"
#include "conductivity.h"
#include "dielec.h"
#include "dynamical.h"
#include "error.h"
#include "fcs_phonon.h"
#include "gruneisen.h"
#include "kpoint.h"
#include "memory.h"
#include "parsephon.h"
#include "phonon_dos.h"
#include "thermodynamics.h"
#include "phonon_velocity.h"
#include "anharmonic_core.h"
#include "symmetry_core.h"
#include "system.h"
#include "write_phonons.h"
#include "mathfunctions.h"
#include "isotope.h"
#include "integration.h"
#include "scph.h"

#ifdef _HDF5

#include "H5Cpp.h"

#endif

using namespace PHON_NS;

Writes::Writes(PHON *phon) : Pointers(phon)
{
    Ry_to_kayser = Hz_to_kayser / time_ry;

    print_ucorr = false;
    print_xsf = false;
    print_anime = false;
    print_msd = false;
    print_zmode = false;
    anime_cellsize[0] = 0;
    anime_cellsize[1] = 0;
    anime_cellsize[2] = 0;
    shift_ucorr[0] = 0;
    shift_ucorr[1] = 0;
    shift_ucorr[2] = 0;
    anime_kpoint[0] = 0.0;
    anime_kpoint[1] = 0.0;
    anime_kpoint[2] = 0.0;
    anime_frames = 20;
    anime_format = "xyz";
    verbosity = 1;
};

Writes::~Writes() {};

void Writes::write_input_vars()
{
    if (verbosity == 0) return;

    unsigned int i;

    std::cout << std::endl;
    std::cout << " Input variables:" << std::endl;
    std::cout << " -----------------------------------------------------------------" << std::endl;
    std::cout << " General:" << std::endl;
    std::cout << "  PREFIX = " << input->job_title << std::endl;
    std::cout << "  MODE = " << phon->mode << std::endl;
    std::cout << "  FCSXML = " << fcs_phonon->file_fcs << std::endl;
    if (fcs_phonon->update_fc2) {
        std::cout << "  FC2XML = " << fcs_phonon->file_fc2 << std::endl;
    }
    std::cout << std::endl;

    std::cout << "  NKD = " << system->nkd << "; KD = ";
    for (i = 0; i < system->nkd; ++i) {
        std::cout << std::setw(4) << system->symbol_kd[i];
    }
    std::cout << std::endl;
    std::cout << "  MASS = ";
    if (system->mass_kd) {
        for (i = 0; i < system->nkd; ++i) {
            std::cout << std::setw(10) << system->mass_kd[i];
        }
    }
    std::cout << std::endl;
    std::cout << "  NSYM = " << symmetry->nsym << "; TOLERANCE = " << symmetry->tolerance;
    std::cout << "; PRINTSYM = " << symmetry->printsymmetry << std::endl;
    // std::cout << "  TREVSYM = " << symmetry->time_reversal_sym << std::endl;

    std::cout << std::endl;

    std::cout << "  NONANALYTIC = " << dynamical->nonanalytic << std::endl;
    if (dynamical->nonanalytic) {
        std::cout << "  BORNINFO = " << dynamical->file_born
                  << "; NA_SIGMA = " << dynamical->na_sigma << std::endl;
    }
    std::cout << std::endl;
    if (writes->nbands >= 0) {
        std::cout << "  NBANDS = " << writes->nbands << std::endl;
    }

    std::cout << "  TMIN = " << system->Tmin
              << "; TMAX = " << system->Tmax
              << "; DT = " << system->dT << std::endl;
    std::cout << "  EMIN = " << dos->emin
              << "; EMAX = " << dos->emax
              << "; DELTA_E = " << dos->delta_e << std::endl;
    std::cout << std::endl;

    std::cout << "  ISMEAR = " << integration->ismear
              << "; EPSILON = " << integration->epsilon << std::endl;
    std::cout << std::endl;
    std::cout << "  CLASSICAL = " << thermodynamics->classical << std::endl;
    std::cout << "  BCONNECT = " << dynamical->band_connection << std::endl;
    std::cout << std::endl;

    if (phon->mode == "RTA") {
        std::cout << "  RESTART = " << conductivity->get_restart_conductivity(3) << std::endl;
        std::cout << "  TRISYM = " << anharmonic_core->use_triplet_symmetry << std::endl;
        std::cout << std::endl;
    } else if (phon->mode == "SCPH") {
        std::cout << " Scph:" << std::endl;
        std::cout << "  KMESH_INTERPOLATE = ";
        for (i = 0; i < 3; ++i) std::cout << std::setw(5) << scph->kmesh_interpolate[i];
        std::cout << std::endl;
        std::cout << "  KMESH_SCPH        = ";
        for (i = 0; i < 3; ++i) std::cout << std::setw(5) << scph->kmesh_scph[i];
        std::cout << std::endl;
        std::cout << "  SELF_OFFDIAG = " << scph->selfenergy_offdiagonal << std::endl;
        std::cout << "  IALGO = " << scph->ialgo << std::endl << std::endl;
        std::cout << "  RESTART_SCPH = " << scph->restart_scph << std::endl;
        std::cout << "  LOWER_TEMP = " << scph->lower_temp << std::endl;
        std::cout << "  WARMSTART = " << scph->warmstart_scph << std::endl << std::endl;
        std::cout << "  TOL_SCPH = " << scph->tolerance_scph << std::endl;
        std::cout << "  MAXITER = " << scph->maxiter << std::endl;
        std::cout << "  MIXALPHA = " << scph->mixalpha << std::endl;
    }
    std::cout << std::endl;

    std::cout << " Kpoint:" << std::endl;
    std::cout << "  KPMODE (1st entry for &kpoint) = "
              << kpoint->kpoint_mode << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;

    std::cout << " Analysis:" << std::endl;
    if (phon->mode == "PHONONS") {
        std::cout << "  PRINTVEL = " << phonon_velocity->print_velocity << std::endl;
        std::cout << "  PRINTVEC = " << dynamical->print_eigenvectors << std::endl;
        std::cout << "  PRINTXSF = " << writes->print_xsf << std::endl;
        std::cout << std::endl;

        if (print_anime) {
            std::cout << "  ANIME = ";
            for (i = 0; i < 3; ++i) std::cout << std::setw(5) << anime_kpoint[i];
            std::cout << std::endl;
            std::cout << "  ANIME_CELL = ";
            for (i = 0; i < 3; ++i) std::cout << std::setw(5) << anime_cellsize[i];
            std::cout << std::endl;
            std::cout << "  ANIME_FORMAT = " << anime_format << std::endl;
            std::cout << std::endl;
        }

        if (kpoint->kpoint_mode == 2) {
            std::cout << "  PDOS = " << dos->projected_dos
                      << "; TDOS = " << dos->two_phonon_dos << std::endl;
            std::cout << "  PRINTMSD = " << print_msd << std::endl;
            std::cout << "  SPS = " << dos->scattering_phase_space << std::endl;
            std::cout << std::endl;
        }
        std::cout << "  GRUNEISEN = " << gruneisen->print_gruneisen << std::endl;
        std::cout << "  NEWFCS = " << gruneisen->print_newfcs;
        if (gruneisen->print_newfcs) {
            std::cout << "; DELTA_A = " << gruneisen->delta_a << std::endl;
            std::cout << "  QUARTIC = " << anharmonic_core->quartic_mode;
        }
        std::cout << std::endl;

    } else if (phon->mode == "RTA") {
        std::cout << "  ISOTOPE = " << isotope->include_isotope << std::endl;
        if (isotope->include_isotope) {
            std::cout << "  ISOFACT = ";
            if (isotope->isotope_factor) {
                for (i = 0; i < system->nkd; ++i) {
                    std::cout << std::scientific
                              << std::setw(13) << isotope->isotope_factor[i];
                }
            }
            std::cout << std::endl;
        }

        std::cout << "  KAPPA_SPEC = " << conductivity->calc_kappa_spec << std::endl;

        //        std::cout << "  KS_INPUT = " << anharmonic_core->ks_input << std::endl;
        //        std::cout << "  QUARTIC = " << anharmonic_core->quartic_mode << std::endl;
        // std::cout << "  REALPART = " << anharmonic_core->calc_realpart << std::endl;
        // std::cout << "  ATOMPROJ = " << anharmonic_core->atom_project_mode << std::endl;
        // std::cout << "  FSTATE_W = " << anharmonic_core->calc_fstate_omega << std::endl;
        //  std::cout << "  FSTATE_K = " << anharmonic_core->calc_fstate_k << std::endl;

    } else if (phon->mode == "SCPH") {
        // Do nothing
    } else {
        exit("write_input_vars", "This cannot happen");
    }

    std::cout << std::endl << std::endl;
    std::cout << " -----------------------------------------------------------------" << std::endl;
    std::cout << std::endl;
}

void Writes::setWriteOptions(const bool print_msd_,
                             const bool print_xsf_,
                             const bool print_anime_,
                             const std::string &anime_format_,
                             const int anime_frames_,
                             const unsigned int anime_cellsize_[3],
                             const double anime_kpoint_[3],
                             const bool print_ucorr_,
                             const int shift_ucorr_[3],
                             const bool print_zmode_)
{
    print_msd = print_msd_;
    print_xsf = print_xsf_;
    print_anime = print_anime_;
    anime_format = anime_format_;
    anime_frames = anime_frames_;
    print_ucorr = print_ucorr_;
    print_zmode = print_zmode_;

    for (auto i = 0; i < 3; ++i) {
        anime_cellsize[i] = anime_cellsize_[i];
        anime_kpoint[i] = anime_kpoint_[i];
        shift_ucorr[i] = shift_ucorr_[i];
    }
}

bool Writes::getPrintMSD() const
{
    return print_msd;
}

bool Writes::getPrintUcorr() const
{
    return print_ucorr;
}

std::array<int, 3> Writes::getShiftUcorr() const
{
    return {shift_ucorr[0], shift_ucorr[1], shift_ucorr[2]};
}

void Writes::print_phonon_energy() const
{
    unsigned int i;
    unsigned int ik, is;
    const auto ns = dynamical->neval;

    const auto kayser_to_THz = 0.0299792458;

    std::cout << std::endl;
    std::cout << " -----------------------------------------------------------------" << std::endl << std::endl;
    std::cout << " Phonon frequencies below:" << std::endl << std::endl;

    if (kpoint->kpoint_mode == 0) {

        auto nk_now = kpoint->kpoint_general->nk;
        auto xk_now = kpoint->kpoint_general->xk;
        auto eval_now = dynamical->dymat_general->get_eigenvalues();

        for (ik = 0; ik < nk_now; ++ik) {
            std::cout << " # k point " << std::setw(5) << ik + 1;
            std::cout << " : (";

            for (i = 0; i < 3; ++i) {
                std::cout << std::fixed << std::setprecision(4)
                          << std::setw(8) << xk_now[ik][i];
                if (i < 2) std::cout << ",";
            }
            std::cout << ")" << std::endl;

            std::cout << "   Mode, Frequency " << std::endl;

            for (is = 0; is < ns; ++is) {
                std::cout << std::setw(7) << is + 1;
                std::cout << std::fixed << std::setprecision(4) << std::setw(12)
                          << in_kayser(eval_now[ik][is]);
                std::cout << " cm^-1  (";
                std::cout << std::fixed << std::setprecision(4) << std::setw(12)
                          << kayser_to_THz * in_kayser(eval_now[ik][is]);
                std::cout << " THz )" << std::endl;
            }
            std::cout << std::endl;
        }

    } else if (kpoint->kpoint_bs) {

        auto nk = kpoint->kpoint_bs->nk;

        for (ik = 0; ik < nk; ++ik) {
            std::cout << " # k point " << std::setw(5) << ik + 1;
            std::cout << " : (";

            for (i = 0; i < 3; ++i) {
                std::cout << std::fixed << std::setprecision(4)
                          << std::setw(8) << kpoint->kpoint_bs->xk[ik][i];
                if (i < 2) std::cout << ",";
            }
            std::cout << ")" << std::endl;

            std::cout << "   Mode, Frequency " << std::endl;

            for (is = 0; is < ns; ++is) {
                std::cout << std::setw(7) << is + 1;
                std::cout << std::fixed << std::setprecision(4) << std::setw(12)
                          << in_kayser(dynamical->dymat_band->get_eigenvalues()[ik][is]);
                std::cout << " cm^-1  (";
                std::cout << std::fixed << std::setprecision(4) << std::setw(12)
                          << kayser_to_THz * in_kayser(dynamical->dymat_band->get_eigenvalues()[ik][is]);
                std::cout << " THz )" << std::endl;
            }
            std::cout << std::endl;
        }

    } else if (kpoint->kpoint_mode == 2) {

        for (ik = 0; ik < dos->kmesh_dos->kpoint_irred_all.size(); ++ik) {

            std::cout << " # Irred. k point" << std::setw(5) << ik + 1;
            std::cout << " : (";

            for (i = 0; i < 3; ++i) {
                std::cout << std::fixed << std::setprecision(4) << std::setw(8)
                          << dos->kmesh_dos->kpoint_irred_all[ik][0].kval[i];
                if (i < 2) std::cout << ",";
            }
            std::cout << ")" << std::endl;

            std::cout << "   Mode, Frequency " << std::endl;

            const auto knum = dos->kmesh_dos->kpoint_irred_all[ik][0].knum;

            for (is = 0; is < ns; ++is) {
                std::cout << std::setw(7) << is + 1;
                std::cout << std::fixed << std::setprecision(4) << std::setw(12)
                          << in_kayser(dos->dymat_dos->get_eigenvalues()[knum][is]);
                std::cout << " cm^-1  (";
                std::cout << std::fixed << std::setprecision(4) << std::setw(12)
                          << kayser_to_THz * in_kayser(dos->dymat_dos->get_eigenvalues()[knum][is]);
                std::cout << " THz )" << std::endl;
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
}

void Writes::write_phonon_info()
{
    if (nbands < 0) {
        nbands = 3 * system->natmin;
    }

    if (print_anime) {
        write_normal_mode_animation(anime_kpoint, anime_cellsize);
    }

    std::cout << std::endl;
    std::cout << " -----------------------------------------------------------------" << std::endl << std::endl;
    std::cout << " The following files are created: " << std::endl;

    if (kpoint->kpoint_mode == 1) {
        write_phonon_bands();
    }

    if (phonon_velocity->print_velocity) {
        if (kpoint->kpoint_bs) {
            write_phonon_vel();
        }
        if (dos->kmesh_dos) {
            write_phonon_vel_all();
        }
    }

    if (dos->flag_dos) {

        if (dos->compute_dos || dos->projected_dos) {
            write_phonon_dos();
        }

        if (dos->two_phonon_dos) {
            write_two_phonon_dos();
        }

        if (dos->scattering_phase_space == 1) {
            write_scattering_phase_space();
        } else if (dos->scattering_phase_space == 2) {
            write_scattering_amplitude();
        }

        write_thermodynamics();
        if (print_msd) write_msd();
        if (print_ucorr) write_disp_correlation();
    }

    if (print_xsf) {
        write_normal_mode_direction();
    }

    if (dynamical->print_eigenvectors) {
        write_eigenvectors();
#ifdef _HDF5
        write_eigenvectors_HDF5();
#endif
    }

    if (dynamical->participation_ratio) {
        write_participation_ratio();
    }

    if (gruneisen->print_gruneisen) {
        write_gruneisen();
    }

    if (dielec->calc_dielectric_constant) {
        write_dielectric_function();
    }

    if (print_anime) {
        if (anime_format == "XSF" || anime_format == "AXSF") {
            std::cout << "  " << std::setw(input->job_title.length() + 12) << std::left
                      << input->job_title + ".anime*.axsf";
            std::cout << " : AXSF files for animate phonon modes" << std::endl;
        } else if (anime_format == "XYZ") {
            std::cout << "  " << std::setw(input->job_title.length() + 12) << std::left
                      << input->job_title + ".anime*.xyz";
            std::cout << " : XYZ files for animate phonon modes" << std::endl;
        }
    }

    if (print_zmode) {
        print_normalmode_borncharge();
    }
}

void Writes::write_phonon_bands() const
{
    std::ofstream ofs_bands;
    auto file_bands = input->job_title + ".bands";

    ofs_bands.open(file_bands.c_str(), std::ios::out);
    if (!ofs_bands)
        exit("write_phonon_bands",
             "cannot open file_bands");

    unsigned int i, j;
    const auto nk = kpoint->kpoint_bs->nk;
    const auto kaxis = kpoint->kpoint_bs->kaxis;
    const auto eval = dynamical->dymat_band->get_eigenvalues();

    auto kcount = 0;

    std::string str_tmp = "NONE";
    std::string str_kpath;
    std::string str_kval;

    for (i = 0; i < kpoint->kpInp.size(); ++i) {
        if (str_tmp != kpoint->kpInp[i].kpelem[0]) {
            str_tmp = kpoint->kpInp[i].kpelem[0];
            str_kpath += " " + str_tmp;

            std::ostringstream ss;
            ss << std::fixed << std::setprecision(6) << kaxis[kcount];
            str_kval += " " + ss.str();
        }
        kcount += std::atoi(kpoint->kpInp[i].kpelem[8].c_str());

        if (str_tmp != kpoint->kpInp[i].kpelem[4]) {
            str_tmp = kpoint->kpInp[i].kpelem[4];
            str_kpath += " " + str_tmp;

            std::ostringstream ss;
            ss << std::fixed << std::setprecision(6) << kaxis[kcount - 1];
            str_kval += " " + ss.str();
        }
    }

    ofs_bands << "# " << str_kpath << std::endl;
    ofs_bands << "#" << str_kval << std::endl;
    ofs_bands << "# k-axis, Eigenvalues [cm^-1]" << std::endl;

    if (dynamical->band_connection == 0) {
        for (i = 0; i < nk; ++i) {
            ofs_bands << std::setw(8) << std::fixed << kaxis[i];
            for (j = 0; j < nbands; ++j) {
                ofs_bands << std::setw(15) << std::scientific << in_kayser(eval[i][j]);
            }
            ofs_bands << std::endl;
        }
    } else {
        for (i = 0; i < nk; ++i) {
            ofs_bands << std::setw(8) << std::fixed << kaxis[i];
            for (j = 0; j < nbands; ++j) {
                ofs_bands << std::setw(15) << std::scientific
                          << in_kayser(eval[i][dynamical->index_bconnect[i][j]]);
            }
            ofs_bands << std::endl;
        }
    }

    ofs_bands.close();

    std::cout << "  " << std::setw(input->job_title.length() + 12) << std::left << file_bands;
    std::cout << " : Phonon band structure" << std::endl;

    if (dynamical->band_connection == 2) {
        std::ofstream ofs_connect;
        auto file_connect = input->job_title + ".connection";

        ofs_connect.open(file_connect.c_str(), std::ios::out);
        if (!ofs_connect)
            exit("write_phonon_bands",
                 "cannot open file_connect");

        ofs_connect << "# " << str_kpath << std::endl;
        ofs_connect << "#" << str_kval << std::endl;
        ofs_connect << "# k-axis, mapping" << std::endl;

        for (i = 0; i < nk; ++i) {
            ofs_connect << std::setw(8) << std::fixed << kaxis[i];
            for (j = 0; j < nbands; ++j) {
                ofs_connect << std::setw(5) << dynamical->index_bconnect[i][j] + 1;
            }
            ofs_connect << std::endl;
        }
        ofs_connect.close();
        std::cout << "  " << std::setw(input->job_title.length() + 12) << std::left << file_connect;
        std::cout << " : Connectivity map information of band dispersion" << std::endl;
    }
}

void Writes::write_phonon_vel() const
{
    std::ofstream ofs_vel;
    auto file_vel = input->job_title + ".phvel";

    ofs_vel.open(file_vel.c_str(), std::ios::out);
    if (!ofs_vel) exit("write_phonon_vel", "cannot open file_vel");

    const auto nk = kpoint->kpoint_bs->nk;
    const auto kaxis = kpoint->kpoint_bs->kaxis;
    const auto Ry_to_SI_vel = Bohr_in_Angstrom * 1.0e-10 / time_ry;

    double **phvel_bs;
    allocate(phvel_bs, nk, dynamical->neval);

    phonon_velocity->get_phonon_group_velocity_bandstructure(kpoint->kpoint_bs,
                                                             system->lavec_p,
                                                             system->rlavec_p,
                                                             fcs_phonon->fc2_ext,
                                                             phvel_bs);

    ofs_vel << "# k-axis, |Velocity| [m / sec]" << std::endl;
    ofs_vel.setf(std::ios::fixed);

    if (dynamical->band_connection == 0) {
        for (auto i = 0; i < nk; ++i) {
            ofs_vel << std::setw(8) << kaxis[i];
            for (auto j = 0; j < nbands; ++j) {
                ofs_vel << std::setw(15)
                        << std::abs(phvel_bs[i][j] * Ry_to_SI_vel);
            }
            ofs_vel << std::endl;
        }
    } else {
        for (auto i = 0; i < nk; ++i) {
            ofs_vel << std::setw(8) << kaxis[i];
            for (auto j = 0; j < nbands; ++j) {
                ofs_vel << std::setw(15)
                        << std::abs(phvel_bs[i][dynamical->index_bconnect[i][j]] * Ry_to_SI_vel);
            }
            ofs_vel << std::endl;
        }
    }

    ofs_vel.close();

    std::cout << "  " << std::setw(input->job_title.length() + 12) << std::left << file_vel;
    std::cout << " : Phonon velocity along given k path" << std::endl;

    deallocate(phvel_bs);
}

void Writes::write_phonon_vel_all() const
{
    std::ofstream ofs_vel;
    auto file_vel = input->job_title + ".phvel_all";

    ofs_vel.open(file_vel.c_str(), std::ios::out);
    if (!ofs_vel) exit("write_phonon_vel_all", "cannot open file_vel_all");

    const auto nk = dos->kmesh_dos->nk;
    const auto nk_irred = dos->kmesh_dos->nk_irred;
    const auto ns = dynamical->neval;
    const auto Ry_to_SI_vel = Bohr_in_Angstrom * 1.0e-10 / time_ry;
    const auto eval = dos->dymat_dos->get_eigenvalues();

    double ***phvel_xyz;
    double **phvel;

    allocate(phvel, nk, ns);
    allocate(phvel_xyz, nk, ns, 3);

    phonon_velocity->get_phonon_group_velocity_mesh(*dos->kmesh_dos,
                                                    system->lavec_p,
                                                    fcs_phonon->fc2_ext,
                                                    false,
                                                    phvel_xyz);
    unsigned int ik, is;
#ifdef _OPENMP
#pragma omp parallel for private(is)
#endif
    for (ik = 0; ik < nk; ++ik) {
        for (is = 0; is < ns; ++is) {
            phvel[ik][is] = std::sqrt(std::pow(phvel_xyz[ik][is][0], 2)
                                      + std::pow(phvel_xyz[ik][is][1], 2)
                                      + std::pow(phvel_xyz[ik][is][2], 2));
        }
    }

    ofs_vel << "# Phonon group velocity at all reducible k points." << std::endl;
    ofs_vel << "# irred. knum, knum, mode num, frequency [cm^-1], "
               "|velocity| [m/sec], velocity_(x,y,z) [m/sec]" << std::endl << std::endl;
    ofs_vel.setf(std::ios::fixed);

    for (unsigned int i = 0; i < nk_irred; ++i) {
        ofs_vel << "# Irreducible k point  : " << std::setw(8) << i + 1;
        ofs_vel << " (" << std::setw(4) << dos->kmesh_dos->kpoint_irred_all[i].size() << ")" << std::endl;

        for (unsigned int j = 0; j < dos->kmesh_dos->kpoint_irred_all[i].size(); ++j) {
            const auto knum = dos->kmesh_dos->kpoint_irred_all[i][j].knum;

            ofs_vel << "## xk =    ";
            for (auto k = 0; k < 3; ++k)
                ofs_vel << std::setw(15) << std::fixed
                        << std::setprecision(10) << dos->kmesh_dos->xk[knum][k];
            ofs_vel << std::endl;

            for (auto k = 0; k < ns; ++k) {
                ofs_vel << std::setw(7) << i + 1;
                ofs_vel << std::setw(8) << knum + 1;
                ofs_vel << std::setw(5) << k + 1;
                ofs_vel << std::setw(10) << std::fixed
                        << std::setprecision(2) << in_kayser(eval[knum][k]);
                ofs_vel << std::setw(10) << std::fixed
                        << std::setprecision(2) << phvel[knum][k] * Ry_to_SI_vel;
                for (auto ii = 0; ii < 3; ++ii) {
                    ofs_vel << std::setw(10) << std::fixed << std::setprecision(2)
                            << phvel_xyz[knum][k][ii] * Ry_to_SI_vel;
                }
                ofs_vel << std::endl;
            }
            ofs_vel << std::endl;
        }

        ofs_vel << std::endl;
    }

    ofs_vel.close();

    std::cout << "  " << std::setw(input->job_title.length() + 12) << std::left << file_vel;
    std::cout << " : Phonon velocity at all k points" << std::endl;

    deallocate(phvel);
    deallocate(phvel_xyz);
}

void Writes::write_phonon_dos() const
{
    int i;
    std::ofstream ofs_dos;
    auto file_dos = input->job_title + ".dos";

    ofs_dos.open(file_dos.c_str(), std::ios::out);
    if (!ofs_dos) exit("write_phonon_dos", "cannot open file_dos");

    ofs_dos << "#";
    for (i = 0; i < system->nkd; ++i) {
        ofs_dos << std::setw(5) << system->symbol_kd[i];
    }
    ofs_dos << std::endl;
    ofs_dos << "#";

    unsigned int *nat_each_kd;
    allocate(nat_each_kd, system->nkd);
    for (i = 0; i < system->nkd; ++i) nat_each_kd[i] = 0;
    for (i = 0; i < system->natmin; ++i) {
        ++nat_each_kd[system->kd[system->map_p2s[i][0]]];
    }
    for (i = 0; i < system->nkd; ++i) {
        ofs_dos << std::setw(5) << nat_each_kd[i];
    }
    ofs_dos << std::endl;
    deallocate(nat_each_kd);

    if (dos->compute_dos) {
        ofs_dos << "# Energy [cm^-1], TOTAL-DOS";
    } else {
        ofs_dos << "# Energy [cm^-1]";
    }
    if (dos->projected_dos) {
        ofs_dos << ", Atom Projected-DOS";
    }
    ofs_dos << std::endl;
    ofs_dos.setf(std::ios::scientific);

    for (i = 0; i < dos->n_energy; ++i) {
        ofs_dos << std::setw(15) << dos->energy_dos[i];
        if (dos->compute_dos) {
            ofs_dos << std::setw(15) << dos->dos_phonon[i];
        }
        if (dos->projected_dos) {
            for (auto iat = 0; iat < system->natmin; ++iat) {
                ofs_dos << std::setw(15) << dos->pdos_phonon[iat][i];
            }
        }
        ofs_dos << std::endl;
    }
    ofs_dos.close();

    std::cout << "  " << std::setw(input->job_title.length() + 12) << std::left << file_dos;

    if (dos->projected_dos & dos->compute_dos) {
        std::cout << " : Phonon DOS and atom projected DOS" << std::endl;
    } else if (dos->projected_dos) {
        std::cout << " : Atom projected phonon DOS" << std::endl;
    } else {
        std::cout << " : Phonon DOS" << std::endl;
    }
}

void Writes::write_two_phonon_dos() const
{
    std::ofstream ofs_tdos;
    auto file_tdos = input->job_title + ".tdos";
    ofs_tdos.open(file_tdos.c_str(), std::ios::out);

    ofs_tdos << "# Two-phonon DOS (TDOS) for all irreducible k points. " << std::endl;
    ofs_tdos << "# Energy [cm^-1], emission delta(e-e1-e2), absorption delta (e-e1+e2)" << std::endl;

    const auto n = dos->n_energy;

    for (auto ik = 0; ik < dos->kmesh_dos->nk_irred; ++ik) {

        ofs_tdos << "# Irred. kpoint : " << std::setw(5) << ik + 1 << std::endl;
        for (auto i = 0; i < n; ++i) {
            ofs_tdos << std::setw(15) << dos->emin + dos->delta_e * static_cast<double>(i);

            for (auto j = 0; j < 2; ++j) ofs_tdos << std::setw(15) << dos->dos2_phonon[ik][i][j];
            ofs_tdos << std::endl;
        }
        ofs_tdos << std::endl;
    }

    ofs_tdos.close();

    std::cout << "  " << std::setw(input->job_title.length() + 12) << std::left << file_tdos;
    std::cout << " : Two-phonon DOS" << std::endl;
}

void Writes::write_scattering_phase_space() const
{
    std::ofstream ofs_sps;

    auto file_sps = input->job_title + ".sps";
    ofs_sps.open(file_sps.c_str(), std::ios::out);

    ofs_sps << "# Total scattering phase space (cm): "
            << std::scientific << dos->total_sps3 << std::endl;
    ofs_sps << "# Mode decomposed scattering phase space are printed below." << std::endl;
    ofs_sps << "# Irred. k, mode, omega (cm^-1), P+ (absorption) (cm), P- (emission) (cm)" << std::endl;

    for (auto ik = 0; ik < dos->kmesh_dos->nk_irred; ++ik) {
        const auto knum = dos->kmesh_dos->kpoint_irred_all[ik][0].knum;

        for (auto is = 0; is < dynamical->neval; ++is) {
            ofs_sps << std::setw(5) << ik + 1;
            ofs_sps << std::setw(5) << is + 1;
            ofs_sps << std::setw(15) << in_kayser(dos->dymat_dos->get_eigenvalues()[knum][is]);
            ofs_sps << std::setw(15) << std::scientific << dos->sps3_mode[ik][is][1];
            ofs_sps << std::setw(15) << std::scientific << dos->sps3_mode[ik][is][0];
            ofs_sps << std::endl;
        }
        ofs_sps << std::endl;
    }

    ofs_sps.close();

    std::cout << "  " << std::setw(input->job_title.length() + 12) << std::left << file_sps;
    std::cout << " : Three-phonon scattering phase space" << std::endl;
}

void Writes::write_scattering_amplitude() const
{
    int i, j;
    unsigned int knum;
    const auto ns = dynamical->neval;

    auto file_w = input->job_title + ".sps_Bose";
    std::ofstream ofs_w;

    const auto Tmin = system->Tmin;
    const auto Tmax = system->Tmax;
    const auto dT = system->dT;
    const auto NT = static_cast<unsigned int>((Tmax - Tmin) / dT) + 1;

    ofs_w.open(file_w.c_str(), std::ios::out);

    ofs_w << "# Scattering phase space with the Bose-Einstein distribution function" << std::endl;
    ofs_w << "# Irreducible kpoints " << std::endl;
    for (i = 0; i < dos->kmesh_dos->kpoint_irred_all.size(); ++i) {
        ofs_w << "#" << std::setw(5) << i + 1;

        knum = dos->kmesh_dos->kpoint_irred_all[i][0].knum;
        for (j = 0; j < 3; ++j) ofs_w << std::setw(15) << dos->kmesh_dos->xk[knum][j];
        ofs_w << std::endl;
    }
    ofs_w << std::endl;
    ofs_w << "# k, mode, frequency (cm^-1), temperature, W+ (absorption) (cm), W- (emission) (cm)"
          << std::endl << std::endl;

    for (i = 0; i < dos->kmesh_dos->kpoint_irred_all.size(); ++i) {

        knum = dos->kmesh_dos->kpoint_irred_all[i][0].knum;

        for (unsigned int is = 0; is < ns; ++is) {

            const auto omega = in_kayser(dos->dymat_dos->get_eigenvalues()[knum][is]);

            for (j = 0; j < NT; ++j) {
                ofs_w << std::setw(5) << i + 1 << std::setw(5) << is + 1 << std::setw(15) << omega;
                ofs_w << std::setw(8) << Tmin + static_cast<double>(j) * dT;
                ofs_w << std::setw(15) << dos->sps3_with_bose[i][is][j][1];
                ofs_w << std::setw(15) << dos->sps3_with_bose[i][is][j][0];
                ofs_w << std::endl;
            }
            ofs_w << std::endl;
        }
        ofs_w << std::endl;
    }

    ofs_w.close();
    std::cout << "  " << std::setw(input->job_title.length() + 12) << std::left << file_w;
    std::cout << " : Three-phonon scattering phase space " << std::endl;
    std::cout << " " << std::setw(input->job_title.length() + 16) << " "
              << "with the Bose distribution function" << std::endl;
}

void Writes::write_normal_mode_direction() const
{
    std::string fname_axsf;

    if (kpoint->kpoint_general && dynamical->dymat_general) {
        fname_axsf = input->job_title + ".axsf";
        write_normal_mode_direction_each(fname_axsf,
                                         kpoint->kpoint_general->nk,
                                         dynamical->dymat_general->get_eigenvectors());
    }

    if (kpoint->kpoint_bs && dynamical->dymat_band) {
        fname_axsf = input->job_title + ".band.axsf";
        write_normal_mode_direction_each(fname_axsf,
                                         kpoint->kpoint_bs->nk,
                                         dynamical->dymat_band->get_eigenvectors());
    }

    if (dos->kmesh_dos && dos->dymat_dos) {
        fname_axsf = input->job_title + ".mesh.axsf";
        write_normal_mode_direction_each(fname_axsf,
                                         dos->kmesh_dos->nk,
                                         dos->dymat_dos->get_eigenvectors());
    }

}

void Writes::write_normal_mode_direction_each(const std::string &fname_axsf,
                                              const unsigned int nk_in,
                                              const std::complex<double> *const *const *evec_in) const
{
    std::ofstream ofs_anime;

    ofs_anime.open(fname_axsf.c_str(), std::ios::out);
    if (!ofs_anime) exit("write_mode_anime", "cannot open fname_axsf");

    ofs_anime.setf(std::ios::scientific);

    unsigned int i, j, k;
    const auto natmin = system->natmin;
    const auto force_factor = 100.0;

    double **xmod;
    std::string *kd_tmp;

    allocate(xmod, natmin, 3);
    allocate(kd_tmp, natmin);

    ofs_anime << "ANIMSTEPS " << nbands * nk_in << std::endl;
    ofs_anime << "CRYSTAL" << std::endl;
    ofs_anime << "PRIMVEC" << std::endl;

    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            ofs_anime << std::setw(15) << system->lavec_p[j][i] * Bohr_in_Angstrom;
        }
        ofs_anime << std::endl;
    }

    for (i = 0; i < natmin; ++i) {
        k = system->map_p2s[i][0];
        for (j = 0; j < 3; ++j) {
            xmod[i][j] = system->xc[k][j];
        }

        for (j = 0; j < 3; ++j) {
            xmod[i][j] *= Bohr_in_Angstrom;
        }
        kd_tmp[i] = system->symbol_kd[system->kd[k]];
    }

    i = 0;

    for (unsigned int ik = 0; ik < nk_in; ++ik) {
        for (unsigned int imode = 0; imode < nbands; ++imode) {
            ofs_anime << "PRIMCOORD " << std::setw(10) << i + 1 << std::endl;
            ofs_anime << std::setw(10) << natmin << std::setw(10) << 1 << std::endl;
            auto norm = 0.0;

            for (j = 0; j < 3 * natmin; ++j) {
                auto evec_tmp = evec_in[ik][imode][j];
                norm += std::pow(evec_tmp.real(), 2) + std::pow(evec_tmp.imag(), 2);
            }

            norm *= force_factor / static_cast<double>(natmin);

            for (j = 0; j < natmin; ++j) {

                const auto m = system->map_p2s[j][0];

                ofs_anime << std::setw(10) << kd_tmp[j];

                for (k = 0; k < 3; ++k) {
                    ofs_anime << std::setw(15) << xmod[j][k];
                }
                for (k = 0; k < 3; ++k) {
                    ofs_anime << std::setw(15)
                              << evec_in[ik][imode][3 * j + k].real()
                                 / (std::sqrt(system->mass[m]) * norm);
                }
                ofs_anime << std::endl;
            }

            ++i;
        }
    }

    deallocate(xmod);
    deallocate(kd_tmp);

    ofs_anime.close();
    std::cout << "  " << std::setw(input->job_title.length() + 12) << std::left << fname_axsf;
    std::cout << " : XcrysDen AXSF file to visualize phonon mode directions" << std::endl;
}

void Writes::write_eigenvectors() const
{
    std::string fname_evec;

    if (kpoint->kpoint_general && dynamical->dymat_general) {
        fname_evec = input->job_title + ".evec";
        write_eigenvectors_each(fname_evec,
                                kpoint->kpoint_general->nk,
                                kpoint->kpoint_general->xk,
                                dynamical->dymat_general->get_eigenvalues(),
                                dynamical->dymat_general->get_eigenvectors());
    }

    if (kpoint->kpoint_bs && dynamical->dymat_band) {
        fname_evec = input->job_title + ".band.evec";
        write_eigenvectors_each(fname_evec,
                                kpoint->kpoint_bs->nk,
                                kpoint->kpoint_bs->xk,
                                dynamical->dymat_band->get_eigenvalues(),
                                dynamical->dymat_band->get_eigenvectors());
    }

    if (dos->kmesh_dos && dos->dymat_dos) {
        fname_evec = input->job_title + ".mesh.evec";
        write_eigenvectors_each(fname_evec,
                                dos->kmesh_dos->nk,
                                dos->kmesh_dos->xk,
                                dos->dymat_dos->get_eigenvalues(),
                                dos->dymat_dos->get_eigenvectors());
    }
}

void Writes::write_eigenvectors_each(const std::string &fname_evec,
                                     const unsigned int nk_in,
                                     const double *const *xk_in,
                                     const double *const *eval_in,
                                     const std::complex<double> *const *const *evec_in) const
{
    unsigned int i, j, k;
    const auto neval = dynamical->neval;
    std::ofstream ofs_evec;

    ofs_evec.open(fname_evec.c_str(), std::ios::out);
    if (!ofs_evec) exit("write_eigenvectors", "cannot open file_evec");
    ofs_evec.setf(std::ios::scientific);

    ofs_evec << "# Lattice vectors of the primitive cell" << std::endl;

    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            ofs_evec << std::setw(15) << system->lavec_p[j][i];
        }
        ofs_evec << std::endl;
    }

    ofs_evec << std::endl;
    ofs_evec << "# Reciprocal lattice vectors of the primitive cell" << std::endl;

    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            ofs_evec << std::setw(15) << system->rlavec_p[i][j];
        }
        ofs_evec << std::endl;
    }

    ofs_evec << std::endl;
    ofs_evec << "# Number of phonon modes: " << std::setw(10) << nbands << std::endl;
    ofs_evec << "# Number of k points : " << std::setw(10) << nk_in << std::endl;
    ofs_evec << "# Number of atomic kinds : " << std::setw(4) << system->nkd << '\n';
    ofs_evec << "# Atomic masses :";
    for (i = 0; i < system->nkd; ++i) {
        ofs_evec << std::setw(15) << system->mass_kd[i];
    }
    ofs_evec << "\n\n";
    ofs_evec << "# Eigenvalues and eigenvectors for each phonon modes below:" << std::endl << std::endl;

    unsigned int **index_bconnect_tmp;
    allocate(index_bconnect_tmp, nk_in, nbands);

    if (dynamical->index_bconnect) {
        for (i = 0; i < nk_in; ++i) {
            for (j = 0; j < nbands; ++j) {
                index_bconnect_tmp[i][j] = dynamical->index_bconnect[i][j];
            }
        }
    } else {
        for (i = 0; i < nk_in; ++i) {
            for (j = 0; j < nbands; ++j) {
                index_bconnect_tmp[i][j] = j;
            }
        }
    }

    for (i = 0; i < nk_in; ++i) {
        ofs_evec << "## kpoint " << std::setw(7) << i + 1 << " : ";
        for (j = 0; j < 3; ++j) {
            ofs_evec << std::setw(15) << xk_in[i][j];
        }
        ofs_evec << std::endl;
        for (j = 0; j < nbands; ++j) {

            k = index_bconnect_tmp[i][j];

            auto omega2 = eval_in[i][k];
            if (omega2 >= 0.0) {
                omega2 = omega2 * omega2;
            } else {
                omega2 = -omega2 * omega2;
            }

            ofs_evec << "### mode " << std::setw(8) << j + 1 << " : ";
            ofs_evec << std::setw(15) << omega2 << std::endl;

            for (unsigned int m = 0; m < neval; ++m) {
                ofs_evec << std::setw(15) << real(evec_in[i][k][m]);
                ofs_evec << std::setw(15) << imag(evec_in[i][k][m]) << std::endl;
            }
            ofs_evec << std::endl;
        }
        ofs_evec << std::endl;
    }
    ofs_evec.close();

    deallocate(index_bconnect_tmp);

    std::cout << "  " << std::setw(input->job_title.length() + 12) << std::left << fname_evec;
    std::cout << " : Eigenvector of all k points" << std::endl;
}

#ifdef _HDF5

void Writes::write_eigenvectors_HDF5() const
{
    std::string fname_evec;

    if (kpoint->kpoint_general && dynamical->dymat_general) {
        fname_evec = input->job_title + ".evec.hdf5";
        write_eigenvectors_each_HDF5(fname_evec,
                                     kpoint->kpoint_general->nk,
                                     kpoint->kpoint_general->xk,
                                     dynamical->dymat_general->get_eigenvalues(),
                                     dynamical->dymat_general->get_eigenvectors(), 0);
    }

    if (kpoint->kpoint_bs && dynamical->dymat_band) {
        fname_evec = input->job_title + ".band.evec.hdf5";
        write_eigenvectors_each_HDF5(fname_evec,
                                     kpoint->kpoint_bs->nk,
                                     kpoint->kpoint_bs->xk,
                                     dynamical->dymat_band->get_eigenvalues(),
                                     dynamical->dymat_band->get_eigenvectors(), 1);
    }

    if (dos->kmesh_dos && dos->dymat_dos) {
        fname_evec = input->job_title + ".mesh.evec.hdf5";
        write_eigenvectors_each_HDF5(fname_evec,
                                     dos->kmesh_dos->nk,
                                     dos->kmesh_dos->xk,
                                     dos->dymat_dos->get_eigenvalues(),
                                     dos->dymat_dos->get_eigenvectors(), 2);
    }
}

void Writes::write_eigenvectors_each_HDF5(const std::string &fname_evec,
                                          const unsigned int nk_in,
                                          const double *const *xk_in,
                                          const double *const *eval_in,
                                          const std::complex<double> *const *const *evec_in,
                                          const unsigned int kpmode_in) const
{
    using namespace H5;

    unsigned int i, j, k;
    const auto neval = dynamical->neval;
    std::ofstream ofs_evec;

    H5File file(fname_evec, H5F_ACC_TRUNC);
    Group group_cell(file.createGroup("/PrimitiveCell"));
    Group group_band(file.createGroup("/Eigenvalues"));
    Group group_kpoint(file.createGroup("/Kpoints"));

    // Write setting information
    hid_t str_datatype = H5Tcopy(H5T_C_S1);
    H5Tset_size(str_datatype, H5T_VARIABLE);
    std::vector<const char *> arr_c_str;
    for (unsigned int ii = 0; ii < system->nkd; ++ii) {
        arr_c_str.push_back(system->symbol_kd[ii].c_str());
    }
    hsize_t str_dim[1]{arr_c_str.size()};
    DataSpace *dataspace = new DataSpace(1, str_dim);
    DataSet *dataset = new DataSet(group_cell.createDataSet("elements",
                                                            str_datatype,
                                                            *dataspace));
    dataset->write(&arr_c_str[0], str_datatype);
    dataset->close();
    dataspace->close();

    std::vector<double> mass_tmp;
    for (i = 0; i < system->nkd; ++i) {
        mass_tmp.push_back(system->mass_kd[i]);
    }
    dataspace = new DataSpace(1, str_dim);
    dataset = new DataSet(group_cell.createDataSet("masses",
                                                   PredType::NATIVE_DOUBLE,
                                                   *dataspace));
    dataset->write(&mass_tmp[0], PredType::NATIVE_DOUBLE);
    dataset->close();
    dataspace->close();


    // Write primitive cell information
    hsize_t dims[2];
    dims[0] = 3;
    dims[1] = 3;
    double lavec_tmp[3][3];
    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            lavec_tmp[i][j] = system->lavec_p[j][i];
        }
    }
    dataspace = new DataSpace(2, dims);
    dataset = new DataSet(group_cell.createDataSet("lattice_vector",
                                                   PredType::NATIVE_DOUBLE,
                                                   *dataspace));
    dataset->write(lavec_tmp, PredType::NATIVE_DOUBLE);
    DataSpace attr_dataspace_str(H5S_SCALAR);
    Attribute myatt_in = dataset->createAttribute("unit",
                                                  str_datatype,
                                                  attr_dataspace_str);
    myatt_in.write(str_datatype, std::string("bohr"));
    myatt_in.close();
    dataset->close();
    dataspace->close();

    dims[0] = system->natmin;
    dims[1] = 3;
    std::vector<double> xfrac_1D(dims[0] * dims[1]);
    hsize_t counter = 0;

    double xtmp[3];
    for (i = 0; i < system->natmin; ++i) {
        for (j = 0; j < 3; ++j) xtmp[j] = system->xr_s[system->map_p2s[i][0]][j];
        rotvec(xtmp, xtmp, system->lavec_s);
        rotvec(xtmp, xtmp, system->rlavec_p);
        for (j = 0; j < 3; ++j) {
            xfrac_1D[counter++] = xtmp[j] / (2.0 * pi);
        }
    }
    dataspace = new DataSpace(2, dims);
    dataset = new DataSet(group_cell.createDataSet("fractional_coordinate",
                                                   PredType::NATIVE_DOUBLE,
                                                   *dataspace));
    dataset->write(&xfrac_1D[0], PredType::NATIVE_DOUBLE);
    dataset->close();
    dataspace->close();

    hsize_t dims2[1];
    dims2[0] = system->natmin;
    dataspace = new DataSpace(1, dims2);
    dataset = new DataSet(group_cell.createDataSet("atomic_kinds",
                                                   PredType::NATIVE_INT,
                                                   *dataspace));
    int kdtmp[dims[0]];
    for (i = 0; i < system->natmin; ++i) {
        k = system->map_p2s[i][0];
        kdtmp[i] = system->kd[k];
    }

    dataset->write(&kdtmp[0], PredType::NATIVE_INT);
    dataset->close();

    // write eigenvalues

    unsigned int **index_bconnect_tmp;
    int band_index_reordered = 0;
    allocate(index_bconnect_tmp, nk_in, nbands);

    if (dynamical->index_bconnect) {
        band_index_reordered = 1;
        for (i = 0; i < nk_in; ++i) {
            for (j = 0; j < nbands; ++j) {
                index_bconnect_tmp[i][j] = dynamical->index_bconnect[i][j];
            }
        }
    } else {
        for (i = 0; i < nk_in; ++i) {
            for (j = 0; j < nbands; ++j) {
                index_bconnect_tmp[i][j] = j;
            }
        }
    }

    // Write band structure information
    dims[0] = nk_in;
    dims[1] = nbands;

    hsize_t dims_evec[4];
    dims_evec[0] = nk_in;
    dims_evec[1] = nbands;
    dims_evec[2] = neval;
    dims_evec[3] = 2;

    double **freq_kayser;
    double ****evec_tmp;
    allocate(freq_kayser, nk_in, nbands);
    allocate(evec_tmp, nk_in, nbands, neval, 2);

    for (i = 0; i < nk_in; ++i) {
        for (j = 0; j < nbands; ++j) {
            k = index_bconnect_tmp[i][j];
            freq_kayser[i][j] = in_kayser(eval_in[i][k]);

            for (unsigned int m = 0; m < neval; ++m) {
                evec_tmp[i][j][m][0] = evec_in[i][k][m].real();
                evec_tmp[i][j][m][1] = evec_in[i][k][m].imag();
            }
        }
    }

    dataspace = new DataSpace(2, dims);
    dataset = new DataSet(group_band.createDataSet("frequencies",
                                                   PredType::NATIVE_DOUBLE,
                                                   *dataspace));
    IntType int_type(PredType::NATIVE_INT);
    DataSpace attr_dataspace_int(H5S_SCALAR);
    myatt_in = dataset->createAttribute("band_index_reordered",
                                        int_type,
                                        attr_dataspace_int);
    myatt_in.write(int_type, &band_index_reordered);
    myatt_in = dataset->createAttribute("unit",
                                        str_datatype,
                                        attr_dataspace_str);
    myatt_in.write(str_datatype, std::string("kayser (cm^-1)"));

    dataset->write(&freq_kayser[0][0], PredType::NATIVE_DOUBLE);
    myatt_in.close();
    dataset->close();
    dataspace->close();
    deallocate(freq_kayser);

    dataspace = new DataSpace(4, dims_evec);
    dataset = new DataSet(group_band.createDataSet("polarization_vectors",
                                                   PredType::NATIVE_DOUBLE,
                                                   *dataspace));
    dataset->write(&evec_tmp[0][0][0][0], PredType::NATIVE_DOUBLE);
    dataset->close();
    dataspace->close();

    deallocate(evec_tmp);

    group_cell.close();
    group_band.close();

    dims[0] = nk_in;
    dims[1] = 3;
    std::vector<double> xk_1D(dims[0] * dims[1]);
    counter = 0;

    for (i = 0; i < nk_in; ++i) {
        for (j = 0; j < 3; ++j) {
            xk_1D[counter++] = xk_in[i][j];
        }
    }
    dataspace = new DataSpace(2, dims);
    dataset = new DataSet(group_kpoint.createDataSet("kpoint_coordinates",
                                                     PredType::NATIVE_DOUBLE,
                                                     *dataspace));

    myatt_in = dataset->createAttribute("kpoint_mode",
                                        int_type,
                                        attr_dataspace_int);
    myatt_in.write(int_type, &kpmode_in);
    dataset->write(&xk_1D[0], PredType::NATIVE_DOUBLE);
    myatt_in.close();
    dataset->close();
    dataspace->close();

    if (kpmode_in == 1 && kpoint->kpoint_bs) {
        const auto kaxis = kpoint->kpoint_bs->kaxis;
        dims2[0] = nk_in;
        dataspace = new DataSpace(1, dims2);
        dataset = new DataSet(group_kpoint.createDataSet("bandstructure_xaxis",
                                                         PredType::NATIVE_DOUBLE, *dataspace));
        dataset->write(&kaxis[0], PredType::NATIVE_DOUBLE);
        dataset->close();
        dataspace->close();
    }

    group_kpoint.close();
}

#endif

double Writes::in_kayser(const double x) const
{
    return x * Ry_to_kayser;
}

void Writes::write_thermodynamics() const
{
    const auto Tmin = system->Tmin;
    const auto Tmax = system->Tmax;
    const auto dT = system->dT;

    const auto NT = static_cast<unsigned int>((Tmax - Tmin) / dT) + 1;

    std::ofstream ofs_thermo;
    auto file_thermo = input->job_title + ".thermo";
    ofs_thermo.open(file_thermo.c_str(), std::ios::out);
    if (!ofs_thermo) exit("write_thermodynamics", "cannot open file_thermo");
    if (thermodynamics->calc_FE_bubble) {
        ofs_thermo << "# The bubble free-energy is also shown." << std::endl;
        ofs_thermo <<
                   "# Temperature [K], Heat capacity / kB, Entropy / kB, Internal energy [Ry], Free energy (QHA) [Ry], Free energy (Bubble) [Ry]"
                   << std::endl;
    } else {
        ofs_thermo <<
                   "# Temperature [K], Heat capacity / kB, Entropy / kB, Internal energy [Ry], Free energy (QHA) [Ry]"
                   << std::endl;
    }

    if (thermodynamics->classical) {
        ofs_thermo << "# CLASSICAL = 1: use classical statistics" << std::endl;
    }

    for (unsigned int i = 0; i < NT; ++i) {
        const auto T = Tmin + dT * static_cast<double>(i);

        const auto heat_capacity = thermodynamics->Cv_tot(T,
                                                          dos->kmesh_dos->nk_irred,
                                                          dynamical->neval,
                                                          dos->kmesh_dos->kpoint_irred_all,
                                                          &dos->kmesh_dos->weight_k[0],
                                                          dos->dymat_dos->get_eigenvalues());

        const auto Svib = thermodynamics->vibrational_entropy(T,
                                                              dos->kmesh_dos->nk_irred,
                                                              dynamical->neval,
                                                              dos->kmesh_dos->kpoint_irred_all,
                                                              &dos->kmesh_dos->weight_k[0],
                                                              dos->dymat_dos->get_eigenvalues());

        const auto Uvib = thermodynamics->internal_energy(T,
                                                          dos->kmesh_dos->nk_irred,
                                                          dynamical->neval,
                                                          dos->kmesh_dos->kpoint_irred_all,
                                                          &dos->kmesh_dos->weight_k[0],
                                                          dos->dymat_dos->get_eigenvalues());

        const auto FE_QHA = thermodynamics->free_energy_QHA(T,
                                                            dos->kmesh_dos->nk_irred,
                                                            dynamical->neval,
                                                            dos->kmesh_dos->kpoint_irred_all,
                                                            &dos->kmesh_dos->weight_k[0],
                                                            dos->dymat_dos->get_eigenvalues());

        ofs_thermo << std::setw(16) << std::fixed << T;
        ofs_thermo << std::setw(18) << std::scientific << heat_capacity / k_Boltzmann;
        ofs_thermo << std::setw(18) << Svib / k_Boltzmann;
        ofs_thermo << std::setw(18) << Uvib;
        ofs_thermo << std::setw(18) << FE_QHA;

        if (thermodynamics->calc_FE_bubble) {
            ofs_thermo << std::setw(18) << thermodynamics->FE_bubble[i];
        }
        ofs_thermo << std::endl;
    }

    ofs_thermo.close();

    std::cout << "  " << std::setw(input->job_title.length() + 12) << std::left << file_thermo;
    std::cout << " : Thermodynamic quantities" << std::endl;
}

void Writes::write_gruneisen()
{
    if (kpoint->kpoint_bs && gruneisen->gruneisen_bs) {
        if (nbands < 0 || nbands > 3 * system->natmin) {
            nbands = 3 * system->natmin;
        }

        std::ofstream ofs_gruneisen;

        auto file_gru = input->job_title + ".gruneisen";
        ofs_gruneisen.open(file_gru.c_str(), std::ios::out);
        if (!ofs_gruneisen) exit("write_gruneisen", "cannot open file_vel");

        const auto nk = kpoint->kpoint_bs->nk;
        const auto kaxis = kpoint->kpoint_bs->kaxis;

        ofs_gruneisen << "# k-axis, gamma" << std::endl;
        ofs_gruneisen.setf(std::ios::fixed);

        if (dynamical->band_connection == 0) {
            for (unsigned int i = 0; i < nk; ++i) {
                ofs_gruneisen << std::setw(8) << kaxis[i];
                for (unsigned int j = 0; j < nbands; ++j) {
                    ofs_gruneisen << std::setw(15) << gruneisen->gruneisen_bs[i][j].real();
                }
                ofs_gruneisen << std::endl;
            }
        } else {
            for (unsigned int i = 0; i < nk; ++i) {
                ofs_gruneisen << std::setw(8) << kaxis[i];
                for (unsigned int j = 0; j < nbands; ++j) {
                    ofs_gruneisen << std::setw(15)
                                  << gruneisen->gruneisen_bs[i][dynamical->index_bconnect[i][j]].real();
                }
                ofs_gruneisen << std::endl;
            }
        }

        ofs_gruneisen.close();

        std::cout << "  " << std::setw(input->job_title.length() + 12) << std::left << file_gru;
        std::cout << " : Gruneisen parameters along given k-path" << std::endl;

    }

    if (dos->kmesh_dos && gruneisen->gruneisen_dos) {

        std::ofstream ofs_gruall;
        auto file_gruall = input->job_title + ".gru_all";
        ofs_gruall.open(file_gruall.c_str(), std::ios::out);
        if (!ofs_gruall) exit("write_gruneisen", "cannot open file_gruall");

        const auto nk = dos->kmesh_dos->nk;
        const auto ns = dynamical->neval;
        const auto xk = dos->kmesh_dos->xk;
        const auto eval = dos->dymat_dos->get_eigenvalues();

        ofs_gruall << "# knum, snum, omega [cm^-1], gruneisen parameter" << std::endl;

        for (unsigned int i = 0; i < nk; ++i) {
            ofs_gruall << "# knum = " << i;
            for (unsigned int k = 0; k < 3; ++k) {
                ofs_gruall << std::setw(15) << xk[i][k];
            }
            ofs_gruall << std::endl;

            for (unsigned int j = 0; j < ns; ++j) {
                ofs_gruall << std::setw(5) << i;
                ofs_gruall << std::setw(5) << j;
                ofs_gruall << std::setw(15) << in_kayser(eval[i][j]);
                ofs_gruall << std::setw(15) << gruneisen->gruneisen_dos[i][j].real();
                ofs_gruall << std::endl;
            }
        }
        ofs_gruall.close();

        std::cout << "  " << std::setw(input->job_title.length() + 12) << std::left << file_gruall;
        std::cout << " : Gruneisen parameters at all k points" << std::endl;
    }
}

void Writes::write_msd() const
{
    // Write room mean square displacement of atoms

    auto file_rmsd = input->job_title + ".msd";
    std::ofstream ofs_rmsd;

    const auto ns = dynamical->neval;

    const auto Tmin = system->Tmin;
    const auto Tmax = system->Tmax;
    const auto dT = system->dT;
    const auto nk = dos->kmesh_dos->nk;
    const auto xk = dos->kmesh_dos->xk;
    const auto eval = dos->dymat_dos->get_eigenvalues();
    const auto evec = dos->dymat_dos->get_eigenvectors();

    ofs_rmsd.open(file_rmsd.c_str(), std::ios::out);
    if (!ofs_rmsd) exit("write_rmsd", "Could not open file_rmsd");

    ofs_rmsd << "# Mean Square Displacements at a function of temperature." << std::endl;
    ofs_rmsd << "# Temperature [K], <(u_{1}^{x})^{2}>, <(u_{1}^{y})^{2}>, <(u_{1}^{z})^{2}>, .... [Angstrom^2]" << std::
    endl;

    const auto NT = static_cast<unsigned int>((Tmax - Tmin) / dT) + 1;

    for (unsigned int i = 0; i < NT; ++i) {

        const auto T = Tmin + static_cast<double>(i) * dT;
        ofs_rmsd << std::setw(15) << T;

        for (unsigned int j = 0; j < ns; ++j) {
            const auto d2_tmp = thermodynamics->disp2_avg(T, j, j, nk, ns, xk, eval, evec);
            ofs_rmsd << std::setw(15) << d2_tmp * std::pow(Bohr_in_Angstrom, 2.0);
        }
        ofs_rmsd << std::endl;
    }
    ofs_rmsd.close();

    std::cout << "  " << std::setw(input->job_title.length() + 12) << std::left << file_rmsd;
    std::cout << " : Mean-square-displacement (MSD)" << std::endl;
}

void Writes::write_scph_msd(double **msd_scph, const int bubble) const
{
    const auto ns = dynamical->neval;
    const auto Tmin = system->Tmin;
    const auto Tmax = system->Tmax;
    const auto dT = system->dT;
    const auto NT = static_cast<unsigned int>((Tmax - Tmin) / dT) + 1;

    std::ofstream ofs_msd;
    std::string file_msd;
    if (bubble == 0) {
        file_msd = input->job_title + ".scph_msd";
    } else if (bubble == 1) {
        file_msd = input->job_title + ".scph+bubble(0)_msd";
    } else if (bubble == 2) {
        file_msd = input->job_title + ".scph+bubble(w)_msd";
    } else if (bubble == 3) {
        file_msd = input->job_title + ".scph+bubble(wQP)_msd";
    }
    ofs_msd.open(file_msd.c_str(), std::ios::out);
    if (!ofs_msd) exit("write_scph_msd", "cannot open file_thermo");
    ofs_msd << "# Mean Square Displacements at a function of temperature." << std::endl;
    ofs_msd << "# Temperature [K], <(u_{1}^{x})^{2}>, <(u_{1}^{y})^{2}>, <(u_{1}^{z})^{2}>, .... [Angstrom^2]" << std::
    endl;

    for (unsigned int iT = 0; iT < NT; ++iT) {
        const auto temp = Tmin + static_cast<double>(iT) * dT;

        ofs_msd << std::setw(15) << temp;
        for (unsigned int i = 0; i < ns; ++i) {
            ofs_msd << std::setw(15) << msd_scph[iT][i] * std::pow(Bohr_in_Angstrom, 2.0);
        }
        ofs_msd << std::endl;
    }

    ofs_msd.close();
    std::cout << "  " << std::setw(input->job_title.length() + 12) << std::left << file_msd;
    if (bubble == 0) {
        std::cout << " : Mean-square-displacement (SCPH level)" << std::endl;
    } else if (bubble == 1) {
        std::cout << " : Mean-square-displacement (SCPH+Bubble(0) level)" << std::endl;
    } else if (bubble == 2) {
        std::cout << " : Mean-square-displacement (SCPH+Bubble(w) level)" << std::endl;
    } else if (bubble == 3) {
        std::cout << " : Mean-square-displacement (SCPH+Bubble(wQP) level)" << std::endl;
    }
}

void Writes::write_disp_correlation() const
{
    if (!dos->kmesh_dos) return;

    auto file_ucorr = input->job_title + ".ucorr";
    std::ofstream ofs;

    const auto ns = dynamical->neval;
    const auto Tmin = system->Tmin;
    const auto Tmax = system->Tmax;
    const auto dT = system->dT;
    const auto NT = static_cast<unsigned int>((Tmax - Tmin) / dT) + 1;

    ofs.open(file_ucorr.c_str(), std::ios::out);
    if (!ofs) exit("write_disp_correlation", "Could not open file_rmsd");

    ofs << "# Displacement-displacement correlation function at various temperatures." << std::endl;
    if (thermodynamics->classical) ofs << "# CLASSICAL = 1: classical statistics is used.\n";

    double shift[3];

    for (auto i = 0; i < 3; ++i) {
        shift[i] = static_cast<double>(shift_ucorr[i]);
    }

    ofs <<
        "# Temperature [K], (atom1,crd1), (atom2,crd2), SHIFT_UCORR, <u_{0,atom1}^{crd1} * u_{L, atom2}^{crd2}> [Angstrom^2]\n";

    for (unsigned int i = 0; i < NT; ++i) {

        const auto T = Tmin + static_cast<double>(i) * dT;

        for (unsigned int j = 0; j < ns; ++j) {
            for (unsigned int k = 0; k < ns; ++k) {

                const auto ucorr = thermodynamics->disp_corrfunc(T, j, k,
                                                                 shift,
                                                                 dos->kmesh_dos->nk,
                                                                 ns,
                                                                 dos->kmesh_dos->xk,
                                                                 dos->dymat_dos->get_eigenvalues(),
                                                                 dos->dymat_dos->get_eigenvectors());

                ofs << std::setw(17) << T;
                ofs << std::setw(11) << j / 3 + 1;
                ofs << std::setw(3) << j % 3 + 1;
                ofs << std::setw(11) << k / 3 + 1;
                ofs << std::setw(3) << k % 3 + 1;
                ofs << std::setw(4) << shift_ucorr[0];
                ofs << std::setw(4) << shift_ucorr[1];
                ofs << std::setw(4) << shift_ucorr[2];
                ofs << std::setw(15) << ucorr * std::pow(Bohr_in_Angstrom, 2.0);
                ofs << '\n';
            }
        }
        ofs << '\n';
    }
    ofs << std::flush;
    ofs.close();

    std::cout << "  " << std::setw(input->job_title.length() + 12) << std::left << file_ucorr;
    std::cout << " : displacement correlation functions" << std::endl;
}

void Writes::write_scph_ucorr(double ***ucorr_scph,
                              const int bubble) const
{
    std::string file_ucorr;
    std::ofstream ofs;

    const auto ns = dynamical->neval;
    const auto Tmin = system->Tmin;
    const auto Tmax = system->Tmax;
    const auto dT = system->dT;
    const auto NT = static_cast<unsigned int>((Tmax - Tmin) / dT) + 1;

    if (bubble == 0) {
        file_ucorr = input->job_title + ".scph_ucorr";
    } else if (bubble == 1) {
        file_ucorr = input->job_title + ".scph+bubble(0)_ucorr";
    } else if (bubble == 2) {
        file_ucorr = input->job_title + ".scph+bubble(w)_ucorr";
    }

    ofs.open(file_ucorr.c_str(), std::ios::out);
    if (!ofs) exit("write_disp_correlation", "Could not open file_rmsd");

    ofs << "# Displacement-displacement correlation function at various temperatures.\n";
    ofs << "# Self-consistent phonon frequencies and eigenvectors are used.\n";
    if (thermodynamics->classical) ofs << "# CLASSICAL = 1: classical statistics is used.\n";

    double shift[3];

    for (auto i = 0; i < 3; ++i) {
        shift[i] = static_cast<double>(shift_ucorr[i]);
    }

    ofs <<
        "# Temperature [K], (atom1,crd1), (atom2,crd2), SHIFT_UCORR, <u_{0,atom1}^{crd1} * u_{L, atom2}^{crd2}> [Angstrom^2]\n";

    for (unsigned int i = 0; i < NT; ++i) {

        const auto T = Tmin + static_cast<double>(i) * dT;

        for (unsigned int j = 0; j < ns; ++j) {
            for (unsigned int k = 0; k < ns; ++k) {

                ofs << std::setw(17) << T;
                ofs << std::setw(11) << j / 3 + 1;
                ofs << std::setw(3) << j % 3 + 1;
                ofs << std::setw(11) << k / 3 + 1;
                ofs << std::setw(3) << k % 3 + 1;
                ofs << std::setw(4) << writes->shift_ucorr[0];
                ofs << std::setw(4) << writes->shift_ucorr[1];
                ofs << std::setw(4) << writes->shift_ucorr[2];
                ofs << std::setw(15) << ucorr_scph[i][j][k] * std::pow(Bohr_in_Angstrom, 2.0);
                ofs << '\n';
            }
        }
        ofs << '\n';
    }
    ofs << std::flush;
    ofs.close();

    std::cout << "  " << std::setw(input->job_title.length() + 12) << std::left << file_ucorr;
    if (bubble == 0) {
        std::cout << " : displacement correlation functions (SCPH level)" << std::endl;
    } else if (bubble == 1) {
        std::cout << " : displacement correlation functions (SCPH+Bubble(0) level)" << std::endl;
    } else if (bubble == 2) {
        std::cout << " : displacement correlation functions (SCPH+Bubble(w) level)" << std::endl;
    }
}

void Writes::write_kappa() const
{
    // Write lattice thermal conductivity

    if (mympi->my_rank == 0) {
        int i, j, k;

        std::string file_kappa;
        std::string file_kappa_3only;

        if (conductivity->fph_rta > 0) {
            file_kappa_3only = input->job_title + ".kl3";
            file_kappa = input->job_title + ".kl4";
        } else {
            file_kappa = input->job_title + ".kl";
        }

        auto file_kappa2 = input->job_title + ".kl_spec";
        auto file_kappa_coherent = input->job_title + ".kl_coherent";

        std::ofstream ofs_kl;

        if (conductivity->fph_rta > 0) {
            ofs_kl.open(file_kappa_3only.c_str(), std::ios::out);
            if (!ofs_kl) exit("write_kappa", "Could not open file_kappa");

            ofs_kl << "# Temperature [K], Thermal Conductivity (xx, xy, xz, yx, yy, yz, zx, zy, zz) [W/mK]"
                   << std::endl;
            ofs_kl << "# three phonon part";

            if (isotope->include_isotope) {
                ofs_kl << "# Isotope effects are included." << std::endl;
            }

            if (conductivity->len_boundary > eps) {
                ofs_kl << "# Size of boundary " << std::scientific << std::setprecision(2)
                       << conductivity->len_boundary * 1e9 << " [nm]" << std::endl;
            }

            for (i = 0; i < conductivity->ntemp; ++i) {
                ofs_kl << std::setw(10) << std::right << std::fixed << std::setprecision(2)
                       << conductivity->temperature[i];
                for (j = 0; j < 3; ++j) {
                    for (k = 0; k < 3; ++k) {
                        ofs_kl << std::setw(15) << std::fixed
                               << std::setprecision(4) << conductivity->kappa_3only[i][j][k];
                    }
                }
                ofs_kl << std::endl;
            }
            ofs_kl.close();
        }

        ofs_kl.open(file_kappa.c_str(), std::ios::out);
        if (!ofs_kl) exit("write_kappa", "Could not open file_kappa");

        ofs_kl << "# Temperature [K], Thermal Conductivity (xx, xy, xz, yx, yy, yz, zx, zy, zz) [W/mK]" << std::endl;

        if (isotope->include_isotope) {
            ofs_kl << "# Isotope effects are included." << std::endl;
        }

        if (conductivity->len_boundary > eps) {
            ofs_kl << "# Size of boundary " << std::scientific << std::setprecision(2)
                   << conductivity->len_boundary * 1e9 << " [nm]" << std::endl;
        }

        for (i = 0; i < conductivity->ntemp; ++i) {
            ofs_kl << std::setw(10) << std::right << std::fixed << std::setprecision(2)
                   << conductivity->temperature[i];
            for (j = 0; j < 3; ++j) {
                for (k = 0; k < 3; ++k) {
                    ofs_kl << std::setw(15) << std::fixed
                           << std::setprecision(4) << conductivity->kappa[i][j][k];
                }
            }
            ofs_kl << std::endl;
        }
        ofs_kl.close();

        if (conductivity->calc_kappa_spec) {

            ofs_kl.open(file_kappa2.c_str(), std::ios::out);
            if (!ofs_kl) exit("write_kappa", "Could not open file_kappa2");

            ofs_kl << "# Temperature [K], Frequency [cm^-1], Thermal Conductivity Spectra (xx, yy, zz) [W/mK * cm]" <<
                   std::endl;

            if (isotope->include_isotope) {
                ofs_kl << "# Isotope effects are included." << std::endl;
            }

            for (i = 0; i < conductivity->ntemp; ++i) {
                for (j = 0; j < dos->n_energy; ++j) {
                    ofs_kl << std::setw(10) << std::right << std::fixed << std::setprecision(2)
                           << conductivity->temperature[i];
                    ofs_kl << std::setw(10) << dos->energy_dos[j];
                    for (k = 0; k < 3; ++k) {
                        ofs_kl << std::setw(15) << std::fixed
                               << std::setprecision(6) << conductivity->kappa_spec[j][i][k];
                    }
                    ofs_kl << std::endl;
                }
                ofs_kl << std::endl;
            }
            ofs_kl.close();
        }

        if (conductivity->calc_coherent) {
            ofs_kl.open(file_kappa_coherent.c_str(), std::ios::out);
            if (!ofs_kl) exit("write_kappa", "Could not open file_kappa_coherent");

            ofs_kl << "# Temperature [K], Coherent part of the lattice thermal Conductivity (xx, yy, zz) [W/mK]" <<
                   std::endl;

            if (isotope->include_isotope) {
                ofs_kl << "# Isotope effects are included." << std::endl;
            }

            for (i = 0; i < conductivity->ntemp; ++i) {
                ofs_kl << std::setw(10) << std::right << std::fixed << std::setprecision(2)
                       << conductivity->temperature[i];
                for (j = 0; j < 3; ++j) {
                    ofs_kl << std::setw(15) << std::fixed
                           << std::setprecision(4) << conductivity->kappa_coherent[i][j][j];
                }
                ofs_kl << std::endl;
            }
            ofs_kl.close();
        }

        std::cout << std::endl;
        std::cout << " -----------------------------------------------------------------" << std::endl << std::endl;
        std::cout << " Lattice thermal conductivity is stored in the file " << file_kappa << std::endl;
        if (conductivity->calc_kappa_spec) {
            std::cout << " Thermal conductivity spectra is stored in the file " << file_kappa2 << std::endl;
        }
        if (conductivity->calc_coherent) {
            std::cout << " Coherent part is stored in the file " << file_kappa_coherent << std::endl;
        }
    }
}

void Writes::write_selfenergy_isotope() const
{
    unsigned int k;
    const auto ns = dynamical->neval;
    const auto eval = dos->dymat_dos->get_eigenvalues();
    const auto gamma_iso = isotope->gamma_isotope;

    if (mympi->my_rank == 0) {
        if (isotope->include_isotope == 2) {

            auto file_iso = input->job_title + ".self_isotope";
            std::ofstream ofs_iso;

            ofs_iso.open(file_iso.c_str(), std::ios::out);
            if (!ofs_iso) exit("write_selfenergy_isotope", "Could not open file_iso");

            ofs_iso << "# Phonon selfenergy due to phonon-isotope scatterings for the irreducible k points." << std::
            endl;
            ofs_iso << "# Irred. knum, mode num, frequency [cm^-1], Gamma_iso [cm^-1]" << std::endl << std::endl;

            for (unsigned int i = 0; i < dos->kmesh_dos->nk_irred; ++i) {
                ofs_iso << "# Irreducible k point  : " << std::setw(8) << i + 1;
                ofs_iso << " (" << std::setw(4) << dos->kmesh_dos->kpoint_irred_all[i].size() << ")" << std::endl;

                const auto knum = dos->kmesh_dos->kpoint_irred_all[i][0].knum;

                ofs_iso << "## xk = " << std::setw(3);
                for (k = 0; k < 3; ++k) ofs_iso << std::setw(15) << dos->kmesh_dos->xk[knum][k];
                ofs_iso << std::endl;

                for (k = 0; k < ns; ++k) {
                    ofs_iso << std::setw(7) << i + 1;
                    ofs_iso << std::setw(5) << k + 1;
                    ofs_iso << std::setw(15) << in_kayser(eval[knum][k]);
                    ofs_iso << std::setw(15) << in_kayser(gamma_iso[i][k]);
                    ofs_iso << std::endl;
                }
                ofs_iso << std::endl;
            }

            std::cout << std::endl;
            std::cout << " ISOTOPE = 2: Phonon selfenergy due to phonon-isotope " << std::endl;
            std::cout << "              scatterings is stored in the file " << file_iso << std::endl;

            ofs_iso.close();
        }
    }
}

void Writes::write_normal_mode_animation(const double xk_in[3],
                                         const unsigned int ncell[3]) const
{
    unsigned int i, j, k;
    unsigned int iband, istep;
    const auto ns = dynamical->neval;
    const auto natmin = system->natmin;
    const auto nsuper = ncell[0] * ncell[1] * ncell[2];
    unsigned int ntmp = nbands;
    unsigned int ndigits = 0;

    double phase_time;
    const auto max_disp_factor = 0.1;
    double lavec_super[3][3];
    double dmod[3];
    double xk[3], kvec[3];

    double *eval, **evec_mag, **evec_theta;
    double **disp_mag, *mass;
    double *phase_cell;
    double ***xmod, **xtmp;

    std::complex<double> **evec;

    std::ofstream ofs_anime;
    std::ostringstream ss;
    std::string file_anime;
    std::string *kd_tmp;

    for (i = 0; i < 3; ++i) {
        xk[i] = xk_in[i];
    }
    std::cout << " -----------------------------------------------------------------" << std::endl << std::endl;
    std::cout << " ANIME-tag is given: Making animation files for the given" << std::endl;
    std::cout << "                     k point ( ";
    std::cout << std::setw(5) << xk[0] << ", "
              << std::setw(5) << xk[1] << ", "
              << std::setw(5) << xk[2] << ")." << std::endl;
    std::cout << " ANIME_CELLSIZE = ";
    std::cout << std::setw(3) << ncell[0]
              << std::setw(3) << ncell[1]
              << std::setw(3) << ncell[2] << std::endl;
    std::cout << " ANIME_FORMAT = " << anime_format << std::endl;

    for (i = 0; i < 3; ++i) dmod[i] = std::fmod(xk[i] * static_cast<double>(ncell[i]), 1.0);

    if (std::sqrt(dmod[0] * dmod[0] + dmod[1] * dmod[1] + dmod[2] * dmod[2]) > eps12) {
        warn("write_normal_mode_animation",
             "The supercell size is not commensurate with given k point.");
    }

    rotvec(kvec, xk, system->rlavec_p, 'T');
    const auto norm = std::sqrt(kvec[0] * kvec[0] + kvec[1] * kvec[1] + kvec[2] * kvec[2]);
    if (norm > eps) {
        for (i = 0; i < 3; ++i) kvec[i] /= norm;
    }

    // Allocation

    allocate(eval, ns);
    allocate(evec, ns, ns);
    allocate(evec_mag, ns, ns);
    allocate(evec_theta, ns, ns);
    allocate(disp_mag, ns, ns);
    allocate(xmod, nsuper, natmin, 3);
    allocate(kd_tmp, natmin);
    allocate(mass, natmin);
    allocate(phase_cell, nsuper);

    // Get eigenvalues and eigenvectors at xk

    dynamical->eval_k(xk, kvec, fcs_phonon->fc2_ext, eval, evec, true);

    for (i = 0; i < ns; ++i) {
        for (j = 0; j < ns; ++j) {
            evec_mag[i][j] = std::abs(evec[i][j]);
            evec_theta[i][j] = std::arg(evec[i][j]);
        }
    }

    // Get fractional coordinates of atoms in a primitive cell

    allocate(xtmp, natmin, 3);

    for (i = 0; i < natmin; ++i) {
        rotvec(xtmp[i], system->xr_s[system->map_p2s[i][0]], system->lavec_s);
        rotvec(xtmp[i], xtmp[i], system->rlavec_p);
        for (j = 0; j < 3; ++j) xtmp[i][j] /= 2.0 * pi;
    }

    // Prepare fractional coordinates of atoms in the supercell
    unsigned int icell = 0;

    for (unsigned int ix = 0; ix < ncell[0]; ++ix) {
        for (unsigned int iy = 0; iy < ncell[1]; ++iy) {
            for (unsigned int iz = 0; iz < ncell[2]; ++iz) {

                phase_cell[icell] = 2.0 * pi * (xk_in[0] * static_cast<double>(ix)
                                                + xk_in[1] * static_cast<double>(iy)
                                                + xk_in[2] * static_cast<double>(iz));

                for (i = 0; i < natmin; ++i) {
                    xmod[icell][i][0] = (xtmp[i][0] + static_cast<double>(ix)) / static_cast<double>(ncell[0]);
                    xmod[icell][i][1] = (xtmp[i][1] + static_cast<double>(iy)) / static_cast<double>(ncell[1]);
                    xmod[icell][i][2] = (xtmp[i][2] + static_cast<double>(iz)) / static_cast<double>(ncell[2]);
                }
                ++icell;
            }
        }
    }

    deallocate(xtmp);

    // Prepare atomic symbols and masses

    for (i = 0; i < natmin; ++i) {
        k = system->map_p2s[i][0];
        kd_tmp[i] = system->symbol_kd[system->kd[k]];
        mass[i] = system->mass[k];
    }

    // Prepare lattice vectors of the supercell

    for (i = 0; i < 3; ++i) {
        lavec_super[i][0] = system->lavec_p[i][0] * ncell[0] * Bohr_in_Angstrom;
        lavec_super[i][1] = system->lavec_p[i][1] * ncell[1] * Bohr_in_Angstrom;
        lavec_super[i][2] = system->lavec_p[i][2] * ncell[2] * Bohr_in_Angstrom;
    }

    // Normalize the magnitude of displacements

    auto mass_min = mass[0];
    for (i = 0; i < natmin; ++i) {
        if (mass[i] < mass_min) mass_min = mass[i];
    }

    for (iband = 0; iband < nbands; ++iband) {
        auto max_disp_mag = 0.0;

        for (j = 0; j < ns; ++j) {
            disp_mag[iband][j] = std::sqrt(mass_min / mass[j / 3]) * evec_mag[iband][j];
        }

        for (j = 0; j < natmin; ++j) {
            auto disp_mag_tmp = 0.0;
            for (k = 0; k < 3; ++k) disp_mag_tmp += std::pow(disp_mag[iband][3 * j + k], 2);
            disp_mag_tmp = std::sqrt(disp_mag_tmp);
            max_disp_mag = std::max(max_disp_mag, disp_mag_tmp);
        }

        for (j = 0; j < ns; ++j) disp_mag[iband][j] *= max_disp_factor / max_disp_mag;
    }

    // Convert atomic positions to Cartesian coordinate
    for (i = 0; i < nsuper; ++i) {
        for (j = 0; j < natmin; ++j) {
            rotvec(xmod[i][j], xmod[i][j], lavec_super);
        }
    }

    while (ntmp > 0) {
        ++ndigits;
        ntmp /= 10;
    }

    if (anime_format == "XSF" || anime_format == "AXSF") {

        // Save animation to AXSF (XcrysDen) files

        for (iband = 0; iband < nbands; ++iband) {

            eval[iband] = dynamical->freq(eval[iband]);
            ss.str("");
            ss.clear();
            ss << std::setw(ndigits) << std::setfill('0') << iband + 1;
            const auto result = ss.str();

            file_anime = input->job_title + ".anime" + result + ".axsf";

            ofs_anime.open(file_anime.c_str(), std::ios::out);
            if (!ofs_anime)
                exit("write_normal_mode_animation",
                     "cannot open file_anime");

            ofs_anime.setf(std::ios::scientific);

            ofs_anime << "ANIMSTEPS " << anime_frames << std::endl;
            ofs_anime << "CRYSTAL" << std::endl;
            ofs_anime << "PRIMVEC" << std::endl;

            for (i = 0; i < 3; ++i) {
                for (j = 0; j < 3; ++j) {
                    ofs_anime << std::setw(15) << lavec_super[j][i];
                }
                ofs_anime << std::endl;
            }

            for (istep = 0; istep < anime_frames; ++istep) {

                phase_time = 2.0 * pi / static_cast<double>(anime_frames) * static_cast<double>(istep);

                ofs_anime << "PRIMCOORD " << std::setw(10) << istep + 1 << std::endl;
                ofs_anime << std::setw(10) << natmin * nsuper << std::setw(10) << 1 << std::endl;

                for (i = 0; i < nsuper; ++i) {
                    for (j = 0; j < natmin; ++j) {

                        ofs_anime << std::setw(10) << kd_tmp[j];

                        for (k = 0; k < 3; ++k) {
                            ofs_anime << std::setw(15)
                                      << xmod[i][j][k]
                                         + disp_mag[iband][3 * j + k]
                                           * std::sin(phase_cell[i] + evec_theta[iband][3 * j + k] + phase_time);
                        }
                        ofs_anime << std::endl;
                    }
                }
            }

            ofs_anime.close();
        }

    } else if (anime_format == "XYZ") {

        // Save animation to XYZ files

        for (iband = 0; iband < nbands; ++iband) {

            eval[iband] = dynamical->freq(eval[iband]);
            ss.str("");
            ss.clear();
            ss << std::setw(ndigits) << std::setfill('0') << iband + 1;
            const auto result = ss.str();

            file_anime = input->job_title + ".anime" + result + ".xyz";

            ofs_anime.open(file_anime.c_str(), std::ios::out);
            if (!ofs_anime)
                exit("write_normal_mode_animation",
                     "cannot open file_anime");

            ofs_anime.setf(std::ios::scientific);

            for (istep = 0; istep < anime_frames; ++istep) {

                phase_time = 2.0 * pi / static_cast<double>(anime_frames) * static_cast<double>(istep);

                ofs_anime.unsetf(std::ios::scientific);

                ofs_anime << natmin * nsuper << std::endl;
                ofs_anime << "Mode " << std::setw(4) << iband + 1 << " at (";
                for (i = 0; i < 3; ++i) ofs_anime << std::setw(8) << xk_in[i];
                ofs_anime << "), Frequency (cm^-1) = " << in_kayser(eval[iband])
                          << ", Time step = " << std::setw(4) << istep + 1 << std::endl;

                ofs_anime.setf(std::ios::scientific);

                for (i = 0; i < nsuper; ++i) {
                    for (j = 0; j < natmin; ++j) {

                        ofs_anime << std::setw(4) << kd_tmp[j];

                        for (k = 0; k < 3; ++k) {
                            ofs_anime << std::setw(15)
                                      << xmod[i][j][k]
                                         + disp_mag[iband][3 * j + k]
                                           * std::sin(phase_cell[i] + evec_theta[iband][3 * j + k] + phase_time);
                        }
                        ofs_anime << std::endl;
                    }
                }
            }
            ofs_anime.close();
        }
    }

    deallocate(xmod);
    deallocate(kd_tmp);
    deallocate(eval);
    deallocate(evec);
    deallocate(phase_cell);
    deallocate(evec_mag);
    deallocate(evec_theta);
    deallocate(disp_mag);
    deallocate(mass);
}

void Writes::print_normalmode_borncharge() const
{

    if (mympi->my_rank == 0) {

        auto zstar_born = dielec->get_zstar_mode();

        const auto ns = dynamical->neval;

        std::string file_zstar = input->job_title + ".Born_mode";
        std::ofstream ofs_zstar;
        ofs_zstar.open(file_zstar.c_str(), std::ios::out);
        if (!ofs_zstar)
            exit("print_normalmode_borncharge",
                 "Cannot open file file_zstar");

        ofs_zstar << "# Born effective charges of each phonon mode at q = (0, 0, 0). Unit is (amu)^{-1/2}\n";
        for (auto is = 0; is < ns; ++is) {
            ofs_zstar << "# Mode " << std::setw(5) << is + 1 << '\n';
            ofs_zstar << "#";
            ofs_zstar << std::setw(14) << 'x';
            ofs_zstar << std::setw(15) << 'y';
            ofs_zstar << std::setw(15) << 'z';
            ofs_zstar << '\n';
            for (auto i = 0; i < 3; ++i) {
                ofs_zstar << std::setw(15) << std::fixed << zstar_born[is][i];
            }
            ofs_zstar << "\n\n";
        }
        ofs_zstar.close();
    }
}

void Writes::write_participation_ratio() const
{
    std::string fname_pr, fname_apr;

    if (kpoint->kpoint_general && dynamical->dymat_general) {
        fname_pr = input->job_title + ".pr";
        fname_apr = input->job_title + ".apr";
        write_participation_ratio_each(fname_pr, fname_apr,
                                       kpoint->kpoint_general->nk,
                                       kpoint->kpoint_general->xk,
                                       dynamical->dymat_general->get_eigenvalues(),
                                       dynamical->dymat_general->get_eigenvectors());
    }

    if (kpoint->kpoint_bs && dynamical->dymat_band) {
        fname_pr = input->job_title + ".band.pr";
        fname_apr = input->job_title + ".band.apr";
        write_participation_ratio_each(fname_pr, fname_apr,
                                       kpoint->kpoint_bs->nk,
                                       kpoint->kpoint_bs->xk,
                                       dynamical->dymat_band->get_eigenvalues(),
                                       dynamical->dymat_band->get_eigenvectors());
    }

    if (dos->kmesh_dos && dos->dymat_dos) {
        fname_pr = input->job_title + ".mesh.pr";
        fname_apr = input->job_title + ".mesh.apr";
        write_participation_ratio_mesh(fname_pr, fname_apr,
                                       dos->kmesh_dos,
                                       dos->dymat_dos->get_eigenvalues(),
                                       dos->dymat_dos->get_eigenvectors());
    }
}

void Writes::write_participation_ratio_each(const std::string &fname_pr,
                                            const std::string &fname_apr,
                                            const unsigned int nk_in,
                                            const double *const *xk_in,
                                            const double *const *eval_in,
                                            const std::complex<double> *const *const *evec_in) const
{
    unsigned int i, j, k;
    const auto neval = dynamical->neval;
    const auto natmin = system->natmin;

    double **participation_ratio = nullptr;
    double ***atomic_participation_ratio = nullptr;

    std::ofstream ofs_pr, ofs_apr;

    ofs_pr.open(fname_pr.c_str(), std::ios::out);
    if (!ofs_pr)
        exit("write_participation_ratio",
             "cannot open file_pr");
    ofs_pr.setf(std::ios::scientific);

    ofs_apr.open(fname_apr.c_str(), std::ios::out);
    if (!ofs_apr)
        exit("write_participation_ratio",
             "cannot open file_apr");
    ofs_apr.setf(std::ios::scientific);

    allocate(participation_ratio, nk_in, neval);
    allocate(atomic_participation_ratio, nk_in, neval, natmin);

    dynamical->calc_participation_ratio_all(nk_in, evec_in,
                                            participation_ratio,
                                            atomic_participation_ratio);

    ofs_pr << "# Participation ratio of each phonon modes at k points" << std::endl;
    ofs_pr << "# kpoint, mode, PR[kpoint][mode]" << std::endl;

    for (i = 0; i < nk_in; ++i) {
        ofs_pr << "#" << std::setw(8) << i + 1;
        ofs_pr << " xk = ";
        for (j = 0; j < 3; ++j) {
            ofs_pr << std::setw(15) << xk_in[i][j];
        }
        ofs_pr << std::endl;
        for (j = 0; j < nbands; ++j) {
            ofs_pr << std::setw(8) << i + 1;
            ofs_pr << std::setw(5) << j + 1;
            ofs_pr << std::setw(15) << participation_ratio[i][j];
            ofs_pr << std::endl;
        }
        ofs_pr << std::endl;
    }
    ofs_pr.close();

    ofs_apr << "# Atomic participation ratio of each phonon modes at k points" << std::endl;
    ofs_apr << "# kpoint, mode, atom, APR[kpoint][mode][atom]" << std::endl;

    for (i = 0; i < nk_in; ++i) {
        ofs_apr << "#" << std::setw(8) << i + 1;
        ofs_apr << " xk = ";
        for (j = 0; j < 3; ++j) {
            ofs_apr << std::setw(15) << xk_in[i][j];
        }
        ofs_apr << std::endl;
        for (j = 0; j < nbands; ++j) {
            for (k = 0; k < natmin; ++k) {
                ofs_apr << std::setw(8) << i + 1;
                ofs_apr << std::setw(5) << j + 1;
                ofs_apr << std::setw(5) << k + 1;
                ofs_apr << std::setw(15) << atomic_participation_ratio[i][j][k];
                ofs_apr << std::endl;
            }
        }
        ofs_apr << std::endl;
    }
    ofs_apr.close();

    deallocate(participation_ratio);
    deallocate(atomic_participation_ratio);

    std::cout << "  " << std::setw(input->job_title.length() + 12) << std::left << fname_pr;
    std::cout << " : Participation ratio for all k points" << std::endl;
    std::cout << "  " << std::setw(input->job_title.length() + 12) << std::left << fname_apr;
    std::cout << " : Atomic participation ratio for all k points" << std::endl;
}

void Writes::write_participation_ratio_mesh(const std::string &fname_pr,
                                            const std::string &fname_apr,
                                            const KpointMeshUniform *kmesh_in,
                                            const double *const *eval_in,
                                            const std::complex<double> *const *const *evec_in) const
{
    unsigned int i, j, k;
    unsigned int knum;
    const auto neval = dynamical->neval;
    const auto natmin = system->natmin;
    const auto nk = kmesh_in->nk;

    double **participation_ratio = nullptr;
    double ***atomic_participation_ratio = nullptr;

    std::ofstream ofs_pr, ofs_apr;

    ofs_pr.open(fname_pr.c_str(), std::ios::out);
    if (!ofs_pr)
        exit("write_participation_ratio",
             "cannot open file_pr");
    ofs_pr.setf(std::ios::scientific);

    ofs_apr.open(fname_apr.c_str(), std::ios::out);
    if (!ofs_apr)
        exit("write_participation_ratio",
             "cannot open file_apr");
    ofs_apr.setf(std::ios::scientific);

    allocate(participation_ratio, nk, neval);
    allocate(atomic_participation_ratio, nk, neval, natmin);

    dynamical->calc_participation_ratio_all(nk, evec_in,
                                            participation_ratio,
                                            atomic_participation_ratio);

    ofs_pr << "# Participation ratio of each phonon modes at k points" << std::endl;
    ofs_pr << "# irred. kpoint, mode, frequency[kpoint][mode] (cm^-1), PR[kpoint][mode]" << std::endl;

    for (i = 0; i < kmesh_in->nk_irred; ++i) {
        knum = kmesh_in->kpoint_irred_all[i][0].knum;
        ofs_pr << "#" << std::setw(8) << i + 1;
        ofs_pr << " xk = ";
        for (j = 0; j < 3; ++j) {
            ofs_pr << std::setw(15) << kmesh_in->xk[knum][j];
        }
        ofs_pr << std::endl;
        for (j = 0; j < nbands; ++j) {
            ofs_pr << std::setw(8) << i + 1;
            ofs_pr << std::setw(5) << j + 1;
            ofs_pr << std::setw(15) << in_kayser(eval_in[knum][j]);
            ofs_pr << std::setw(15) << participation_ratio[knum][j];
            ofs_pr << std::endl;
        }
        ofs_pr << std::endl;
    }
    ofs_pr.close();

    ofs_apr << "# Atomic participation ratio of each phonon modes at k points" << std::endl;
    ofs_apr << "# irred. kpoint, mode, atom, frequency[kpoint][mode] (cm^-1), APR[kpoint][mode][atom]" << std::endl;

    for (i = 0; i < kmesh_in->nk_irred; ++i) {
        knum = kmesh_in->kpoint_irred_all[i][0].knum;

        ofs_apr << "#" << std::setw(8) << i + 1;
        ofs_apr << " xk = ";
        for (j = 0; j < 3; ++j) {
            ofs_apr << std::setw(15) << kmesh_in->xk[knum][j];
        }
        ofs_apr << std::endl;
        for (j = 0; j < nbands; ++j) {
            for (k = 0; k < natmin; ++k) {
                ofs_apr << std::setw(8) << i + 1;
                ofs_apr << std::setw(5) << j + 1;
                ofs_apr << std::setw(5) << k + 1;
                ofs_apr << std::setw(15) << in_kayser(eval_in[knum][j]);
                ofs_apr << std::setw(15) << atomic_participation_ratio[knum][j][k];
                ofs_apr << std::endl;
            }
        }
        ofs_apr << std::endl;
    }
    ofs_apr.close();

    deallocate(participation_ratio);
    deallocate(atomic_participation_ratio);

    std::cout << "  " << std::setw(input->job_title.length() + 12) << std::left << fname_pr;
    std::cout << " : Participation ratio for all k points" << std::endl;
    std::cout << "  " << std::setw(input->job_title.length() + 12) << std::left << fname_apr;
    std::cout << " : Atomic participation ratio for all k points" << std::endl;
}

void Writes::write_dielectric_function() const
{
    std::ofstream ofs_dielec;
    auto file_dielec = input->job_title + ".dielec";

    ofs_dielec.open(file_dielec.c_str(), std::ios::out);
    if (!ofs_dielec) exit("write_phonon_vel", "cannot open file_vel");

    unsigned int nomega;
    auto omega_grid = dielec->get_omega_grid(nomega);
    auto dielecfunc = dielec->get_dielectric_func();

    ofs_dielec << "# Real part of dielectric function (phonon part only)\n";
    ofs_dielec << "# Frequency (cm^-1), xx, yy, zz,   xy, xz, yx, yz, zx, zy\n";
    for (auto iomega = 0; iomega < nomega; ++iomega) {
        ofs_dielec << std::setw(10) << omega_grid[iomega];
        for (auto i = 0; i < 3; ++i) {
            ofs_dielec << std::setw(15) << dielecfunc[iomega][i][i];
        }
        for (auto i = 0; i < 3; ++i) {
            for (auto j = 0; j < 3; ++j) {
                if (i == j) continue;
                ofs_dielec << std::setw(15) << dielecfunc[iomega][i][j];
            }
        }
        ofs_dielec << '\n';
    }
    ofs_dielec << std::endl;
    ofs_dielec.close();

    std::cout << "  " << std::setw(input->job_title.length() + 12) << std::left << file_dielec;
    std::cout << " : Frequency-dependent dielectric function" << std::endl;
}

void Writes::write_scph_energy(const unsigned int nk_in,
                               const double *const *const *eval_in,
                               const int bubble) const
{
    const auto ns = dynamical->neval;
    const auto Tmin = system->Tmin;
    const auto Tmax = system->Tmax;
    const auto dT = system->dT;
    const auto NT = static_cast<unsigned int>((Tmax - Tmin) / dT) + 1;

    std::ofstream ofs_energy;
    std::string file_energy;

    if (bubble == 0) {
        file_energy = input->job_title + ".scph_eval";
    } else if (bubble == 1) {
        file_energy = input->job_title + ".scph+bubble(0)_eval";
    } else if (bubble == 2) {
        file_energy = input->job_title + ".scph+bubble(w)_eval";
    }

    ofs_energy.open(file_energy.c_str(), std::ios::out);
    if (!ofs_energy) exit("write_scph_energy", "cannot open file_energy");

    ofs_energy << "# K point, mode, Temperature [K], Eigenvalues [cm^-1]" << std::endl;

    for (unsigned int ik = 0; ik < nk_in; ++ik) {
        for (unsigned int is = 0; is < ns; ++is) {
            for (unsigned int iT = 0; iT < NT; ++iT) {
                const auto temp = Tmin + static_cast<double>(iT) * dT;

                ofs_energy << std::setw(5) << ik + 1;
                ofs_energy << std::setw(5) << is + 1;
                ofs_energy << std::setw(8) << temp;
                ofs_energy << std::setw(15) << writes->in_kayser(eval_in[iT][ik][is]);
                ofs_energy << std::endl;
            }
            ofs_energy << std::endl;
        }
        ofs_energy << std::endl;
    }

    ofs_energy.close();
}

void Writes::write_scph_bands(const unsigned int nk_in,
                              const double *kaxis_in,
                              const double *const *const *eval,
                              const int bubble) const
{
    std::ofstream ofs_bands;
    std::string file_bands;

    if (bubble == 0) {
        file_bands = input->job_title + ".scph_bands";
    } else if (bubble == 1) {
        file_bands = input->job_title + ".scph+bubble(0)_bands";
    } else if (bubble == 2) {
        file_bands = input->job_title + ".scph+bubble(w)_bands";
    } else if (bubble == 3) {
        file_bands = input->job_title + ".scph+bubble(wQP)_bands";
    }

    ofs_bands.open(file_bands.c_str(), std::ios::out);
    if (!ofs_bands) exit("write_scph_bands", "cannot open file_bands");

    unsigned int i;
    const auto Tmin = system->Tmin;
    const auto Tmax = system->Tmax;
    const auto dT = system->dT;
    const auto NT = static_cast<unsigned int>((Tmax - Tmin) / dT) + 1;
    const auto ns = dynamical->neval;
    auto kcount = 0;

    std::string str_tmp = "NONE";
    std::string str_kpath;
    std::string str_kval;

    for (i = 0; i < kpoint->kpInp.size(); ++i) {
        if (str_tmp != kpoint->kpInp[i].kpelem[0]) {
            str_tmp = kpoint->kpInp[i].kpelem[0];
            str_kpath += " " + str_tmp;

            std::ostringstream ss;
            ss << std::fixed << std::setprecision(6) << kaxis_in[kcount];
            str_kval += " " + ss.str();
        }
        kcount += std::atoi(kpoint->kpInp[i].kpelem[8].c_str());

        if (str_tmp != kpoint->kpInp[i].kpelem[4]) {
            str_tmp = kpoint->kpInp[i].kpelem[4];
            str_kpath += " " + str_tmp;

            std::ostringstream ss;
            ss << std::fixed << std::setprecision(6) << kaxis_in[kcount - 1];
            str_kval += " " + ss.str();
        }
    }

    ofs_bands << "# " << str_kpath << std::endl;
    ofs_bands << "#" << str_kval << std::endl;
    ofs_bands << "# Temperature [K], k-axis, Eigenvalues [cm^-1]" << std::endl;

    for (unsigned int iT = 0; iT < NT; ++iT) {
        const auto temp = Tmin + static_cast<double>(iT) * dT;

        for (i = 0; i < nk_in; ++i) {
            ofs_bands << std::setw(15) << std::fixed << temp;
            ofs_bands << std::setw(15) << std::fixed << kaxis_in[i];
            for (unsigned int j = 0; j < ns; ++j) {
                ofs_bands << std::setw(15) << std::scientific << writes->in_kayser(eval[iT][i][j]);
            }
            ofs_bands << std::endl;
        }
        ofs_bands << std::endl;
    }

    ofs_bands.close();
    std::cout << "  " << std::setw(input->job_title.length() + 12) << std::left << file_bands;
    if (bubble == 0) {
        std::cout << " : SCPH band structure" << std::endl;
    } else if (bubble == 1) {
        std::cout << " : SCPH+Bubble(0) band structure" << std::endl;
    } else if (bubble == 2) {
        std::cout << " : SCPH+Bubble(w) band structure" << std::endl;
    } else if (bubble == 3) {
        std::cout << " : SCPH+Bubble(wQP) band structure" << std::endl;
    }
}

void Writes::write_scph_dos(double **dos_scph, const int bubble) const
{
    unsigned int iT;
    const auto Tmin = system->Tmin;
    const auto Tmax = system->Tmax;
    const auto dT = system->dT;
    const auto NT = static_cast<unsigned int>((Tmax - Tmin) / dT) + 1;

    std::ofstream ofs_dos;
    std::string file_dos;

    if (bubble == 0) {
        file_dos = input->job_title + ".scph_dos";
    } else if (bubble == 1) {
        file_dos = input->job_title + ".scph+bubble(0)_dos";
    } else if (bubble == 2) {
        file_dos = input->job_title + ".scph+bubble(w)_dos";
    }

    ofs_dos.open(file_dos.c_str(), std::ios::out);
    if (!ofs_dos) exit("write_scph_dos", "cannot open file_dos");

    ofs_dos << "# ";

    for (iT = 0; iT < NT; ++iT) {
        ofs_dos << std::setw(15) << Tmin + static_cast<double>(iT) * dT;
    }
    ofs_dos << std::endl;

    for (unsigned int j = 0; j < dos->n_energy; ++j) {
        ofs_dos << std::setw(15) << dos->energy_dos[j];

        for (iT = 0; iT < NT; ++iT) {
            ofs_dos << std::setw(15) << dos_scph[iT][j];
        }
        ofs_dos << std::endl;
    }

    ofs_dos << std::endl;
    ofs_dos.close();
    std::cout << "  " << std::setw(input->job_title.length() + 12) << std::left << file_dos;
    if (bubble == 0) {
        std::cout << " : SCPH DOS" << std::endl;
    } else if (bubble == 1) {
        std::cout << " : SCPH+Bubble(0) DOS" << std::endl;
    } else if (bubble == 2) {
        std::cout << " : SCPH+Bubble(w) DOS" << std::endl;
    }
}

void Writes::write_scph_thermodynamics(double *heat_capacity,
                                       double *heat_capacity_correction,
                                       double *FE_QHA,
                                       double *dFE_scph) const
{
    const auto Tmin = system->Tmin;
    const auto Tmax = system->Tmax;
    const auto dT = system->dT;
    const auto NT = static_cast<unsigned int>((Tmax - Tmin) / dT) + 1;

    bool print_anharmonic_correction_Cv = false;

    if (heat_capacity_correction) {
        print_anharmonic_correction_Cv = true;
    }

    std::ofstream ofs_thermo;
    auto file_thermo = input->job_title + ".scph_thermo";
    ofs_thermo.open(file_thermo.c_str(), std::ios::out);
    if (!ofs_thermo)
        exit("write_scph_thermodynamics",
             "cannot open file_thermo");

    if (thermodynamics->calc_FE_bubble) {
        ofs_thermo << "# The bubble free-energy calculated on top of the SCPH wavefunction is also shown." << std::endl;
        ofs_thermo <<
                   "# Temperature [K], Cv [in kB unit], F_{vib} (QHA term) [Ry], F_{vib} (SCPH correction) [Ry], F_{vib} (Bubble correction) [Ry]"
                   << std::endl;
    } else {
        ofs_thermo << "# Temperature [K], Cv [in kB unit], F_{vib} (QHA term) [Ry], F_{vib} (SCPH correction) [Ry]"
                   << std::endl;
    }

    if (thermodynamics->classical) {
        ofs_thermo << "# CLASSICAL = 1: Use classical limit." << std::endl;
    }

    for (unsigned int iT = 0; iT < NT; ++iT) {

        const auto temp = Tmin + static_cast<double>(iT) * dT;

        ofs_thermo << std::setw(16) << std::fixed << temp;
        ofs_thermo << std::setw(18) << std::scientific << heat_capacity[iT] / k_Boltzmann;
        if (print_anharmonic_correction_Cv) {
            ofs_thermo << std::setw(18) << std::scientific << heat_capacity_correction[iT] / k_Boltzmann;
        }
        ofs_thermo << std::setw(18) << FE_QHA[iT];
        ofs_thermo << std::setw(18) << dFE_scph[iT];

        if (thermodynamics->calc_FE_bubble) {
            ofs_thermo << std::setw(18) << thermodynamics->FE_bubble[iT];
        }
        ofs_thermo << std::endl;
    }

    ofs_thermo.close();
    std::cout << "  " << std::setw(input->job_title.length() + 12) << std::left << file_thermo;
    std::cout << " : SCPH heat capcaity and free energy" << std::endl;
}

void Writes::write_scph_dielec(double ****dielec_scph) const
{
    const auto Tmin = system->Tmin;
    const auto Tmax = system->Tmax;
    const auto dT = system->dT;
    const auto NT = static_cast<unsigned int>((Tmax - Tmin) / dT) + 1;

    std::ofstream ofs_dielec;
    auto file_dielec = input->job_title + ".scph_dielec";

    ofs_dielec.open(file_dielec.c_str(), std::ios::out);
    if (!ofs_dielec) exit("write_scph_dielec", "cannot open PREFIX.scph_dielec");

    unsigned int nomega;
    auto omega_grid = dielec->get_omega_grid(nomega);

    ofs_dielec << "# Real part of dielectric function (phonon part only)\n";
    ofs_dielec << "# Temperature (K), Frequency (cm^-1), xx, yy, zz\n";

    for (unsigned int iT = 0; iT < NT; ++iT) {

        const auto temp = Tmin + static_cast<double>(iT) * dT;

        for (auto iomega = 0; iomega < nomega; ++iomega) {
            ofs_dielec << std::setw(16) << std::fixed << temp;
            ofs_dielec << std::setw(15) << std::scientific << omega_grid[iomega];
            for (auto i = 0; i < 3; ++i) {
                ofs_dielec << std::setw(15) << dielec_scph[iT][iomega][i][i];
            }
            ofs_dielec << '\n';
        }
        ofs_dielec << std::endl;
    }

    ofs_dielec.close();

    std::cout << "  " << std::setw(input->job_title.length() + 12) << std::left << file_dielec;
    std::cout << " : SCPH frequency-dependent dielectric function" << std::endl;
}

unsigned int Writes::getVerbosity() const
{
    return verbosity;
}

void Writes::setVerbosity(unsigned int verbosity_in)
{
    verbosity = verbosity_in;
}
