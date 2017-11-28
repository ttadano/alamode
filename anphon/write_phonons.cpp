/*
write_phonons.cpp

Copyright (c) 2014, 2015, 2016 Terumasa Tadano

This file is distributed under the terms of the MIT license.
Please see the file 'LICENCE.txt' in the root directory 
or http://opensource.org/licenses/mit-license.php for information.
*/

#include "mpi_common.h"
#include <fstream>
#include <iomanip>
#include <sys/stat.h>
#include "constants.h"
#include "conductivity.h"
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
#include "relaxation.h"
#include "symmetry_core.h"
#include "system.h"
#include "write_phonons.h"
#include "mathfunctions.h"
#include "isotope.h"
#include "integration.h"
#include "scph.h"
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/lexical_cast.hpp>

using namespace PHON_NS;

Writes::Writes(PHON *phon): Pointers(phon)
{
    Ry_to_kayser = Hz_to_kayser / time_ry;
};

Writes::~Writes()
{
};

void Writes::write_input_vars()
{
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
    for (i = 0; i < system->nkd; ++i) {
        std::cout << std::setw(10) << system->mass_kd[i];
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
        std::cout << "  RESTART = " << phon->restart_flag << std::endl;
        std::cout << "  TRISYM = " << relaxation->use_triplet_symmetry << std::endl;
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
            std::cout << "  PRINTMSD = " << writes->print_msd << std::endl;
            std::cout << "  SPS = " << dos->scattering_phase_space << std::endl;
            std::cout << std::endl;
        }
        std::cout << "  GRUNEISEN = " << gruneisen->print_gruneisen << std::endl;
        std::cout << "  NEWFCS = " << gruneisen->print_newfcs;
        if (gruneisen->print_newfcs) {
            std::cout << "; DELTA_A = " << gruneisen->delta_a << std::endl;
            std::cout << "  QUARTIC = " << relaxation->quartic_mode;
        }
        std::cout << std::endl;

    } else if (phon->mode == "RTA") {
        std::cout << "  ISOTOPE = " << isotope->include_isotope << std::endl;
        if (isotope->include_isotope) {
            std::cout << "  ISOFACT = ";
            for (i = 0; i < system->nkd; ++i) {
                std::cout << std::scientific
                    << std::setw(13) << isotope->isotope_factor[i];
            }
            std::cout << std::endl;
        }

        std::cout << "  KAPPA_SPEC = " << conductivity->calc_kappa_spec << std::endl;

        //        std::cout << "  KS_INPUT = " << relaxation->ks_input << std::endl;
        //        std::cout << "  QUARTIC = " << relaxation->quartic_mode << std::endl;
        // std::cout << "  REALPART = " << relaxation->calc_realpart << std::endl;
        // std::cout << "  ATOMPROJ = " << relaxation->atom_project_mode << std::endl;
        // std::cout << "  FSTATE_W = " << relaxation->calc_fstate_omega << std::endl;
        //  std::cout << "  FSTATE_K = " << relaxation->calc_fstate_k << std::endl;

    } else if (phon->mode == "SCPH") {


    } else {
        error->exit("write_input_vars", "This cannot happen");
    }

    std::cout << std::endl << std::endl;
    std::cout << " -----------------------------------------------------------------" << std::endl;
    std::cout << std::endl;
}

void Writes::setup_result_io()
{
    if (mympi->my_rank == 0) {

        if (phon->restart_flag) {

            std::cout << " RESTART = 1 : Restart from the interrupted run." << std::endl;
            std::cout << "               Phonon lifetimes will be load from file " << file_result << std::endl;
            std::cout << "               and check the consistency of the computatioal settings." << std::endl;

            // Restart
            fs_result.open(file_result.c_str(), std::ios::in | std::ios::out);
            if (!fs_result) {
                error->exit("setup_result_io",
                            "Could not open file_result");
            }

            // Check the consistency

            std::string line_tmp, str_tmp;
            int natmin_tmp, nkd_tmp;
            int nk_tmp[3], nksym_tmp;
            int ismear, is_classical;
            double epsilon_tmp, T1, T2, delta_T;


            bool found_tag;

            found_tag = false;
            while (fs_result >> line_tmp) {
                if (line_tmp == "#SYSTEM") {
                    found_tag = true;
                    break;
                }
            }
            if (!found_tag)
                error->exit("setup_result_io",
                            "Could not find #SYSTEM tag");

            fs_result >> natmin_tmp >> nkd_tmp;

            if (!(natmin_tmp == system->natmin && nkd_tmp == system->nkd)) {
                error->exit("setup_result_io",
                            "SYSTEM information is not consistent");
            }

            found_tag = false;
            while (fs_result >> line_tmp) {
                if (line_tmp == "#KPOINT") {
                    found_tag = true;
                    break;
                }
            }
            if (!found_tag)
                error->exit("setup_result_io",
                            "Could not find #KPOINT tag");

            fs_result >> nk_tmp[0] >> nk_tmp[1] >> nk_tmp[2];
            fs_result >> nksym_tmp;

            if (!(kpoint->nkx == nk_tmp[0] &&
                kpoint->nky == nk_tmp[1] &&
                kpoint->nkz == nk_tmp[2] &&
                kpoint->nk_reduced == nksym_tmp)) {
                error->exit("setup_result_io",
                            "KPOINT information is not consistent");
            }

            found_tag = false;
            while (fs_result >> line_tmp) {
                if (line_tmp == "#CLASSICAL") {
                    found_tag = true;
                    break;
                }
            }
            if (!found_tag) {
                std::cout << " Could not find the #CLASSICAL tag in the restart file." << std::endl;
                std::cout << " CLASSIACAL = 0 is assumed." << std::endl;
                is_classical = 0;
            } else {
                fs_result >> is_classical;
            }
            if (static_cast<bool>(is_classical) != thermodynamics->classical) {
                error->warn("setup_result_io",
                            "CLASSICAL val is not consistent");
            }

            found_tag = false;
            while (fs_result >> line_tmp) {
                if (line_tmp == "#FCSXML") {
                    found_tag = true;
                    break;
                }
            }
            if (!found_tag)
                error->exit("setup_result_io",
                            "Could not find #FCSXML tag");

            fs_result >> str_tmp;
            if (str_tmp != fcs_phonon->file_fcs) {
                error->warn("setup_result_io",
                            "FCSXML is not consistent");
            }

            found_tag = false;
            while (fs_result >> line_tmp) {
                if (line_tmp == "#SMEARING") {
                    found_tag = true;
                    break;
                }
            }
            if (!found_tag)
                error->exit("setup_result_io",
                            "Could not find #SMEARING tag");

            fs_result >> ismear;
            fs_result >> epsilon_tmp;

            if (ismear != integration->ismear) {
                error->exit("setup_result_io",
                            "Smearing method is not consistent");
            }
            if ((ismear != -1) && (std::abs(epsilon_tmp - integration->epsilon * Ry_to_kayser) >= eps4)) {
                std::cout << "epsilon from file : " << std::setw(15)
                    << std::setprecision(10) << epsilon_tmp * Ry_to_kayser << std::endl;
                std::cout << "epsilon from input: " << std::setw(15)
                    << std::setprecision(10) << integration->epsilon * Ry_to_kayser << std::endl;
                error->exit("setup_result_io",
                            "Smearing width is not consistent");
            }

            found_tag = false;
            while (fs_result >> line_tmp) {
                if (line_tmp == "#TEMPERATURE") {
                    found_tag = true;
                    break;
                }
            }
            if (!found_tag)
                error->exit("setup_result_io",
                            "Could not find #TEMPERATURE tag");

            fs_result >> T1 >> T2 >> delta_T;

            if (!(T1 == system->Tmin &&
                T2 == system->Tmax &&
                delta_T == system->dT)) {
                error->exit("setup_result_io",
                            "Temperature information is not consistent");
            }

        } else {
            // From scratch
            fs_result.open(file_result.c_str(), std::ios::out);
            if (!fs_result) {
                error->exit("setup_result_io",
                            "Could not open file_result");
            }

            fs_result << "## General information" << std::endl;
            fs_result << "#SYSTEM" << std::endl;
            fs_result << system->natmin << " " << system->nkd << std::endl;
            fs_result << system->volume_p << std::endl;
            fs_result << "#END SYSTEM" << std::endl;

            fs_result << "#KPOINT" << std::endl;
            fs_result << kpoint->nkx << " " << kpoint->nky << " " << kpoint->nkz << std::endl;
            fs_result << kpoint->nk_reduced << std::endl;

            for (int i = 0; i < kpoint->nk_reduced; ++i) {
                fs_result << std::setw(6) << i + 1 << ":";
                for (int j = 0; j < 3; ++j) {
                    fs_result << std::setw(15)
                        << std::scientific << kpoint->kpoint_irred_all[i][0].kval[j];
                }
                fs_result << std::setw(12)
                    << std::fixed << kpoint->weight_k[i] << std::endl;
            }
            fs_result.unsetf(std::ios::fixed);

            fs_result << "#END KPOINT" << std::endl;

            fs_result << "#CLASSICAL" << std::endl;
            fs_result << thermodynamics->classical << std::endl;
            fs_result << "#END CLASSICAL" << std::endl;

            fs_result << "#FCSXML" << std::endl;
            fs_result << fcs_phonon->file_fcs << std::endl;
            fs_result << "#END  FCSXML" << std::endl;

            fs_result << "#SMEARING" << std::endl;
            fs_result << integration->ismear << std::endl;
            fs_result << integration->epsilon * Ry_to_kayser << std::endl;
            fs_result << "#END SMEARING" << std::endl;

            fs_result << "#TEMPERATURE" << std::endl;
            fs_result << system->Tmin << " " << system->Tmax << " " << system->dT << std::endl;
            fs_result << "#END TEMPERATURE" << std::endl;

            fs_result << "##END General information" << std::endl;
        }
    }
}

void Writes::print_phonon_energy()
{
    unsigned int i;
    unsigned int ik, is;
    unsigned int nk = kpoint->nk;
    unsigned int ns = dynamical->neval;
    unsigned int knum;

    double kayser_to_THz = 0.0299792458;

    std::cout << std::endl;
    std::cout << " -----------------------------------------------------------------" << std::endl << std::endl;
    std::cout << " Phonon frequencies below:" << std::endl << std::endl;

    if (kpoint->kpoint_mode == 0 || kpoint->kpoint_mode == 1) {

        for (ik = 0; ik < nk; ++ik) {
            std::cout << " # k point " << std::setw(5) << ik + 1;
            std::cout << " : (";

            for (i = 0; i < 3; ++i) {
                std::cout << std::fixed << std::setprecision(4)
                    << std::setw(8) << kpoint->xk[ik][i];
                if (i < 2) std::cout << ",";
            }
            std::cout << ")" << std::endl;

            std::cout << "   Mode, Frequency " << std::endl;

            for (is = 0; is < ns; ++is) {
                std::cout << std::setw(7) << is + 1;
                std::cout << std::fixed << std::setprecision(4) << std::setw(12)
                    << in_kayser(dynamical->eval_phonon[ik][is]);
                std::cout << " cm^-1  (";
                std::cout << std::fixed << std::setprecision(4) << std::setw(12)
                    << kayser_to_THz * in_kayser(dynamical->eval_phonon[ik][is]);
                std::cout << " THz )" << std::endl;
            }
            std::cout << std::endl;
        }

    } else if (kpoint->kpoint_mode == 2) {

        for (ik = 0; ik < kpoint->kpoint_irred_all.size(); ++ik) {

            std::cout << " # Irred. k point" << std::setw(5) << ik + 1;
            std::cout << " : (";

            for (i = 0; i < 3; ++i) {
                std::cout << std::fixed << std::setprecision(4) << std::setw(8)
                    << kpoint->kpoint_irred_all[ik][0].kval[i];
                if (i < 2) std::cout << ",";
            }
            std::cout << ")" << std::endl;

            std::cout << "   Mode, Frequency " << std::endl;

            knum = kpoint->kpoint_irred_all[ik][0].knum;

            for (is = 0; is < ns; ++is) {
                std::cout << std::setw(7) << is + 1;
                std::cout << std::fixed << std::setprecision(4) << std::setw(12)
                    << in_kayser(dynamical->eval_phonon[knum][is]);
                std::cout << " cm^-1  (";
                std::cout << std::fixed << std::setprecision(4) << std::setw(12)
                    << kayser_to_THz * in_kayser(dynamical->eval_phonon[knum][is]);
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
        if (kpoint->kpoint_mode == 1) {
            write_phonon_vel();
        } else if (kpoint->kpoint_mode == 2) {
            write_phonon_vel_all();
        }
    }

    if (dos->flag_dos) {
        write_phonon_dos();

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
    }

    if (print_xsf) {
        write_normal_mode_direction();
    }

    if (dynamical->print_eigenvectors) {
        write_eigenvectors();
    }

    if (dynamical->participation_ratio) {
        write_participation_ratio();
    }

    if (gruneisen->print_gruneisen) {
        write_gruneisen();
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
}

void Writes::write_phonon_bands()
{
    std::ofstream ofs_bands;
    std::string file_bands = input->job_title + ".bands";

    ofs_bands.open(file_bands.c_str(), std::ios::out);
    if (!ofs_bands)
        error->exit("write_phonon_bands",
                    "cannot open file_bands");

    unsigned int i, j;

    unsigned int nk = kpoint->nk;

    double *kaxis = kpoint->kaxis;
    double **eval = dynamical->eval_phonon;

    int kcount = 0;

    std::string str_tmp = "NONE";
    std::string str_kpath = "";
    std::string str_kval = "";

    for (i = 0; i < kpoint->kpInp.size(); ++i) {
        if (str_tmp != kpoint->kpInp[i].kpelem[0]) {
            str_tmp = kpoint->kpInp[i].kpelem[0];
            str_kpath += " " + str_tmp;

            std::ostringstream ss;
            ss << std::fixed << std::setprecision(6) << kpoint->kaxis[kcount];
            str_kval += " " + ss.str();
        }
        kcount += std::atoi(kpoint->kpInp[i].kpelem[8].c_str());

        if (str_tmp != kpoint->kpInp[i].kpelem[4]) {
            str_tmp = kpoint->kpInp[i].kpelem[4];
            str_kpath += " " + str_tmp;

            std::ostringstream ss;
            ss << std::fixed << std::setprecision(6) << kpoint->kaxis[kcount - 1];
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
        std::string file_connect = input->job_title + ".connection";

        ofs_connect.open(file_connect.c_str(), std::ios::out);
        if (!ofs_connect)
            error->exit("write_phonon_bands",
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

void Writes::write_phonon_vel()
{
    std::ofstream ofs_vel;
    std::string file_vel = input->job_title + ".phvel";

    ofs_vel.open(file_vel.c_str(), std::ios::out);
    if (!ofs_vel) error->exit("write_phonon_vel", "cannot open file_vel");

    unsigned int i, j;
    unsigned int nk = kpoint->nk;

    double *kaxis = kpoint->kaxis;
    double Ry_to_SI_vel = Bohr_in_Angstrom * 1.0e-10 / time_ry;

    ofs_vel << "# k-axis, |Velocity| [m / sec]" << std::endl;
    ofs_vel.setf(std::ios::fixed);

    for (i = 0; i < nk; ++i) {
        ofs_vel << std::setw(8) << kaxis[i];
        for (j = 0; j < nbands; ++j) {
            ofs_vel << std::setw(15)
                << std::abs(phonon_velocity->phvel[i][j] * Ry_to_SI_vel);
        }
        ofs_vel << std::endl;
    }

    ofs_vel.close();

    std::cout << "  " << std::setw(input->job_title.length() + 12) << std::left << file_vel;
    std::cout << " : Phonon velocity along given k path" << std::endl;
}

void Writes::write_phonon_vel_all()
{
    std::ofstream ofs_vel;
    std::string file_vel = input->job_title + ".phvel_all";

    ofs_vel.open(file_vel.c_str(), std::ios::out);
    if (!ofs_vel) error->exit("write_phonon_vel_all", "cannot open file_vel_all");

    unsigned int i, j, k, ii;
    unsigned int knum;
    unsigned int nk = kpoint->nk;
    unsigned int ns = dynamical->neval;

    double Ry_to_SI_vel = Bohr_in_Angstrom * 1.0e-10 / time_ry;
    double **eval = dynamical->eval_phonon;

    ofs_vel << "# Phonon group velocity at all reducible k points." << std::endl;
    ofs_vel << "# irred. knum, knum, mode num, frequency [cm^-1], "
        "|velocity| [m/sec], velocity_(x,y,z) [m/sec]" << std::endl << std::endl;
    ofs_vel.setf(std::ios::fixed);

    for (i = 0; i < kpoint->nk_reduced; ++i) {
        ofs_vel << "# Irreducible k point  : " << std::setw(8) << i + 1;
        ofs_vel << " (" << std::setw(4) << kpoint->kpoint_irred_all[i].size() << ")" << std::endl;

        for (j = 0; j < kpoint->kpoint_irred_all[i].size(); ++j) {
            knum = kpoint->kpoint_irred_all[i][j].knum;

            ofs_vel << "## xk =    ";
            for (k = 0; k < 3; ++k)
                ofs_vel << std::setw(15) << std::fixed
                    << std::setprecision(10) << kpoint->xk[knum][k];
            ofs_vel << std::endl;

            for (k = 0; k < ns; ++k) {
                ofs_vel << std::setw(7) << i + 1;
                ofs_vel << std::setw(8) << knum + 1;
                ofs_vel << std::setw(5) << k + 1;
                ofs_vel << std::setw(10) << std::fixed
                    << std::setprecision(2) << in_kayser(eval[knum][k]);
                ofs_vel << std::setw(10) << std::fixed
                    << std::setprecision(2) << phonon_velocity->phvel[knum][k] * Ry_to_SI_vel;
                for (ii = 0; ii < 3; ++ii) {
                    ofs_vel << std::setw(10) << std::fixed << std::setprecision(2)
                        << phonon_velocity->phvel_xyz[knum][k][ii] * Ry_to_SI_vel;
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
}

void Writes::write_phonon_dos()
{
    int i, iat;
    std::ofstream ofs_dos;
    std::string file_dos = input->job_title + ".dos";

    ofs_dos.open(file_dos.c_str(), std::ios::out);
    if (!ofs_dos) error->exit("write_phonon_dos", "cannot open file_dos");

    ofs_dos << "#";
    for (i = 0; i < system->nkd; ++i) {
        ofs_dos << std::setw(5) << system->symbol_kd[i];
    }
    ofs_dos << std::endl;
    ofs_dos << "#";

    unsigned int *nat_each_kd;
    memory->allocate(nat_each_kd, system->nkd);
    for (i = 0; i < system->nkd; ++i) nat_each_kd[i] = 0;
    for (i = 0; i < system->natmin; ++i) {
        ++nat_each_kd[system->kd[system->map_p2s[i][0]]];
    }
    for (i = 0; i < system->nkd; ++i) {
        ofs_dos << std::setw(5) << nat_each_kd[i];
    }
    ofs_dos << std::endl;
    memory->deallocate(nat_each_kd);

    ofs_dos << "# Energy [cm^-1], TOTAL-DOS";
    if (dos->projected_dos) {
        ofs_dos << ", Atom Projected-DOS";
    }
    ofs_dos << std::endl;
    ofs_dos.setf(std::ios::scientific);

    for (i = 0; i < dos->n_energy; ++i) {
        ofs_dos << std::setw(15) << dos->energy_dos[i]
            << std::setw(15) << dos->dos_phonon[i];
        if (dos->projected_dos) {
            for (iat = 0; iat < system->natmin; ++iat) {
                ofs_dos << std::setw(15) << dos->pdos_phonon[iat][i];
            }
        }
        ofs_dos << std::endl;
    }
    ofs_dos.close();


    std::cout << "  " << std::setw(input->job_title.length() + 12) << std::left << file_dos;

    if (dos->projected_dos) {
        std::cout << " : Phonon DOS and atom projected DOS" << std::endl;
    } else {
        std::cout << " : Phonon DOS" << std::endl;
    }
}

void Writes::write_two_phonon_dos()
{
    int i, j;
    int ik;

    std::string file_tdos;
    std::ofstream ofs_tdos;

    file_tdos = input->job_title + ".tdos";
    ofs_tdos.open(file_tdos.c_str(), std::ios::out);

    ofs_tdos << "# Two-phonon DOS (TDOS) for all irreducible k points. " << std::endl;
    ofs_tdos << "# Energy [cm^-1], emission delta(e-e1-e2), absorption delta (e-e1+e2)" << std::endl;

    int n = dos->n_energy;

    for (ik = 0; ik < kpoint->nk_reduced; ++ik) {

        ofs_tdos << "# Irred. kpoint : " << std::setw(5) << ik + 1 << std::endl;
        for (i = 0; i < n; ++i) {
            ofs_tdos << std::setw(15) << dos->emin + dos->delta_e * static_cast<double>(i);

            for (j = 0; j < 2; ++j) ofs_tdos << std::setw(15) << dos->dos2_phonon[ik][i][j];
            ofs_tdos << std::endl;
        }
        ofs_tdos << std::endl;
    }

    ofs_tdos.close();

    std::cout << "  " << std::setw(input->job_title.length() + 12) << std::left << file_tdos;
    std::cout << " : Two-phonon DOS" << std::endl;
}

void Writes::write_scattering_phase_space()
{
    int ik, is;
    unsigned int knum;

    std::string file_sps;
    std::ofstream ofs_sps;

    file_sps = input->job_title + ".sps";
    ofs_sps.open(file_sps.c_str(), std::ios::out);

    ofs_sps << "# Total scattering phase space (cm): "
        << std::scientific << dos->total_sps3 << std::endl;
    ofs_sps << "# Mode decomposed scattering phase space are printed below." << std::endl;
    ofs_sps << "# Irred. k, mode, omega (cm^-1), P+ (absorption) (cm), P- (emission) (cm)" << std::endl;

    for (ik = 0; ik < kpoint->nk_reduced; ++ik) {
        knum = kpoint->kpoint_irred_all[ik][0].knum;

        for (is = 0; is < dynamical->neval; ++is) {
            ofs_sps << std::setw(5) << ik + 1;
            ofs_sps << std::setw(5) << is + 1;
            ofs_sps << std::setw(15) << in_kayser(dynamical->eval_phonon[knum][is]);
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


void Writes::write_scattering_amplitude()
{
    int i, j;
    unsigned int knum;
    unsigned int is;
    unsigned int ns = dynamical->neval;

    std::string file_w = input->job_title + ".sps_Bose";
    std::ofstream ofs_w;

    double Tmin = system->Tmin;
    double Tmax = system->Tmax;
    double dT = system->dT;
    unsigned int NT = static_cast<unsigned int>((Tmax - Tmin) / dT) + 1;

    double omega;

    ofs_w.open(file_w.c_str(), std::ios::out);

    ofs_w << "# Scattering phase space with the Bose-Einstein distribution function" << std::endl;
    ofs_w << "# Irreducible kpoints " << std::endl;
    for (i = 0; i < kpoint->kpoint_irred_all.size(); ++i) {
        ofs_w << "#" << std::setw(5) << i + 1;

        knum = kpoint->kpoint_irred_all[i][0].knum;
        for (j = 0; j < 3; ++j) ofs_w << std::setw(15) << kpoint->xk[knum][j];
        ofs_w << std::endl;
    }
    ofs_w << std::endl;
    ofs_w << "# k, mode, frequency (cm^-1), temperature, W+ (absorption) (cm), W- (emission) (cm)"
        << std::endl << std::endl;

    for (i = 0; i < kpoint->kpoint_irred_all.size(); ++i) {

        knum = kpoint->kpoint_irred_all[i][0].knum;

        for (is = 0; is < ns; ++is) {

            omega = in_kayser(dynamical->eval_phonon[knum][is]);

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

void Writes::write_normal_mode_direction()
{
    std::ofstream ofs_anime;
    std::string file_anime = input->job_title + ".axsf";

    ofs_anime.open(file_anime.c_str(), std::ios::out);
    if (!ofs_anime) error->exit("write_mode_anime", "cannot open file_anime");

    ofs_anime.setf(std::ios::scientific);

    unsigned int i, j, k;
    unsigned int natmin = system->natmin;
    unsigned int nk = kpoint->nk;

    double force_factor = 100.0;

    double **xmod;
    std::string *kd_tmp;

    memory->allocate(xmod, natmin, 3);
    memory->allocate(kd_tmp, natmin);

    ofs_anime << "ANIMSTEPS " << nbands * nk << std::endl;
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

    unsigned int ik, imode;
    double norm;
    std::complex<double> evec_tmp;
    unsigned int m;
    i = 0;

    for (ik = 0; ik < nk; ++ik) {
        for (imode = 0; imode < nbands; ++imode) {
            ofs_anime << "PRIMCOORD " << std::setw(10) << i + 1 << std::endl;
            ofs_anime << std::setw(10) << natmin << std::setw(10) << 1 << std::endl;
            norm = 0.0;

            for (j = 0; j < 3 * natmin; ++j) {
                evec_tmp = dynamical->evec_phonon[ik][imode][j];
                norm += std::pow(evec_tmp.real(), 2) + std::pow(evec_tmp.imag(), 2);
            }

            norm *= force_factor / static_cast<double>(natmin);

            for (j = 0; j < natmin; ++j) {

                m = system->map_p2s[j][0];

                ofs_anime << std::setw(10) << kd_tmp[j];

                for (k = 0; k < 3; ++k) {
                    ofs_anime << std::setw(15) << xmod[j][k];
                }
                for (k = 0; k < 3; ++k) {
                    ofs_anime << std::setw(15)
                        << dynamical->evec_phonon[ik][imode][3 * j + k].real()
                        / (std::sqrt(system->mass[m]) * norm);
                }
                ofs_anime << std::endl;
            }

            ++i;
        }
    }

    memory->deallocate(xmod);
    memory->deallocate(kd_tmp);

    ofs_anime.close();
    std::cout << "  " << std::setw(input->job_title.length() + 12) << std::left << file_anime;
    std::cout << " : XcrysDen AXSF file to visualize phonon mode directions" << std::endl;
}

void Writes::write_eigenvectors()
{
    unsigned int i, j, k;
    unsigned int nk = kpoint->nk;
    unsigned int neval = dynamical->neval;
    double omega2;
    std::ofstream ofs_evec;
    std::string file_evec = input->job_title + ".evec";

    ofs_evec.open(file_evec.c_str(), std::ios::out);
    if (!ofs_evec) error->exit("write_eigenvectors", "cannot open file_evec");
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
    ofs_evec << "# Number of k points : " << std::setw(10) << nk << std::endl << std::endl;
    ofs_evec << "# Eigenvalues and eigenvectors for each phonon modes below:" << std::endl << std::endl;

    for (i = 0; i < nk; ++i) {
        ofs_evec << "## kpoint " << std::setw(7) << i + 1 << " : ";
        for (j = 0; j < 3; ++j) {
            ofs_evec << std::setw(15) << kpoint->xk[i][j];
        }
        ofs_evec << std::endl;
        for (j = 0; j < nbands; ++j) {
            omega2 = dynamical->eval_phonon[i][j];
            if (omega2 >= 0.0) {
                omega2 = omega2 * omega2;
            } else {
                omega2 = - omega2 * omega2;
            }

            ofs_evec << "### mode " << std::setw(8) << j + 1 << " : ";
            ofs_evec << std::setw(15) << omega2 << std::endl;

            for (k = 0; k < neval; ++k) {
                ofs_evec << std::setw(15) << real(dynamical->evec_phonon[i][j][k]);
                ofs_evec << std::setw(15) << imag(dynamical->evec_phonon[i][j][k]) << std::endl;
            }
            ofs_evec << std::endl;
        }
        ofs_evec << std::endl;
    }
    ofs_evec.close();

    std::cout << "  " << std::setw(input->job_title.length() + 12) << std::left << file_evec;
    std::cout << " : Eigenvector of all k points" << std::endl;
}

double Writes::in_kayser(const double x)
{
    return x * Ry_to_kayser;
}

void Writes::write_thermodynamics()
{
    unsigned int i, NT;
    double Tmin = system->Tmin;
    double Tmax = system->Tmax;
    double dT = system->dT;

    double T;
    std::string file_thermo;

    NT = static_cast<unsigned int>((Tmax - Tmin) / dT) + 1;

    std::ofstream ofs_thermo;
    file_thermo = input->job_title + ".thermo";
    ofs_thermo.open(file_thermo.c_str(), std::ios::out);
    if (!ofs_thermo) error->exit("write_thermodynamics", "cannot open file_cv");
    ofs_thermo << "# Temperature [K], Heat capacity / kB, Entropy / kB, Internal energy [Ry], Free energy [Ry]" << std::endl;


    //     for (i = 0; i <= NT; ++i) {
    //    
    //         T = Tmin + dT * static_cast<double>(i);
    //         std::cout << " T = " << std::setw(15) << T;
    //         TD = 1000.0;
    //         phonon_thermodynamics->Debye_T(T, TD);
    //         std::cout << "TD = " << TD << std::endl;
    //     }


    for (i = 0; i < NT; ++i) {
        T = Tmin + dT * static_cast<double>(i);

        ofs_thermo << std::setw(16) << std::fixed << T;
        ofs_thermo << std::setw(18) << std::scientific << thermodynamics->Cv_tot(T) / k_Boltzmann;
        ofs_thermo << std::setw(18) << thermodynamics->vibrational_entropy(T) / k_Boltzmann;
        ofs_thermo << std::setw(18) << thermodynamics->internal_energy(T);
        ofs_thermo << std::setw(18) << thermodynamics->free_energy(T) << std::endl;
    }

    ofs_thermo.close();

    std::cout << "  " << std::setw(input->job_title.length() + 12) << std::left << file_thermo;
    std::cout << " : Thermodynamic quantities" << std::endl;
}

void Writes::write_gruneisen()
{
    if (kpoint->kpoint_mode == 1) {
        if (nbands < 0 || nbands > 3 * system->natmin) {
            nbands = 3 * system->natmin;
        }

        std::ofstream ofs_gruneisen;

        std::string file_gru = input->job_title + ".gruneisen";
        ofs_gruneisen.open(file_gru.c_str(), std::ios::out);
        if (!ofs_gruneisen) error->exit("write_gruneisen", "cannot open file_vel");

        unsigned int i, j;
        unsigned int nk = kpoint->nk;

        double *kaxis = kpoint->kaxis;

        ofs_gruneisen << "# k-axis, gamma" << std::endl;
        ofs_gruneisen.setf(std::ios::fixed);

        for (i = 0; i < nk; ++i) {
            ofs_gruneisen << std::setw(8) << kaxis[i];
            for (j = 0; j < nbands; ++j) {
                ofs_gruneisen << std::setw(15) << gruneisen->gruneisen[i][j].real();
            }
            ofs_gruneisen << std::endl;
        }

        ofs_gruneisen.close();

        std::cout << "  " << std::setw(input->job_title.length() + 12) << std::left << file_gru;
        std::cout << " : Gruneisen parameters along given k-path" << std::endl;

    } else {

        std::ofstream ofs_gruall;
        std::string file_gruall;
        file_gruall = input->job_title + ".gru_all";
        ofs_gruall.open(file_gruall.c_str(), std::ios::out);
        if (!ofs_gruall) error->exit("write_gruneisen", "cannot open file_gruall");

        unsigned int i, j, k;
        unsigned int nk = kpoint->nk;
        unsigned int ns = dynamical->neval;

        ofs_gruall << "# knum, snum, omega [cm^-1], gruneisen parameter" << std::endl;

        for (i = 0; i < nk; ++i) {
            ofs_gruall << "# knum = " << i;
            for (k = 0; k < 3; ++k) {
                ofs_gruall << std::setw(15) << kpoint->xk[i][k];
            }
            ofs_gruall << std::endl;

            for (j = 0; j < ns; ++j) {
                ofs_gruall << std::setw(5) << i;
                ofs_gruall << std::setw(5) << j;
                ofs_gruall << std::setw(15) << in_kayser(dynamical->eval_phonon[i][j]);
                ofs_gruall << std::setw(15) << gruneisen->gruneisen[i][j].real();
                ofs_gruall << std::endl;
            }
        }

        ofs_gruall.close();


        std::cout << "  " << std::setw(input->job_title.length() + 12) << std::left << file_gruall;
        std::cout << " : Gruneisen parameters at all k points" << std::endl;
    }
}

void Writes::write_msd()
{
    // Write room mean square displacement of atoms

    std::string file_rmsd = input->job_title + ".msd";
    std::ofstream ofs_rmsd;

    unsigned int i, j;
    unsigned int NT;
    unsigned int ns = dynamical->neval;

    double Tmin = system->Tmin;
    double Tmax = system->Tmax;
    double dT = system->dT;
    double T, d2_tmp;

    ofs_rmsd.open(file_rmsd.c_str(), std::ios::out);
    if (!ofs_rmsd) error->exit("write_rmsd", "Could not open file_rmsd");

    ofs_rmsd << "# Mean Square Displacements at a function of temperature." << std::endl;
    ofs_rmsd << "# Temperature [K], <(u_{1}^{x})^{2}>, <(u_{1}^{y})^{2}>, <(u_{1}^{z})^{2}>, .... [Angstrom^2]" << std::endl;

    NT = static_cast<unsigned int>((Tmax - Tmin) / dT) + 1;

    for (i = 0; i < NT; ++i) {

        T = Tmin + static_cast<double>(i) * dT;
        ofs_rmsd << std::setw(15) << T;

        for (j = 0; j < ns; ++j) {
            d2_tmp = thermodynamics->disp2_avg(T, j, j);
            ofs_rmsd << std::setw(15) << d2_tmp * std::pow(Bohr_in_Angstrom, 2.0);
        }
        ofs_rmsd << std::endl;
    }
    ofs_rmsd.close();

    std::cout << "  " << std::setw(input->job_title.length() + 12) << std::left << file_rmsd;
    std::cout << " : Mean-square-displacement (MSD)" << std::endl;
}

void Writes::write_kappa()
{
    // Write lattice thermal conductivity

    if (mympi->my_rank == 0) {
        int i, j, k;

        std::string file_kappa = input->job_title + ".kl";
        std::string file_kappa2 = input->job_title + ".kl_spec";

        std::ofstream ofs_kl;

        ofs_kl.open(file_kappa.c_str(), std::ios::out);
        if (!ofs_kl) error->exit("write_kappa", "Could not open file_kappa");

        ofs_kl << "# Temperature [K], Thermal Conductivity (xx, xy, xz, yx, yy, yz, zx, zy, zz) [W/mK]" << std::endl;

        if (isotope->include_isotope) {
            ofs_kl << "# Isotope effects are included." << std::endl;
        }

        for (i = 0; i < conductivity->ntemp; ++i) {
            ofs_kl << std::setw(10) << std::right << std::fixed << std::setprecision(2)
                << conductivity->Temperature[i];
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
            if (!ofs_kl) error->exit("write_kappa", "Could not open file_kappa2");

            ofs_kl << "# Temperature [K], Frequency [cm^-1], Thermal Conductivity Spectra (xx, yy, zz) [W/mK * cm]" << std::endl;

            if (isotope->include_isotope) {
                ofs_kl << "# Isotope effects are included." << std::endl;
            }

            for (i = 0; i < conductivity->ntemp; ++i) {
                for (j = 0; j < dos->n_energy; ++j) {
                    ofs_kl << std::setw(10) << std::right << std::fixed << std::setprecision(2)
                        << conductivity->Temperature[i];
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

        std::cout << std::endl;
        std::cout << " -----------------------------------------------------------------" << std::endl << std::endl;
        std::cout << " Lattice thermal conductivity is stored in the file " << file_kappa << std::endl;
        if (conductivity->calc_kappa_spec) {
            std::cout << " Thermal conductivity spectra is stored in the file " << file_kappa2 << std::endl;
        }
    }
}

void Writes::write_selfenergy_isotope()
{
    unsigned int i, k;
    unsigned int ns = dynamical->neval;
    unsigned int knum;
    double **eval = dynamical->eval_phonon;
    double **gamma_iso = isotope->gamma_isotope;

    if (mympi->my_rank == 0) {
        if (isotope->include_isotope == 2) {

            std::string file_iso = input->job_title + ".self_isotope";
            std::ofstream ofs_iso;

            ofs_iso.open(file_iso.c_str(), std::ios::out);
            if (!ofs_iso) error->exit("write_selfenergy_isotope", "Could not open file_iso");

            ofs_iso << "# Phonon selfenergy due to phonon-isotope scatterings for the irreducible k points." << std::endl;
            ofs_iso << "# Irred. knum, mode num, frequency [cm^-1], Gamma_iso [cm^-1]" << std::endl << std::endl;

            for (i = 0; i < kpoint->nk_reduced; ++i) {
                ofs_iso << "# Irreducible k point  : " << std::setw(8) << i + 1;
                ofs_iso << " (" << std::setw(4) << kpoint->kpoint_irred_all[i].size() << ")" << std::endl;

                knum = kpoint->kpoint_irred_all[i][0].knum;

                ofs_iso << "## xk = " << std::setw(3);
                for (k = 0; k < 3; ++k) ofs_iso << std::setw(15) << kpoint->xk[knum][k];
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
                                         const unsigned int ncell[3])
{
    unsigned int i, j, k;
    unsigned int iband, istep;
    unsigned int ns = dynamical->neval;
    unsigned int natmin = system->natmin;
    unsigned int nsuper = ncell[0] * ncell[1] * ncell[2];
    unsigned int ix, iy, iz;
    unsigned int icell;
    unsigned int nsteps = 20;
    unsigned int ntmp = nbands;
    unsigned int ndigits = 0;

    double phase_time;
    double max_disp_factor = 0.1;
    double lavec_super[3][3];
    double norm, dmod[3];
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
        error->warn("write_normal_mode_animation",
                    "The supercell size is not commensurate with given k point.");
    }

    rotvec(kvec, xk, system->rlavec_p, 'T');
    norm = std::sqrt(kvec[0] * kvec[0] + kvec[1] * kvec[1] + kvec[2] * kvec[2]);
    if (norm > eps) {
        for (i = 0; i < 3; ++i) kvec[i] /= norm;
    }

    // Allocation

    memory->allocate(eval, ns);
    memory->allocate(evec, ns, ns);
    memory->allocate(evec_mag, ns, ns);
    memory->allocate(evec_theta, ns, ns);
    memory->allocate(disp_mag, ns, ns);
    memory->allocate(xmod, nsuper, natmin, 3);
    memory->allocate(kd_tmp, natmin);
    memory->allocate(mass, natmin);
    memory->allocate(phase_cell, nsuper);

    // Get eigenvalues and eigenvectors at xk

    dynamical->eval_k(xk, kvec, fcs_phonon->fc2_ext, eval, evec, true);

    for (i = 0; i < ns; ++i) {
        for (j = 0; j < ns; ++j) {
            evec_mag[i][j] = std::abs(evec[i][j]);
            evec_theta[i][j] = std::arg(evec[i][j]);
        }
    }

    // Get fractional coordinates of atoms in a primitive cell

    memory->allocate(xtmp, natmin, 3);

    for (i = 0; i < natmin; ++i) {
        rotvec(xtmp[i], system->xr_s[system->map_p2s[i][0]], system->lavec_s);
        rotvec(xtmp[i], xtmp[i], system->rlavec_p);
        for (j = 0; j < 3; ++j) xtmp[i][j] /= 2.0 * pi;
    }

    // Prepare fractional coordinates of atoms in the supercell
    icell = 0;

    for (ix = 0; ix < ncell[0]; ++ix) {
        for (iy = 0; iy < ncell[1]; ++iy) {
            for (iz = 0; iz < ncell[2]; ++iz) {

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

    memory->deallocate(xtmp);

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

    double mass_min = mass[0];

    for (i = 0; i < natmin; ++i) {
        if (mass[i] < mass_min) mass_min = mass[i];
    }

    double max_disp_mag, disp_mag_tmp;

    for (i = 0; i < ns; ++i) {

        max_disp_mag = 0.0;

        for (j = 0; j < ns; ++j) {
            disp_mag[i][j] = std::sqrt(mass_min / mass[j / 3]) * evec_mag[i][j];
        }

        for (j = 0; j < natmin; ++j) {
            disp_mag_tmp = 0.0;
            for (k = 0; k < 3; ++k) disp_mag_tmp += std::pow(disp_mag[i][3 * j + k], 2);
            disp_mag_tmp = std::sqrt(disp_mag_tmp);
            max_disp_mag = std::max(max_disp_mag, disp_mag_tmp);
        }

        for (j = 0; j < ns; ++j) disp_mag[i][j] *= max_disp_factor / max_disp_mag;
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
            std::string result = ss.str();

            file_anime = input->job_title + ".anime" + result + ".axsf";

            ofs_anime.open(file_anime.c_str(), std::ios::out);
            if (!ofs_anime)
                error->exit("write_normal_mode_animation",
                            "cannot open file_anime");

            ofs_anime.unsetf(std::ios::scientific);

            ofs_anime << "# Animation of a phonon mode" << std::endl;
            ofs_anime << "# K point : ";
            for (i = 0; i < 3; ++i) ofs_anime << std::setw(8) << xk_in[i];
            ofs_anime << std::endl;
            ofs_anime << "# Mode index : " << std::setw(4) << iband + 1 << std::endl;
            ofs_anime << "# Frequency (cm^-1) : " << in_kayser(eval[iband]) << std::endl;
            ofs_anime << std::endl;

            ofs_anime.setf(std::ios::scientific);

            ofs_anime << "ANIMSTEPS " << nsteps << std::endl;
            ofs_anime << "CRYSTAL" << std::endl;
            ofs_anime << "PRIMVEC" << std::endl;

            for (i = 0; i < 3; ++i) {
                for (j = 0; j < 3; ++j) {
                    ofs_anime << std::setw(15) << lavec_super[j][i];
                }
                ofs_anime << std::endl;
            }

            for (istep = 0; istep < nsteps; ++istep) {

                phase_time = 4.0 * pi / static_cast<double>(nsteps) * static_cast<double>(istep);

                ofs_anime << "PRIMCOORD " << std::setw(10) << istep + 1 << std::endl;
                ofs_anime << std::setw(10) << natmin * nsuper << std::setw(10) << 1 << std::endl;

                for (i = 0; i < nsuper; ++i) {
                    for (j = 0; j < natmin; ++j) {

                        ofs_anime << std::setw(10) << kd_tmp[j];

                        for (k = 0; k < 3; ++k) {
                            ofs_anime << std::setw(15)
                                << xmod[i][j][k]
                                + disp_mag[j][k]
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
            std::string result = ss.str();

            file_anime = input->job_title + ".anime" + result + ".xyz";

            ofs_anime.open(file_anime.c_str(), std::ios::out);
            if (!ofs_anime)
                error->exit("write_normal_mode_animation",
                            "cannot open file_anime");

            ofs_anime.setf(std::ios::scientific);

            for (istep = 0; istep < nsteps; ++istep) {

                phase_time = 4.0 * pi / static_cast<double>(nsteps) * static_cast<double>(istep);

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
                                + disp_mag[j][k]
                                * std::sin(phase_cell[i] + evec_theta[iband][3 * j + k] + phase_time);
                        }
                        ofs_anime << std::endl;
                    }
                }
            }

            ofs_anime.close();
        }

    }

    memory->deallocate(xmod);
    memory->deallocate(kd_tmp);
    memory->deallocate(eval);
    memory->deallocate(evec);
    memory->deallocate(phase_cell);
    memory->deallocate(evec_mag);
    memory->deallocate(evec_theta);
    memory->deallocate(disp_mag);
    memory->deallocate(mass);
}


void Writes::write_participation_ratio()
{
    unsigned int i, j, k;
    unsigned int knum;
    unsigned int nk = kpoint->nk;
    unsigned int neval = dynamical->neval;
    unsigned int natmin = system->natmin;

    double **participation_ratio;
    double ***atomic_participation_ratio;

    std::ofstream ofs_pr, ofs_apr;
    std::string file_pr = input->job_title + ".pr";
    std::string file_apr = input->job_title + ".apr";

    ofs_pr.open(file_pr.c_str(), std::ios::out);
    if (!ofs_pr)
        error->exit("write_participation_ratio",
                    "cannot open file_pr");

    ofs_pr.setf(std::ios::scientific);

    ofs_apr.open(file_apr.c_str(), std::ios::out);
    if (!ofs_apr)
        error->exit("write_participation_ratio",
                    "cannot open file_apr");

    ofs_apr.setf(std::ios::scientific);

    memory->allocate(participation_ratio, nk, neval);
    memory->allocate(atomic_participation_ratio, nk, neval, natmin);

    dynamical->calc_participation_ratio_all(dynamical->evec_phonon,
                                            participation_ratio,
                                            atomic_participation_ratio);

    ofs_pr << "# Participation ratio of each phonon modes at k points" << std::endl;

    if (kpoint->kpoint_mode == 0 || kpoint->kpoint_mode == 1) {

        ofs_pr << "# kpoint, mode, PR[kpoint][mode]" << std::endl;

        for (i = 0; i < nk; ++i) {
            ofs_pr << "#" << std::setw(8) << i + 1;
            ofs_pr << " xk = ";
            for (j = 0; j < 3; ++j) {
                ofs_pr << std::setw(15) << kpoint->xk[i][j];
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

    } else if (kpoint->kpoint_mode == 2) {

        ofs_pr << "# irred. kpoint, mode, frequency[kpoint][mode] (cm^-1), PR[kpoint][mode]" << std::endl;

        for (i = 0; i < kpoint->nk_reduced; ++i) {
            knum = kpoint->kpoint_irred_all[i][0].knum;
            ofs_pr << "#" << std::setw(8) << i + 1;
            ofs_pr << " xk = ";
            for (j = 0; j < 3; ++j) {
                ofs_pr << std::setw(15) << kpoint->xk[knum][j];
            }
            ofs_pr << std::endl;
            for (j = 0; j < nbands; ++j) {
                ofs_pr << std::setw(8) << i + 1;
                ofs_pr << std::setw(5) << j + 1;
                ofs_pr << std::setw(15) << in_kayser(dynamical->eval_phonon[knum][j]);
                ofs_pr << std::setw(15) << participation_ratio[knum][j];
                ofs_pr << std::endl;
            }
            ofs_pr << std::endl;
        }
    }

    ofs_pr.close();

    ofs_apr << "# Atomic participation ratio of each phonon modes at k points" << std::endl;

    if ((kpoint->kpoint_mode == 0) || (kpoint->kpoint_mode == 1)) {

        ofs_apr << "# kpoint, mode, atom, APR[kpoint][mode][atom]" << std::endl;

        for (i = 0; i < nk; ++i) {
            ofs_apr << "#" << std::setw(8) << i + 1;
            ofs_apr << " xk = ";
            for (j = 0; j < 3; ++j) {
                ofs_apr << std::setw(15) << kpoint->xk[i][j];
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
    } else if (kpoint->kpoint_mode == 2) {

        ofs_apr << "# irred. kpoint, mode, atom, frequency[kpoint][mode] (cm^-1), APR[kpoint][mode][atom]" << std::endl;

        for (i = 0; i < kpoint->nk_reduced; ++i) {
            knum = kpoint->kpoint_irred_all[i][0].knum;

            ofs_apr << "#" << std::setw(8) << i + 1;
            ofs_apr << " xk = ";
            for (j = 0; j < 3; ++j) {
                ofs_apr << std::setw(15) << kpoint->xk[knum][j];
            }
            ofs_apr << std::endl;
            for (j = 0; j < nbands; ++j) {
                for (k = 0; k < natmin; ++k) {
                    ofs_apr << std::setw(8) << i + 1;
                    ofs_apr << std::setw(5) << j + 1;
                    ofs_apr << std::setw(5) << k + 1;
                    ofs_apr << std::setw(15) << in_kayser(dynamical->eval_phonon[knum][j]);
                    ofs_apr << std::setw(15) << atomic_participation_ratio[knum][j][k];
                    ofs_apr << std::endl;
                }
            }
            ofs_apr << std::endl;
        }
    }


    memory->deallocate(participation_ratio);
    memory->deallocate(atomic_participation_ratio);

    std::cout << "  " << std::setw(input->job_title.length() + 12) << std::left << file_pr;
    std::cout << " : Participation ratio for all k points" << std::endl;
    std::cout << "  " << std::setw(input->job_title.length() + 12) << std::left << file_apr;
    std::cout << " : Atomic participation ratio for all k points" << std::endl;
}
