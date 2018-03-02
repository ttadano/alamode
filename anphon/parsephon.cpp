/*
 parsephon.cpp

 Copyright (c) 2014, 2015, 2016 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#include "mpi_common.h"
#include "parsephon.h"
#include "conductivity.h"
#include "dynamical.h"
#include "error.h"
#include "ewald.h"
#include "fcs_phonon.h"
#include "gruneisen.h"
#include "integration.h"
#include "isotope.h"
#include "kpoint.h"
#include "memory.h"
#include "phonon_dos.h"
#include "phonon_velocity.h"
#include "relaxation.h"
#include "scph.h"
#include "symmetry_core.h"
#include "system.h"
#include "thermodynamics.h"
#include "write_phonons.h"
#include <sys/stat.h>
#include <sstream>
#include <istream>
#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <string>
#include <algorithm>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

using namespace PHON_NS;

Input::Input(PHON *phon): Pointers(phon)
{
    from_stdin = false;
    job_title = "";
}

Input::~Input()
{
    if (ifs_input.is_open()) ifs_input.close();
}

void Input::parce_input(int narg, char **arg)
{
    if (narg == 1) {

        from_stdin = true;

    } else {

        from_stdin = false;

        ifs_input.open(arg[1], std::ios::in);
        if (!ifs_input) {
            std::cout << "No such file or directory: " << arg[1] << std::endl;
            exit(EXIT_FAILURE);
        }
    }

    if (!locate_tag("&general"))
        error->exit("parse_input",
                    "&general entry not found in the input file");
    parse_general_vars();

    if (!locate_tag("&cell"))
        error->exit("parse_input",
                    "&cell entry not found in the input file");
    parse_cell_parameter();

    bool use_defaults_for_analysis = !locate_tag("&analysis");
    parse_analysis_vars(use_defaults_for_analysis);

    if (!locate_tag("&kpoint"))
        error->exit("parse_input",
                    "&kpoint entry not found in the input file");
    parse_kpoints();

    if (phon->mode == "SCPH") {
        if (!locate_tag("&scph"))
            error->exit("parse_input",
                        "&scph entry not found in the input file");
        parse_scph_vars();
    }
}


void Input::parse_general_vars()
{
    // Read input parameters in the &general-field.

    int i;
    int nkd;
    struct stat st;
    std::string str_tmp;
    std::vector<std::string> input_list{
        "PREFIX", "MODE", "NSYM", "TOLERANCE", "PRINTSYM", "FCSXML", "FC2XML",
        "TMIN", "TMAX", "DT", "NBANDS", "NONANALYTIC", "BORNINFO", "NA_SIGMA",
        "ISMEAR", "EPSILON", "EMIN", "EMAX", "DELTA_E", "RESTART", "TREVSYM",
        "NKD", "KD", "MASS", "TRISYM", "PREC_EWALD", "CLASSICAL", "BCONNECT", "BORNSYM"
    };

    std::vector<std::string> no_defaults{"PREFIX", "MODE", "FCSXML", "NKD", "KD"};
    std::vector<std::string> kdname_v, masskd_v;
    std::map<std::string, std::string> general_var_dict;

    double *masskd = nullptr;
    std::string *kdname = nullptr;

    if (from_stdin) {
        std::cin.ignore();
    } else {
        ifs_input.ignore();
    }

    get_var_dict(input_list, general_var_dict);

    for (auto &no_default : no_defaults) {
        if (general_var_dict.find(no_default) == general_var_dict.end()) {
            error->exit("parse_general_vars",
                        "The following variable is not found in &general input region: ",
                        no_default.c_str());
        }
    }

    std::string prefix = general_var_dict["PREFIX"];
    std::string mode = general_var_dict["MODE"];
    std::string file_result = prefix + ".result";

    std::transform(mode.begin(), mode.end(), mode.begin(), toupper);

    std::string fcsinfo = general_var_dict["FCSXML"];
    assign_val(nkd, "NKD", general_var_dict);

    split_str_by_space(general_var_dict["KD"], kdname_v);

    if (kdname_v.size() != nkd) {
        error->exit("parse_general_vars",
                    "The number of entries for KD is inconsistent with NKD");
    } else {
        memory->allocate(kdname, nkd);
        for (i = 0; i < nkd; ++i) {
            kdname[i] = kdname_v[i];
        }
    }

    if (!general_var_dict["MASS"].empty()) {
        split_str_by_space(general_var_dict["MASS"], masskd_v);

        if (masskd_v.size() != nkd) {
            error->exit("parse_general_vars",
                        "The number of entries for MASS is inconsistent with NKD");
        } else {
            memory->allocate(masskd, nkd);
            for (i = 0; i < nkd; ++i) {
                masskd[i] = my_cast<double>(masskd_v[i]);
            }
        }
    }


    // Default values

    double Tmin = 0.0;
    double Tmax = 1000.0;
    double dT = 10.0;

    double emin = 0.0;
    double emax = 1000.0;
    double delta_e = 10.0;

    unsigned int nonanalytic = 0;
    int nsym = 0;
    double tolerance = 1.0e-6;
    bool printsymmetry = false;
    bool sym_time_reversal = false;
    bool use_triplet_symmetry = true;
    bool classical = false;
    unsigned int band_connection = 0;
    unsigned int bornsym = 0;

    double prec_ewald = 1.0e-12;

    // if file_result exists in the current directory, 
    // restart mode will be automatically turned on.
    bool restart = stat(file_result.c_str(), &st) == 0;

    int nbands = -1;
    std::string borninfo = "";
    std::string fc2info = "";

    int ismear = -1;
    double epsilon = 10.0;
    double na_sigma = 0.1;

    // Assign given values

    assign_val(Tmin, "TMIN", general_var_dict);
    assign_val(Tmax, "TMAX", general_var_dict);
    assign_val(dT, "DT", general_var_dict);

    assign_val(emin, "EMIN", general_var_dict);
    assign_val(emax, "EMAX", general_var_dict);
    assign_val(delta_e, "DELTA_E", general_var_dict);

    assign_val(nsym, "NSYM", general_var_dict);
    assign_val(sym_time_reversal, "TREVSYM", general_var_dict);
    assign_val(tolerance, "TOLERANCE", general_var_dict);
    assign_val(printsymmetry, "PRINTSYM", general_var_dict);

    assign_val(nonanalytic, "NONANALYTIC", general_var_dict);
    assign_val(restart, "RESTART", general_var_dict);

    assign_val(nbands, "NBANDS", general_var_dict);
    assign_val(borninfo, "BORNINFO", general_var_dict);
    assign_val(fc2info, "FC2XML", general_var_dict);

    assign_val(ismear, "ISMEAR", general_var_dict);
    assign_val(epsilon, "EPSILON", general_var_dict);
    assign_val(na_sigma, "NA_SIGMA", general_var_dict);
    assign_val(classical, "CLASSICAL", general_var_dict);
    assign_val(band_connection, "BCONNECT", general_var_dict);
    assign_val(use_triplet_symmetry, "TRISYM", general_var_dict);
    assign_val(bornsym, "BORNSYM", general_var_dict);

    if (band_connection > 2) {
        error->exit("parse_general_vars", "BCONNECT-tag can take 0, 1, or 2.");
    }

    if (nonanalytic == 3) {
        assign_val(prec_ewald, "PREC_EWALD", general_var_dict);
        if (prec_ewald <= 0.0 || prec_ewald >= 1.0) {
            error->exit("parse_general_vars",
                        "PREC_EWALD should be a small positive value.");
        }
        ewald->is_longrange = true;
        ewald->file_longrange = boost::lexical_cast<std::string>(general_var_dict["BORNINFO"]);
        ewald->prec_ewald = prec_ewald;
        ewald->rate_ab = 6.0 / pi;

    } else {
        ewald->is_longrange = false;
    }

    if (nonanalytic > 3) {
        error->exit("parse_general_vars",
                    "NONANALYTIC-tag can take 0, 1, 2, or 3.");
    }
    if (nonanalytic == 3) {
        if (mode == "SCPH") {
            error->exit("parse_general_vars",
                        "Sorry. NONANALYTIC=3 is not supported for MODE = SCPH.");
        }
    }

    // Copy the values to appropriate classes.

    job_title = prefix;
    writes->file_result = file_result;
    phon->mode = mode;
    phon->restart_flag = restart;
    symmetry->nsym = nsym;
    symmetry->tolerance = tolerance;
    symmetry->printsymmetry = printsymmetry;
    symmetry->time_reversal_sym = sym_time_reversal;

    system->Tmin = Tmin;
    system->Tmax = Tmax;
    system->dT = dT;
    system->nkd = nkd;

    memory->allocate(system->symbol_kd, nkd);
    for (i = 0; i < nkd; ++i) {
        system->symbol_kd[i] = kdname[i];
    }
    if (kdname) {
        memory->deallocate(kdname);
    }

    if (!general_var_dict["MASS"].empty()) {
        memory->allocate(system->mass_kd, nkd);
        for (i = 0; i < nkd; ++i) {
            system->mass_kd[i] = masskd[i];
        }
    }
    if (masskd) {
        memory->deallocate(masskd);
    }


    dos->emax = emax;
    dos->emin = emin;
    dos->delta_e = delta_e;

    dynamical->nonanalytic = nonanalytic;
    dynamical->na_sigma = na_sigma;
    dynamical->symmetrize_borncharge = bornsym;
    writes->nbands = nbands;
    dynamical->file_born = borninfo;
    dynamical->band_connection = band_connection;
    integration->epsilon = epsilon;
    fcs_phonon->file_fcs = fcsinfo;
    fcs_phonon->file_fc2 = fc2info;
    fcs_phonon->update_fc2 = !fc2info.empty();
    thermodynamics->classical = classical;
    integration->ismear = ismear;
    relaxation->use_triplet_symmetry = use_triplet_symmetry;

    general_var_dict.clear();
}

void Input::parse_scph_vars()
{
    // Read input parameters in the &scph-field.

    int i;

    struct stat st;
    std::string str_tmp;
    std::vector<std::string> input_list{
        "KMESH_SCPH", "KMESH_INTERPOLATE", "MIXALPHA", "MAXITER",
        "RESTART_SCPH", "IALGO", "SELF_OFFDIAG", "TOL_SCPH",
        "LOWER_TEMP", "WARMSTART"
    };
    std::vector<std::string> no_defaults{"KMESH_SCPH", "KMESH_INTERPOLATE"};
    std::vector<int> kmesh_v, kmesh_interpolate_v;
    std::map<std::string, std::string> scph_var_dict;

    if (from_stdin) {
        std::cin.ignore();
    } else {
        ifs_input.ignore();
    }

    get_var_dict(input_list, scph_var_dict);

    for (auto &no_default : no_defaults) {
        if (scph_var_dict.find(no_default) == scph_var_dict.end()) {
            error->exit("parse_general_vars",
                        "The following variable is not found in &scph input region: ",
                        no_default.c_str());
        }
    }

    std::string file_dymat = this->job_title + ".scph_dymat";

    // Default values

    double tolerance_scph = 1.0e-10;
    unsigned int maxiter = 1000;
    double mixalpha = 0.1;
    bool selenergy_offdiagonal = true;
    unsigned int ialgo_scph = 0;
    bool lower_temp = true;
    bool warm_start = true;

    // if file_dymat exists in the current directory, 
    // restart mode will be automatically turned on for SCPH calculations.

    bool restart_scph = stat(file_dymat.c_str(), &st) == 0;

    // Assign given values

    assign_val(restart_scph, "RESTART_SCPH", scph_var_dict);
    assign_val(maxiter, "MAXITER", scph_var_dict);
    assign_val(mixalpha, "MIXALPHA", scph_var_dict);
    assign_val(selenergy_offdiagonal, "SELF_OFFDIAG", scph_var_dict);
    assign_val(ialgo_scph, "IALGO", scph_var_dict);
    assign_val(tolerance_scph, "TOL_SCPH", scph_var_dict);
    assign_val(lower_temp, "LOWER_TEMP", scph_var_dict);
    assign_val(warm_start, "WARMSTART", scph_var_dict);

    str_tmp = scph_var_dict["KMESH_SCPH"];

    if (!str_tmp.empty()) {

        std::istringstream is(str_tmp);

        while (1) {
            str_tmp.clear();
            is >> str_tmp;
            if (str_tmp.empty()) {
                break;
            }
            kmesh_v.push_back(my_cast<unsigned int>(str_tmp));
        }

        if (kmesh_v.size() != 3) {
            error->exit("parse_general_vars",
                        "The number of entries for KMESH_SCPH has to be 3.");
        }
    } else {
        error->exit("parse_general_vars",
                    "Please specify KMESH_SCPH for mode = SCPH");
    }

    str_tmp = scph_var_dict["KMESH_INTERPOLATE"];
    if (!str_tmp.empty()) {

        std::istringstream is(str_tmp);

        while (1) {
            str_tmp.clear();
            is >> str_tmp;
            if (str_tmp.empty()) {
                break;
            }
            kmesh_interpolate_v.push_back(my_cast<unsigned int>(str_tmp));
        }

        if (kmesh_interpolate_v.size() != 3) {
            error->exit("parse_general_vars",
                        "The number of entries for KMESH_INTERPOLATE has to be 3.");
        }
    } else {
        error->exit("parse_general_vars",
                    "Please specify KMESH_INTERPOLATE for mode = SCPH");
    }

    // Copy the values to appropriate classes.

    for (i = 0; i < 3; ++i) {
        scph->kmesh_scph[i] = kmesh_v[i];
        scph->kmesh_interpolate[i] = kmesh_interpolate_v[i];
    }
    scph->mixalpha = mixalpha;
    scph->maxiter = maxiter;
    scph->restart_scph = restart_scph;
    scph->selfenergy_offdiagonal = selenergy_offdiagonal;
    scph->ialgo = ialgo_scph;
    scph->tolerance_scph = tolerance_scph;
    scph->lower_temp = lower_temp;
    scph->warmstart_scph = warm_start;

    kmesh_v.clear();
    kmesh_interpolate_v.clear();

    scph_var_dict.clear();
}


void Input::parse_analysis_vars(const bool use_default_values)
{
    // Read input parameters in the &analysis field.
    int i;

    std::vector<std::string> input_list{
        "PRINTEVEC", "PRINTXSF", "PRINTVEL", "QUARTIC", "KS_INPUT",
        "ATOMPROJ", "REALPART", "ISOTOPE", "ISOFACT",
        "FSTATE_W", "FSTATE_K", "PRIMTMSD", "DOS", "PDOS", "TDOS",
        "GRUNEISEN", "NEWFCS", "DELTA_A", "ANIME", "ANIME_CELLSIZE",
        "ANIME_FORMAT", "SPS", "PRINTV3", "PRINTPR", "FC2_EWALD",
        "KAPPA_SPEC", "SELF_W"
    };

    unsigned int cellsize[3];

    double *isotope_factor = nullptr;
    std::string ks_input, anime_format;
    std::map<std::string, std::string> analysis_var_dict;
    std::vector<std::string> isofact_v, anime_kpoint, anime_cellsize;

    // Default values

    bool print_xsf = false;
    bool print_anime = false;

    bool print_vel = false;
    bool print_evec = false;
    bool print_msd = false;

    bool compute_dos = true;
    bool projected_dos = false;
    bool two_phonon_dos = false;
    int scattering_phase_space = 0;
    bool print_gruneisen = false;
    bool print_newfcs = false;
    bool print_V3 = false;
    bool participation_ratio = false;

    double delta_a = 0.001;

    int quartic_mode = 0;
    bool ks_analyze_mode = false;
    bool atom_project_mode = false;
    bool calc_realpart = false;
    int include_isotope = 0;
    bool fstate_omega = false;
    bool fstate_k = false;
    bool bubble_omega = false;

    int calculate_kappa_spec = 0;

    bool print_fc2_ewald = false;
    bool print_self_consistent_fc2 = false;


    // Assign values to variables

    if (!use_default_values) {
        get_var_dict(input_list, analysis_var_dict);

        assign_val(print_vel, "PRINTVEL", analysis_var_dict);
        assign_val(print_evec, "PRINTEVEC", analysis_var_dict);
        assign_val(print_msd, "PRINTMSD", analysis_var_dict);

        assign_val(compute_dos, "DOS", analysis_var_dict);
        assign_val(projected_dos, "PDOS", analysis_var_dict);
        assign_val(two_phonon_dos, "TDOS", analysis_var_dict);
        assign_val(scattering_phase_space, "SPS", analysis_var_dict);
        assign_val(print_gruneisen, "GRUNEISEN", analysis_var_dict);
        assign_val(print_newfcs, "NEWFCS", analysis_var_dict);
        assign_val(delta_a, "DELTA_A", analysis_var_dict);

        assign_val(quartic_mode, "QUARTIC", analysis_var_dict);
        assign_val(atom_project_mode, "ATOMPROJ", analysis_var_dict);
        assign_val(calc_realpart, "REALPART", analysis_var_dict);
        assign_val(include_isotope, "ISOTOPE", analysis_var_dict);
        assign_val(fstate_omega, "FSTATE_W", analysis_var_dict);
        assign_val(fstate_k, "FSTATE_K", analysis_var_dict);
        assign_val(ks_input, "KS_INPUT", analysis_var_dict);
        assign_val(calculate_kappa_spec, "KAPPA_SPEC", analysis_var_dict);
        assign_val(bubble_omega, "SELF_W", analysis_var_dict);

        assign_val(print_xsf, "PRINTXSF", analysis_var_dict);
        assign_val(print_V3, "PRINTV3", analysis_var_dict);
        assign_val(participation_ratio, "PRINTPR", analysis_var_dict);
        assign_val(print_fc2_ewald, "FC2_EWALD", analysis_var_dict);

        if (analysis_var_dict.find("ANIME") == analysis_var_dict.end()) {
            print_anime = false;
        } else {
            print_anime = true;
        }
    }

    if (include_isotope) {

        if (!analysis_var_dict["ISOFACT"].empty()) {
            split_str_by_space(analysis_var_dict["ISOFACT"], isofact_v);

            if (isofact_v.size() != system->nkd) {
                error->exit("parse_analysis_vars",
                    "The number of entries for ISOFACT is inconsistent with NKD");
            }
            else {
                memory->allocate(isotope_factor, system->nkd);
                for (i = 0; i < system->nkd; ++i) {
                    isotope_factor[i] = my_cast<double>(isofact_v[i]);
                }
            }
        }
        
    }

    if (print_anime) {
        split_str_by_space(analysis_var_dict["ANIME"], anime_kpoint);

        if (anime_kpoint.size() != 3) {
            error->exit("parse_analysis_vars",
                        "The number of entries for ANIME should be 3.");
        }

        split_str_by_space(analysis_var_dict["ANIME_CELLSIZE"], anime_cellsize);

        if (anime_cellsize.size() != 3) {
            error->exit("parse_analysis_vars",
                        "The number of entries for ANIME_CELLSIZE should be 3.");
        }

        for (i = 0; i < 3; ++i) {
            try {
                cellsize[i] = boost::lexical_cast<unsigned int>(anime_cellsize[i]);
            }
            catch (std::exception &e) {
                std::cout << e.what() << std::endl;
                error->exit("parse_analysis_vars",
                            "ANIME_CELLSIZE must be a set of positive integers.");
            }
            if (cellsize[i] < 1) {
                error->exit("parse_analysis_vars",
                            "Please give positive integers in ANIME_CELLSIZE.");
            }
        }

        assign_val(anime_format, "ANIME_FORMAT", analysis_var_dict);
        std::transform(anime_format.begin(), anime_format.end(),
                       anime_format.begin(), toupper);

        if (anime_format.empty()) anime_format = "XYZ";

        if (anime_format != "XSF" && anime_format != "AXSF" && anime_format != "XYZ") {
            error->exit("parse_analysis_vars", "Invalid ANIME_FORMAT");
        }
    }

    // Copy the values to appropriate classes

    phonon_velocity->print_velocity = print_vel;
    dynamical->print_eigenvectors = print_evec;
    dynamical->participation_ratio = participation_ratio;
    writes->print_xsf = print_xsf;
    writes->print_anime = print_anime;

    if (print_anime) {
        for (i = 0; i < 3; ++i) {
            writes->anime_kpoint[i] = my_cast<double>(anime_kpoint[i]);
            writes->anime_cellsize[i] = cellsize[i];
        }
        writes->anime_format = anime_format;
    }

    writes->print_msd = print_msd;

    dos->compute_dos = compute_dos;
    dos->projected_dos = projected_dos;
    dos->two_phonon_dos = two_phonon_dos;
    dos->scattering_phase_space = scattering_phase_space;

    conductivity->calc_kappa_spec = calculate_kappa_spec;
    relaxation->quartic_mode = quartic_mode;
    relaxation->atom_project_mode = atom_project_mode;
    relaxation->calc_realpart = calc_realpart;
    relaxation->calc_fstate_omega = fstate_omega;
    relaxation->calc_fstate_k = fstate_k;
    relaxation->print_V3 = print_V3;
    relaxation->spectral_func = bubble_omega;
    isotope->include_isotope = include_isotope;
    relaxation->ks_input = ks_input;

    gruneisen->print_gruneisen = print_gruneisen;
    gruneisen->print_newfcs = print_newfcs;
    gruneisen->delta_a = delta_a;

    ewald->print_fc2_ewald = print_fc2_ewald;

    if (include_isotope) {
        if (!analysis_var_dict["ISOFACT"].empty()) {
            memory->allocate(isotope->isotope_factor, system->nkd);
            for (i = 0; i < system->nkd; ++i) {
                isotope->isotope_factor[i] = isotope_factor[i];
            }
        }
    }
    if (isotope_factor) {
        memory->deallocate(isotope_factor);
    }

    if (phon->mode == "SCPH") {
        scph->print_self_consistent_fc2 = print_self_consistent_fc2;
    }

    analysis_var_dict.clear();
}

void Input::parse_cell_parameter()
{
    // Read the cell parameter

    int i, j;
    double a;
    double lavec_tmp[3][3];
    std::string line;
    std::string line_wo_comment, line_tmp;
    std::vector<std::string> line_vec, line_split;
    std::string::size_type pos_first_comment_tag;

    line_vec.clear();

    if (from_stdin) {

        while (std::getline(std::cin, line)) {

            // Ignore comment region
            pos_first_comment_tag = line.find_first_of('#');

            if (pos_first_comment_tag == std::string::npos) {
                line_wo_comment = line;
            } else {
                line_wo_comment = line.substr(0, pos_first_comment_tag);
            }

            boost::trim_if(line_wo_comment, boost::is_any_of("\t\r\n "));

            if (line_wo_comment.empty()) continue;
            if (is_endof_entry(line_wo_comment)) break;

            line_vec.push_back(line_wo_comment);
        }

    } else {
        while (std::getline(ifs_input, line)) {

            // Ignore comment region
            pos_first_comment_tag = line.find_first_of('#');

            if (pos_first_comment_tag == std::string::npos) {
                line_wo_comment = line;
            } else {
                line_wo_comment = line.substr(0, pos_first_comment_tag);
            }

            boost::trim_if(line_wo_comment, boost::is_any_of("\t\r\n "));

            if (line_wo_comment.empty()) continue;
            if (is_endof_entry(line_wo_comment)) break;

            line_vec.push_back(line_wo_comment);
        }
    }

    if (line_vec.size() != 4) {
        error->exit("parse_cell_parameter",
                    "Too few or too much lines for the &cell field.\n \
                                            The number of valid lines for the &cell field should be 4.");
    }

    for (i = 0; i < 4; ++i) {

        line = line_vec[i];
        boost::split(line_split, line,
                     boost::is_any_of("\t "), boost::token_compress_on);

        if (i == 0) {
            // Lattice factor a
            if (line_split.size() == 1) {
                a = boost::lexical_cast<double>(line_split[0]);
            } else {
                error->exit("parse_cell_parameter",
                            "Unacceptable format for &cell field.");
            }

        } else {
            // Lattice vectors a1, a2, a3
            if (line_split.size() == 3) {
                for (j = 0; j < 3; ++j) {
                    lavec_tmp[j][i - 1] = boost::lexical_cast<double>(line_split[j]);
                }
            } else {
                error->exit("parse_cell_parameter",
                            "Unacceptable format for &cell field.");
            }
        }
    }

    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            system->lavec_p[i][j] = a * lavec_tmp[i][j];
        }
    }
}


void Input::parse_kpoints()
{
    // Read the settings in the &kpoint field.

    int i, kpmode;
    std::string line, line_wo_comment, str_tmp;
    std::vector<std::string> kpelem, line_vec;
    std::string::size_type pos_first_comment_tag;

    line_vec.clear();

    if (from_stdin) {

        while (std::getline(std::cin, line)) {

            // Ignore comment region
            pos_first_comment_tag = line.find_first_of('#');

            if (pos_first_comment_tag == std::string::npos) {
                line_wo_comment = line;
            } else {
                line_wo_comment = line.substr(0, pos_first_comment_tag);
            }

            boost::trim_if(line_wo_comment, boost::is_any_of("\t\r\n "));

            if (line_wo_comment.empty()) continue;
            if (is_endof_entry(line_wo_comment)) break;

            line_vec.push_back(line_wo_comment);
        }

    } else {
        while (std::getline(ifs_input, line)) {

            // Ignore comment region
            pos_first_comment_tag = line.find_first_of('#');

            if (pos_first_comment_tag == std::string::npos) {
                line_wo_comment = line;
            } else {
                line_wo_comment = line.substr(0, pos_first_comment_tag);
            }

            boost::trim_if(line_wo_comment, boost::is_any_of("\t\r\n "));

            if (line_wo_comment.empty()) continue;
            if (is_endof_entry(line_wo_comment)) break;

            line_vec.push_back(line_wo_comment);
        }
    }

    for (i = 0; i < line_vec.size(); ++i) {
        line = line_vec[i];
        boost::split(kpelem, line, boost::is_any_of("\t "), boost::token_compress_on);

        if (i == 0) {
            // kpmode 
            if (kpelem.size() == 1) {
                try {
                    kpmode = boost::lexical_cast<int>(kpelem[0]);
                }
                catch (std::exception &e) {
                    std::cout << e.what() << std::endl;
                    error->exit("parse_kpoints",
                                "KPMODE must be an integer. [0, 1, or 2]");
                }

                if (!(kpmode >= 0 && kpmode <= 3)) {
                    error->exit("parse_kpoints",
                                "KPMODE must be 0, 1, or 2.");
                }

            } else {
                error->exit("parse_kpoints",
                            "Unacceptable format for the &kpoint field.");
            }

        } else {
            // Read each entry of kpoint

            if (kpmode == 0 && kpelem.size() != 3) {
                error->exit("parse_kpoints",
                            "The number of columns must be 3 when KPMODE = 0");
            }
            if (kpmode == 1 && kpelem.size() != 9) {
                error->exit("parse_kpoints",
                            "The number of columns must be 9 when KPMODE = 1");
            }
            if (kpmode == 2 && kpelem.size() != 3) {
                error->exit("parse_kpoints",
                            "The number of columns must be 3 when KPMODE = 2");
            }
            if (kpmode == 3 && kpelem.size() != 8) {
                error->exit("parse_kpoints",
                            "The number of columns must be 8 when KPMODE = 3");
            }

            kpoint->kpInp.push_back(kpelem);
        }
    }

    kpoint->kpoint_mode = kpmode;
}


int Input::locate_tag(std::string key)
{
    int ret = 0;
    std::string line, line2;

    if (from_stdin) {

        // The following two lines do nothing when MPI version is executed.
        // I don't know why this happens.
        std::cin.clear();
        std::cin.seekg(0, std::ios_base::beg);

        while (std::cin >> line) {
#ifdef _USE_BOOST
            boost::to_lower(line);
            boost::trim(line);
#else
            std::transform(line.begin(), line.end(), line.begin(), tolower);
            line2 = line;
            line = trim(line2);
#endif
            if (line == key) {
                ret = 1;
                break;
            }
        }
        return ret;

    } else {

        ifs_input.clear();
        ifs_input.seekg(0, std::ios_base::beg);

        while (ifs_input >> line) {
#ifdef _USE_BOOST
            boost::to_lower(line);
            boost::trim(line);
#else
            std::transform(line.begin(), line.end(), line.begin(), tolower);
            line2 = line;
            line = trim(line2);
#endif
            if (line == key) {
                ret = 1;
                break;
            }
        }
        return ret;
    }
}

void Input::get_var_dict(const std::vector<std::string> &input_list,
                         std::map<std::string, std::string> &var_dict)
{
    std::string line, key, val;
    std::string line_wo_comment, line_tmp;
    std::string::size_type pos_first_comment_tag;
    std::vector<std::string> str_entry, str_varval;
    std::set<std::string> keyword_set;

    for (const auto &it : input_list) {
        keyword_set.insert(it);
    }

    var_dict.clear();

    if (from_stdin) {

        while (std::getline(std::cin, line)) {

            // Ignore comment region
            pos_first_comment_tag = line.find_first_of('#');

            if (pos_first_comment_tag == std::string::npos) {
                line_wo_comment = line;
            } else {
                line_wo_comment = line.substr(0, pos_first_comment_tag);
            }
#ifdef _USE_BOOST
            boost::trim_left(line_wo_comment);
#else
            line_tmp = line_wo_comment;
            line_wo_comment = ltrim(line_tmp);
#endif
            if (line_wo_comment.empty()) continue;
            if (is_endof_entry(line_wo_comment)) break;

            // Split the input line by ';'
#ifdef _USE_BOOST
            boost::split(str_entry, line_wo_comment, boost::is_any_of(";"));
#else
            str_entry = my_split(line_wo_comment, ';');
#endif


            for (auto it = str_entry.begin(); it != str_entry.end(); ++it) {

                // Split the input entry by '='
#ifdef _USE_BOOST
                std::string str_tmp = boost::trim_copy(*it);
#else
                std::string str_tmp = trim((*it));
#endif
                if (!str_tmp.empty()) {
#ifdef _USE_BOOST
                    boost::split(str_varval, str_tmp, boost::is_any_of("="));
#else
                    str_varval = my_split(str_tmp, '=');
#endif
                    if (str_varval.size() != 2) {
                        error->exit("get_var_dict",
                                    "Unacceptable format");
                    }
#ifdef _USE_BOOST
                    key = boost::to_upper_copy(boost::trim_copy(str_varval[0]));
                    val = boost::trim_copy(str_varval[1]);
#else
                    key = trim(str_varval[0]);
                    std::transform(key.begin(), key.end(), key.begin(), toupper);
                    val = trim(str_varval[1]);
#endif
                    if (keyword_set.find(key) == keyword_set.end()) {
                        std::cout << "Could not recognize the variable " << key << std::endl;
                        error->exit("get_var_dict",
                                    "Invalid variable found");
                    }

                    if (var_dict.find(key) != var_dict.end()) {
                        std::cout << "Variable " << key
                            << " appears twice in the input file." << std::endl;
                        error->exit("get_var_dict",
                                    "Redundant input parameter");
                    }

                    // If everything is OK, add the variable and the corresponding value
                    // to the dictionary.

                    var_dict.insert(std::map<std::string, std::string>::value_type(key, val));
                }
            }
        }
    } else {

        while (std::getline(ifs_input, line)) {

            // Ignore comment region
            pos_first_comment_tag = line.find_first_of('#');

            if (pos_first_comment_tag == std::string::npos) {
                line_wo_comment = line;
            } else {
                line_wo_comment = line.substr(0, pos_first_comment_tag);
            }
#ifdef _USE_BOOST
            boost::trim_left(line_wo_comment);
#else
            line_tmp = line_wo_comment;
            line_wo_comment = ltrim(line_tmp);
#endif
            if (line_wo_comment.empty()) continue;
            if (is_endof_entry(line_wo_comment)) break;

            // Split the input line by ';'
#ifdef _USE_BOOST
            boost::split(str_entry, line_wo_comment, boost::is_any_of(";"));
#else
            str_entry = my_split(line_wo_comment, ';');
#endif
            for (auto &it : str_entry) {

                // Split the input entry by '='

#ifdef _USE_BOOST
                std::string str_tmp = boost::trim_copy(*it);
#else
                std::string str_tmp = trim(it);
#endif
                if (!str_tmp.empty()) {
#ifdef _USE_BOOST
                    boost::split(str_varval, str_tmp, boost::is_any_of("="));
#else
                    str_varval = my_split(str_tmp, '=');
#endif

                    if (str_varval.size() != 2) {
                        error->exit("get_var_dict",
                                    "Unacceptable format");
                    }

#ifdef _USE_BOOST
                    key = boost::to_upper_copy(boost::trim_copy(str_varval[0]));
                    val = boost::trim_copy(str_varval[1]);
#else
                    key = trim(str_varval[0]);
                    std::transform(key.begin(), key.end(), key.begin(), toupper);
                    val = trim(str_varval[1]);
#endif

                    if (keyword_set.find(key) == keyword_set.end()) {
                        std::cout << "Could not recognize the variable "
                            << key << std::endl;
                        error->exit("get_var_dict",
                                    "Invalid variable found");
                    }

                    if (var_dict.find(key) != var_dict.end()) {
                        std::cout << "Variable " << key
                            << " appears twice in the input file." << std::endl;
                        error->exit("get_var_dict",
                                    "Redundant input parameter");
                    }

                    // If everything is OK, add the variable and the corresponding value
                    // to the dictionary.

                    var_dict.insert(std::map<std::string, std::string>::value_type(key, val));
                }
            }
        }

    }
    keyword_set.clear();
}


bool Input::is_endof_entry(const std::string str)
{
    return str[0] == '/';
}

void Input::split_str_by_space(const std::string str,
                               std::vector<std::string> &str_vec)
{
    std::string str_tmp;
    std::istringstream is(str);

    str_vec.clear();

    while (1) {
        str_tmp.clear();
        is >> str_tmp;
        if (str_tmp.empty()) {
            break;
        }
        str_vec.push_back(str_tmp);
    }
    str_tmp.clear();
}

template <typename T>
void Input::assign_val(T &val,
                       const std::string key,
                       std::map<std::string, std::string> dict)
{
    // Assign a value to the variable "key" using the boost::lexica_cast.

    if (!dict[key].empty()) {
        try {
            val = boost::lexical_cast<T>(dict[key]);
        }
        catch (std::exception &e) {
            std::string str_tmp;
            std::cout << e.what() << std::endl;
            str_tmp = "Invalid entry for the " + key + " tag.\n";
            str_tmp += " Please check the input value.";
            error->exit("assign_val", str_tmp.c_str());
        }
    }
}


template <typename T_to, typename T_from>
T_to Input::my_cast(T_from const &x)
{
    std::stringstream ss;
    T_to ret;

    ss << x;
    ss >> ret;

    return ret;
}
