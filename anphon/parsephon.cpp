/*
 parsephon.cpp

 Copyright (c) 2014, 2015, 2016 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#include "parsephon.h"
#include "conductivity.h"
#include "dielec.h"
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
#include "anharmonic_core.h"
#include "mode_analysis.h"
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

Input::Input(PHON *phon) : Pointers(phon)
{
    from_stdin = false;
    job_title = "";
}

Input::~Input()
{
    if (ifs_input.is_open()) ifs_input.close();
}

void Input::parce_input(int narg,
                        char **arg)
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

    const auto use_defaults_for_analysis = !locate_tag("&analysis");
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
    const std::vector<std::string> input_list{
            "PREFIX", "MODE", "NSYM", "TOLERANCE", "PRINTSYM", "FCSXML", "FC2XML",
            "TMIN", "TMAX", "DT", "NBANDS", "NONANALYTIC", "BORNINFO", "NA_SIGMA",
            "ISMEAR", "EPSILON", "EMIN", "EMAX", "DELTA_E", "RESTART", "TREVSYM",
            "NKD", "KD", "MASS", "TRISYM", "PREC_EWALD", "CLASSICAL", "BCONNECT", "BORNSYM",
            "VERBOSITY"
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

    const auto prefix = general_var_dict["PREFIX"];
    auto mode = general_var_dict["MODE"];
    auto file_result = prefix + ".result";

    std::transform(mode.begin(), mode.end(), mode.begin(), toupper);

    const auto fcsinfo = general_var_dict["FCSXML"];
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

    auto Tmin = 0.0;
    auto Tmax = 1000.0;
    auto dT = 10.0;

    auto emin = 0.0;
    auto emax = 1000.0;
    auto delta_e = 10.0;

    unsigned int nonanalytic = 0;
    auto nsym = 0;
    auto tolerance = 1.0e-6;
    auto printsymmetry = false;
    auto sym_time_reversal = false;
    auto use_triplet_symmetry = true;
    auto classical = false;
    unsigned int band_connection = 0;
    unsigned int bornsym = 0;
    unsigned int verbosity = 1;

    auto prec_ewald = 1.0e-12;

    // if file_result exists in the current directory, 
    // restart mode will be automatically turned on.
    auto restart = stat(file_result.c_str(), &st) == 0;

    auto nbands = -1;
    std::string borninfo;
    std::string fc2info;

    auto ismear = -1;
    auto epsilon = 10.0;
    auto na_sigma = 0.1;

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
    assign_val(verbosity, "VERBOSITY", general_var_dict);

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
    if (nonanalytic && borninfo == "") {
        error->exit("parse_general_vars",
                    "BORNINFO must be specified when NONANALYTIC > 0.");
    }
    // if (nonanalytic == 3) {
    //     if (mode == "SCPH") {
    //         error->exit("parse_general_vars",
    //                     "Sorry. NONANALYTIC=3 is not supported for MODE = SCPH.");
    //     }
    // }

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
    anharmonic_core->use_triplet_symmetry = use_triplet_symmetry;

    general_var_dict.clear();
}

void Input::parse_scph_vars()
{
    // Read input parameters in the &scph-field.

    struct stat st{};
    const std::vector<std::string> input_list{
            "KMESH_SCPH", "KMESH_INTERPOLATE", "MIXALPHA", "MAXITER",
            "RESTART_SCPH", "IALGO", "SELF_OFFDIAG", "TOL_SCPH",
            "LOWER_TEMP", "WARMSTART", "BUBBLE"
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

    auto file_dymat = this->job_title + ".scph_dymat";

    // Default values

    auto tolerance_scph = 1.0e-10;
    unsigned int maxiter = 1000;
    auto mixalpha = 0.1;
    auto selenergy_offdiagonal = true;
    unsigned int ialgo_scph = 0;
    auto lower_temp = true;
    auto warm_start = true;
    unsigned int bubble = 0;

    // if file_dymat exists in the current directory, 
    // restart mode will be automatically turned on for SCPH calculations.

    auto restart_scph = stat(file_dymat.c_str(), &st) == 0;

    // Assign given values

    assign_val(restart_scph, "RESTART_SCPH", scph_var_dict);
    assign_val(maxiter, "MAXITER", scph_var_dict);
    assign_val(mixalpha, "MIXALPHA", scph_var_dict);
    assign_val(selenergy_offdiagonal, "SELF_OFFDIAG", scph_var_dict);
    assign_val(ialgo_scph, "IALGO", scph_var_dict);
    assign_val(tolerance_scph, "TOL_SCPH", scph_var_dict);
    assign_val(lower_temp, "LOWER_TEMP", scph_var_dict);
    assign_val(warm_start, "WARMSTART", scph_var_dict);
    assign_val(bubble, "BUBBLE", scph_var_dict);

    auto str_tmp = scph_var_dict["KMESH_SCPH"];

    if (!str_tmp.empty()) {

        std::istringstream is(str_tmp);

        while (true) {
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

        while (true) {
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

    for (auto i = 0; i < 3; ++i) {
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
    scph->bubble = bubble;

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
            "REALPART", "ISOTOPE", "ISOFACT",
            "FSTATE_W", "FSTATE_K", "PRINTMSD", "DOS", "PDOS", "TDOS",
            "GRUNEISEN", "NEWFCS", "DELTA_A", "ANIME", "ANIME_CELLSIZE",
            "ANIME_FORMAT", "ANIME_FRAMES", "SPS", "PRINTV3", "PRINTPR",
            "FC2_EWALD", "KAPPA_SPEC", "SELF_W", "UCORR", "SHIFT_UCORR",
        "KAPPA_COHERENT",
        "DIELEC", "SELF_ENERGY", "PRINTV4", "ZMODE", "PROJECTION_AXES"
    };

#ifdef _FE_BUBBLE
    input_list.push_back("FE_BUBBLE");
#endif

    unsigned int cellsize[3] = {1, 1, 1};
    int shift_ucorr[3] = {0, 0, 0};
    double anime_kpoint_double[3] = {0., 0., 0.};

    double *isotope_factor = nullptr;
    std::string ks_input, anime_format;
    std::map<std::string, std::string> analysis_var_dict;
    std::vector<std::string> isofact_v, anime_kpoint, anime_cellsize;
    std::vector<std::string> projection_axes;

    // Default values

    auto print_xsf = false;
    auto print_anime = false;

    auto print_vel = false;
    auto print_evec = false;
    auto print_msd = false;
    auto print_ucorr = false;
    auto anime_frames = 20;

    bool compute_dos = true;
    bool projected_dos = false;
    bool two_phonon_dos = false;
    int scattering_phase_space = 0;
    bool print_gruneisen = false;
    bool print_newfcs = false;
    int print_V3 = 0;
    int print_V4 = 0;
    bool participation_ratio = false;

    auto delta_a = 0.001;

    auto quartic_mode = 0;
    auto calc_realpart = false;
    auto include_isotope = 0;
    auto fstate_omega = false;
    auto fstate_k = false;
    auto bubble_omega = false;

    auto calculate_kappa_spec = 0;

    auto print_fc2_ewald = false;
    auto print_self_consistent_fc2 = false;
    auto calc_FE_bubble = false;
    auto calc_coherent = 0;

    auto calculate_dielectric_constant = 0;
    auto calc_selfenergy = 0;
    auto print_zmode = false;

    auto do_projection = false;

    // Assign values to variables

    if (!use_default_values) {
        get_var_dict(input_list, analysis_var_dict);

        assign_val(print_vel, "PRINTVEL", analysis_var_dict);
        assign_val(print_evec, "PRINTEVEC", analysis_var_dict);
        assign_val(print_msd, "PRINTMSD", analysis_var_dict);
        assign_val(print_ucorr, "UCORR", analysis_var_dict);

        assign_val(compute_dos, "DOS", analysis_var_dict);
        assign_val(projected_dos, "PDOS", analysis_var_dict);
        assign_val(two_phonon_dos, "TDOS", analysis_var_dict);
        assign_val(scattering_phase_space, "SPS", analysis_var_dict);
        assign_val(print_gruneisen, "GRUNEISEN", analysis_var_dict);
        assign_val(print_newfcs, "NEWFCS", analysis_var_dict);
        assign_val(delta_a, "DELTA_A", analysis_var_dict);

        assign_val(quartic_mode, "QUARTIC", analysis_var_dict);
        assign_val(calc_realpart, "REALPART", analysis_var_dict);
        assign_val(include_isotope, "ISOTOPE", analysis_var_dict);
        assign_val(fstate_omega, "FSTATE_W", analysis_var_dict);
        assign_val(fstate_k, "FSTATE_K", analysis_var_dict);
        assign_val(ks_input, "KS_INPUT", analysis_var_dict);
        assign_val(calculate_kappa_spec, "KAPPA_SPEC", analysis_var_dict);
        assign_val(calc_coherent, "KAPPA_COHERENT", analysis_var_dict);
        assign_val(bubble_omega, "SELF_W", analysis_var_dict);
        assign_val(calc_selfenergy, "SELF_ENERGY", analysis_var_dict);

        assign_val(print_xsf, "PRINTXSF", analysis_var_dict);
        assign_val(print_V3, "PRINTV3", analysis_var_dict);
        assign_val(print_V4, "PRINTV4", analysis_var_dict);
        assign_val(participation_ratio, "PRINTPR", analysis_var_dict);
        assign_val(print_fc2_ewald, "FC2_EWALD", analysis_var_dict);
        assign_val(calculate_dielectric_constant, "DIELEC", analysis_var_dict);
        assign_val(print_zmode, "ZMODE", analysis_var_dict);
#ifdef _FE_BUBBLE
        assign_val(calc_FE_bubble, "FE_BUBBLE", analysis_var_dict);
#endif

        if (analysis_var_dict.find("ANIME") == analysis_var_dict.end()) {
            print_anime = false;
        } else {
            print_anime = true;
        }

        if (analysis_var_dict.find("PROJECTION_AXES") != analysis_var_dict.end()) {
            do_projection = true;
        }
    }

    if (include_isotope) {

        if (!analysis_var_dict["ISOFACT"].empty()) {
            split_str_by_space(analysis_var_dict["ISOFACT"], isofact_v);

            if (isofact_v.size() != system->nkd) {
                error->exit("parse_analysis_vars",
                            "The number of entries for ISOFACT is inconsistent with NKD");
            } else {
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
        for (i = 0; i < 3; ++i) {
            anime_kpoint_double[i] = my_cast<double>(anime_kpoint[i]);
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
                            "Please give positive integers for ANIME_CELLSIZE.");
            }
        }

        assign_val(anime_format, "ANIME_FORMAT", analysis_var_dict);
        std::transform(anime_format.begin(), anime_format.end(),
                       anime_format.begin(), toupper);

        if (anime_format.empty()) anime_format = "XYZ";

        if (anime_format != "XSF" && anime_format != "AXSF" && anime_format != "XYZ") {
            error->exit("parse_analysis_vars", "Invalid ANIME_FORMAT");
        }

        assign_val(anime_frames, "ANIME_FRAMES", analysis_var_dict);
    }

    if (print_ucorr) {
        std::string str_shift_ucorr;
        std::vector<std::string> list_shift_ucorr;
        assign_val(str_shift_ucorr, "SHIFT_UCORR", analysis_var_dict);

        if (!str_shift_ucorr.empty()) {
            split_str_by_space(str_shift_ucorr, list_shift_ucorr);
            if (list_shift_ucorr.size() != 3) {
                error->exit("parse_analysis_vars",
                            "The number of entries for SHIFT_UCORR must be 3.");
            }

            for (i = 0; i < 3; ++i) {
                try {
                    shift_ucorr[i] = boost::lexical_cast<int>(list_shift_ucorr[i]);
                }
                catch (std::exception &e) {
                    std::cout << e.what() << std::endl;
                    error->exit("parse_analysis_vars",
                                "SHIFT_UCORR must be an array of integers.");
                }
            }
        }
    }

    if (do_projection) {
        std::string str_projection_axes;
        assign_val(str_projection_axes, "PROJECTION_AXES", analysis_var_dict);
        std::vector<double> direction(3);
        std::vector<std::vector<double>> projection_directions;
        if (!str_projection_axes.empty()) {
            std::vector<std::string> str_projection_each, str_vec;
            boost::split(str_projection_each, str_projection_axes, boost::is_any_of(","));

            if (str_projection_each.size() > 2) {
                error->warn("parse_analysis_vars",
                            "Too many entries for PROJECTION_AXES. Only the first two will be used.");
            }

            for (i = 0; i < str_projection_each.size(); ++i) {
                split_str_by_space(str_projection_each[i], str_vec);
                if (str_vec.size() != 3) {
                    error->exit("parse_analysis_vars",
                                "The number of entries for each vector in PROJECTION_AXES must be 3.");
                }
                for (auto j = 0; j < 3; ++j) {
                    try {
                        direction[j] = boost::lexical_cast<double>(str_vec[j]);
                    } catch (std::exception &e) {
                        std::cout << e.what() << std::endl;
                        error->exit("parse_analysis_vars",
                                    "subset of PROJECTION_AXES must be an array of doubles.");
                    }
                }
                projection_directions.push_back(direction);
            }
        }

        dynamical->set_projection_directions(projection_directions);
    }

    // Copy the values to appropriate classes

    phonon_velocity->print_velocity = print_vel;
    dynamical->print_eigenvectors = print_evec;
    dynamical->participation_ratio = participation_ratio;
//    writes->print_xsf = print_xsf;
//    writes->print_anime = print_anime;
//    writes->print_ucorr = print_ucorr;

    writes->setWriteOptions(print_msd,
                            print_xsf,
                            print_anime,
                            anime_format,
                            anime_frames,
                            cellsize,
                            anime_kpoint_double,
                            print_ucorr,
                            shift_ucorr);

//    if (print_anime) {
//        for (i = 0; i < 3; ++i) {
//            writes->anime_kpoint[i] = my_cast<double>(anime_kpoint[i]);
//            writes->anime_cellsize[i] = cellsize[i];
//        }
//        writes->anime_format = anime_format;
//    }
//
//    writes->print_msd = print_msd;

    dos->compute_dos = compute_dos;
    dos->projected_dos = projected_dos;
    dos->two_phonon_dos = two_phonon_dos;
    dos->scattering_phase_space = scattering_phase_space;

    conductivity->calc_kappa_spec = calculate_kappa_spec;
    conductivity->calc_coherent = calc_coherent;
    anharmonic_core->quartic_mode = quartic_mode;
    dielec->calc_dielectric_constant = calculate_dielectric_constant;

    mode_analysis->ks_input = ks_input;
    mode_analysis->calc_realpart = calc_realpart;
    mode_analysis->calc_fstate_omega = fstate_omega;
    mode_analysis->calc_fstate_k = fstate_k;
    mode_analysis->print_V3 = print_V3;
    mode_analysis->print_V4 = print_V4;
    mode_analysis->spectral_func = bubble_omega;
    mode_analysis->calc_selfenergy = calc_selfenergy;
    isotope->include_isotope = include_isotope;

    gruneisen->print_gruneisen = print_gruneisen;
    gruneisen->print_newfcs = print_newfcs;
    gruneisen->delta_a = delta_a;
    thermodynamics->calc_FE_bubble = calc_FE_bubble;

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

            trim_if(line_wo_comment, boost::is_any_of("\t\r\n "));

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

            trim_if(line_wo_comment, boost::is_any_of("\t\r\n "));

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
        split(line_split, line,
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

    int kpmode;
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

            trim_if(line_wo_comment, boost::is_any_of("\t\r\n "));

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

            trim_if(line_wo_comment, boost::is_any_of("\t\r\n "));

            if (line_wo_comment.empty()) continue;
            if (is_endof_entry(line_wo_comment)) break;

            line_vec.push_back(line_wo_comment);
        }
    }

    for (int i = 0; i < line_vec.size(); ++i) {
        line = line_vec[i];
        split(kpelem, line, boost::is_any_of("\t "), boost::token_compress_on);

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

            kpoint->kpInp.emplace_back(kpelem);
        }
    }

    kpoint->kpoint_mode = kpmode;
}


int Input::locate_tag(const std::string &key)
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

    }
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
                std::string str_tmp = trim(*it);
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


bool Input::is_endof_entry(const std::string &str) const
{
    return str[0] == '/';
}

void Input::split_str_by_space(const std::string &str,
                               std::vector<std::string> &str_vec) const
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

template<typename T>
void Input::assign_val(T &val,
                       const std::string &key,
                       std::map<std::string, std::string> dict)
{
    // Assign a value to the variable "key" using the boost::lexica_cast.

    if (!dict[key].empty()) {
        try {
            val = boost::lexical_cast<T>(dict[key]);
        }
        catch (std::exception &e) {
            std::cout << e.what() << std::endl;
            std::string str_tmp = "Invalid entry for the " + key + " tag.\n";
            str_tmp += " Please check the input value.";
            error->exit("assign_val", str_tmp.c_str());
        }
    }
}


template<typename T_to, typename T_from>
T_to Input::my_cast(T_from const &x)
{
    std::stringstream ss;
    T_to ret;

    ss << x;
    ss >> ret;

    return ret;
}
