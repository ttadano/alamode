/*
 parsephon.cpp

 Copyright (c) 2014 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#include "mpi_common.h"
#include "parsephon.h"
#include "error.h"
#include "gruneisen.h"
#include "system.h"
#include "kpoint.h"
#include "fcs_phonon.h"
#include "dynamical.h"
#include "write_phonons.h"
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include "relaxation.h"
#include "conductivity.h"
#include "symmetry_core.h"
#include <istream>
#include <map>
#include <vector>
#include "phonon_dos.h"
#include <sstream>
#include <sys/stat.h>
#include "memory.h"
#include "isotope.h"
#include "phonon_velocity.h"
#include "integration.h"

#ifdef _USE_BOOST
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#endif

using namespace PHON_NS;

Input::Input(PHON *phon): Pointers(phon) {}

Input::~Input() {
    if (!from_stdin && mympi->my_rank == 0) ifs_input.close();
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

    if (!locate_tag("&general")) error->exit("parse_input", "&general entry not found in the input file");
    parse_general_vars();
   
    if (!locate_tag("&cell")) error->exit("parse_input", "&cell entry not found in the input file");
    parse_cell_parameter();
    
    bool use_defaults_for_analysis;
    if (!locate_tag("&analysis")) {
        use_defaults_for_analysis = true;
    } else {
        use_defaults_for_analysis = false;
    }
    parse_analysis_vars(use_defaults_for_analysis);

    if (!locate_tag("&kpoint")) error->exit("parse_input", "&kpoint entry not found in the input file");
    parse_kpoints();

}

void Input::parse_general_vars()
{
    int i;
    int nsym, celldim[3], nbands, ismear, nkd;
    double *masskd;
    double Tmin, Tmax, dT, na_sigma, epsilon;
    double emin, emax, delta_e, delta_a;
    double tolerance;
    bool printsymmetry;
    bool nonanalytic, restart;
    bool sym_time_reversal, use_triplet_symmetry;

    struct stat st;
    std::string prefix, mode, fcsinfo;
    std::string borninfo, file_result;
    std::string *kdname;
    std::string str_tmp;
    std::string str_allowed_list = "PREFIX MODE NSYM TOLERANCE PRINTSYM CELLDIM FCSXML TMIN TMAX DT \
                                   NBANDS NONANALYTIC BORNINFO NA_SIGMA ISMEAR EPSILON EMIN EMAX DELTA_E \
                                   DELTA_A RESTART TREVSYM NKD KD MASS TRISYM";
    std::string str_no_defaults = "PREFIX MODE FCSXML NKD KD MASS";
    std::vector<std::string> no_defaults, celldim_v;
    std::vector<std::string> kdname_v, masskd_v;
    std::map<std::string, std::string> general_var_dict;

    if (from_stdin) {
        std::cin.ignore();
    } else {
        ifs_input.ignore();
    }

    get_var_dict(str_allowed_list, general_var_dict);
#if _USE_BOOST
    boost::split(no_defaults, str_no_defaults, boost::is_space());
#else 
    no_defaults = my_split(str_no_defaults, ' ');
#endif

    for (std::vector<std::string>::iterator it = no_defaults.begin(); it != no_defaults.end(); ++it){
        if (general_var_dict.find(*it) == general_var_dict.end()) {
            error->exit("parse_general_vars", "The following variable is not found in &general input region: ", (*it).c_str());
        }
    }

    prefix = general_var_dict["PREFIX"];
    mode = general_var_dict["MODE"];

    file_result = prefix + ".result";

#ifdef _USE_BOOST
    boost::to_upper(mode);
    nsym = boost::lexical_cast<int>(general_var_dict["NSYM"]);
#else
    std::transform (mode.begin(), mode.end(), mode.begin(), toupper);
    nsym = my_cast<int>(general_var_dict["NSYM"]);
#endif

    fcsinfo = general_var_dict["FCSXML"];
    assign_val(nkd, "NKD", general_var_dict);

    split_str_by_space(general_var_dict["KD"], kdname_v);

    if (kdname_v.size() != nkd) {
        error->exit("parse_general_vars", "The number of entries for KD is inconsistent with NKD");
    } else {
        memory->allocate(kdname, nkd);
        for (i = 0; i < nkd; ++i){
            kdname[i] = kdname_v[i];
        }
    }

    split_str_by_space(general_var_dict["MASS"], masskd_v);

    if (masskd_v.size() != nkd) {
        error->exit("parse_general_vars", "The number of entries for MASS is inconsistent with NKD");
    } else {
        memory->allocate(masskd, nkd);
        for (i = 0; i < nkd; ++i) {
            masskd[i] = my_cast<double>(masskd_v[i]);
        }
    }

    // Default values

    Tmin = 0.0;
    Tmax = 1000.0;
    dT = 10.0;

    emin = 0.0;
    emax = 1000.0;
    delta_e = 10.0;

    nonanalytic = false;
    nsym = 0;
    tolerance = 1.0e-8;
    printsymmetry = false;
    sym_time_reversal = false;
    use_triplet_symmetry = true;

    // if file_result exists in the current directory, restart mode will be automatically turned on.

    if (stat(file_result.c_str(), &st) == 0) {
        restart = true;
    } else {
        restart = false;
    }

    nbands = -1;
    borninfo = "";

    ismear = -1;
    epsilon = 10.0;
    na_sigma = 0.1;

    delta_a = 0.001;

    for (i = 0; i < 3; ++i) celldim[i] = 0;

    // Assign given values

    assign_val(Tmin, "TMIN", general_var_dict);
    assign_val(Tmax, "TMAX", general_var_dict);
    assign_val(dT, "DT", general_var_dict);

    assign_val(emin, "EMIN", general_var_dict);
    assign_val(emax, "EMAX", general_var_dict);
    assign_val(delta_e, "DELTA_E", general_var_dict);

    assign_val(sym_time_reversal, "TREVSYM", general_var_dict);
    assign_val(tolerance, "TOLERANCE", general_var_dict);
    assign_val(printsymmetry, "PRINTSYM", general_var_dict);

    assign_val(nonanalytic, "NONANALYTIC", general_var_dict);
    assign_val(restart, "RESTART", general_var_dict);

    assign_val(nbands, "NBANDS", general_var_dict);
    assign_val(borninfo, "BORNINFO", general_var_dict);

    assign_val(ismear, "ISMEAR", general_var_dict);
    assign_val(epsilon, "EPSILON", general_var_dict);
    assign_val(na_sigma, "NA_SIGMA", general_var_dict);

    assign_val(delta_a, "DELTA_A", general_var_dict);

    assign_val(use_triplet_symmetry, "TRISYM", general_var_dict);


    str_tmp = general_var_dict["CELLDIM"];

    if (!str_tmp.empty()) {

        std::istringstream is(str_tmp);

        while (1) {
            str_tmp.clear();
            is >> str_tmp;
            if (str_tmp.empty()) {
                break;
            }
            celldim_v.push_back(str_tmp);
        }

        if (celldim_v.size() != 3) {
            error->exit("parse_general_vars", "The number of entries for CELLDIM has to be 3.");
        }

        for (i = 0; i < 3; ++i) {
#ifdef _USE_BOOST
            celldim[i] = boost::lexical_cast<int>(celldim_v[i]);
#else
            celldim[i] = my_cast<int>(celldim_v[i]);
#endif
        }
    }

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

    memory->allocate(system->mass_kd, nkd);
    memory->allocate(system->symbol_kd, nkd);
    for (i = 0; i < nkd; ++i) {
        system->mass_kd[i] = masskd[i];
        system->symbol_kd[i] = kdname[i];
    }
    memory->deallocate(masskd);
    memory->deallocate(kdname);

    dos->emax = emax;
    dos->emin = emin;
    dos->delta_e = delta_e;


    for (i = 0; i < 3; ++i) {
        system->cell_dimension[i] = celldim[i];
    }

    dynamical->nonanalytic = nonanalytic;
    dynamical->na_sigma = na_sigma;
    writes->nbands = nbands;
    dynamical->file_born = borninfo;
    integration->epsilon = epsilon;
    fcs_phonon->file_fcs = fcsinfo;
    gruneisen->delta_a = delta_a;
    integration->ismear = ismear;
    relaxation->use_triplet_symmetry = use_triplet_symmetry;

    general_var_dict.clear();
}

void Input::parse_analysis_vars(const bool use_default_values)
{
    int i;

    std::string str_allowed_list = "LCLASSICAL PRINTEVEC PRINTXSF PRINTVEL QUARTIC KS_INPUT ATOMPROJ REALPART \
                                   ISOTOPE ISOFACT FSTATE_W FSTATE_K PRINTMSD PDOS TDOS GRUNEISEN NEWFCS";

    bool include_isotope;
    bool fstate_omega, fstate_k;
    bool lclassical;
    bool quartic_mode, ks_analyze_mode, atom_project_mode, calc_realpart;
    bool print_vel, print_evec, print_xsf, print_msd;
    bool projected_dos, print_gruneisen, print_newfcs;
    bool two_phonon_dos;

    double *isotope_factor;
    std::string ks_input;
    std::map<std::string, std::string> analysis_var_dict;
    std::vector<std::string> isofact_v;

    print_xsf = false;
    print_vel = false;
    print_evec = false;
    print_msd = false;

    projected_dos = false;
    two_phonon_dos = false;
    print_gruneisen = false;
    print_newfcs = false;

    lclassical = false;
    quartic_mode = false;
    ks_analyze_mode = false;
    atom_project_mode = false;
    calc_realpart = false;
    include_isotope = false;
    fstate_omega = false;
    fstate_k = false;

    if (!use_default_values) {
        get_var_dict(str_allowed_list, analysis_var_dict);

        assign_val(print_xsf, "PRINTXSF", analysis_var_dict);
        assign_val(print_vel, "PRINTVEL", analysis_var_dict);
        assign_val(print_evec, "PRINTEVEC", analysis_var_dict);
        assign_val(print_msd, "PRINTMSD", analysis_var_dict);

        assign_val(projected_dos, "PDOS", analysis_var_dict);
        assign_val(two_phonon_dos, "TDOS", analysis_var_dict);
        assign_val(print_gruneisen, "GRUNEISEN", analysis_var_dict);
        assign_val(print_newfcs, "NEWFCS", analysis_var_dict);

        assign_val(lclassical, "LCLASSICAL", analysis_var_dict);
        assign_val(quartic_mode, "QUARTIC", analysis_var_dict);
        assign_val(atom_project_mode, "ATOMPROJ", analysis_var_dict);
        assign_val(calc_realpart, "REALPART", analysis_var_dict);
        assign_val(include_isotope, "ISOTOPE", analysis_var_dict);
        assign_val(fstate_omega, "FSTATE_W", analysis_var_dict);
        assign_val(fstate_k, "FSTATE_K", analysis_var_dict);
        assign_val(ks_input, "KS_INPUT", analysis_var_dict);
    }

    if (include_isotope) {
        split_str_by_space(analysis_var_dict["ISOFACT"], isofact_v);

        if (isofact_v.size() != system->nkd) {
            error->exit("parse_general_vars", "The number of entries for ISOFACT is inconsistent with NKD");
        } else {
            memory->allocate(isotope_factor, system->nkd);
            for (i = 0; i < system->nkd; ++i){
                isotope_factor[i] = my_cast<double>(isofact_v[i]);
            }
        }
    }

    phonon_velocity->print_velocity = print_vel;
    dynamical->print_eigenvectors = print_evec;
    writes->writeanime = print_xsf;
    writes->print_msd = print_msd;

    dos->projected_dos = projected_dos;
    dos->two_phonon_dos = two_phonon_dos;

    conductivity->use_classical_Cv = lclassical;
    relaxation->quartic_mode = quartic_mode;
    relaxation->atom_project_mode = atom_project_mode;
    relaxation->calc_realpart = calc_realpart;
    relaxation->calc_fstate_omega = fstate_omega;
    relaxation->calc_fstate_k = fstate_k;
    isotope->include_isotope = include_isotope;
    relaxation->ks_input = ks_input;

    gruneisen->print_gruneisen = print_gruneisen;
    gruneisen->print_newfcs = print_newfcs;


    if (include_isotope) {
        memory->allocate(isotope->isotope_factor, system->nkd);
        for (i = 0; i < system->nkd; ++i) {
            isotope->isotope_factor[i] = isotope_factor[i];
        }
        memory->deallocate(isotope_factor);
    }

    analysis_var_dict.clear();
}

void Input::parse_cell_parameter()
{

    int i, j;
    double a;
    double lavec_tmp[3][3];

    if (from_stdin) {
        std::cin >> a;

        std::cin >> lavec_tmp[0][0] >> lavec_tmp[1][0] >> lavec_tmp[2][0]; // a1
        std::cin >> lavec_tmp[0][1] >> lavec_tmp[1][1] >> lavec_tmp[2][1]; // a2
        std::cin >> lavec_tmp[0][2] >> lavec_tmp[1][2] >> lavec_tmp[2][2]; // a3
    } else {
        ifs_input >> a;

        ifs_input >> lavec_tmp[0][0] >> lavec_tmp[1][0] >> lavec_tmp[2][0]; // a1
        ifs_input >> lavec_tmp[0][1] >> lavec_tmp[1][1] >> lavec_tmp[2][1]; // a2
        ifs_input >> lavec_tmp[0][2] >> lavec_tmp[1][2] >> lavec_tmp[2][2]; // a3
    }

    for (i = 0; i < 3; ++i) {
        for (j = 0; j < 3; ++j) {
            system->lavec_p[i][j] = a * lavec_tmp[i][j];
        }
    }
}

void Input::parse_kpoints() {

    int kpmode;
    std::string line, str_tmp;
    std::vector<std::string> kpelem;

    if (from_stdin) {

        std::cin >> kpmode;

        if (!(kpmode >= 0 && kpmode <= 3)) {
            error->exit("parse_kpoints", "Invalid KPMODE");
        }

        std::cin.ignore();

        while (std::getline(std::cin, line)) {

            if (is_endof_entry(line)) {
                break;
            }
            kpelem.clear();

            std::istringstream is(line);

            while (1) {
                str_tmp.clear();
                is >> str_tmp;
                if (str_tmp.empty()) {
                    break;
                }
                kpelem.push_back(str_tmp);
            }

            if (kpmode == 0 && kpelem.size() != 3) {
                error->exit("parse_kpoints", "The number of entries has to be 3 when KPMODE = 0");
            }
            if (kpmode == 1 && kpelem.size() != 9) {
                error->exit("parse_kpoints", "The number of entries has to be 9 when KPMODE = 1");
            }
            if (kpmode == 2 && kpelem.size() != 3) {
                error->exit("parse_kpoints", "The number of entries has to be 3 when KPMODE = 2");
            }
            if (kpmode == 3 && kpelem.size() != 8) {
                error->exit("parse_kpoints", "The number of entries has to be 8 when KPMODE = 3");
            }

            kpoint->kpInp.push_back(kpelem);
        }
    } else {

        ifs_input >> kpmode;

        if (!(kpmode >= 0 && kpmode <= 3)) {
            error->exit("parse_kpoints", "Invalid KPMODE");
        } 

        if (!relaxation->calc_fstate_k && kpmode == 3) {
            error->exit("parse_kpoints", "KPMODE = 3 is valid only when FSTATE_K is true.");
        }

        ifs_input.ignore();

        while (std::getline(ifs_input, line)) {

            if (is_endof_entry(line)) {
                break;
            }
            kpelem.clear();

            std::istringstream is(line);

            while (1) {
                str_tmp.clear();
                is >> str_tmp;
                if (str_tmp.empty()) {
                    break;
                }
                kpelem.push_back(str_tmp);
            }

            if (kpmode == 0 && kpelem.size() != 3) {
                error->exit("parse_kpoints", "The number of entries has to be 3 when KPMODE = 0");
            }
            if (kpmode == 1 && kpelem.size() != 9) {
                error->exit("parse_kpoints", "The number of entries has to be 9 when KPMODE = 1");
            }
            if (kpmode == 2 && kpelem.size() != 3) {
                error->exit("parse_kpoints", "The number of entries has to be 3 when KPMODE = 2");
            }
            if (kpmode == 3 && kpelem.size() != 8) {
                error->exit("parse_kpoints", "The number of entries has to be 8 when KPMODE = 3");
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
            if (line == key){
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
            if (line == key){
                ret = 1;
                break;
            }
        }
        return ret; 
    }

}

void Input::get_var_dict(const std::string keywords, std::map<std::string, std::string> &var_dict)
{

    std::string line, key, val;
    std::string line_wo_comment, line_tmp;
    std::string::size_type pos_first_comment_tag;
    std::vector<std::string> str_entry, str_varval;


    std::set<std::string> keyword_set;
#ifdef _USE_BOOST
    boost::split(keyword_set, keywords, boost::is_space());
#else
    std::vector<std::string> strvec_tmp;
    strvec_tmp = my_split(keywords, ' ');
    for (std::vector<std::string>::iterator it = strvec_tmp.begin(); it != strvec_tmp.end(); ++it) {
        keyword_set.insert(*it);
    }
    strvec_tmp.clear();
#endif

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

            //	std::cout << line_wo_comment << std::endl;

            // Split the input line by ';'
#ifdef _USE_BOOST
            boost::split(str_entry, line_wo_comment, boost::is_any_of(";"));
#else
            str_entry = my_split(line_wo_comment, ';');
#endif


            for (std::vector<std::string>::iterator it = str_entry.begin(); it != str_entry.end(); ++it) {

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
                        error->exit("get_var_dict", "Unacceptable format");
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
                        error->exit("get_var_dict", "Invalid variable found");
                    }

                    if (var_dict.find(key) != var_dict.end()) {
                        std::cout << "Variable " << key << " appears twice in the input file." << std::endl;
                        error->exit("get_var_dict", "Redundant input parameter");
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

            //	std::cout << line_wo_comment << std::endl;

            // Split the input line by ';'
#ifdef _USE_BOOST
            boost::split(str_entry, line_wo_comment, boost::is_any_of(";"));
#else
            str_entry = my_split(line_wo_comment, ';');
#endif
            for (std::vector<std::string>::iterator it = str_entry.begin(); it != str_entry.end(); ++it) {

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
                        error->exit("get_var_dict", "Unacceptable format");
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
                        error->exit("get_var_dict", "Invalid variable found");
                    }

                    if (var_dict.find(key) != var_dict.end()) {
                        std::cout << "Variable " << key << " appears twice in the input file." << std::endl;
                        error->exit("get_var_dict", "Redundant input parameter");
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


bool Input::is_endof_entry(std::string str)
{

    if (str[0] == '/') {
        return true;
    } else {
        return false;
    }
}

void Input::split_str_by_space(const std::string str, std::vector<std::string> &str_vec)
{

    std::string str_tmp;
    std::istringstream is(str);

    str_vec.clear();

    while(1) {
        str_tmp.clear();
        is >> str_tmp;
        if (str_tmp.empty()) {
            break;
        }
        str_vec.push_back(str_tmp);
    }
    str_tmp.clear();
}

template<typename T> void Input::assign_val(T &val, std::string key, std::map<std::string, std::string> dict)
{

    if (!dict[key].empty()) {
#ifdef _USE_BOOST
        val = boost::lexical_cast<T>(dict[key]);
#else
        val = my_cast<T>(dict[key]);
#endif
    }
}

template<typename T_to, typename T_from> T_to Input::my_cast(T_from const &x)
{
    std::stringstream ss;
    T_to ret;

    ss << x;
    ss >> ret;

    return ret;
}
