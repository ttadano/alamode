/*
 input.cpp

 Copyright (c) 2014, 2015, 2016 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#include "input.h"
#include <iostream>
#include <iomanip>
#include <string>
#include "memory.h"
#include "files.h"
#include "interaction.h"
#include "system.h"
#include "symmetry.h"
#include "error.h"
#include "fitting.h"
#include "constraint.h"
#include "writes.h"
#include "patterndisp.h"
#include <algorithm>
#include <map>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

using namespace ALM_NS;

Input::Input(ALM *alm, int narg, char **arg): Pointers(alm)
{
}

Input::~Input()
{
}

void Input::parse_input(int narg, char **arg)
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

    if (!locate_tag("&general")) {
        error->exit("parse_input",
                    "&general entry not found in the input file");
    }
    parse_general_vars();

    if (!locate_tag("&cell")) {
        error->exit("parse_input",
                    "&cell entry not found in the input file");
    }
    parse_cell_parameter();

    if (!locate_tag("&interaction")) {
        error->exit("parse_input",
                    "&interaction entry not found in the input file");
    }
    parse_interaction_vars();

    if (!locate_tag("&cutoff")) {
        error->exit("parse_input",
                    "&cutoff entry not found in the input file");
    }
    parse_cutoff_radii();

    if (alm->mode == "fitting") {
        if (!locate_tag("&fitting")) {
            error->exit("parse_input",
                        "&fitting entry not found in the input file");
        }
        parse_fitting_vars();
    }

    if (!locate_tag("&position")) {
        error->exit("parse_input",
                    "&position entry not found in the input file");
    }
    parse_atomic_positions();
}

void Input::parse_general_vars()
{
    int i, j;
    std::string prefix, mode, str_tmp, str_disp_basis;
    int nat, nkd, nsym;
    int is_printsymmetry, is_periodic[3];
    int icount, ncount;
    bool trim_dispsign_for_evenfunc;
    bool lspin;
    bool print_hessian;
    int noncollinear, trevsym;
    std::string *kdname;
    double **magmom, magmag;
    double tolerance;
    double tolerance_constraint;

    std::vector<std::string> kdname_v, periodic_v, magmom_v, str_split;
    std::string str_allowed_list = "PREFIX MODE NAT NKD NSYM KD PERIODIC PRINTSYM TOLERANCE DBASIS TRIMEVEN\
                                   MAGMOM NONCOLLINEAR TREVSYM HESSIAN TOL_CONST";
    std::string str_no_defaults = "PREFIX MODE NAT NKD KD";
    std::vector<std::string> no_defaults;
    std::map<std::string, std::string> general_var_dict;

    if (from_stdin) {
        std::cin.ignore();
    } else {
        ifs_input.ignore();
    }

    get_var_dict(str_allowed_list, general_var_dict);

    boost::split(no_defaults, str_no_defaults, boost::is_space());

    for (std::vector<std::string>::iterator it = no_defaults.begin();
         it != no_defaults.end(); ++it) {
        if (general_var_dict.find(*it) == general_var_dict.end()) {
            error->exit("parse_general_vars",
                        "The following variable is not found in &general input region: ",
                        (*it).c_str());
        }
    }

    prefix = general_var_dict["PREFIX"];
    mode = general_var_dict["MODE"];

    std::transform(mode.begin(), mode.end(), mode.begin(), tolower);
    if (mode != "fitting" && mode != "suggest") {
        error->exit("parse_general_vars", "Invalid MODE variable");
    }

    assign_val(nat, "NAT", general_var_dict);
    assign_val(nkd, "NKD", general_var_dict);

    if (general_var_dict["NSYM"].empty()) {
        nsym = 0;
    } else {
        assign_val(nsym, "NSYM", general_var_dict);
    }

    if (general_var_dict["PRINTSYM"].empty()) {
        is_printsymmetry = 0;
    } else {
        assign_val(is_printsymmetry, "PRINTSYM", general_var_dict);
    }

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

    split_str_by_space(general_var_dict["PERIODIC"], periodic_v);

    if (periodic_v.size() == 0) {
        for (i = 0; i < 3; ++i) {
            is_periodic[i] = 1;
        }
    } else if (periodic_v.size() == 3) {
        for (i = 0; i < 3; ++i) {
            try {
                is_periodic[i] = boost::lexical_cast<int>(periodic_v[i]);
            }
            catch (std::exception &e) {
                std::cout << e.what() << std::endl;
                error->exit("parse_general_vars",
                            "The PERIODIC tag must be a set of integers.");
            }
        }
    } else {
        error->exit("parse_general_vars",
                    "Invalid number of entries for PERIODIC");
    }

    if (general_var_dict["TOLERANCE"].empty()) {
        tolerance = 1.0e-6;
    } else {
        assign_val(tolerance, "TOLERANCE", general_var_dict);
    }

    if (general_var_dict["TOL_CONST"].empty()) {
        tolerance_constraint = eps6;
    } else {
        assign_val(tolerance_constraint, "TOL_CONST", general_var_dict);
    }

    // Convert MAGMOM input to array
    memory->allocate(magmom, nat, 3);
    lspin = false;

    for (i = 0; i < nat; ++i) {
        for (j = 0; j < 3; ++j) {
            magmom[i][j] = 0.0;
        }
    }

    if (general_var_dict["NONCOLLINEAR"].empty()) {
        noncollinear = 0;
    } else {
        assign_val(noncollinear, "NONCOLLINEAR", general_var_dict);
    }
    if (general_var_dict["TREVSYM"].empty()) {
        trevsym = 1;
    } else {
        assign_val(trevsym, "TREVSYM", general_var_dict);
    }
    if (general_var_dict["HESSIAN"].empty()) {
        print_hessian = false;
    } else {
        assign_val(print_hessian, "HESSIAN", general_var_dict);
    }

    if (!general_var_dict["MAGMOM"].empty()) {
        lspin = true;

        if (noncollinear) {
            icount = 0;
            split_str_by_space(general_var_dict["MAGMOM"], magmom_v);
            for (std::vector<std::string>::const_iterator it = magmom_v.begin();
                 it != magmom_v.end(); ++it) {
                if ((*it).find("*") != std::string::npos) {
                    error->exit("parse_general_vars",
                                "Wild card '*' is not supported when NONCOLLINEAR = 1.");
                } else {
                    magmag = boost::lexical_cast<double>((*it));
                    if (icount / 3 >= nat) {
                        error->exit("parse_general_vars", "Too many entries for MAGMOM.");
                    }
                    magmom[icount / 3][icount % 3] = magmag;
                    ++icount;
                }
            }

            if (icount != 3 * nat) {
                error->exit("parse_general_vars",
                            "Number of entries for MAGMOM must be 3*NAT when NONCOLLINEAR = 1.");
            }
        } else {
            icount = 0;
            split_str_by_space(general_var_dict["MAGMOM"], magmom_v);
            for (std::vector<std::string>::const_iterator it = magmom_v.begin();
                 it != magmom_v.end(); ++it) {

                if ((*it).find("*") != std::string::npos) {
                    if ((*it) == "*") {
                        error->exit("parse_general_vars",
                                    "Please place '*' without space for the MAGMOM-tag.");
                    }
                    boost::split(str_split, (*it), boost::is_any_of("*"));
                    if (str_split.size() != 2) {
                        error->exit("parse_general_vars",
                                    "Invalid format for the MAGMOM-tag.");
                    } else {
                        if (str_split[0].empty() || str_split[1].empty()) {
                            error->exit("parse_general_vars",
                                        "Please place '*' without space for the MAGMOM-tag.");
                        }
                        try {
                            magmag = boost::lexical_cast<double>(str_split[1]);
                            ncount = static_cast<int>(boost::lexical_cast<double>(str_split[0]));
                        }
                        catch (std::exception &e) {
                            error->exit("parse_general_vars", "Bad format for MAGMOM.");
                        }

                        for (i = icount; i < icount + ncount; ++i) {
                            magmom[i][2] = magmag;
                        }
                        icount += ncount;
                    }

                } else {
                    magmag = boost::lexical_cast<double>((*it));
                    magmom[icount++][2] = magmag;
                }
            }
            if (icount != nat) {
                error->exit("parse_general_vars",
                            "Number of entries for MAGMOM must be NAT.");
            }
        }
    }


    if (mode == "suggest") {
        if (general_var_dict["DBASIS"].empty()) {
            str_disp_basis = "Cart";
        } else {
            str_disp_basis = general_var_dict["DBASIS"];
        }
        std::transform(str_disp_basis.begin(), str_disp_basis.end(),
                       str_disp_basis.begin(), toupper);
        if ((str_disp_basis[0] != 'C') && (str_disp_basis[0] != 'F')) {
            error->exit("parse_general_vars", "Invalid DBASIS");
        }

        if (general_var_dict["TRIMEVEN"].empty()) {
            trim_dispsign_for_evenfunc = true;
        } else {
            assign_val(trim_dispsign_for_evenfunc,
                       "TRIMEVEN", general_var_dict);
        }

    }

    files->job_title = prefix;
    alm->mode = mode;
    system->nat = nat;
    system->nkd = nkd;
    symmetry->nsym = nsym;
    symmetry->is_printsymmetry = is_printsymmetry;
    symmetry->tolerance = tolerance;
    this->str_magmom = general_var_dict["MAGMOM"];
    memory->allocate(system->kdname, nkd);
    memory->allocate(system->magmom, nat, 3);

    for (i = 0; i < nkd; ++i) {
        system->kdname[i] = kdname[i];
    }
    for (i = 0; i < 3; ++i) {
        interaction->is_periodic[i] = is_periodic[i];
    }
    for (i = 0; i < nat; ++i) {
        for (j = 0; j < 3; ++j) {
            system->magmom[i][j] = magmom[i][j];
        }
    }
    system->lspin = lspin;
    system->noncollinear = noncollinear;
    symmetry->trev_sym_mag = trevsym;
    writes->print_hessian = print_hessian;
    constraint->tolerance_constraint = tolerance_constraint;

    if (mode == "suggest") {
        displace->disp_basis = str_disp_basis;
        displace->trim_dispsign_for_evenfunc = trim_dispsign_for_evenfunc;
    }

    memory->deallocate(kdname);
    memory->deallocate(magmom);

    kdname_v.clear();
    periodic_v.clear();
    no_defaults.clear();
    general_var_dict.clear();
}

void Input::parse_cell_parameter()
{
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

            boost::trim_if(line_wo_comment, boost::is_any_of("\t\n\r "));

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

            boost::trim_if(line_wo_comment, boost::is_any_of("\t\n\r "));

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
        boost::split(line_split, line, boost::is_any_of("\t "), boost::token_compress_on);

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
            system->lavec[i][j] = a * lavec_tmp[i][j];
        }
    }
    line_vec.clear();
    line_split.clear();
}

void Input::parse_interaction_vars()
{
    int i;
    int maxorder;
    int *nbody_include;

    std::vector<std::string> nbody_v;
    std::string str_allowed_list = "NORDER NBODY";
    std::string str_no_defaults = "NORDER";
    std::vector<std::string> no_defaults;
    std::map<std::string, std::string> interaction_var_dict;


    if (from_stdin) {
        std::cin.ignore();
    } else {
        ifs_input.ignore();
    }

    get_var_dict(str_allowed_list, interaction_var_dict);

    boost::split(no_defaults, str_no_defaults, boost::is_space());

    for (std::vector<std::string>::iterator it = no_defaults.begin();
         it != no_defaults.end(); ++it) {
        if (interaction_var_dict.find(*it) == interaction_var_dict.end()) {
            error->exit("parse_interaction_vars",
                        "The following variable is not found in &interaction input region: ",
                        (*it).c_str());
        }
    }

    assign_val(maxorder, "NORDER", interaction_var_dict);
    if (maxorder < 1)
        error->exit("parse_interaction_vars",
                    "maxorder has to be a positive integer");

    memory->allocate(nbody_include, maxorder);

    boost::split(nbody_v, interaction_var_dict["NBODY"], boost::is_space());

    if (nbody_v[0].empty()) {
        for (i = 0; i < maxorder; ++i) {
            nbody_include[i] = i + 2;
        }
    } else if (nbody_v.size() == maxorder) {
        for (i = 0; i < maxorder; ++i) {
            try {
                nbody_include[i] = boost::lexical_cast<int>(nbody_v[i]);
            }
            catch (std::exception &e) {
                std::cout << e.what() << std::endl;
                error->exit("parse_interaction_vars",
                            "NBODY must be an integer.");
            }
        }
    } else {
        error->exit("parse_interaction_vars",
                    "The number of entry of NBODY has to be equal to NORDER");
    }

    if (nbody_include[0] != 2) {
        error->warn("parce_input",
                    "Harmonic interaction is always 2 body (except on-site 1 body)");
    }

    interaction->maxorder = maxorder;
    memory->allocate(interaction->nbody_include, maxorder);

    for (i = 0; i < maxorder; ++i) {
        interaction->nbody_include[i] = nbody_include[i];
    }

    memory->deallocate(nbody_include);
    nbody_v.clear();
    no_defaults.clear();
}

void Input::parse_cutoff_radii()
{
    std::string line, line_wo_comment;
    std::string::size_type pos_first_comment_tag;
    std::vector<std::string> str_cutoff;

    if (from_stdin) {
        std::cin.ignore();
    } else {
        ifs_input.ignore();
    }

    str_cutoff.clear();

    if (from_stdin) {

        while (std::getline(std::cin, line)) {

            pos_first_comment_tag = line.find_first_of('#');

            if (pos_first_comment_tag == std::string::npos) {
                line_wo_comment = line;
            } else {
                line_wo_comment = line.substr(0, pos_first_comment_tag);
            }

            boost::trim_left(line_wo_comment);
            if (line_wo_comment.empty()) continue;
            if (is_endof_entry(line_wo_comment)) break;

            str_cutoff.push_back(line_wo_comment);
        }
    } else {

        while (std::getline(ifs_input, line)) {

            pos_first_comment_tag = line.find_first_of('#');

            if (pos_first_comment_tag == std::string::npos) {
                line_wo_comment = line;
            } else {
                line_wo_comment = line.substr(0, pos_first_comment_tag);
            }

            boost::trim_left(line_wo_comment);
            if (line_wo_comment.empty()) continue;
            if (is_endof_entry(line_wo_comment)) break;

            str_cutoff.push_back(line_wo_comment);
        }

    }

    int i, j, k;
    int nkd = system->nkd;
    int maxorder = interaction->maxorder;
    int ikd, jkd;
    int order;
    std::vector<std::string> cutoff_line;
    std::set<std::string> element_allowed;
    std::vector<std::string> str_pair;
    std::map<std::string, int> kd_map;

    double cutoff_tmp;
    double ***rcs;
    bool ***undefined_cutoff;

    memory->allocate(undefined_cutoff, maxorder, nkd, nkd);

    for (order = 0; order < maxorder; ++order) {
        for (i = 0; i < nkd; ++i) {
            for (j = 0; j < nkd; ++j) {
                undefined_cutoff[order][i][j] = true;
            }
        }
    }

    memory->allocate(rcs, maxorder, nkd, nkd);

    element_allowed.clear();

    for (i = 0; i < nkd; ++i) {
        element_allowed.insert(system->kdname[i]);
        kd_map.insert(std::map<std::string, int>::value_type(system->kdname[i], i));
    }

    element_allowed.insert("*");
    kd_map.insert(std::map<std::string, int>::value_type("*", -1));

    for (std::vector<std::string>::const_iterator it = str_cutoff.begin();
         it != str_cutoff.end(); ++it) {

        split_str_by_space(*it, cutoff_line);

        if (cutoff_line.size() < maxorder + 1) {
            error->exit("parse_cutoff_radii",
                        "Invalid format for &cutoff entry");
        }

        boost::split(str_pair, cutoff_line[0], boost::is_any_of("-"));

        if (str_pair.size() != 2) {
            error->exit("parse_cutoff_radii2",
                        "Invalid format for &cutoff entry");
        }

        for (i = 0; i < 2; ++i) {
            if (element_allowed.find(str_pair[i]) == element_allowed.end()) {
                error->exit("parse_cutoff_radii2",
                            "Invalid format for &cutoff entry");
            }
        }

        ikd = kd_map[str_pair[0]];
        jkd = kd_map[str_pair[1]];

        for (order = 0; order < maxorder; ++order) {
            // Accept any strings starting with 'N' or 'n' as 'None'
            if ((cutoff_line[order + 1][0] == 'N') || (cutoff_line[order + 1][0] == 'n')) {
                // Minus value for cutoff radius.
                // This is a flag for neglecting cutoff radius
                cutoff_tmp = -1.0;
            } else {
                cutoff_tmp = boost::lexical_cast<double>(cutoff_line[order + 1]);
            }

            if (ikd == -1 && jkd == -1) {
                for (i = 0; i < nkd; ++i) {
                    for (j = 0; j < nkd; ++j) {
                        rcs[order][i][j] = cutoff_tmp;
                        undefined_cutoff[order][i][j] = false;
                    }
                }
            } else if (ikd == -1) {
                for (i = 0; i < nkd; ++i) {
                    rcs[order][i][jkd] = cutoff_tmp;
                    rcs[order][jkd][i] = cutoff_tmp;
                    undefined_cutoff[order][i][jkd] = false;
                    undefined_cutoff[order][jkd][i] = false;
                }
            } else if (jkd == -1) {
                for (j = 0; j < nkd; ++j) {
                    rcs[order][j][ikd] = cutoff_tmp;
                    rcs[order][ikd][j] = cutoff_tmp;
                    undefined_cutoff[order][j][ikd] = false;
                    undefined_cutoff[order][ikd][j] = false;
                }
            } else {
                rcs[order][ikd][jkd] = cutoff_tmp;
                rcs[order][jkd][ikd] = cutoff_tmp;
                undefined_cutoff[order][ikd][jkd] = false;
                undefined_cutoff[order][jkd][ikd] = false;
            }
        }
    }
    element_allowed.clear();
    str_cutoff.clear();

    for (order = 0; order < maxorder; ++order) {
        for (j = 0; j < nkd; ++j) {
            for (k = 0; k < nkd; ++k) {
                if (undefined_cutoff[order][j][k]) {
                    std::cout << " Cutoff radius for " << std::setw(3)
                        << order + 2 << "th-order terms" << std::endl;
                    std::cout << " are not defined between elements " << std::setw(3) << j + 1
                        << " and " << std::setw(3) << k + 1 << std::endl;
                    error->exit("parse_cutoff_radii", "Incomplete cutoff radii");
                }
            }
        }
    }
    memory->deallocate(undefined_cutoff);

    memory->allocate(interaction->rcs, maxorder, nkd, nkd);

    for (i = 0; i < maxorder; ++i) {
        for (j = 0; j < nkd; ++j) {
            for (k = 0; k < nkd; ++k) {
                interaction->rcs[i][j][k] = rcs[i][j][k];
            }
        }
    }

    memory->deallocate(rcs);
}

void Input::parse_fitting_vars()
{
    int ndata, nstart, nend, nskip, nboot;
    std::string dfile, ffile;
    int multiply_data, constraint_flag;
    std::string rotation_axis;
    std::string fc2_file, fc3_file;

    std::string str_allowed_list = "NDATA NSTART NEND NSKIP NBOOT DFILE FFILE MULTDAT ICONST ROTAXIS FC2XML FC3XML";
    std::string str_no_defaults = "NDATA DFILE FFILE";
    std::vector<std::string> no_defaults;

    std::map<std::string, std::string> fitting_var_dict;

    bool fix_harmonic, fix_cubic;

    if (from_stdin) {
        std::cin.ignore();
    } else {
        ifs_input.ignore();
    }

    get_var_dict(str_allowed_list, fitting_var_dict);

    boost::split(no_defaults, str_no_defaults, boost::is_space());

    for (std::vector<std::string>::iterator it = no_defaults.begin();
         it != no_defaults.end(); ++it) {
        if (fitting_var_dict.find(*it) == fitting_var_dict.end()) {
            error->exit("parse_fitting_vars",
                        "The following variable is not found in &fitting input region: ",
                        (*it).c_str());
        }
    }

    assign_val(ndata, "NDATA", fitting_var_dict);

    if (fitting_var_dict["NSTART"].empty()) {
        nstart = 1;
    } else {
        assign_val(nstart, "NSTART", fitting_var_dict);
    }
    if (fitting_var_dict["NEND"].empty()) {
        nend = ndata;
    } else {
        assign_val(nend, "NEND", fitting_var_dict);
    }
    if (fitting_var_dict["NSKIP"].empty()) {
        nskip = 0;
    } else {
        assign_val(nskip, "NSKIP", fitting_var_dict);
    }

    if (ndata <= 0 || nstart <= 0 || nend <= 0
        || nstart > ndata || nend > ndata || nstart > nend) {
        error->exit("parce_fitting_vars",
                    "ndata, nstart, nend are not consistent with each other");
    }

    if (nskip < -1)
        error->exit("parce_fitting_vars",
                    "nskip has to be larger than -2.");
    if (nskip == -1) {

        if (fitting_var_dict["NBOOT"].empty()) {
            error->exit("parse_fitting_vars",
                        "NBOOT has to be given when NSKIP=-1");
        } else {
            assign_val(nboot, "NBOOT", fitting_var_dict);
        }

        if (nboot <= 0)
            error->exit("parce_input",
                        "nboot has to be a positive integer");
    } else {
        nboot = 0;
    }

    dfile = fitting_var_dict["DFILE"];
    ffile = fitting_var_dict["FFILE"];


    if (fitting_var_dict["MULTDAT"].empty()) {
        multiply_data = 1;
    } else {
        assign_val(multiply_data, "MULTDAT", fitting_var_dict);
    }

    if (fitting_var_dict["ICONST"].empty()) {
        constraint_flag = 1;
    } else {
        assign_val(constraint_flag, "ICONST", fitting_var_dict);
    }

    fc2_file = fitting_var_dict["FC2XML"];
    if (fc2_file.empty()) {
        fix_harmonic = false;
    } else {
        fix_harmonic = true;
    }

    fc3_file = fitting_var_dict["FC3XML"];
    if (fc3_file.empty()) {
        fix_cubic = false;
    } else {
        fix_cubic = true;
    }

    if (constraint_flag % 10 >= 2) {
        rotation_axis = fitting_var_dict["ROTAXIS"];
        if (rotation_axis.empty()) {
            error->exit("parse_fitting_vars",
                        "ROTAXIS has to be given when ICONST=2 or 3");
        }
    }

    system->ndata = ndata;
    system->nstart = nstart;
    system->nend = nend;
    system->nskip = nskip;

    fitting->nboot = nboot;
    files->file_disp = dfile;
    files->file_force = ffile;
    symmetry->multiply_data = multiply_data;
    constraint->constraint_mode = constraint_flag;
    constraint->rotation_axis = rotation_axis;
    constraint->fc2_file = fc2_file;
    constraint->fix_harmonic = fix_harmonic;
    constraint->fc3_file = fc3_file;
    constraint->fix_cubic = fix_cubic;

    fitting_var_dict.clear();
}

void Input::parse_atomic_positions()
{
    int i, j;
    std::string line, line_wo_comment;
    std::string str_tmp;
    std::string::size_type pos_first_comment_tag;
    std::vector<std::string> str_v, pos_line;
    double **xeq;
    int *kd;

    int nat = system->nat;

    if (from_stdin) {
        std::cin.ignore();
    } else {
        ifs_input.ignore();
    }

    str_v.clear();

    if (from_stdin) {

        while (std::getline(std::cin, line)) {

            pos_first_comment_tag = line.find_first_of('#');

            if (pos_first_comment_tag == std::string::npos) {
                line_wo_comment = line;
            } else {
                line_wo_comment = line.substr(0, pos_first_comment_tag);
            }

            boost::trim_if(line_wo_comment, boost::is_any_of("\t\n\r "));
            //            boost::trim_left(line_wo_comment);
            if (line_wo_comment.empty()) continue;
            if (is_endof_entry(line_wo_comment)) break;

            str_v.push_back(line_wo_comment);
        }

    } else {

        while (std::getline(ifs_input, line)) {

            pos_first_comment_tag = line.find_first_of('#');

            if (pos_first_comment_tag == std::string::npos) {
                line_wo_comment = line;
            } else {
                line_wo_comment = line.substr(0, pos_first_comment_tag);
            }

            boost::trim_if(line_wo_comment, boost::is_any_of("\t\n\r "));
            //            boost::trim_left(line_wo_comment);
            if (line_wo_comment.empty()) continue;
            if (is_endof_entry(line_wo_comment)) break;

            str_v.push_back(line_wo_comment);
        }
    }


    if (str_v.size() != nat) {
        error->exit("parse_atomic_positions",
                    "The number of entries for atomic positions should be NAT");
    }

    memory->allocate(xeq, nat, 3);
    memory->allocate(kd, nat);


    for (i = 0; i < nat; ++i) {

        split_str_by_space(str_v[i], pos_line);

        if (pos_line.size() == 4) {
            try {
                kd[i] = boost::lexical_cast<int>(pos_line[0]);
            }
            catch (std::exception &e) {
                std::cout << e.what() << std::endl;
                error->exit("parse_atomic_positions",
                            "Invalid entry for the &position field at line ",
                            i + 1);
            }

            for (j = 0; j < 3; ++j) {
                xeq[i][j] = boost::lexical_cast<double>(pos_line[j + 1]);
            }

        } else {
            error->exit("parse_atomic_positions",
                        "Bad format for &position region");
        }
    }


    memory->allocate(system->xcoord, nat, 3);
    memory->allocate(system->kd, nat);

    for (i = 0; i < nat; ++i) {

        system->kd[i] = kd[i];

        for (j = 0; j < 3; ++j) {
            system->xcoord[i][j] = xeq[i][j];
        }
    }

    memory->deallocate(xeq);
    memory->deallocate(kd);
    pos_line.clear();
    str_v.clear();
}

void Input::get_var_dict(const std::string keywords,
                         std::map<std::string, std::string> &var_dict)
{
    std::string line, key, val;
    std::string line_wo_comment;
    std::string::size_type pos_first_comment_tag;
    std::vector<std::string> str_entry, str_varval;

    std::set<std::string> keyword_set;

    boost::split(keyword_set, keywords, boost::is_space());

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

            boost::trim_left(line_wo_comment);
            if (line_wo_comment.empty()) continue;
            if (is_endof_entry(line_wo_comment)) break;

            //	std::cout << line_wo_comment << std::endl;

            // Split the input line by ';'

            boost::split(str_entry, line_wo_comment, boost::is_any_of(";"));

            for (std::vector<std::string>::iterator it = str_entry.begin();
                 it != str_entry.end(); ++it) {

                // Split the input entry by '='

                std::string str_tmp = boost::trim_copy(*it);

                if (!str_tmp.empty()) {

                    boost::split(str_varval, str_tmp, boost::is_any_of("="));

                    if (str_varval.size() != 2) {
                        error->exit("get_var_dict", "Unacceptable format");
                    }

                    key = boost::to_upper_copy(boost::trim_copy(str_varval[0]));
                    val = boost::trim_copy(str_varval[1]);

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

            boost::trim_left(line_wo_comment);
            if (line_wo_comment.empty()) continue;
            if (is_endof_entry(line_wo_comment)) break;

            // Split the input line by ';'

            boost::split(str_entry, line_wo_comment, boost::is_any_of(";"));

            for (std::vector<std::string>::iterator it = str_entry.begin();
                 it != str_entry.end(); ++it) {

                // Split the input entry by '='

                std::string str_tmp = boost::trim_copy(*it);

                if (!str_tmp.empty()) {

                    boost::split(str_varval, str_tmp, boost::is_any_of("="));

                    if (str_varval.size() != 2) {
                        error->exit("get_var_dict", "Unacceptable format");
                    }

                    key = boost::to_upper_copy(boost::trim_copy(str_varval[0]));
                    val = boost::trim_copy(str_varval[1]);

                    if (keyword_set.find(key) == keyword_set.end()) {
                        std::cout << "Could not recognize the variable "
                            << key << std::endl;
                        error->exit("get_var_dict", "Invalid variable found");
                    }

                    if (var_dict.find(key) != var_dict.end()) {
                        std::cout << "Variable " << key
                            << " appears twice in the input file." << std::endl;
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


int Input::locate_tag(std::string key)
{
    int ret = 0;
    std::string line;

    if (from_stdin) {

        std::cin.clear();
        std::cin.seekg(0, std::ios_base::beg);

        while (std::cin >> line) {
            boost::to_lower(line);
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
            boost::to_lower(line);
            if (line == key) {
                ret = 1;
                break;
            }
        }
        return ret;
    }
}

bool Input::is_endof_entry(std::string str)
{
    if (str[0] == '/') {
        return true;
    } else {
        return false;
    }
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
