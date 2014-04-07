/*
 input.cpp

 Copyright (c) 2014 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#include "input.h"
#include <iostream>
#include <string>
#include "memory.h"
#include "files.h"
#include "interaction.h"
#include "system.h"
#include "symmetry.h"
#include "error.h"
#include "fitting.h"
#include "constraint.h"
#include "fcs.h"
#include "ewald.h"
#include "patterndisp.h"
#include <algorithm>
#include <map>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

using namespace ALM_NS;

Input::Input(ALM *alm, int narg, char **arg): Pointers(alm) {}

Input::~Input() {}

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

    if (!locate_tag("&general")) {
        error->exit("parse_input", "&general entry not found in the input file");
    }
    parse_general_vars();

    if (!locate_tag("&cell")) {
        error->exit("parse_input", "&cell entry not found in the input file");
    }
    parse_cell_parameter();

    if (!locate_tag("&interaction")) {
        error->exit("parse_input", "&interaction entry not found in the input file");
    }
    parse_interaction_vars();

    if (!locate_tag("&cutoff")) {
        error->exit("parse_input", "&cutoff entry not found in the input file");
    }
    parse_cutoff_radii();


    if (alm->mode == "fitting") {
        if (!locate_tag("&fitting")) {
            error->exit("parse_input", "&fitting entry not found in the input file");
        }
        parse_fitting_vars();
    }

    if (!locate_tag("&position")) {
        error->exit("parse_input", "&position entry not found in the input file");
    }
    parse_atomic_positions();
}

void Input::parse_general_vars(){

    int i;
    std::string prefix, mode, str_tmp, str_disp_basis;
    int nat, nkd, nsym;
    int is_printsymmetry;
    bool is_periodic[3];
    std::string *kdname;
    double tolerance;

    std::vector<std::string> kdname_v, periodic_v;
    std::string str_allowed_list = "PREFIX MODE NAT NKD NSYM KD PERIODIC PRINTSYM TOLERANCE DBASIS";
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

    for (std::vector<std::string>::iterator it = no_defaults.begin(); it != no_defaults.end(); ++it){
        if (general_var_dict.find(*it) == general_var_dict.end()) {
            error->exit("parse_general_vars", "The following variable is not found in &general input region: ", (*it).c_str());
        }
    }

    prefix = general_var_dict["PREFIX"];
    mode = general_var_dict["MODE"];

    std::transform(mode.begin(), mode.end(), mode.begin(), tolower);
    if (mode != "fitting" && mode != "suggest") {
        error->exit("parse_general_vars", "Invalid MODE variable");
    }

    nat = boost::lexical_cast<int>(general_var_dict["NAT"]);
    nkd = boost::lexical_cast<int>(general_var_dict["NKD"]);


    if (general_var_dict["NSYM"].empty()) {
        nsym = 0;
    } else {
        nsym= boost::lexical_cast<int>(general_var_dict["NSYM"]);
    }

    if (general_var_dict["PRINTSYM"].empty()) {
        is_printsymmetry = 0;
    } else {
        is_printsymmetry = boost::lexical_cast<int>(general_var_dict["PRINTSYM"]);
    }

    split_str_by_space(general_var_dict["KD"], kdname_v);

    if (kdname_v.size() != nkd) {
        error->exit("parse_general_vars", "The number of entries for KD is inconsistent with NKD");
    } else {
        memory->allocate(kdname, nkd);
        for (i = 0; i < nkd; ++i){
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
            is_periodic[i] = boost::lexical_cast<int>(periodic_v[i]);
        }
    } else {
        error->exit("parse_general_vars", "Invalid number of entries for PERIODIC");
    }

    if (general_var_dict["TOLERANCE"].empty()) {
        tolerance = 1.0e-8;
    } else {
        tolerance = boost::lexical_cast<double>(general_var_dict["TOLERANCE"]);
    }

    if (mode == "suggest") {
        if (general_var_dict["DBASIS"].empty()) {
            str_disp_basis = "Cart";
        } else {
            str_disp_basis = general_var_dict["DBASIS"];
        }
        std::transform(str_disp_basis.begin(), str_disp_basis.end(), str_disp_basis.begin(), toupper);
        if (str_disp_basis[0] != 'C' && str_disp_basis[0] != 'F') {
            error->exit("parse_general_vars", "Invalid DBASIS");
        }
    }

    files->job_title = prefix;
    alm->mode = mode;
    system->nat = nat;
    system->nkd = nkd;
    symmetry->nsym = nsym;
    symmetry->is_printsymmetry = is_printsymmetry;
    symmetry->tolerance = tolerance;

    memory->allocate(system->kdname, nkd);

    for (i = 0; i < nkd; ++i){
        system->kdname[i] = kdname[i];
    }
    for (i = 0; i < 3; ++i) {
        interaction->is_periodic[i] = is_periodic[i];
    }
   

    if (mode == "suggest") displace->disp_basis = str_disp_basis;

    memory->deallocate(kdname);

    kdname_v.clear();
    periodic_v.clear();
    no_defaults.clear();
    general_var_dict.clear();
}

void Input::parse_cell_parameter() {

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
            system->lavec[i][j] = a * lavec_tmp[i][j];
        }
    }
}

void Input::parse_interaction_vars() {

    int i, interaction_type;
    int maxorder;
    bool is_longrange;
    int *nbody_include;
    std::string file_longrange;

    std::vector<std::string> nbody_v;
    std::string str_allowed_list = "NORDER NBODY INTERTYPE ILONG FLONG";
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

    for (std::vector<std::string>::iterator it = no_defaults.begin(); it != no_defaults.end(); ++it){
        if (interaction_var_dict.find(*it) == interaction_var_dict.end()) {
            error->exit("parse_interaction_vars", "The following variable is not found in &interaction input region: ", (*it).c_str());
        }
    }

    maxorder = boost::lexical_cast<int>(interaction_var_dict["NORDER"]);

    if (maxorder < 1) error->exit("parse_interaction_vars", "maxorder has to be a positive integer");

    if (interaction_var_dict["INTERTYPE"].empty()) {
        interaction_type = 0;
    } else {
        interaction_type = boost::lexical_cast<int>(interaction_var_dict["INTERTYPE"]);
    }
    if (interaction_type < 0 || interaction_type > 3) error->exit("parse_general_vars", "INTERTYPE should be 0, 1, 2 or 3.");


    memory->allocate(nbody_include, maxorder);

    boost::split(nbody_v, interaction_var_dict["NBODY"], boost::is_space());

    if (nbody_v[0].empty()) {
        for (i = 0; i < maxorder; ++i) {
            nbody_include[i] = i + 2;
        }
    } else if (nbody_v.size() == maxorder) {
        for (i = 0; i < maxorder; ++i) {
            nbody_include[i] = boost::lexical_cast<int>(nbody_v[i]);
        }
    } else {
        error->exit("parse_interaction_vars", "The number of entry of NBODY has to be equal to NORDER");
    }

    if (nbody_include[0] != 2) {
        error->warn("parce_input", "Harmonic interaction is always 2 body (except on-site 1 body)");
    }


    if (interaction_var_dict["ILONG"].empty()) {
        is_longrange = false;
    } else {
        is_longrange = boost::lexical_cast<int>(interaction_var_dict["ILONG"]);
    }

    if (is_longrange) {
        if (interaction_var_dict["FLONG"].empty()) {
            error->exit("parse_interaction_vars", "FLONG is necessary when ILONG = 1");
        } else {
            file_longrange = interaction_var_dict["FLONG"];
        }
    } else {
        file_longrange = "";
    }

    interaction->maxorder = maxorder;
    interaction->interaction_type = interaction_type;
    memory->allocate(interaction->nbody_include, maxorder);

    for (i = 0; i < maxorder; ++i){
        interaction->nbody_include[i] = nbody_include[i];    
    }

    memory->deallocate(nbody_include);
    nbody_v.clear();
    no_defaults.clear();

}

void Input::parse_cutoff_radii() {

    int i, j, k;
    double ***rcs;
    int nkd = system->nkd;
    int maxorder = interaction->maxorder;

    if (from_stdin) {
        std::cin.ignore();
    } else {
        ifs_input.ignore();
    }

    memory->allocate(rcs, maxorder, nkd, nkd);

    for (i = 0; i < maxorder; ++i) {
        for (j = 0; j < nkd; ++j) { 
            for (k = 0; k < nkd; ++k) {

                if (from_stdin) {
                    std::cin >> rcs[i][j][k];
                } else {
                    ifs_input >> rcs[i][j][k];
                }

            }
        }
        for (j = 0; j < nkd; ++j) {
            for (k = j + 1; k < nkd; ++k){
                if (rcs[i][j][k] != rcs[i][k][j]) error->exit("input", "Inconsistent cutoff radius rcs for order =", i + 2);
            }
        }
    }

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

void Input::parse_fitting_vars() {

    int ndata, nstart, nend, nskip, nboot;
    std::string dfile, ffile;
    int multiply_data, constraint_flag;
    std::string rotation_axis;
    std::string refsys_file, fc2_file;

    std::string str_allowed_list = "NDATA NSTART NEND NSKIP NBOOT DFILE FFILE MULTDAT ICONST ROTAXIS REFINFO FC2INFO";
    std::string str_no_defaults = "NDATA DFILE FFILE";
    std::vector<std::string> no_defaults;

    std::map<std::string, std::string> fitting_var_dict;


    if (from_stdin) {
        std::cin.ignore();
    } else {
        ifs_input.ignore();
    }

    get_var_dict(str_allowed_list, fitting_var_dict);

    boost::split(no_defaults, str_no_defaults, boost::is_space());

    for (std::vector<std::string>::iterator it = no_defaults.begin(); it != no_defaults.end(); ++it){
        if (fitting_var_dict.find(*it) == fitting_var_dict.end()) {
            error->exit("parse_fitting_vars", "The following variable is not found in &fitting input region: ", (*it).c_str());
        }
    }

    ndata = boost::lexical_cast<int>(fitting_var_dict["NDATA"]);

    if (fitting_var_dict["NSTART"].empty()) {
        nstart = 1;
    } else {
        nstart= boost::lexical_cast<int>(fitting_var_dict["NSTART"]);
    }
    if (fitting_var_dict["NEND"].empty()) {
        nend = ndata;
    } else {
        nend = boost::lexical_cast<int>(fitting_var_dict["NEND"]);
    }
    if (fitting_var_dict["NSKIP"].empty()) {
        nskip = 0;
    } else {
        nskip = boost::lexical_cast<int>(fitting_var_dict["NSKIP"]);
    }

    if(ndata <= 0 || nstart <= 0 || nend <= 0 
        || nstart > ndata || nend > ndata || nstart > nend) {
            error->exit("parce_fitting_vars", "ndata, nstart, nend are not consistent with each other");
    }


    if(nskip < -1) error->exit("parce_fitting_vars", "nskip has to be larger than -2.");
    if(nskip == -1) {

        if (fitting_var_dict["NBOOT"].empty()) {
            error->exit("parse_fitting_vars", "NBOOT has to be given when NSKIP=-1");
        } else {
            nboot = boost::lexical_cast<int>(fitting_var_dict["NBOOT"]);
        }

        if(nboot <= 0) error->exit("parce_input", "nboot has to be a positive integer");
    } else {
        nboot = 0;
    }

    dfile = fitting_var_dict["DFILE"];
    ffile = fitting_var_dict["FFILE"];


    if (fitting_var_dict["MULTDAT"].empty()) {
        multiply_data = 1;
    } else {
        multiply_data = boost::lexical_cast<int>(fitting_var_dict["MULTDAT"]);
    }

    if (multiply_data == 3) {
        refsys_file = fitting_var_dict["REFINFO"];
        if (refsys_file.empty()) {
            error->exit("parse_fitting_vars", "REFINFO tag has to be given when MULTDAT=3");
        } 
    }

    if (fitting_var_dict["ICONST"].empty()) {
        constraint_flag = 1;
    } else {
        constraint_flag = boost::lexical_cast<int>(fitting_var_dict["ICONST"]);
    }

    if(constraint_flag == 2 || constraint_flag == 4 || constraint_flag == 6) {
        fc2_file = fitting_var_dict["FC2INFO"];
        if (fc2_file.empty()) {
            error->exit("parse_fitting_vars", "FC2INFO has to be given when ICONST=2, 4, 6");
        }
    }

    if(constraint_flag >= 3) {
        rotation_axis = fitting_var_dict["ROTAXIS"];
        if (rotation_axis.empty()) {
            error->exit("parse_fitting_vars", "ROTAXIS has to be given when ICONST>=3");
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
    symmetry->refsys_file = refsys_file;

    fitting_var_dict.clear();
}

void Input::parse_atomic_positions() {

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

        while(std::getline(std::cin, line)) {

            pos_first_comment_tag = line.find_first_of('#');

            if (pos_first_comment_tag == std::string::npos) {
                line_wo_comment = line;
            } else {
                line_wo_comment = line.substr(0, pos_first_comment_tag);
            }


            boost::trim_left(line_wo_comment);
            if (line_wo_comment.empty()) continue;
            if (is_endof_entry(line_wo_comment)) break;

            str_v.push_back(line_wo_comment);
        }

    } else {

        while(std::getline(ifs_input, line)) {

            pos_first_comment_tag = line.find_first_of('#');

            if (pos_first_comment_tag == std::string::npos) {
                line_wo_comment = line;
            } else {
                line_wo_comment = line.substr(0, pos_first_comment_tag);
            }


            boost::trim_left(line_wo_comment);
            if (line_wo_comment.empty()) continue;
            if (is_endof_entry(line_wo_comment)) break;

            str_v.push_back(line_wo_comment);
        }
    }
    

    if (str_v.size() != nat) {
        error->exit("parse_atomic_positions", "The number of entries for atomic positions should be NAT");
    }

    memory->allocate(xeq, nat, 3);
    memory->allocate(kd, nat);


    for (i = 0; i < nat; ++i) {

        split_str_by_space(str_v[i], pos_line);

        if (pos_line.size() == 4) {
            kd[i] = boost::lexical_cast<int>(pos_line[0]);

            for (j = 0; j < 3; ++j) {
                xeq[i][j] = boost::lexical_cast<double>(pos_line[j + 1]);
            }

        } else {
            error->exit("parse_atomic_positions", "Bad format for &position region");
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

void Input::get_var_dict(const std::string keywords, std::map<std::string, std::string> &var_dict) {

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

            for (std::vector<std::string>::iterator it = str_entry.begin(); it != str_entry.end(); ++it) {

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

            //	std::cout << line_wo_comment << std::endl;

            // Split the input line by ';'

            boost::split(str_entry, line_wo_comment, boost::is_any_of(";"));

            for (std::vector<std::string>::iterator it = str_entry.begin(); it != str_entry.end(); ++it) {

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
    }

    keyword_set.clear();
}


int Input::locate_tag(std::string key){

    int ret = 0;
    std::string line;

    if (from_stdin) {
     
        std::cin.clear();
        std::cin.seekg(0, std::ios_base::beg);

        while (std::cin >> line) {
            boost::to_lower(line);
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
            boost::to_lower(line);
            if (line == key) {
                ret = 1;
                break;
            }
        }
        return ret;
    }
 
}

bool Input::is_endof_entry(std::string str) {

    if (str[0] == '/') {
        return true;
    } else {
        return false;
    }
}

void Input::split_str_by_space(const std::string str, std::vector<std::string> &str_vec) {

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
