/*
 input_parser.cpp

 Copyright (c) 2014, 2015, 2016 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory
 or http://opensource.org/licenses/mit-license.php for information.
*/

#include "input_parser.h"
#include "alm.h"
#include "error.h"
#include "files.h"
#include "optimize.h"
#include "input_setter.h"
#include "memory.h"
#include <iostream>
#include <iomanip>
#include <string>
#include <algorithm>
#include <map>
#include <set>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

using namespace ALM_NS;

InputParser::InputParser()
{
    input_setter = new InputSetter();
}

InputParser::~InputParser()
{
    delete input_setter;
}

void InputParser::run(ALM *alm,
                      const int narg,
                      const char *const *arg)
{
    if (narg == 1) {

        from_stdin = true;

    } else {

        from_stdin = false;

        ifs_input.open(arg[1], std::ios::in);
        if (!ifs_input) {
            std::cout << "No such file or directory: " << arg[1] << std::endl;
            std::exit(EXIT_FAILURE);
        }
    }

    parse_input(alm);
}

std::string InputParser::get_run_mode() const
{
    return mode;
}

void InputParser::parse_displacement_and_force_files(std::vector<std::vector<double>> &u,
                                                     std::vector<std::vector<double>> &f,
                                                     DispForceFile &datfile_in) const
{
    int nrequired;

    if (datfile_in.ndata == 0) {
        nrequired = -1;
    } else {
        // Total number of data entries (displacement + force)
        nrequired = 6 * nat * datfile_in.ndata;
    }

    std::vector<double> value_arr;

    std::string line;
    double val;
    auto nline_u = 0;

    // Open the target file and copy the data to 1D temporary vector
    std::ifstream ifs_data;
    ifs_data.open(datfile_in.filename.c_str(), std::ios::in);
    if (!ifs_data) exit("openfiles", "cannot open DFSET file");
    auto reach_end = false;
    while (std::getline(ifs_data >> std::ws, line)) {
        if (line[0] != '#') {

            std::istringstream iss(line);

            while (iss >> val) {
                value_arr.push_back(val);
                ++nline_u;

                if (nline_u == nrequired) {
                    reach_end = true;
                    break;
                }
            }
        }

        if (reach_end) break;
    }
    ifs_data.close();

    // Check if the length of the vector is correct.
    // Also, estimate ndata if it is not set. 
    const auto n_entries = value_arr.size();

    if (nrequired == -1) {
        if (n_entries % (6 * nat) == 0) {
            datfile_in.ndata = n_entries / (6 * nat);
        } else {
            exit("parse_displacement_and_force_files",
                 "The number of lines in DFSET is indivisible by NAT");
        }
    } else {
        if (n_entries < nrequired) {
            exit("parse_displacement_and_force_files",
                 "The number of lines in DFSET is too small for the given NDATA = ",
                 datfile_in.ndata);
        }
    }

    if (datfile_in.nstart == 0) datfile_in.nstart = 1;
    if (datfile_in.nend == 0) datfile_in.nend = datfile_in.ndata;

    // Copy the data into 2D array
    const auto ndata_used = datfile_in.nend - datfile_in.nstart
                            + 1 - datfile_in.skip_e + datfile_in.skip_s;

    u.resize(ndata_used, std::vector<double>(3 * nat));
    f.resize(ndata_used, std::vector<double>(3 * nat));

    auto idata = 0;
    for (size_t i = 0; i < datfile_in.ndata; ++i) {
        if (i < datfile_in.nstart - 1) continue;
        if (i >= datfile_in.skip_s - 1 && i < datfile_in.skip_e - 1) continue; // When skip_s == skip_e, skip nothing.
        if (i > datfile_in.nend - 1) break;

        for (auto j = 0; j < nat; ++j) {
            for (auto k = 0; k < 3; ++k) {
                u[idata][3 * j + k] = value_arr[6 * nat * i + 6 * j + k];
                f[idata][3 * j + k] = value_arr[6 * nat * i + 6 * j + k + 3];
            }
        }
        ++idata;
    }
    value_arr.clear();
}

void InputParser::parse_input(ALM *alm)
{
    // The order of calling methods in this method is important.
    // Since following methods rely on variables those already
    // parsed.
    // Those below are set as the private class variables. See input_parser.h.
    //  std::string *kdname;
    //  std::string mode;
    //  int maxorder;
    //  int nat;
    //  int nkd;

    if (!locate_tag("&general")) {
        exit("parse_input",
             "&general entry not found in the input file");
    }

    // kdname is allocated in this method.
    parse_general_vars(alm);

    if (!locate_tag("&cell")) {
        exit("parse_input",
             "&cell entry not found in the input file");
    }
    parse_cell_parameter();

    if (!locate_tag("&position")) {
        exit("parse_input",
             "&position entry not found in the input file");
    }
    parse_atomic_positions();
    input_setter->set_geometric_structure(alm);

    if (!locate_tag("&interaction")) {
        exit("parse_input",
             "&interaction entry not found in the input file");
    }
    parse_interaction_vars();

    if (!locate_tag("&cutoff")) {
        exit("parse_input",
             "&cutoff entry not found in the input file");
    }
    parse_cutoff_radii();
    input_setter->define(alm);

    if (mode == "optimize") {
        if (!locate_tag("&optimize")) {
            exit("parse_input",
                 "&optimize entry not found in the input file");
        }
        parse_optimize_vars(alm);
    }

    deallocate(kdname);
}

void InputParser::parse_general_vars(ALM *alm)
{
    size_t i;
    std::string str_tmp, str_disp_basis, basis_force_constant;
    int printsymmetry, is_periodic[3];
    size_t icount, ncount;
    auto trim_dispsign_for_evenfunc = true;
    bool print_hessian;
    int noncollinear, trevsym;
    double **magmom, magmag{0.0};
    double tolerance;
    double tolerance_constraint;
    int verbosity, nmaxsave;

    std::vector<std::string> kdname_v, periodic_v, magmom_v, str_split;
    const std::vector<std::string> input_list{
            "PREFIX", "MODE", "NAT", "NKD", "KD", "PERIODIC", "PRINTSYM", "TOLERANCE",
            "DBASIS", "TRIMEVEN", "VERBOSITY",
            "MAGMOM", "NONCOLLINEAR", "TREVSYM", "HESSIAN", "TOL_CONST", "FCSYM_BASIS",
            "NMAXSAVE"
    };
    std::vector<std::string> no_defaults{"PREFIX", "MODE", "NAT", "NKD", "KD"};
    std::map<std::string, std::string> general_var_dict;

    if (from_stdin) {
        std::cin.ignore();
    } else {
        ifs_input.ignore();
    }

    get_var_dict(input_list, general_var_dict);

    for (const auto &it : no_defaults) {
        if (general_var_dict.find(it) == general_var_dict.end()) {
            exit("parse_general_vars",
                 "The following variable is not found in &general input region: ",
                 it.c_str());
        }
    }

    const auto prefix = general_var_dict["PREFIX"];
    mode = general_var_dict["MODE"];

    std::transform(mode.begin(), mode.end(), mode.begin(), tolower);
    if (mode != "suggest" && mode != "optimize" && mode != "opt") {
        exit("parse_general_vars", "Invalid MODE variable");
    }
    if (mode == "opt") mode = "optimize";

    assign_val(nat, "NAT", general_var_dict);
    assign_val(nkd, "NKD", general_var_dict);

    if (general_var_dict["VERBOSITY"].empty()) {
        verbosity = 1;
    } else {
        assign_val(verbosity, "VERBOSITY", general_var_dict);
    }

    if (general_var_dict["PRINTSYM"].empty()) {
        printsymmetry = 0;
    } else {
        assign_val(printsymmetry, "PRINTSYM", general_var_dict);
    }

    split_str_by_space(general_var_dict["KD"], kdname_v);

    if (kdname_v.size() != nkd) {
        exit("parse_general_vars",
             "The number of entries for KD is inconsistent with NKD");
    } else {
        allocate(kdname, nkd);
        for (i = 0; i < nkd; ++i) {
            kdname[i] = kdname_v[i];
        }
    }

    split_str_by_space(general_var_dict["PERIODIC"], periodic_v);

    if (periodic_v.empty()) {
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
                exit("parse_general_vars",
                     "The PERIODIC tag must be a set of integers.");
            }
        }
    } else {
        exit("parse_general_vars",
             "Invalid number of entries for PERIODIC");
    }

    if (general_var_dict["TOLERANCE"].empty()) {
        tolerance = 1.0e-3;
    } else {
        assign_val(tolerance, "TOLERANCE", general_var_dict);
    }
    if (general_var_dict["TOL_CONST"].empty()) {
        tolerance_constraint = 1.0e-6;
    } else {
        assign_val(tolerance_constraint, "TOL_CONST", general_var_dict);
    }

    if (general_var_dict["FCSYM_BASIS"].empty()) {
        basis_force_constant = "Lattice";
    } else {
        basis_force_constant = general_var_dict["FCSYM_BASIS"];
        boost::to_lower(basis_force_constant);

        if (basis_force_constant[0] != 'c' && basis_force_constant[0] != 'l') {
            exit("parse_general_vars", "Invalid FC_BASIS.",
                 basis_force_constant.c_str());
        }
    }

    if (general_var_dict["NMAXSAVE"].empty()) {
        nmaxsave = 5;
    } else {
        assign_val(nmaxsave, "NMAXSAVE", general_var_dict);
    }

    // Convert MAGMOM input to array
    allocate(magmom, nat, 3);
    auto lspin = false;

    for (i = 0; i < nat; ++i) {
        for (size_t j = 0; j < 3; ++j) {
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
            for (const auto & it : magmom_v) {
                if (it.find('*') != std::string::npos) {
                    exit("parse_general_vars",
                         "Wild card '*' is not supported when NONCOLLINEAR = 1.");
                } else {
                    magmag = boost::lexical_cast<double>(it);
                    if (icount / 3 >= nat) {
                        exit("parse_general_vars", "Too many entries for MAGMOM.");
                    }
                    magmom[icount / 3][icount % 3] = magmag;
                    ++icount;
                }
            }

            if (icount != 3 * nat) {
                exit("parse_general_vars",
                     "Number of entries for MAGMOM must be 3*NAT when NONCOLLINEAR = 1.");
            }
        } else {
            icount = 0;
            split_str_by_space(general_var_dict["MAGMOM"], magmom_v);
            for (const auto & it : magmom_v) {

                if (it.find('*') != std::string::npos) {
                    if (it == "*") {
                        exit("parse_general_vars",
                             "Please place '*' without space for the MAGMOM-tag.");
                    }
                    boost::split(str_split, it, boost::is_any_of("*"));
                    if (str_split.size() != 2) {
                        exit("parse_general_vars",
                             "Invalid format for the MAGMOM-tag.");
                    } else {
                        if (str_split[0].empty() || str_split[1].empty()) {
                            exit("parse_general_vars",
                                 "Please place '*' without space for the MAGMOM-tag.");
                        }
                        ncount = 0;
                        try {
                            magmag = boost::lexical_cast<double>(str_split[1]);
                            ncount = static_cast<size_t>(boost::lexical_cast<double>(str_split[0]));
                        }
                        catch (std::exception) {
                            exit("parse_general_vars", "Bad format for MAGMOM.");
                        }

                        for (i = icount; i < icount + ncount; ++i) {
                            magmom[i][2] = magmag;
                        }
                        icount += ncount;
                    }

                } else {
                    magmag = boost::lexical_cast<double>(it);
                    if (icount == nat) {
                        icount = 0;
                        break;
                    }
                    magmom[icount++][2] = magmag;
                }
            }
            if (icount != nat) {
                exit("parse_general_vars",
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
            exit("parse_general_vars", "Invalid DBASIS");
        }

        if (!general_var_dict["TRIMEVEN"].empty()) {
            assign_val(trim_dispsign_for_evenfunc,
                       "TRIMEVEN", general_var_dict);
        }
    }

    input_setter->set_general_vars(alm,
                                   prefix,
                                   mode,
                                   verbosity,
                                   str_disp_basis,
                                   general_var_dict["MAGMOM"],
                                   nat,
                                   nkd,
                                   printsymmetry,
                                   is_periodic,
                                   trim_dispsign_for_evenfunc,
                                   lspin,
                                   print_hessian,
                                   noncollinear,
                                   trevsym,
                                   kdname,
                                   magmom,
                                   tolerance,
                                   tolerance_constraint,
                                   basis_force_constant,
                                   nmaxsave);

    allocate(magmom, nat, 3);

    kdname_v.clear();
    periodic_v.clear();
    no_defaults.clear();
    general_var_dict.clear();
}

void InputParser::parse_cell_parameter()
{
    auto a = 0.0;
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
        exit("parse_cell_parameter",
             "Too few or too much lines for the &cell field.\n"
             "The number of valid lines for the &cell field should be 4.");
    }

    for (auto i = 0; i < 4; ++i) {

        line = line_vec[i];
        boost::split(line_split, line, boost::is_any_of("\t "), boost::token_compress_on);

        if (i == 0) {
            // Lattice factor a
            if (line_split.size() == 1) {
                a = boost::lexical_cast<double>(line_split[0]);
            } else {
                exit("parse_cell_parameter",
                     "Unacceptable format for &cell field.");
            }

        } else {
            // Lattice vectors a1, a2, a3
            if (line_split.size() == 3) {
                for (auto j = 0; j < 3; ++j) {
                    lavec_tmp[j][i - 1] = boost::lexical_cast<double>(line_split[j]);
                }
            } else {
                exit("parse_cell_parameter",
                     "Unacceptable format for &cell field.");
            }
        }
    }

    input_setter->set_cell_parameter(a, lavec_tmp);
}


void InputParser::parse_interaction_vars()
{
    int i;
    int *nbody_include;

    std::vector<std::string> nbody_v;
    const std::vector<std::string> input_list{"NORDER", "NBODY"};
    std::vector<std::string> no_defaults{"NORDER"};
    std::map<std::string, std::string> interaction_var_dict;


    if (from_stdin) {
        std::cin.ignore();
    } else {
        ifs_input.ignore();
    }

    get_var_dict(input_list, interaction_var_dict);

    for (const auto &it : no_defaults) {
        if (interaction_var_dict.find(it) == interaction_var_dict.end()) {
            exit("parse_interaction_vars",
                 "The following variable is not found in &interaction input region: ",
                 it.c_str());
        }
    }

    assign_val(maxorder, "NORDER", interaction_var_dict);
    if (maxorder < 1)
        exit("parse_interaction_vars",
             "maxorder has to be a positive integer");

    allocate(nbody_include, maxorder);

    boost::split(nbody_v, interaction_var_dict["NBODY"], boost::is_space());

    if (nbody_v[0].empty()) {
        for (i = 0; i < maxorder; ++i) {
            nbody_include[i] = i + 2;
        }
    } else if (nbody_v.size() == static_cast<size_t>(maxorder)) {
        for (i = 0; i < maxorder; ++i) {
            try {
                nbody_include[i] = boost::lexical_cast<int>(nbody_v[i]);
            }
            catch (std::exception &e) {
                std::cout << e.what() << std::endl;
                exit("parse_interaction_vars",
                     "NBODY must be an integer.");
            }
        }
    } else {
        exit("parse_interaction_vars",
             "The number of entry of NBODY has to be equal to NORDER");
    }

    if (nbody_include[0] != 2) {
        warn("parse_interaction_vars",
             "Harmonic interaction is always 2 body (except on-site 1 body)");
    }

    input_setter->set_interaction_vars(maxorder, nbody_include);

    deallocate(nbody_include);

    nbody_v.clear();
    no_defaults.clear();
}


void InputParser::parse_optimize_vars(ALM *alm)
{
    // This method must not be called before setting NAT by parse_general_vars.

    int constraint_flag;
    auto flag_sparse = 0;
    std::string rotation_axis;

    OptimizerControl optcontrol;
    std::vector<std::vector<double>> u_tmp1, f_tmp1;
    std::vector<std::vector<double>> u_tmp2, f_tmp2;

    const std::vector<std::string> input_list{
            "LMODEL", "SPARSE", "SPARSESOLVER",
            "ICONST", "ROTAXIS", "FC2XML", "FC3XML",
            "NDATA", "NSTART", "NEND", "SKIP", "DFSET",
            "NDATA_CV", "NSTART_CV", "NEND_CV", "DFSET_CV",
            "L1_RATIO", "STANDARDIZE", "ENET_DNORM",
            "L1_ALPHA", "CV_MAXALPHA", "CV_MINALPHA", "CV_NALPHA",
            "CV", "MAXITER", "CONV_TOL", "NWRITE", "SOLUTION_PATH", "DEBIAS_OLS"
    };

    std::map<std::string, std::string> optimize_var_dict;

    if (from_stdin) {
        std::cin.ignore();
    } else {
        ifs_input.ignore();
    }

    get_var_dict(input_list, optimize_var_dict);

    if (!optimize_var_dict["LMODEL"].empty()) {
        auto str_LMODEL = optimize_var_dict["LMODEL"];
        boost::to_lower(str_LMODEL);

        if (str_LMODEL == "ols" || str_LMODEL == "ls"
            || str_LMODEL == "least-squares" || str_LMODEL == "1") {
            optcontrol.linear_model = 1;
        } else if (str_LMODEL == "enet" || str_LMODEL == "elastic-net"
                   || str_LMODEL == "2") {
            optcontrol.linear_model = 2;
        } else if (str_LMODEL == "adaptive-lasso" || str_LMODEL == "3") {
            optcontrol.linear_model = 3;
        } else {
            exit("parse_optimize_vars", "Invalid OPTIMIZER-tag");
        }
    }

    if (!optimize_var_dict["SPARSE"].empty()) {
        assign_val(flag_sparse, "SPARSE", optimize_var_dict);
        optcontrol.use_sparse_solver = flag_sparse;
    }
    if (!optimize_var_dict["SPARSESOLVER"].empty()) {
        std::string str_sparsesolver;
        assign_val(str_sparsesolver, "SPARSESOLVER", optimize_var_dict);
        const auto str_lower = boost::algorithm::to_lower_copy(str_sparsesolver);

        if (str_lower != "simplicialldlt"
            && str_lower != "sparseqr"
            && str_lower != "conjugategradient"
            && str_lower != "leastsquaresconjugategradient"
            && str_lower != "bicgstab") {
            exit("parse_optimize_vars", "Unsupported SPARSESOLVER :", str_sparsesolver.c_str());
        }
        optcontrol.sparsesolver = str_sparsesolver;
    }

    if (!optimize_var_dict["ENET_DNORM"].empty()) {
        optcontrol.displacement_normalization_factor
                = boost::lexical_cast<double>(optimize_var_dict["ENET_DNORM"]);
    }
    if (!optimize_var_dict["L1_ALPHA"].empty()) {
        optcontrol.l1_alpha = boost::lexical_cast<double>(optimize_var_dict["L1_ALPHA"]);
    }
    if (!optimize_var_dict["CV_MINALPHA"].empty()) {
        optcontrol.l1_alpha_min = boost::lexical_cast<double>(optimize_var_dict["CV_MINALPHA"]);
    }
    if (!optimize_var_dict["CV_MAXALPHA"].empty()) {
        optcontrol.l1_alpha_max = boost::lexical_cast<double>(optimize_var_dict["CV_MAXALPHA"]);
    }
    if (!optimize_var_dict["CV_NALPHA"].empty()) {
        optcontrol.num_l1_alpha = boost::lexical_cast<int>(optimize_var_dict["CV_NALPHA"]);
    }
    if (!optimize_var_dict["CONV_TOL"].empty()) {
        optcontrol.tolerance_iteration = boost::lexical_cast<double>(optimize_var_dict["CONV_TOL"]);
    }
    if (!optimize_var_dict["MAXITER"].empty()) {
        optcontrol.maxnum_iteration = boost::lexical_cast<int>(optimize_var_dict["MAXITER"]);
    }
    if (!optimize_var_dict["CV"].empty()) {
        optcontrol.cross_validation = boost::lexical_cast<int>(optimize_var_dict["CV"]);
        if (optcontrol.cross_validation == 1) {
            exit("parse_optimize_vars", "CV must not be 1.");
        }
    }
    if (!optimize_var_dict["NWRITE"].empty()) {
        optcontrol.output_frequency = boost::lexical_cast<int>(optimize_var_dict["NWRITE"]);
    }
    if (!optimize_var_dict["STANDARDIZE"].empty()) {
        optcontrol.standardize = boost::lexical_cast<int>(optimize_var_dict["STANDARDIZE"]);
    }
    if (!optimize_var_dict["SOLUTION_PATH"].empty()) {
        optcontrol.save_solution_path = boost::lexical_cast<int>(optimize_var_dict["SOLUTION_PATH"]);
    }
    if (!optimize_var_dict["DEBIAS_OLS"].empty()) {
        optcontrol.debiase_after_l1opt = boost::lexical_cast<int>(optimize_var_dict["DEBIAS_OLS"]);
    }
    if (!optimize_var_dict["L1_RATIO"].empty()) {
        optcontrol.l1_ratio = boost::lexical_cast<double>(optimize_var_dict["L1_RATIO"]);
    }


    DispForceFile datfile_train;

    if (optimize_var_dict["DFSET"].empty()) {
        exit("parse_optimize_vars", "DFSET tag must be given.");
    }
    datfile_train.filename = optimize_var_dict["DFSET"];

    if (!optimize_var_dict["NDATA"].empty()) {
        assign_val(datfile_train.ndata, "NDATA", optimize_var_dict);
    }
    if (!optimize_var_dict["NSTART"].empty()) {
        assign_val(datfile_train.nstart, "NSTART", optimize_var_dict);
    }
    if (!optimize_var_dict["NEND"].empty()) {
        assign_val(datfile_train.nend, "NEND", optimize_var_dict);
    }

    if (!optimize_var_dict["SKIP"].empty()) {
        std::vector<std::string> str_entry;
        std::string str_skip;
        assign_val(str_skip, "SKIP", optimize_var_dict);
        boost::split(str_entry, str_skip, boost::is_any_of("-"));
        if (str_entry.size() == 1) {
            datfile_train.skip_s = boost::lexical_cast<int>(str_entry[0]);
            datfile_train.skip_e = datfile_train.skip_s + 1;
        } else if (str_entry.size() == 2) {
            datfile_train.skip_s = boost::lexical_cast<int>(str_entry[0]);
            datfile_train.skip_e = boost::lexical_cast<int>(str_entry[1]) + 1;
        } else {
            exit("parse_optimize_vars", "Invalid format for the SKIP-tag.");
        }
    }

    if (!is_data_range_consistent(datfile_train)) {
        exit("parse_optimize_vars",
             "NDATA, NSTART, NEND and SKIP tags are inconsistent.");
    }

    // Parse u_tmp1 and f_tmp1 from DFSET and set ndata if it's not.
    parse_displacement_and_force_files(u_tmp1,
                                       f_tmp1,
                                       datfile_train);

    // Check consistency again
    if (!is_data_range_consistent(datfile_train)) {
        exit("parse_optimize_vars",
             "NDATA, NSTART, NEND and SKIP tags are inconsistent.");
    }

    auto datfile_validation = datfile_train;
    datfile_validation.skip_s = 0;
    datfile_validation.skip_e = 0;

    if (!optimize_var_dict["NDATA_CV"].empty()) {
        datfile_validation.ndata = boost::lexical_cast<int>(optimize_var_dict["NDATA_CV"]);
    }

    if (!optimize_var_dict["NSTART_CV"].empty()) {
        datfile_validation.nstart = boost::lexical_cast<int>(optimize_var_dict["NSTART_CV"]);
    }
    if (!optimize_var_dict["NEND_CV"].empty()) {
        datfile_validation.nend = boost::lexical_cast<int>(optimize_var_dict["NEND_CV"]);
    }

    if (!optimize_var_dict["DFSET_CV"].empty()) {
        datfile_validation.filename = optimize_var_dict["DFSET_CV"];
    }

    if (!is_data_range_consistent(datfile_validation)) {
        exit("parse_optimize_vars",
             "NDATA_CV, NSTART_CV and NEND_CV tags are inconsistent.");
    }

    if (optcontrol.cross_validation == -1) {
        parse_displacement_and_force_files(u_tmp2,
                                           f_tmp2,
                                           datfile_validation);

        if (!is_data_range_consistent(datfile_validation)) {
            exit("parse_optimize_vars",
                 "NDATA_CV, NSTART_CV and NEND_CV tags are inconsistent.");
        }
    }


    input_setter->set_optimize_vars(alm,
                                    u_tmp1, f_tmp1,
                                    u_tmp2, f_tmp2,
                                    optcontrol);

    input_setter->set_file_vars(alm,
                                datfile_train,
                                datfile_validation);


    if (optimize_var_dict["ICONST"].empty()) {
        constraint_flag = 11;
    } else {
        assign_val(constraint_flag, "ICONST", optimize_var_dict);
    }

    auto fc2_file = optimize_var_dict["FC2XML"];
    auto fc3_file = optimize_var_dict["FC3XML"];
    const auto fix_harmonic = !fc2_file.empty();
    const auto fix_cubic = !fc3_file.empty();

    if (constraint_flag % 10 >= 2) {
        rotation_axis = optimize_var_dict["ROTAXIS"];
        if (rotation_axis.empty()) {
            exit("parse_optimize_vars",
                 "ROTAXIS has to be given when ICONST=2 or 3");
        }
    }

    input_setter->set_constraint_vars(alm,
                                      constraint_flag,
                                      rotation_axis,
                                      fc2_file,
                                      fc3_file,
                                      fix_harmonic,
                                      fix_cubic);

    optimize_var_dict.clear();
}

void InputParser::parse_atomic_positions()
{
    std::string line, line_wo_comment;
    std::string str_tmp;
    std::string::size_type pos_first_comment_tag;
    std::vector<std::string> str_v, pos_line;
    double (*xeq)[3];
    int *kd;

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
            if (line_wo_comment.empty()) continue;
            if (is_endof_entry(line_wo_comment)) break;

            str_v.push_back(line_wo_comment);
        }
    }


    if (str_v.size() != nat) {
        exit("parse_atomic_positions",
             "The number of entries for atomic positions should be NAT");
    }

    allocate(xeq, nat);
    allocate(kd, nat);


    for (size_t i = 0; i < nat; ++i) {

        split_str_by_space(str_v[i], pos_line);

        if (pos_line.size() == 4) {
            try {
                kd[i] = boost::lexical_cast<int>(pos_line[0]);
            }
            catch (std::exception &e) {
                std::cout << e.what() << std::endl;
                exit("parse_atomic_positions",
                     "Invalid entry for the &position field at line ",
                     i + 1);
            }

            for (auto j = 0; j < 3; ++j) {
                xeq[i][j] = boost::lexical_cast<double>(pos_line[j + 1]);
            }

        } else {
            exit("parse_atomic_positions",
                 "Bad format for &position region");
        }
    }


    input_setter->set_atomic_positions(nat, kd, xeq);

    deallocate(xeq);
    deallocate(kd);
    pos_line.clear();
    str_v.clear();
}

void InputParser::parse_cutoff_radii()
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

    size_t i, j, k;
    int order;
    std::vector<std::string> cutoff_line;
    std::set<std::string> element_allowed;
    std::vector<std::string> str_pair;
    std::map<std::string, int> kd_map;

    double cutoff_tmp;
    double ***cutoff_radii_tmp;
    bool ***undefined_cutoff;

    allocate(undefined_cutoff, maxorder, nkd, nkd);

    for (order = 0; order < maxorder; ++order) {
        for (i = 0; i < nkd; ++i) {
            for (j = 0; j < nkd; ++j) {
                undefined_cutoff[order][i][j] = true;
            }
        }
    }

    allocate(cutoff_radii_tmp, maxorder, nkd, nkd);

    element_allowed.clear();

    for (i = 0; i < nkd; ++i) {
        element_allowed.insert(kdname[i]);
        kd_map.insert(std::map<std::string, int>::value_type(kdname[i], i));
    }

    element_allowed.insert("*");
    kd_map.insert(std::map<std::string, int>::value_type("*", -1));

    for (const auto &it : str_cutoff) {

        split_str_by_space(it, cutoff_line);

        if (cutoff_line.size() < maxorder + 1) {
            exit("parse_cutoff_radii",
                 "Invalid format for &cutoff entry");
        }

        boost::split(str_pair, cutoff_line[0], boost::is_any_of("-"));

        if (str_pair.size() != 2) {
            exit("parse_cutoff_radii2",
                 "Invalid format for &cutoff entry");
        }

        for (i = 0; i < 2; ++i) {
            if (element_allowed.find(str_pair[i]) == element_allowed.end()) {
                exit("parse_cutoff_radii2",
                     "Invalid format for &cutoff entry");
            }
        }

        const auto ikd = kd_map[str_pair[0]];
        const auto jkd = kd_map[str_pair[1]];

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
                        cutoff_radii_tmp[order][i][j] = cutoff_tmp;
                        undefined_cutoff[order][i][j] = false;
                    }
                }
            } else if (ikd == -1) {
                for (i = 0; i < nkd; ++i) {
                    cutoff_radii_tmp[order][i][jkd] = cutoff_tmp;
                    cutoff_radii_tmp[order][jkd][i] = cutoff_tmp;
                    undefined_cutoff[order][i][jkd] = false;
                    undefined_cutoff[order][jkd][i] = false;
                }
            } else if (jkd == -1) {
                for (j = 0; j < nkd; ++j) {
                    cutoff_radii_tmp[order][j][ikd] = cutoff_tmp;
                    cutoff_radii_tmp[order][ikd][j] = cutoff_tmp;
                    undefined_cutoff[order][j][ikd] = false;
                    undefined_cutoff[order][ikd][j] = false;
                }
            } else {
                cutoff_radii_tmp[order][ikd][jkd] = cutoff_tmp;
                cutoff_radii_tmp[order][jkd][ikd] = cutoff_tmp;
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
                    exit("parse_cutoff_radii", "Incomplete cutoff radii");
                }
            }
        }
    }
    deallocate(undefined_cutoff);

    std::vector<double> cutoff_information_flatten;
    cutoff_information_flatten.resize(maxorder * nkd * nkd);

    i = 0;
    for (order = 0; order < maxorder; ++order) {
        for (j = 0; j < nkd; ++j) {
            for (k = 0; k < nkd; ++k) {
                cutoff_information_flatten[i++] = cutoff_radii_tmp[order][j][k];
            }
        }
    }
    deallocate(cutoff_radii_tmp);

    input_setter->set_cutoff_radii(maxorder,
                                   nkd,
                                   cutoff_information_flatten);
}

void InputParser::get_var_dict(const std::vector<std::string> &input_list,
                               std::map<std::string, std::string> &var_dict)
{
    std::string line, key, val;
    std::string line_wo_comment;
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

            boost::trim_left(line_wo_comment);
            if (line_wo_comment.empty()) continue;
            if (is_endof_entry(line_wo_comment)) break;

            // Split the input line by ';'

            boost::split(str_entry, line_wo_comment, boost::is_any_of(";"));

            for (auto &it : str_entry) {

                // Split the input entry by '='

                auto str_tmp = boost::trim_copy(it);

                if (!str_tmp.empty()) {

                    boost::split(str_varval, str_tmp, boost::is_any_of("="));

                    if (str_varval.size() != 2) {
                        exit("get_var_dict", "Unacceptable format");
                    }

                    key = boost::to_upper_copy(boost::trim_copy(str_varval[0]));
                    val = boost::trim_copy(str_varval[1]);

                    if (keyword_set.find(key) == keyword_set.end()) {
                        std::cout << "Could not recognize the variable " << key << std::endl;
                        exit("get_var_dict", "Invalid variable found");
                    }

                    if (var_dict.find(key) != var_dict.end()) {
                        std::cout << "Variable " << key << " appears twice in the input file." << std::endl;
                        exit("get_var_dict", "Redundant input parameter");
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

            for (auto &it : str_entry) {

                // Split the input entry by '='

                auto str_tmp = boost::trim_copy(it);

                if (!str_tmp.empty()) {

                    boost::split(str_varval, str_tmp, boost::is_any_of("="));

                    if (str_varval.size() != 2) {
                        exit("get_var_dict", "Unacceptable format");
                    }

                    key = boost::to_upper_copy(boost::trim_copy(str_varval[0]));
                    val = boost::trim_copy(str_varval[1]);

                    if (keyword_set.find(key) == keyword_set.end()) {
                        std::cout << "Could not recognize the variable "
                                  << key << std::endl;
                        exit("get_var_dict", "Invalid variable found");
                    }

                    if (var_dict.find(key) != var_dict.end()) {
                        std::cout << "Variable " << key
                                  << " appears twice in the input file." << std::endl;
                        exit("get_var_dict", "Redundant input parameter");
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

bool InputParser::is_data_range_consistent(const DispForceFile &datfile_in) const
{
    const auto ndata = datfile_in.ndata;
    const auto nstart = datfile_in.nstart;
    const auto nend = datfile_in.nend;
    const auto skip_s = datfile_in.skip_s;
    const auto skip_e = datfile_in.skip_e;

    if (nstart > 0 && skip_s > 0 && skip_s < nstart) return false;
    if (nend > 0 && skip_e > 0 && (skip_e - 1) > nend) return false;
    if (nstart > 0 && nend > 0 && nstart > nend) return false;
    if (skip_s > 0 && skip_e > 0 && skip_s >= skip_e) return false;

    if (ndata > 0) {
        if (nstart > 0 && nstart > ndata) return false;
        if (nend > 0 && nend > ndata) return false;
        if (skip_s > 0 && skip_s > ndata) return false;
        if (skip_e > 0 && (skip_e - 1) > ndata) return false;
    }

    return true;
}


int InputParser::locate_tag(const std::string key)
{
    auto ret = 0;
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

    }

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

bool InputParser::is_endof_entry(const std::string str)
{
    return str[0] == '/';
}

void InputParser::split_str_by_space(const std::string str,
                                     std::vector<std::string> &str_vec)
{
    std::string str_tmp;
    std::istringstream is(str);

    str_vec.clear();

    while (true) {
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
void InputParser::assign_val(T &val,
                             const std::string& key,
                             std::map<std::string, std::string> dict)
{
    // Assign a value to the variable "key" using the boost::lexica_cast.

    if (!dict[key].empty()) {
        try {
            val = boost::lexical_cast<T>(dict[key]);
        }
        catch (std::exception &e) {
            std::cout << e.what() << std::endl;
            auto str_tmp = "Invalid entry for the " + key + " tag.\n";
            str_tmp += " Please check the input value.";
            exit("assign_val", str_tmp.c_str());
        }
    }
}
