/*
 input_parser.h

 Copyright (c) 2014, 2015, 2016 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory
 or http://opensource.org/licenses/mit-license.php for information.
*/

#pragma once

#include "alm.h"
#include "input_setter.h"

#include <fstream>
#include <string>
#include <map>
#include <vector>


namespace ALM_NS {
    class InputParser {
    public:
        InputParser();

        ~InputParser();

        void run(ALM *alm,
                 const int narg,
                 const char *const *arg);

        std::string str_magmom;

        std::string get_run_mode() const;

    private:
        std::ifstream ifs_input;
        bool from_stdin{};
        std::string *kdname{};
        std::string mode;
        int maxorder{};
        size_t nat{};
        size_t nkd{};
        InputSetter *input_setter;

        void parse_input(ALM *alm);

        void parse_general_vars(ALM *alm);

        void parse_cell_parameter();

        void parse_atomic_positions();

        void parse_interaction_vars();

        void parse_cutoff_radii();

        void parse_optimize_vars(ALM *alm);

        int locate_tag(const std::string);

        static void split_str_by_space(const std::string,
                                std::vector<std::string> &) ;

        static bool is_endof_entry(const std::string) ;

        void get_var_dict(const std::vector<std::string> &,
                          std::map<std::string,
                                  std::string> &);

        bool is_data_range_consistent(const DispForceFile &datfile_in) const;

        template<typename T>
        void assign_val(T &,
                        const std::string&,
                        std::map<std::string, std::string>);

        void parse_displacement_and_force_files(std::vector<std::vector<double>> &u,
                                                std::vector<std::vector<double>> &f,
                                                DispForceFile &datfile_in) const;
    };
}
