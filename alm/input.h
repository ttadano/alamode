/*
 input.h

 Copyright (c) 2014 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#pragma once

#include "pointers.h"
#include <fstream>
#include <string>
#include <map>
#include <vector>

namespace ALM_NS {
    class Input: protected Pointers {
    public:
        Input(class ALM *, int, char **);
        ~Input();
        void parce_input(int, char **);

    private:
        std::ifstream ifs_input;
        bool from_stdin;

        int locate_tag(std::string);
        void split_str_by_space(const std::string, std::vector<std::string>&);
        void parse_general_vars();
        void parse_cell_parameter();
        void parse_interaction_vars();
        void parse_cutoff_radii();
        void parse_cutoff_radii2();
        void parse_fitting_vars();
        void parse_atomic_positions();
        bool is_endof_entry(std::string);
        void get_var_dict(const std::string, std::map<std::string, std::string>&);
    };
}
