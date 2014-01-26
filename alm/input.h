#pragma once

#include "pointers.h"
#include <string>
#include <map>
#include <vector>

namespace ALM_NS {
    class Input: protected Pointers {
    public:
        Input(class ALM *, int, char **);
        ~Input();
        void parce_input();

	private:
		int locate_tag(std::string);
		void split_str_by_space(const std::string, std::vector<std::string>&);
		void parse_general_vars();
		void parse_cell_parameter();
		void parse_interaction_vars();
		void parse_cutoff_radii();
		void parse_fitting_vars();
		void parse_atomic_positions();
		bool is_endof_entry(std::string);
		void get_var_dict(const std::string, std::map<std::string, std::string>&);
    };
}
