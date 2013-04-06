#pragma once

#include "pointers.h"
#include <string>
#include <map>
#include <vector>

namespace PHON_NS {
    class Input: protected Pointers {
    public:
        Input(class PHON *, int, char **);
        ~Input();
        void parce_input();

        std::string job_title;

    private:
        void read_input_phonons();
        void read_input_boltzmann();
        void read_input_gruneisen();

		int locate_tag(std::string);
		void parse_general_vars();
		void parse_cell_parameter();
		void parse_kpoints();
		void get_var_dict(std::string, std::map<std::string, std::string> &);
		void split_str_by_space(const std::string, std::vector<std::string>&);
		bool is_endof_entry(std::string);
		template<typename T> void assign_val(T&, std::string, std::map<std::string, std::string>);
    };
}
