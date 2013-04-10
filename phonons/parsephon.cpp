#include "mpi_common.h"
#include "parsephon.h"
#include "error.h"
#include "gruneisen.h"
#include "system.h"
#include "kpoint.h"
#include "fcs_phonon.h"
#include "dynamical.h"
#include "mpi_common.h"
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

	if (!locate_tag("&general") && mympi->my_rank == 0) {
		error->exit("parse_input", "&general entry not found in the input file");
	}
	if (mympi->my_rank == 0) {
		parse_general_vars();
	}

	if (!locate_tag("&cell") && mympi->my_rank == 0) {
		error->exit("parse_input", "&cell entry not found in the input file");
	}
	if (mympi->my_rank == 0) {
		parse_cell_parameter();
	}

	if (!locate_tag("&kpoint")) {
		error->exit("parse_input", "&kpoint entry not found in the input file");
	}
	parse_kpoints();
}

void Input::parse_general_vars() {

	int i;
	std::string prefix, mode, fcsinfo, borninfo;
	int nsym, nnp, celldim[3], nbands, ismear;
	double Tmin, Tmax, dT, na_sigma, epsilon;
	double emin, emax, delta_e, delta_a;
	bool eigenvector, printxsf, nonanalytic, lclassical;

	std::string str_tmp;
	std::string str_allowed_list = "PREFIX MODE NSYM NNP CELLDIM FCSINFO TMIN TMAX DT EIGENVECTOR PRINTXSF NBANDS NONANALYTIC BORNINFO NA_SIGMA LCLASSICAL ISMEAR EPSILON EMIN EMAX DELTA_E DELTA_A";
	std::string str_no_defaults = "PREFIX MODE NSYM NNP FCSINFO";
	std::vector<std::string> no_defaults, celldim_v;
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
	
#ifdef _USE_BOOST
	boost::to_lower(mode);
	nsym = boost::lexical_cast<int>(general_var_dict["NSYM"]);
	nnp = boost::lexical_cast<int>(general_var_dict["NNP"]);
#else
	std::transform (mode.begin(), mode.end(), mode.begin(), tolower);
	nsym = my_cast<int>(general_var_dict["NSYM"]);
	nnp = my_cast<int>(general_var_dict["NNP"]);
#endif

	fcsinfo = general_var_dict["FCSINFO"];

	// Default values

	Tmin = 0.0;
	Tmax = 1000.0;
	dT = 1.0;

	emin = 0.0;
	emax = 1000.0;
	delta_e = 1.0;

	eigenvector = false;
	printxsf = false;
	nonanalytic = false;
	lclassical = false;

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

	assign_val(printxsf, "PRINTXSF", general_var_dict);

	if (printxsf || mode == "boltzmann") {
		eigenvector = true;
	} else {
		eigenvector = false;
	}

	assign_val(eigenvector, "EIGENVECTOR", general_var_dict);

	assign_val(nonanalytic, "NONANALYTIC", general_var_dict);
	assign_val(lclassical, "LCLASSICAL", general_var_dict);

	assign_val(nbands, "NBANDS", general_var_dict);
	assign_val(borninfo, "BORNINFO", general_var_dict);

	assign_val(ismear, "ISMEAR", general_var_dict);
	assign_val(epsilon, "EPSILON", general_var_dict);
	assign_val(na_sigma, "NA_SIGMA", general_var_dict);

	assign_val(delta_a, "DELTA_A", general_var_dict);


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

	if ((printxsf || mode == "boltzmann") & !eigenvector) {
		error->warn("parse_general_vars", "EIGENVECTOR is automatically changed to 1");
		eigenvector = true;
	}


	job_title = prefix;
	phon->mode = mode;
	symmetry->nsym = nsym;
	symmetry->nnp = nnp;

	system->Tmin = Tmin;
	system->Tmax = Tmax;
	system->dT = dT;

	dos->emax = emax;
	dos->emin = emin;
	dos->delta_e = delta_e;

	for (i = 0; i < 3; ++i) {
		system->cell_dimension[i] = celldim[i];
	}

	dynamical->eigenvectors = eigenvector;
	dynamical->nonanalytic = nonanalytic;
	dynamical->na_sigma = na_sigma;
	writes->writeanime = printxsf;
	writes->nbands = nbands;
	dynamical->file_born = borninfo;

	relaxation->epsilon = epsilon;
	fcs_phonon->file_fcs = fcsinfo;
	conductivity->use_classical_Cv = lclassical;

	gruneisen->delta_a = delta_a;
	relaxation->ksum_mode = ismear;

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

		if (!(kpmode >= 0 && kpmode <= 2)) {
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

			kpoint->kpInp.push_back(kpelem);
		}
	} else {

		ifs_input >> kpmode;

		if (!(kpmode >= 0 && kpmode <= 2)) {
			error->exit("parse_kpoints", "Invalid KPMODE");
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

			kpoint->kpInp.push_back(kpelem);
		}

	}

	kpoint->kpoint_mode = kpmode;
}

/*
void Input::read_input_phonons()
{
using namespace std;

string file_fcs;
double lavec[3][3];
int kpoint_mode;
bool eigenvectors, writeanime, nonanalytic;
int nbands;
unsigned int cell_dimension[3];

string file_born;
double na_sigma;
double Tmin, Tmax, dT;

unsigned int i, j;
unsigned int nsym, nnp;

cin >> nsym >> nnp;
for (i = 0; i < 3; ++i){
cin >> lavec[0][i] >> lavec[1][i] >> lavec[2][i];
}
// This is a tentative modification. 
// When cell_dimension = 1 or 2, algorithm will be slightly changed.

cin >> cell_dimension[0] >> cell_dimension[1] >> cell_dimension[2];
cin >> file_fcs;

cin >> eigenvectors >> writeanime >> nonanalytic;
cin >> nbands;
cin >> Tmin >> Tmax >> dT;
cin >> kpoint_mode;
if(nonanalytic) cin >> file_born >> na_sigma;

// distribute input parameters to each class

symmetry->nsym = nsym;
symmetry->nnp = nnp;

for (i = 0; i < 3; ++i){
for(j = 0; j < 3; ++j){
system->lavec_p[i][j] = lavec[i][j];
}
system->cell_dimension[i] = cell_dimension[i];
}

dynamical->eigenvectors = eigenvectors;
writes->writeanime = writeanime;
dynamical->nonanalytic = nonanalytic;
writes->nbands = nbands;

if(nonanalytic) {
dynamical->file_born = file_born;
dynamical->na_sigma = na_sigma;
}

fcs_phonon->file_fcs = file_fcs;
kpoint->kpoint_mode = kpoint_mode;

if (Tmin > Tmax) error->exit("read_input_phonon", "Tmin is larger than Tmax");

system->Tmin = Tmin;
system->Tmax = Tmax;
system->dT = dT;

}

void Input::read_input_boltzmann()
{
using namespace std;

string file_fcs;
double lavec[3][3];
double epsilon;
double Tmin, Tmax, dT;

unsigned int i, j;
unsigned int nsym, nnp;
unsigned int cell_dimension[3];
int ksum_mode;

int use_classical_Cv;

cin >> nsym >> nnp;
for (i = 0; i < 3; ++i) {
cin >> lavec[0][i] >> lavec[1][i] >> lavec[2][i];
}
for (i = 0; i < 3; ++i) {
cin >> cell_dimension[i];
}

cin >> file_fcs;
cin >> ksum_mode >> epsilon;
cin >> use_classical_Cv;
cin >> Tmin >> Tmax >> dT;

if (Tmin > Tmax) {
error->exit("read_input_boltzmann", "Tmin is bigger than Tmax");
}

symmetry->nsym = nsym;
symmetry->nnp = nnp;

for (i = 0; i < 3; ++i){
for (j = 0; j < 3; ++j){
system->lavec_p[i][j] = lavec[i][j];
}
system->cell_dimension[i] = cell_dimension[i];
}
fcs_phonon->file_fcs = file_fcs;
relaxation->ksum_mode = ksum_mode;
relaxation->epsilon = epsilon;

conductivity->use_classical_Cv = use_classical_Cv;

system->Tmin = Tmin;
system->Tmax = Tmax;
system->dT = dT;
}

void Input::read_input_gruneisen()
{
using namespace std;

string file_fcs;
double lavec[3][3];
int kpoint_mode;
bool eigenvectors, writeanime, nonanalytic;
int nbands;
unsigned int cell_dimension[3];

string file_born;
double na_sigma;
double delta_a;

unsigned int i, j;
unsigned int nsym, nnp;

cin >> nsym >> nnp;
for (i = 0; i < 3; ++i){
cin >> lavec[0][i] >> lavec[1][i] >> lavec[2][i];
}
// This is a tentative modification. 
// When cell_dimension = 1 or 2, algorithm will be slightly changed.

cin >> cell_dimension[0] >> cell_dimension[1] >> cell_dimension[2];
cin >> file_fcs;

cin >> eigenvectors >> writeanime >> nonanalytic;
cin >> nbands;
cin >> delta_a;
cin >> kpoint_mode;
if(nonanalytic) cin >> file_born >> na_sigma;

// distribute input parameters to each class

symmetry->nsym = nsym;
symmetry->nnp = nnp;

for (i = 0; i < 3; ++i){
for(j = 0; j < 3; ++j){
system->lavec_p[i][j] = lavec[i][j];
}
system->cell_dimension[i] = cell_dimension[i];
}

dynamical->eigenvectors = eigenvectors;
writes->writeanime = writeanime;
dynamical->nonanalytic = nonanalytic;
writes->nbands = nbands;

if(nonanalytic) {
dynamical->file_born = file_born;
dynamical->na_sigma = na_sigma;
}

fcs_phonon->file_fcs = file_fcs;
kpoint->kpoint_mode = kpoint_mode;

// Fractional change in cell parameters
gruneisen->delta_a = delta_a;
}

*/

int Input::locate_tag(std::string key){

	int ret = 0;
	std::string line, line2;

	if (from_stdin) {

		// An error detected in the following two lines when MPI is used.
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

		// An error detected in the following two lines when MPI is used.
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


void Input::get_var_dict(const std::string keywords, std::map<std::string, std::string> &var_dict) {

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

template<typename T> void Input::assign_val(T &val, std::string key, std::map<std::string, std::string> dict) {

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
