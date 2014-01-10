#pragma once

#include "pointers.h"
#include <complex>
#include <vector>
#include <string>

namespace PHON_NS {

	class KsList {
	public:
		std::vector<int> ks;

		KsList();

		KsList(const KsList &a)
		{
			for(std::vector<int>::const_iterator p = a.ks.begin(); p != a.ks.end(); ++p){
				ks.push_back(*p);
			}
		}

		KsList(const int n, int *ks_in) {
			for (int i = 0; i < n; ++i) {
				ks.push_back(ks_in[i]);
			}
		}
	};

    inline bool operator<(const KsList a, const KsList b){
			return std::lexicographical_compare(a.ks.begin(), a.ks.end(), b.ks.begin(), b.ks.end());
	}

	class KsListGroup {
	public:
		std::vector<KsList> group;

		KsListGroup();
		KsListGroup(const std::vector<KsList> &a) {
			for (std::vector<KsList>::const_iterator it = a.begin(); it != a.end(); ++it) {
				group.push_back(*it);
			}
		}
	};

	class Relaxation: protected Pointers {
	public:
		Relaxation(class PHON *);
		~Relaxation();

		void setup_relaxation();
		void finish_relaxation();
		void setup_mode_analysis();
		void compute_mode_tau();
		void calc_damping(const unsigned int, double *, const double, const unsigned int, const unsigned int, double *);
		void calc_damping2(const unsigned int, double *, const double, const unsigned int, const unsigned int, double *);
		void calc_damping_atom(const unsigned int, double *, const double, const unsigned int, const unsigned int, double ***);
		void calc_damping_tetra(const unsigned int, double *, const double, const unsigned int, const unsigned int, double *);	 
		void calc_damping_tetra_atom(const unsigned int, double *, const double, const unsigned int, const unsigned int, double ***); 
		void calc_realpart_V4(const unsigned int, double *, const double, const unsigned int, const unsigned int, double *);
		double epsilon;
		int ksum_mode;
		bool quartic_mode;
		bool ks_analyze_mode;
		bool atom_project_mode;
		bool calc_realpart;
		std::string ks_input;

		double delta_lorentz(const double);
		double delta_gauss(const double);

		std::complex<double> V3(const unsigned int [3]);
		std::complex<double> V4(const unsigned int [4]);

	private:
		unsigned int nk, ns, nks;
		std::vector<unsigned int> kslist;
		std::complex<double> im;

		double **e_tmp, **f_tmp;
		double ***vec_for_v3, *invmass_for_v3;
		double ***vec_for_v4, *invmass_for_v4;
		int **evec_index;
		int **evec_index4;

	    std::vector<KsListGroup> *pair_uniq;

		void gensym_kpairs();
		void gen_pair_uniq();

		int knum_sym(const int, const int);
	};
}
