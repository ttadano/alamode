#pragma once

#include "pointers.h"
#include <complex>
#include <vector>
#include <string>

namespace PHON_NS {

	class Relaxation: protected Pointers {
    public:
        Relaxation(class PHON *);
        ~Relaxation();

        void setup_relaxation();
        void finish_relaxation();
        void setup_mode_analysis();
        void compute_mode_tau();
        void calc_damping(const unsigned int, double *, const double, const unsigned int, const unsigned int, double *);
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
        struct StructKS {
        public:
            unsigned int ks1, ks2, ks3;
        };
        unsigned int nk, ns, nks;
        std::vector<unsigned int> kslist;
        std::complex<double> im;

        double **e_tmp, **f_tmp;
        double ***vec_for_v3, *invmass_for_v3;
        double ***vec_for_v4, *invmass_for_v4;
        int **evec_index;
        int **evec_index4;
    };
}
