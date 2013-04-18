#pragma once

#include "pointers.h"
#include <vector>
#include <set>

namespace PHON_NS {

    class Conductivity: protected Pointers {
    public:
        Conductivity(class PHON *);
        ~Conductivity();
        void setup_kl();
        void finish_kl();
        void calc_kl();
		void prepare_restart();
		void calc_kl2();

        int use_classical_Cv;
        double **tau;

    private:
        double ***vel;
        void calc_kl_at_T(const double, double [3][3]);
        void calc_kl_mpi(const unsigned int, unsigned int **, double *, 
            unsigned int *, unsigned int **, const unsigned int, double *, double ***);
        unsigned int nk, ns;
		unsigned int ntemp;
		int nshift_restart;
		double *Temperature;
		double *tau_l;
		std::vector<int> vks, vks_l, vks_done;
		std::set<int> vks_job;
    };
}
