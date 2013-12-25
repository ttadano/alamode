#pragma once

#include "pointers.h"
#include <vector>
#include <set>

namespace PHON_NS {

    class Conductivity: protected Pointers {
    public:
        Conductivity(class PHON *);
        ~Conductivity();
        void setup_kappa();
		void prepare_restart();
		void calc_anharmonic_tau();
		void compute_kappa();
		void finish_kappa();

        int use_classical_Cv;
		unsigned int ntemp;
        double **tau;
		double ***kappa;
		double *Temperature;

    private:
        double ***vel;
        unsigned int nk, ns;
		int nshift_restart;
		double *tau_l;
		std::vector<int> vks, vks_l, vks_done;
		std::set<int> vks_job;
    };
}
