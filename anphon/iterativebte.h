/*
 conductivity.h

 Copyright (c) 2014, 2015, 2016 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/


/*
implementation of the iterative bte
we can add additional parameter
    max_cycle
    convergence_criteria
we need to calculate and store: 
    Phi1, phi2, phi3 <- calculate once?
    phi2 should have the same size as V123 (q1,v1,q2,v2,v3) which is quite large
    maybe phi2, phi3 can be calculated for each iteration?
    v123 <- ? is it too big to store?

    => let's at least try it on silicon which is manageable to store

1) we need to calculate and store V3, phi2 and phi3, as well as the three phonon relaxation
   time
*/

#pragma once

#include "pointers.h"
#include "kpoint.h"
#include <vector>
#include <set>
#include <complex>

namespace PHON_NS {
    class Iterativebte : protected Pointers {
    public:
        Iterativebte(class PHON *);

        ~Iterativebte();

        void setup_iterative();  // initialize variables

        void do_iterativebte();  // wrapper

        void write_kappa();

        bool do_iterative;
        double *Temperature;
        unsigned int ntemp;
        
        int max_cycle;
        double convergence_criteria;  // dF(i+1) - dF(i) < cc
        
        double ***kappa;
        bool use_triplet_symmetry;
        bool sym_permutation;

    private:

        void set_default_variables();

        void deallocate_variables();

        int kplength;
        int nktot, nklocal, ns, ns2;

        // calculated at equilibrium
        double ***v3;
        double ***vel;

        double ***delta1;
        double ***delta2;
        double ***delta3;

        // linear response to deltaT
        double ***dFold;
        double ***dFnew;

        std::vector<KsListGroup> localnk_triplets;

        std::vector<int> num_unipair;
        std::vector<int> start_unipair;
        std::vector<int> nk_l, nk_job;

        void iterative_solver(); // calculate kappa iteratively

        void calc_kappa(int, double*** &, double** &);   // calculate kappa with off equilibrium part

        //void calc_phi_smear(int, double*** &, double*** &, double*** &);
        
        void calc_delta_smear(double*** &, double*** &, double*** &);

        void calc_delta_tetra(double*** &, double*** &, double*** &);

        void get_triplets();        // set up all triplets

        void setup_v3();            // calculate and store the scattering matrix elements
        
        void calc_Q_from_phi1(double** &, double** &);  // calculate Q

        void average_Q(double** &);

        void calc_boson(int, double** &, double** &);

        bool check_convergence(double** &, double** &); // check if convergence cirteria is meet

        void average_W_at_k(int, double** &);
        //void write_result_gamma(unsigned int,unsigned int,double ***,double **) const;

        //void average_self_energy_at_degenerate_point(int, int, double **) const;
        void write_result();
        void write_Q_dF(int , double** &, double*** &);

    };
}
