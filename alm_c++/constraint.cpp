#include "constraint.h"
#include "interaction.h"
#include "memory.h"
#include <Eigen/Dense>
#include "fcs.h"
#include "symmetry.h"
#include "system.h"
#include "combination.h"
#include "constants.h"
#include "error.h"

using namespace ALM_NS;

Constraint::Constraint(ALM *alm) : Pointers(alm){};
Constraint::~Constraint() {
    if (constraint_mode != 0) {
        memory->deallocate(const_mat);
        memory->deallocate(const_rhs);
    }
};

void Constraint::setup(){

    if (constraint_mode != 0){

        // Generate constraint matrix

        int i;
        int maxorder = interaction->maxorder;
        int Pmax, order;
        int N = 0;

        for(i = 0; i < maxorder; ++i){
            N += fcs->ndup[i].size();
        }

        memory->allocate(const_translation, maxorder);
        translational_invariance();

        memory->allocate(const_rotation_self, maxorder);
        memory->allocate(const_rotation_cross, maxorder);
        rotational_invariance();

        Pmax = 0;
        for (order = 0; order < maxorder; ++order){
            //   Pmax += const_translation[order].size() + const_rotation[order].size();
            Pmax += const_translation[order].size();
        }
        if(constraint_mode == 2){
            //    Pmax -= const_translation[0].size() + const_rotation[0].size();
            Pmax -= const_translation[0].size();
            Pmax += fcs->ndup[0].size();
        }
        memory->allocate(const_mat, Pmax, N);
        memory->allocate(const_rhs, Pmax);

        calc_constraint_matrix(N, P);
        std::cout << "Total number of constraints: " << P << std::endl;

        memory->deallocate(const_translation);
        memory->deallocate(const_rotation_self);
        memory->deallocate(const_rotation_cross);
    }
}

void Constraint::calc_constraint_matrix(const int N, int &P){

    int i, j;
    int maxorder = interaction->maxorder;
    int order;
    int icol, irow;
    int nrow, ncol;
    int *nrank, *nparam;
    double *arr_tmp;
    std::vector<ConstraintClass> *const_vec;

    memory->allocate(nrank, maxorder);
    memory->allocate(nparam, maxorder);
    memory->allocate(const_vec, maxorder);

    std::cout << "Removing redundant constraints ...";

    using namespace Eigen;

    for(order = 0; order < maxorder; ++order){

        const_vec[order].clear();

        nrow = fcs->ndup[order].size();
        ncol = const_translation[order].size();

        nparam[order] = nrow;

        memory->allocate(arr_tmp, nrow);

        MatrixXd mat_tmp(nrow, ncol);
        icol = 0;

        for (std::set<ConstraintClass>::iterator p = const_translation[order].begin(); p != const_translation[order].end(); ++p){
            ConstraintClass const_now = *p;
            for (i = 0; i < nrow; ++i){
                mat_tmp(i, icol) = const_now.w_const[i];
            }
            ++icol;
        }

        FullPivLU<MatrixXd> lu_decomp(mat_tmp);
        nrank[order] = lu_decomp.rank();
        MatrixXd c_reduced = lu_decomp.image(mat_tmp);

        for(icol = 0; icol < nrank[order]; ++icol){
            for(irow = 0; irow < nrow; ++irow){
                arr_tmp[irow] = c_reduced(irow, icol);
            }

            const_vec[order].push_back(ConstraintClass(nrow, arr_tmp));
        }
        memory->deallocate(arr_tmp);
    }

    std::cout << " done." << std::endl << std::endl;

    std::cout << "Rank of the constraint matrices for each order ..." << std::endl;
    P = 0;
    for (order = 0; order < maxorder; ++order){
        P += nrank[order];
        std::cout << std::setw(9) << interaction->str_order[order] << ": " << const_vec[order].size() << std::endl;
    }
    std::cout << std::endl;

    if(constraint_mode == 2) {
        std::cout << "Harmonic Force Constants will be fixed to the values in the given reference file: " << fc2_file << std::endl;
        std::cout << "Constraint Matrix for Harmonic fcs will be updated." << std::endl << std::endl;
        P =  P - nrank[0] + nparam[0];
    }

    for(i = 0; i < P; ++i){
        for(j = 0; j < N; ++j){
            const_mat[i][j] = 0.0;
        }
        const_rhs[i] = 0.0;
    }

    int minorder= 0;

    irow = 0;
    icol = 0;

    if(constraint_mode == 2){
        std::ifstream ifs_fc2;
        ifs_fc2.open(fc2_file.c_str(), std::ios::in);
        if(!ifs_fc2) error->exit("calc_constraint_matrix", "cannot open file fc2_file");

        bool is_found = false;

        int nparam_harmonic;
        std::string str_tmp;
        while(!ifs_fc2.eof())
        {
            std::getline(ifs_fc2, str_tmp);

            if(str_tmp == "##HARMONIC FORCE CONSTANTS")
            {
                ifs_fc2 >> nparam_harmonic;
                if(nparam_harmonic != nparam[0]) error->exit("calc_constraint_matrix", "Number of fc2 not the same");

                is_found = true;

                for (i = 0; i < nparam[0]; ++i){
                    const_mat[i][i] = 1.0;
                    ifs_fc2 >> const_rhs[i];
                }
                break;
            }
        }
        minorder = 1;
        irow += nparam[0];
        icol += nparam[0];
        ifs_fc2.close();
        if(!is_found) error->exit("calc_constraint_matrix", "HARMONIC FORCE CONSTANTS flag not found in the fc2_file");
    }

    for(order = minorder; order < maxorder; ++order){

        for(std::vector<ConstraintClass>::iterator it = const_vec[order].begin(); it != const_vec[order].end(); ++it){
            ConstraintClass const_now = *it;

            for(i = 0; i < nparam[order]; ++i){
                const_mat[irow][i + icol] = const_now.w_const[i];
            }
            ++irow;
        }
        icol += nparam[order];
    }

    memory->deallocate(nrank);
    memory->deallocate(const_vec);
}

void Constraint::translational_invariance()
{
    // create constraint matrix arising from translational invariance.
    // might be a little tricky. (hard to follow)

    int i, j;
    int iat, jat, icrd, jcrd;
    int order;
    int maxorder = interaction->maxorder;

    int *ind;
    int *intarr, *intarr_copy;
    int **xyzcomponent;

    int nxyz;
    int natmin = symmetry->natmin;
    int nat = system->nat;

    double *arr_constraint;

    std::vector<int> intlist;
    std::set<FcProperty> list_found;
    std::set<FcProperty>::iterator iter_found;

    std::cout << "Start generating constraints for translational invariance ..." << std::endl;

    memory->allocate(ind, maxorder + 1);

    for (order = 0; order < maxorder; ++order){

        std::cout << std::setw(8) << interaction->str_order[order] << " ...";

        const_translation[order].clear();
        int nparams = fcs->ndup[order].size();

        if(nparams == 0) {
            std::cout << "No parameters! ... skipped."<< std::endl;
            continue;
        }

        // make interaction list

        list_found.clear();
        for(std::vector<FcProperty>::iterator p = fcs->fc_set[order].begin(); p != fcs->fc_set[order].end(); ++p){
            FcProperty list_tmp = *p; //using copy constructor
            for (i = 0; i < order + 2; ++i){
                ind[i] = list_tmp.elems[i];
            }
            if(list_found.find(FcProperty(order + 2, list_tmp.coef, ind, list_tmp.mother)) != list_found.end()) {
                error->exit("translational invariance", "Duplicate interaction list found");
            }
            list_found.insert(FcProperty(order + 2, list_tmp.coef, ind, list_tmp.mother));
        }

        // generate xyz component for each order

        nxyz = static_cast<int>(std::pow(static_cast<double>(3), order + 1));
        memory->allocate(xyzcomponent, nxyz, order + 1);
        fcs->get_xyzcomponent(order + 1, xyzcomponent);

        memory->allocate(arr_constraint, nparams);
        memory->allocate(intarr     , order + 2);
        memory->allocate(intarr_copy, order + 2);

        for(i = 0; i < natmin; ++i){

            iat = symmetry->map_p2s[i][0];

            // generate atom pairs for each order

            if(order == 0){
                for(icrd = 0; icrd < 3; ++icrd){

                    intarr[0] = 3 * iat + icrd;

                    for (jcrd = 0; jcrd < 3; ++jcrd){

                        // Reset the temporary array for next constraint
                        for (j = 0; j < nparams; ++j) arr_constraint[j] = 0.0;

                        for (jat = 0; jat < 3 * nat; jat += 3){
                            intarr[1] = jat + jcrd;

                            iter_found = list_found.find(FcProperty(order + 2, 1.0, intarr, 1));

                            // corresponding fcs found.
                            if(iter_found != list_found.end()){
                                FcProperty arrtmp = *iter_found;
                                arr_constraint[arrtmp.mother] += arrtmp.coef;
                            }
                        }
                        if(!is_allzero(nparams,arr_constraint)){
                            // add to constraint list
                            const_translation[order].insert(ConstraintClass(nparams, arr_constraint));
                        }
                    }
                }
            } else {
                for(j = 0; j < interaction->ninter[i][order]; ++j){
                    intlist.push_back(interaction->intpairs[i][order][j]);
                }
                std::sort(intlist.begin(), intlist.end());

                CombinationWithRepetition<int> g(intlist.begin(), intlist.end(), order);
                do {
                    std::vector<int> data = g.now();

                    intarr[0] = iat;
                    for (unsigned int isize = 0; isize < data.size(); ++isize){
                        intarr[isize + 1] = data[isize];
                    }

                    if(!interaction->is_incutoff(order + 1, intarr)) continue;

                    for(int ixyz = 0; ixyz < nxyz; ++ixyz){
                        for(int jcrd = 0; jcrd < 3; ++jcrd){

                            // Reset the temporary array for next constraint
                            for (j = 0; j < nparams; ++j) arr_constraint[j] = 0.0;

                            for(int jat = 0; jat < 3 * nat; jat += 3){
                                intarr[order + 1] = jat / 3;

                                if(!interaction->is_incutoff(order + 2, intarr)) continue;

                                for(j = 0; j < order + 1; ++j)  intarr_copy[j] = 3 * intarr[j] + xyzcomponent[ixyz][j];
                                intarr_copy[order + 1] = jat + jcrd;

                                fcs->sort_tail(order + 2, intarr_copy);

                                iter_found = list_found.find(FcProperty(order + 2, 1.0, intarr_copy, 1));
                                if(iter_found != list_found.end()){
                                    FcProperty arrtmp = *iter_found;
                                    arr_constraint[arrtmp.mother] += arrtmp.coef;                                
                                } 

                            }
                            if(!is_allzero(nparams,arr_constraint)){
                                const_translation[order].insert(ConstraintClass(nparams, arr_constraint));
                            }
                        }
                    }

                } while(g.next());
                intlist.clear();
            }
        }

        memory->deallocate(xyzcomponent);
        memory->deallocate(arr_constraint);
        memory->deallocate(intarr);
        memory->deallocate(intarr_copy);

        remove_redundant_rows(nparams, const_translation[order]);
        std::cout << " done." << std::endl;
    }
    memory->deallocate(ind);

    std::cout << "Finished !" << std::endl << std::endl;
    for(order = 0;  order < maxorder; ++order){
        std::cout << "Number of Translational Constraints for" << std::setw(9) << interaction->str_order[order] << " : " << const_translation[order].size() << std::endl;
    }
    std::cout << std::endl;
}

void Constraint::rotational_invariance()
{

    // A first implemention of complicated constraints
    // 2012.4.26

    std::cout << "Start generating constraint matrix for rotational invariance..." << std::endl;

#ifdef _DEBUG
    std::ofstream ofs_constraint;
    ofs_constraint.open("CONSTRAINT", std::ios::out);
#endif

    int i, j;
    int iat, jat;
    int icrd, jcrd;
    int order;
    int maxorder = interaction->maxorder;
    int natmin = symmetry->natmin;
    int nat = system->nat;
    int mu, nu;
    int ixyz, nxyz, nxyz2;
    int mu_lambda, lambda;
    int levi_factor;

    int *ind;
    int **xyzcomponent, **xyzcomponent2;
    int *nparams, nparam_sub;
    int *interaction_index, *interaction_atom;
    int *interaction_tmp;

    double *arr_constraint;
    double *arr_constraint_self;

    std::vector<int> interaction_list, interaction_list_old, interaction_list_now;

    std::set<FcProperty> list_found;
    std::set<FcProperty> list_found_last;
    std::set<FcProperty>::iterator iter_found;

    CombinationWithRepetition<int> g;

    memory->allocate(ind, maxorder + 1);
    memory->allocate(nparams, maxorder);

    for (order = 0; order < maxorder; ++order) {
        const_rotation_self[order].clear();
        const_rotation_cross[order].clear();
    }

    for (order = 0; order < maxorder; ++order){

        nparams[order] = fcs->ndup[order].size();

        if (order == 0) {
            std::cout << "Constraints between " << std::setw(8) << "1st-order IFCs (which are zero) and " 
                << std::setw(8) << interaction->str_order[order] << " ..." << std::endl;
            nparam_sub = nparams[order];
        } else {
            std::cout << "Constraints between " << std::setw(8) << interaction->str_order[order - 1] << " and "
                << std::setw(8) << interaction->str_order[order] << " ..." << std::endl;
            nparam_sub = nparams[order] + nparams[order - 1];
        }

        memory->allocate(arr_constraint, nparam_sub);
        memory->allocate(arr_constraint_self, nparams[order]);
        memory->allocate(interaction_atom , order + 2);
        memory->allocate(interaction_index, order + 2);
        memory->allocate(interaction_tmp  , order + 2);

        if (order > 0) {
            list_found_last = list_found;
            nxyz = static_cast<int>(pow(static_cast<double>(3), order));
            memory->allocate(xyzcomponent, nxyz, order);
            fcs->get_xyzcomponent(order, xyzcomponent);
        }

        list_found.clear();

        for (std::vector<FcProperty>::iterator p = fcs->fc_set[order].begin(); p != fcs->fc_set[order].end(); ++p){
            FcProperty list_tmp = *p; // using copy constructor
            for (i = 0; i < order + 2; ++i){
                ind[i] = list_tmp.elems[i];
            }
            list_found.insert(FcProperty(order + 2, list_tmp.coef, ind, list_tmp.mother));
        }

        for (i = 0; i < natmin; ++i){

            iat = symmetry->map_p2s[i][0];

            interaction_atom[0] = iat;

            if (order == 0){

                interaction_list_now.clear();
                for (j = 0; j < interaction->ninter[i][order]; ++j){
                    interaction_list_now.push_back(interaction->intpairs[i][order][j]);
                }
                std::sort(interaction_list_now.begin(), interaction_list_now.end());

                // special treatment for harmonic force constants

                for (icrd = 0; icrd < 3; ++icrd){

                    interaction_index[0] = 3 * iat + icrd;

                    for (mu = 0; mu < 3; ++mu){
                        for (nu = 0; nu < 3; ++nu){

                            if (mu == nu) continue;

                            // Clear history

                            for (j = 0; j < nparam_sub; ++j) arr_constraint[j] = 0.0;

                            for (std::vector<int>::iterator iter_list = interaction_list_now.begin(); iter_list != interaction_list_now.end(); ++iter_list){

                                jat = *iter_list;
                                interaction_index[1] = 3 * jat + mu;
                                iter_found = list_found.find(FcProperty(order + 2, 1.0, interaction_index, 1));
                                if(iter_found != list_found.end()){
                                    FcProperty arrtmp = *iter_found;              
                                    arr_constraint[arrtmp.mother] += arrtmp.coef * interaction->minvec[i][jat][nu];
                                }

                                // Exchange mu <--> nu and repeat again. 
                                // Note that the sign is inverted (+ --> -) in the summation

                                interaction_index[1] = 3 * jat + nu;
                                iter_found = list_found.find(FcProperty(order + 2, 1.0, interaction_index, 1));
                                if(iter_found != list_found.end()){
                                    FcProperty arrtmp = *iter_found;                        
                                    arr_constraint[arrtmp.mother] -= arrtmp.coef * interaction->minvec[i][jat][mu];                             
                                }
                            }

                            if(!is_allzero(nparam_sub,arr_constraint)){
                                // Add to constraint list
                                const_rotation_self[order].insert(ConstraintClass(nparam_sub, arr_constraint));
                            }

                        } // nu
                    }  // mu
                }
            } else {

                // constraint between force constants of different orders

                interaction_list_old.clear();
                interaction_list_now.clear();

                for (j = 0; j < interaction->ninter[i][order]; ++j){
                    interaction_list_now.push_back(interaction->intpairs[i][order][j]);
                }
                for (j = 0; j < interaction->ninter[i][order - 1]; ++j){
                    interaction_list_old.push_back(interaction->intpairs[i][order - 1][j]);
                }

                std::sort(interaction_list_now.begin(), interaction_list_now.end());
                std::sort(interaction_list_old.begin(), interaction_list_old.end());

                for (icrd = 0; icrd < 3; ++icrd){

                    interaction_index[0] = 3 * iat + icrd;

                    CombinationWithRepetition<int> g_now(interaction_list_now.begin(), interaction_list_now.end(), order);
                    CombinationWithRepetition<int> g_old(interaction_list_old.begin(), interaction_list_old.end(), order);

                    // m    -th order --> (m-1)-th order
                    // (m-1)-th order -->     m-th order
                    // 2-different directions to find all constraints

                    for (unsigned int direction = 0; direction < 2; ++direction){

                        if(direction == 0) {
                            g = g_now;
                            interaction_list = interaction_list_now;
                        } else {
                            g = g_old;
                            interaction_list = interaction_list_old;
                        }

                        // loop for the interacting pairs

                        do {
                            std::vector<int> data = g.now();

                            for (unsigned int idata = 0; idata < data.size(); ++idata)  interaction_atom[idata + 1] = data[idata];

                            for (ixyz = 0; ixyz < nxyz; ++ixyz){

                                for (j = 0; j < order; ++j) interaction_index[j + 1] = 3 * interaction_atom[j + 1] + xyzcomponent[ixyz][j];

                                for (mu = 0; mu < 3; ++mu){
                                    for (nu = 0; nu < 3; ++nu){

                                        if (mu == nu) continue;

                                        // Search for a new constraint below

                                        for (j = 0; j < nparam_sub; ++j) arr_constraint[j] = 0.0;
#ifdef _DEBUG
                                        ofs_constraint << "-------------------------------------------------------------" << std::endl;
                                        ofs_constraint << "New m_1 a_1 ... m_N a_N, mu_1 ... mu_N, mu, nu" << std::endl;
                                        for (j = 0; j < order + 1; ++j){
                                            ofs_constraint << std::setw(10) << fcs->easyvizint(interaction_index[j]);
                                        }
                                        ofs_constraint << ", mu = " << mu << ", nu = " << nu << std::endl;
                                        ofs_constraint << "Entering loop for m_N+1, a_N+1" << std::endl;
                                        ofs_constraint << "********************************" << std::endl;
#endif

                                        // loop for m_{N+1}, a_{N+1}
                                        for (std::vector<int>::iterator iter_list = interaction_list.begin(); iter_list != interaction_list.end(); ++iter_list){
                                            jat = *iter_list;

                                            interaction_atom[order + 1] = jat;
                                            if(!interaction->is_incutoff(order + 2, interaction_atom)) continue;

                                            // mu, nu

                                            interaction_index[order + 1] = 3 * jat + mu;
                                            for (j = 0; j < order + 2; ++j) interaction_tmp[j] = interaction_index[j];
#ifdef _DEBUG
                                            ofs_constraint << "m_N+1 a_N+1 = " << std::setw(5) << fcs->easyvizint(interaction_index[order + 1]);
#endif
                                            fcs->sort_tail(order + 2, interaction_tmp);

                                            iter_found = list_found.find(FcProperty(order + 2, 1.0, interaction_tmp, 1));
                                            if(iter_found != list_found.end()){
                                                FcProperty arrtmp = *iter_found;
#ifdef _DEBUG
                                                ofs_constraint << " --> Index = " << std::setw(10) << nparams[order - 1] + arrtmp.mother << ", Coefficient = " 
                                                    /* << std::setw(15) << arrtmp.coef * (system->x_cartesian[jat][nu] - system->x_cartesian[iat][nu]) << std::endl;
                                                    arr_constraint[nparams[order - 1] + arrtmp.mother] += arrtmp.coef * system->x_cartesian[jat][nu]; */
                                                    << std::setw(15) << arrtmp.coef * interaction->minvec[i][jat][nu] << std::endl;
#endif
                                                arr_constraint[nparams[order - 1] + arrtmp.mother] += arrtmp.coef * interaction->minvec[i][jat][nu];


                                            } else {
#ifdef _DEBUG
                                                ofs_constraint <<" --> Not Found!" << std::endl;
#endif
                                            }

                                            // Exchange mu <--> nu and repeat again.

                                            interaction_index[order + 1] = 3 * jat + nu;
                                            for (j = 0; j < order + 2; ++j) interaction_tmp[j] = interaction_index[j];
#ifdef _DEBUG
                                            ofs_constraint << "m_N+1 a_N+1 = " << std::setw(5) << fcs->easyvizint(interaction_index[order + 1]);
#endif
                                            fcs->sort_tail(order + 2, interaction_tmp);

                                            iter_found = list_found.find(FcProperty(order + 2, 1.0, interaction_tmp, 1));
                                            if(iter_found != list_found.end()){
                                                FcProperty arrtmp = *iter_found;
#ifdef _DEBUG
                                                ofs_constraint << " --> Index = " << std::setw(10) << nparams[order - 1] + arrtmp.mother << ", Coefficient = " 
                                                    /*  << std::setw(15) << -arrtmp.coef * (system->x_cartesian[jat][mu] - system->x_cartesian[iat][mu]) << std::endl;
                                                    arr_constraint[nparams[order - 1] + arrtmp.mother] -= arrtmp.coef * system->x_cartesian[jat][mu];*/
                                                    << std::setw(15) << arrtmp.coef * interaction->minvec[i][jat][mu] << std::endl;
#endif
                                                arr_constraint[nparams[order - 1] + arrtmp.mother] -= arrtmp.coef * interaction->minvec[i][jat][mu];
                                            } else {
#ifdef _DEBUG
                                                ofs_constraint << " --> Not Found!" << std::endl;
#endif
                                            }
                                        }
#ifdef _DEBUG
                                        ofs_constraint << "********************************" << std::endl;
                                        ofs_constraint << "Loop for lambda start" << std::endl;
#endif

                                        for (lambda = 0; lambda < order + 1; ++lambda){

#ifdef _DEBUG
                                            ofs_constraint << "lambda = " << lambda << std::endl;
#endif

                                            mu_lambda = interaction_index[lambda] % 3;

                                            for (jcrd = 0; jcrd < 3; ++jcrd){

                                                for (j = 0; j < order + 1; ++j) interaction_tmp[j] = interaction_index[j];

                                                interaction_tmp[lambda] = 3 * interaction_atom[lambda] + jcrd;

                                                levi_factor = 0;

                                                for (j = 0; j < 3; ++j){
                                                    levi_factor += levi_civita(j, mu, nu)*levi_civita(j, mu_lambda, jcrd); 
                                                }

                                                /*ofs_constraint << "#!" << "(" << std::setw(5) << mu << "," << std::setw(5) << mu_lambda << ")"
                                                << "(" << std::setw(5) << nu << "," << std::setw(5) << jcrd << ") - "
                                                << "(" << std::setw(5) << nu << "," << std::setw(5) << mu_lambda << ")"
                                                << "(" << std::setw(5) << mu << "," << std::setw(5) << jcrd << ")" << " = " << levi_factor << std::endl; */

                                                if(levi_factor == 0) continue;
#ifdef _DEBUG
                                                ofs_constraint << "mu_lambda = " << jcrd << ":";
                                                for (j = 0; j < order + 1; ++j){
                                                    ofs_constraint << std::setw(10) << fcs->easyvizint(interaction_tmp[j]);
                                                }
#endif

                                                fcs->sort_tail(order + 1, interaction_tmp);

                                                iter_found = list_found_last.find(FcProperty(order + 1, 1.0, interaction_tmp, 1));
                                                if(iter_found != list_found_last.end()){
                                                    FcProperty arrtmp = *iter_found;
#ifdef _DEBUG
                                                    ofs_constraint << " --> Index = " << std::setw(10) << arrtmp.mother << ", Coefficient = " 
                                                        << std::setw(15) << arrtmp.coef * static_cast<double>(levi_factor) << std::endl;
#endif
                                                    arr_constraint[arrtmp.mother] += arrtmp.coef * static_cast<double>(levi_factor);                                
                                                } else {
#ifdef _DEBUG
                                                    ofs_constraint << " --> Not Found!" << std::endl;
#endif
                                                }
                                            }
                                        }

                                        if(!is_allzero(nparam_sub,arr_constraint)){

                                            // A Candidate for another constraint found !
                                            // Add to the appropriate set

                                            if(is_allzero(nparam_sub, arr_constraint, nparams[order - 1])){
                                                const_rotation_self[order - 1].insert(ConstraintClass(nparams[order - 1], arr_constraint));                                      
                                            } else if (is_allzero(nparams[order - 1], arr_constraint)) {
                                                const_rotation_self[order].insert(ConstraintClass(nparam_sub, arr_constraint, nparams[order - 1]));
                                            } else {
                                                const_rotation_cross[order].insert(ConstraintClass(nparam_sub, arr_constraint)); 
                                            }
#ifdef _DEBUG
                                            for(j = 0; j < nparam_sub; ++j){
                                                if(std::abs(arr_constraint[j]) > eps10) {
                                                    ofs_constraint << "index = " << std::setw(5) <<  j << ": " 
                                                        << std::setw(5) << arr_constraint[j] << std::endl;
                                                }
                                            }
#endif
                                        }

                                    } // nu
                                } // mu

                            } // ixyz

                        } while(g.next());

                    } // direction
                } // icrd
            }

            // Additional constraint for the last order.s
            // All IFCs over maxorder-th order are neglected.

            if (order == maxorder - 1) {

                nxyz2 = static_cast<int>(pow(static_cast<double>(3), order + 1));
                memory->allocate(xyzcomponent2, nxyz2, order + 1);
                fcs->get_xyzcomponent(order + 1, xyzcomponent2);

                for (icrd = 0; icrd < 3; ++icrd){

                    interaction_index[0] = 3 * interaction_atom[0] + icrd;

                    CombinationWithRepetition<int> g_now(interaction_list_now.begin(), interaction_list_now.end(), order + 1);
                    do {

                        std::vector<int> data = g_now.now();

                        for (unsigned int idata = 0; idata < data.size(); ++idata)  interaction_atom[idata + 1] = data[idata];

                        for (ixyz = 0; ixyz < nxyz2; ++ixyz){



                            for (j = 0; j < order + 1; ++j) interaction_index[j + 1] = 3 * interaction_atom[j + 1] + xyzcomponent2[ixyz][j];
#ifdef _DEBUG
                            for (j = 0; j < order + 2; ++j){
                                ofs_constraint << "#$$ i_1 .. i_N = " << std::setw(5) << fcs->easyvizint(interaction_index[j]);
                            }
                            ofs_constraint << std::endl;
#endif

                            for (mu = 0; mu < 3; ++mu){
                                for (nu = 0; nu < 3; ++nu){

                                    //   ofs_constraint << "#$$ mu = " << mu << " nu = " << nu << std::endl;

                                    if (mu == nu) continue;

                                    for (j = 0; j < nparams[order]; ++j) arr_constraint_self[j] = 0.0;

                                    for (lambda = 0; lambda < order + 2; ++lambda){

                                        mu_lambda = interaction_index[lambda] % 3;

                                        for (jcrd = 0; jcrd < 3; ++jcrd){

                                            for (j = 0; j < order + 2; ++j) interaction_tmp[j] = interaction_index[j];

                                            interaction_tmp[lambda] = 3 * interaction_atom[lambda] + jcrd;

                                            levi_factor = 0;
                                            for (j = 0; j < 3; ++j){
                                                levi_factor += levi_civita(j, mu, nu)*levi_civita(j, mu_lambda, jcrd); 
                                            }

                                            if(levi_factor == 0) continue;

                                            fcs->sort_tail(order + 2, interaction_tmp);

                                            iter_found = list_found.find(FcProperty(order + 2, 1.0, interaction_tmp, 1));
                                            if(iter_found != list_found.end()){
                                                FcProperty arrtmp = *iter_found;
                                                arr_constraint_self[arrtmp.mother] += arrtmp.coef * static_cast<double>(levi_factor);                                
                                            }
                                        } // jcrd
                                    } // lambda

                                    if(!is_allzero(nparams[order], arr_constraint_self)){
                                        const_rotation_self[order].insert(ConstraintClass(nparams[order], arr_constraint_self));         
                                    }

                                } // nu
                            } // mu

                        } // ixyz


                    } while(g_now.next());

                } // icrd

                memory->deallocate(xyzcomponent2);
            }
        } // iat

        std::cout << " done" << std::endl;

        if (order > 0) {
            memory->deallocate(xyzcomponent);
        }
        memory->deallocate(arr_constraint);
        memory->deallocate(arr_constraint_self);
        memory->deallocate(interaction_tmp);
        memory->deallocate(interaction_index);
        memory->deallocate(interaction_atom);
    } // order

    for (order = 0; order < maxorder; ++order) {
        remove_redundant_rows(nparam_sub, const_rotation_cross[order]);
        remove_redundant_rows(nparams[order], const_rotation_self[order]);

        /*   if (order == 0) {
        nparam_sub = nparams[order];
        } else {
        nparam_sub = nparams[order] + nparams[order - 1];
        }
        std::cout << const_rotation_cross[order].size() << std::endl;
        for(std::set<ConstraintClass>::iterator p = const_rotation_cross[order].begin(); p != const_rotation_cross[order].end(); ++p){
        ConstraintClass const_tmp = *p;
        for (j = 0; j < nparam_sub; ++j){
        std::cout << std::setw(15) << std::scientific << const_tmp.w_const[j];
        }
        std::cout << std::endl;
        }

        std::cout << const_rotation_self[order].size() << std::endl;
        for(std::set<ConstraintClass>::iterator p = const_rotation_self[order].begin(); p != const_rotation_self[order].end(); ++p){
        ConstraintClass const_tmp = *p;
        for (j = 0; j < nparams[order]; ++j){
        std::cout << std::setw(15) << std::scientific << const_tmp.w_const[j];
        }
        std::cout << std::endl;
        }
        */
    }
    memory->deallocate(ind);
    memory->deallocate(nparams);
}

void Constraint::remove_redundant_rows(const int n, std::set<ConstraintClass> &Constraint_Set)
{
    using namespace Eigen;

    int nrow = n;
    int ncol = Constraint_Set.size();
    double *arr_tmp;

    if(ncol > 0) {
        memory->allocate(arr_tmp, nrow);
        MatrixXd mat_tmp(nrow, ncol);

        int icol = 0;

        for (std::set<ConstraintClass>::iterator p = Constraint_Set.begin(); p != Constraint_Set.end(); ++p){
            ConstraintClass const_now = *p;
            for (int i = 0; i < nrow; ++i){
                mat_tmp(i, icol) = const_now.w_const[i];
            }
            ++icol;
        }

        FullPivLU<MatrixXd> lu_decomp(mat_tmp);
        int nrank = lu_decomp.rank();
        MatrixXd c_reduced = lu_decomp.image(mat_tmp);

        Constraint_Set.clear();

        for(icol = 0; icol < nrank; ++icol){
            for(int irow = 0; irow < nrow; ++irow){
                arr_tmp[irow] = c_reduced(irow, icol);
            }

            Constraint_Set.insert(ConstraintClass(nrow, arr_tmp));
        }

        memory->deallocate(arr_tmp);
    }
}

int Constraint::levi_civita(const int i, const int j, const int k)
{
    int epsilon = (j - i) * (k - i) * (k - j) / 2;
    return epsilon;
}


bool Constraint::is_allzero(const int n, const double *arr, const int nshift){

    for(int i = nshift; i < n; ++i){
        if(std::abs(arr[i]) > eps10) {
            return false;
        }
    }
    return true;
}