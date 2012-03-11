#ifndef ALM_FITTING_HEADER
#define ALM_FITTING_HEADER

#include "pointers.h"
#include <vector>
#include <set>
#include <string>

namespace ALM_NS {

    class Constraint {
    public:
        std::vector<double> w_const;

        Constraint();
        Constraint(const Constraint &a){
            for(std::vector<double>::const_iterator p = a.w_const.begin(); p != a.w_const.end(); ++p){
                w_const.push_back(*p);
            }
        }
        Constraint(const int n, const double *arr){
            for(int i = 0; i < n; ++i){
                w_const.push_back(arr[i]);
            }
        }
    };

    inline bool operator<(const Constraint a, const Constraint b){
        return std::lexicographical_compare(a.w_const.begin(), a.w_const.end(), b.w_const.begin(), b.w_const.end());
    }

    class Fitting: protected Pointers {
    public:
        Fitting(class ALM *);
        ~Fitting();

        void fitmain();
       // int rank(const int, const int, double **);
        int rank(int, int, double *);
        // int getRankEigen(const int, const int,const  int);
        int constraint;
        std::string fc2_file;

        double *params;

    private:

        int inprim_index(const int);
        void wrtfcs(const double *);
        void fit_without_constraints(int, int, int);
        void fit_with_constraints(int, int, int, int);
        void fit_consecutively(int, int, const int, const int,
            const int, const int, const int, const int);
        void calc_matrix_elements(const int, const int, const int, 
            const int, const int, const int, const int);
        double gamma(const int, const int *);
        int factorial(const int);

        void translational_invariance();
        void rotational_invariance();
        void calc_constraint_matrix(const int, int &);
        bool is_allzero(const int, const double *);

        double **amat;
        double *fsum;

        double **const_mat;
        double *const_rhs;
        std::set<Constraint> *const_translation;
    };

    extern "C" {
        
        void dgelss_(int *m, int *n, int *nrhs, double *a, int *lda,	
        double *b, int *ldb, double *s, double *rcond, int *rank,
        double *work, int *lwork, int *info);
       
        void dgglse_(int *m, int *n, int *p, double *a, int *lda,
        double *b, int *ldb, double *c, double *d, double *x,
        double *work, int *lwork, int *info);

        void dgesdd_(const char *jobz, int *m, int *n, double *a, int *lda,
        double *s, double *u, int *ldu, double *vt, int *ldvt, double *work,
        int *lwork, int *iwork, int *info);

        void dgeqrf_(int *m, int *n, double *a, double *tau,
            double *work, int *lwork, int *info);
    }

}
#endif