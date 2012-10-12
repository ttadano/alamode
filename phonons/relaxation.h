#pragma once

#include "pointers.h"
#include <complex>
#include <vector>

namespace PHON_NS {

    class ReciprocalVs {
    public:
        std::vector<unsigned int> ks;
        std::complex<double> v;

        ReciprocalVs();
        ReciprocalVs(const ReciprocalVs &obj){
            v = obj.v;
            for (std::vector<unsigned int>::const_iterator it = obj.ks.begin(); it != obj.ks.end(); ++it){
                ks.push_back(*it);
            }
        }

        ReciprocalVs(const std::complex<double> val, const std::vector<unsigned int> vec)
        {
            v = val;
            for (std::vector<unsigned int>::const_iterator it = vec.begin(); it != vec.end(); ++it){
                ks.push_back(*it);
            }
        }

        ReciprocalVs(const std::complex<double> val, const unsigned int *vec, const unsigned int n)
        {
            v = val;
            for (unsigned int i = 0; i < n; ++i){
                ks.push_back(vec[i]);
            }
        }
    };

    class Relaxation: protected Pointers {
    public:
        Relaxation(class PHON *);
        ~Relaxation();

        void setup_relaxation();
        void finish_relaxation();
        void calc_ReciprocalV();
        void calc_selfenergy();
        void calc_selfenergy_at_T(const double);
        void calc_damping(const unsigned int, double *, const double, const unsigned int, const unsigned int, double *);
        void calc_damping_tetra(const unsigned int, double *, const double, const unsigned int, const unsigned int, double *);
        std::complex<double> selfenergy(const double, const double, const unsigned int, const unsigned int);
        std::complex<double> selfenergy2(const double, const double, const unsigned int, const unsigned int);
        std::complex<double> *self_E;
        double self_tetra(const double, const double, const unsigned int, const unsigned int);
        double epsilon;
        int ksum_mode;

     private:
        struct StructKS {
        public:
            unsigned int ks1, ks2, ks3;
        };
        unsigned int nk, ns, nks;
        std::vector<ReciprocalVs> *V;
        void modify_eigenvectors();
        double delta_lorentz(const double);
        std::complex<double> im;
        std::complex<double> V3(const unsigned int, const unsigned int, const unsigned int);
        inline std::complex<double> V3new(const unsigned int [3]);
        inline std::complex<double> V3new2(const unsigned int [3]);
        double **vec_s;
        double ***relvec;
        double *invsqrt_mass_p;
        double mat_convert[3][3];
        double **e_tmp, **f_tmp;
        double ***fc3;
        double ***vec_for_v3, *invmass_for_v3;
        std::complex<double> ***cexp_phase;
    };
}
