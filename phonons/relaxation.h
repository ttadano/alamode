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

        double *tau;
        std::complex<double> *self_E;

        void test_delta(const double);

    private:
        std::vector<ReciprocalVs> *V;
        std::complex<double> delta_lorentz(const double);
        std::complex<double> im;
        std::complex<double> V3(const unsigned int, const unsigned int, const unsigned int);
        double **vec_s;
        double *mass_p;
        double epsilon;
        double mat_convert[3][3];
        double freq2(const double);
    };
}
