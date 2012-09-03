#pragma once

#include "pointers.h"
#include <string>
#include <vector>
#include "../alm_c++/constants.h"

namespace PHON_NS {

    class KpointList {
    public:
        std::vector<double> kval;
        unsigned int knum;

        KpointList(){};
        KpointList(const KpointList &obj){
            knum = obj.knum;
            for (std::vector<double>::const_iterator it = obj.kval.begin(); it != obj.kval.end(); ++it){
                kval.push_back(*it);
            }
        }

        KpointList(const unsigned int knum_in, const std::vector<double> vec)
        {
            knum = knum_in;
            for (std::vector<double>::const_iterator it = vec.begin(); it != vec.end(); ++it){
                kval.push_back(*it);
            }
        }

    };

    inline bool operator<(const KpointList a, const KpointList b){
        return std::lexicographical_compare(a.kval.begin(), a.kval.end(), b.kval.begin(), b.kval.end());
    }

    inline bool operator==(const KpointList a, const KpointList b){
        double tmp = 0.0;
        for (unsigned int i = 0; i < 3; ++i){
            tmp += std::pow(a.kval[i] - b.kval[i], 2);
        }
        return std::sqrt(tmp) < eps;
    }

    class Kpoint: protected Pointers {
    public:
        Kpoint(class PHON *);
        ~Kpoint();

        void kpoint_setups();
        void reduce_kpoints();

        int kpoint_mode;
        unsigned int nkx, nky, nkz;

        unsigned int npath, nk;

        unsigned int *knum_minus;

        double **xk;
        double **kpoint_direction;
        double *kaxis;

        std::vector<KpointList> kpIBZ;
        std::vector<unsigned int> nk_equiv;
        std::vector<double> weight_k;

    private:
        void gen_kpoints_band();
        void gen_kmesh();
        void gen_nkminus();
       
        std::string **kp_symbol;
        double ***kp_bound;
        unsigned int *nkp;
    };
}
