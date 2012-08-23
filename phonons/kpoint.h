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

        bool flag = true;
        std::vector<double>::const_iterator first1, first2, last1, last2;

        first1 = a.kval.begin();
        first2 = b.kval.begin();
        last1 = a.kval.end();
        last2 = b.kval.end();

        while (first1 != last1)
        {
            if (first2 == last2 || *first2 < *first1) { 
                return false;
            } else if (*first1 < *first2) {
                return true;
            }
            ++first1;
            ++first2;
        }
        return first2 != last2;

/*        for (unsigned int i = 0; i < a.kval.size(); ++i){
            flag = flag && (a.kval[i] < b.kval[i]) && (std::abs(a.kval[i] - b.kval[i]) > eps12);
        }
        return flag;
*/
    }
    
    class Kpoint: protected Pointers {
    public:
        Kpoint(class PHON *);
        ~Kpoint();

      void kpoint_setups();

      int kpoint_mode;
      unsigned int nkx, nky, nkz;
      
      unsigned int npath, nk;

      unsigned int *knum_minus;

      double **xk;
      double **kpoint_direction;
      double *kaxis;

    private:
        void gen_kpoints_band();
        void gen_kmesh();
        void gen_nkminus();
        std::string **kp_symbol;
        double ***kp_bound;
        unsigned int *nkp;
    };
}
