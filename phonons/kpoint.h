#pragma once

#include "pointers.h"
#include <iostream>
#include <string>
#include <vector>
#include <set>
#include "constants.h"

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

	class KpointInp {
	public:
		std::vector<std::string> kpelem;

		KpointInp(){};

		KpointInp(const std::vector<std::string> &obj) {
			for (std::vector<std::string>::const_iterator it = obj.begin(); it != obj.end(); ++it) {
				kpelem.push_back(*it);
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

	class KpointPlane {
	public:
		double k[3];
		int n[2];

		KpointPlane(){};

		KpointPlane(double *xk_in, int *n_in) {
			for (int i = 0; i < 3; ++i) k[i] = xk_in[i];
			for (int i = 0; i < 2; ++i) n[i] = n_in[i];
		}
	};

    class Kpoint: protected Pointers {
    public:
        Kpoint(class PHON *);
        ~Kpoint();

        void kpoint_setups(std::string);

        int kpoint_mode;
        unsigned int nkx, nky, nkz;
        unsigned int npath, nk;
        unsigned int *knum_minus;

        double **xk;
        double **kpoint_direction;
        double *kaxis;

        std::vector<KpointList> kpIBZ;
		std::vector<KpointInp> kpInp;
        std::vector<unsigned int> nk_equiv;
        std::vector<double> weight_k;
        std::set<unsigned int> kpset_uniq;

		unsigned int nplanes;
		std::vector<KpointPlane> *kp_planes;

        unsigned int **k_reduced, *nk_equiv_arr;
        unsigned int nk_reduced, nequiv_max;

        int get_knum(const double, const double, const double);
		void gen_kmesh(bool, unsigned int [3], double **, std::vector<unsigned int> &, std::vector<KpointList> &);

    private:
        void gen_kpoints_band();
		void reduce_kpoints(double **, unsigned int [3], std::vector<unsigned int> &, std::vector<KpointList> &);
        void gen_nkminus();
		void gen_kpoints_plane(std::vector<KpointInp>, std::vector<KpointPlane> *);
        bool in_first_BZ(double *);

        std::string **kp_symbol;
        double ***kp_bound;
        unsigned int *nkp;
    };
}
