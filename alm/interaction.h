#pragma once

#include <string>
#include <vector>
#include <set>
#include "listcomparison.h"
#include "pointers.h"
#include "constants.h"


namespace ALM_NS {
    class InteractionCluster {
    public:
        std::vector<double> x;

        InteractionCluster();
        InteractionCluster(const double *arr){
            for (int i = 0; i < 3; ++i){
                x.push_back(arr[i]);
            }
        }
    };

    inline bool operator<(const InteractionCluster a, const InteractionCluster b){
        return std::lexicographical_compare(a.x.begin(), a.x.end(), b.x.begin(), b.x.end());
    }

	class DistInfo {
	public:
		int cell;
		double dist;
		double relvec[3];

		DistInfo();
		DistInfo(const int n, const double d, const double x[3]) {
			cell = n;
			dist = d;
			for (int i = 0; i < 3; ++i) relvec[i] = x[i];
		}

		DistInfo(const DistInfo &obj) {
			cell = obj.cell;
			dist = obj.dist;
			for (int i = 0; i < 3; ++i) relvec[i] = obj.relvec[i];
		}
	};

	inline bool operator<(const DistInfo a, const DistInfo b) {
		return a.dist < b.dist;
	}

	class DistList {
	public:
		int atom;
		double dist;

		DistList();
		DistList(int n, double dist_tmp) {
			atom = n;
			dist = dist_tmp;
		}
	};
	inline bool operator<(const DistList a, const DistList b) {
		if (std::abs(a.dist - b.dist) > eps8) {
			return a.dist < b.dist;
		} else {
			return a.atom < b.atom;
		}
	}

    class Interaction: protected Pointers {
    public:
        Interaction(class ALM *);
        ~Interaction();
        void init();

		int *nbody_include;

        double ***rcs;
        double **distlist;
		std::vector<DistInfo> **mindist_pairs;
		std::set<IntList> *pairs;


        double distance(double *, double *);

        bool is_periodic[3];

        std::string *str_order;

        int nneib;
        int maxorder;
		int interaction_type;
        double ***xcrd;

        int ***intpairs;
        int **ninter;
        double ****relvec;

        double ***minvec;

        bool is_incutoff(int, int *);

        template <typename T>
        T maxval(int n, T *arr)
        {
            T tmp;
            tmp = arr[0];

            for (int i = 0; i < n; i++) {
                tmp = std::max<T>(tmp, arr[i]);
            }
            return tmp;
        }

        template <typename T>
        T maxval(int n1, int n2, T **arr)
        {
            T tmp;
            tmp = arr[0][0];

            for (int i = 0; i < n1; i++) {
                for (int j = 0; j < n2; j++){
                    tmp = std::max<T>(tmp, arr[i][j]);
                } 
            }
            return tmp;
        }

        template <typename T>
        T maxval(int n1, int n2, int n3, T ***arr)
        {
            T tmp;
            tmp = arr[0][0][0];

            for (int i = 0; i < n1; i++) {
                for (int j = 0; j < n2; j++){
                    for (int k = 0; k < n3; k++){
                        tmp = std::max<T>(tmp, arr[i][j][k]);
                    } 
                }
            }
            return tmp;
        }

        template <typename T>
        void insort(int n, T *arr)
        {
            int i, j;
            T tmp;

            for (i = 1; i < n; ++i){
                tmp = arr[i];
                for (j = i - 1; j >= 0 && arr[j] > tmp; --j){
                    arr[j + 1] = arr[j];
                }
                arr[j + 1] = tmp;
            }
        }

    private:
        int nsize[3];
        void calc_distlist(int, double **);
        void search_interactions();
        void set_ordername();
        void calc_minvec();
		int nbody(const int, const int *);

		std::vector<DistInfo> **distall;
		std::set<IntList> *interacting_atom_pairs;
    };
}
