#pragma once

#include "pointers.h"
#include <string>
#include <vector>

#ifdef _USE_EIGEN
#include <Eigen/Core>
#endif

namespace PHON_NS {

    class SymmetryOperation {
    public:
        std::vector<int> symop;

        SymmetryOperation();

        // Declaration construction
        
        SymmetryOperation(const SymmetryOperation &a)
        {
            for(std::vector<int>::const_iterator p = a.symop.begin(); p != a.symop.end(); ++p){
                symop.push_back(*p);
            }
        }

        SymmetryOperation(const int rot[3][3], const int trans[3])
        {
            for (int i = 0; i < 3; ++i){
                for (int j = 0; j < 3; ++j){
                    symop.push_back(rot[i][j]);
                }
            }
            for (int i = 0; i < 3; ++i){
                symop.push_back(trans[i]);
            }
        }
    };

    inline bool operator<(const SymmetryOperation a, const SymmetryOperation b){
        return std::lexicographical_compare(a.symop.begin(), a.symop.end(), b.symop.begin(), b.symop.end());
    }

	class SymmetryOperationWithMapping {
	public:
		std::vector<int> symop;
		std::vector<unsigned int> mapping;
		double shift[3];

		SymmetryOperationWithMapping();

		SymmetryOperationWithMapping(const SymmetryOperationWithMapping &a) 
		{
			for (std::vector<int>::const_iterator p = a.symop.begin(); p != a.symop.end(); ++p) {
				symop.push_back(*p);
			}
			for (std::vector<unsigned int>::const_iterator p = a.mapping.begin(); p != a.mapping.end(); ++p) {
				mapping.push_back(*p);
			}
			for (unsigned int i = 0; i < 3; ++i) {
				shift[i] = a.shift[i];
			}
		}

		SymmetryOperationWithMapping(const int S[3][3], unsigned int *mapping_info, const unsigned int n, const double shift_in[3])
		{
			unsigned int i, j;

			for (i = 0; i < 3; ++i) {
				for (j = 0; j < 3; ++j) {
					symop.push_back(S[i][j]);
				}
			}
			for (i = 0; i < n; ++i) {
				mapping.push_back(mapping_info[i]);
			}
			for (i = 0; i < 3; ++i) {
				shift[i] = shift_in[i];
			}
		} 
	};

	inline bool operator<(const SymmetryOperationWithMapping a, const SymmetryOperationWithMapping b) {
		return std::lexicographical_compare(a.symop.begin(), a.symop.end(), b.symop.begin(), b.symop.end());
	}

    class Symmetry: protected Pointers {
    public:
        Symmetry(class PHON *);
        ~Symmetry();

        unsigned int nsym, nnp;
        bool symmetry_flag, time_reversal_sym;
        std::string file_sym;
        std::vector<SymmetryOperation> SymmList;
		std::vector<SymmetryOperationWithMapping> SymmListWithMap;
        void setup_symmetry();
        void gensym(unsigned int, unsigned int&, unsigned int, double[3][3], double[3][3], double **, unsigned int *);
		

    private:
        void findsym(unsigned int, unsigned int, unsigned int *, double [3][3], double [3][3], double **);
		void gensym_withmap(double **, unsigned int *);

#ifdef _USE_EIGEN
        bool is_ortho(Eigen::Matrix3d, Eigen::Matrix3d, Eigen::Matrix3d);
        bool is_invariant(Eigen::Matrix3d, unsigned int, unsigned int *, double **, int[3], unsigned int);
#else
		bool is_ortho(double [3][3], double [3][3], double [3][3]);
		bool is_invariant(double [3][3], unsigned int, unsigned int *, double **, int [3], unsigned int);
#endif
    };
}
