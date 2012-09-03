#pragma once

#include "pointers.h"
#include <string>
#include <vector>
#include <Eigen/Core>

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

    class Symmetry: protected Pointers {
    public:
        Symmetry(class PHON *);
        ~Symmetry();

        unsigned int nsym, nnp;
        std::string file_sym;
        std::vector<SymmetryOperation> SymmList;
        void setup_symmetry();
        void gensym(unsigned int, unsigned int&, unsigned int, double[3][3], double[3][3], double **, unsigned int *);
       

    private:
        void findsym(unsigned int, unsigned int, unsigned int *, double [3][3], double [3][3], double **);
        bool is_ortho(Eigen::Matrix3d, Eigen::Matrix3d, Eigen::Matrix3d);
        bool is_invariant(Eigen::Matrix3d, unsigned int, unsigned int *, double **, int[3], unsigned int);
    };
}