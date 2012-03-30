#ifndef ALM_SYMMETRY_HEADER
#define ALM_SYMMETRY_HEADER

#include "pointers.h"
#include <Eigen/Core>
#include <string>
#include <fstream>
#include <vector>

namespace ALM_NS {
    
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
        Symmetry(class ALM *);
        ~Symmetry();

        void init();
        void gensym(int, int&, int, double[3][3], double[3][3], 
            double **, int *);
     /*   void findsym(int, int, int *, double[3][3], double[3][3],
            double **, int &, int ***, int **);*/
        void findsym(int, int, int *, double[3][3], double[3][3],
            double **);

        int nsym, nnp;
        int ntran, natmin;

        bool multiply_data;

        int ***symrel_int;
        int *symnum_tran;
        double ***symrel;
        double **tnons;

        int **map_sym;
        int **map_p2s;
        class Maps {
        public:
            int atom_num;
            int tran_num;
        };
        Maps *map_s2p;

        void genmaps(int, double **, int **, int **, class Symmetry::Maps *);


        std::string file_sym;
        std::ofstream ofs_sym;
        std::ifstream ifs_sym;

    private:
        bool is_ortho(Eigen::Matrix3d, Eigen::Matrix3d, Eigen::Matrix3d);
        bool is_invariant(Eigen::Matrix3d, int, int*, double **, int[3], int);
        void symop_in_cart(double [3][3], double[3][3]);
        void pure_translations();
        void data_multiplier(int, int);

        std::vector<SymmetryOperation> SymmList;
    };

    
}
#endif
