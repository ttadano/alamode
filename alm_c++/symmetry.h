#ifndef ALM_SYMMETRY_HEADER
#define ALM_SYMMETRY_HEADER

#include "pointers.h"
#include <Eigen/core>
#include <string>
#include <fstream>

namespace ALM_NS {
    class Symmetry: protected Pointers {
    public:
        Symmetry(class ALM *);
        ~Symmetry();

        void init();
        void gensym(int, int&, int, double[3][3], double[3][3], 
            double **, int *, int ***, double **);
        void findsym(int, int, int *, double[3][3], double[3][3],
            double **, int &, int ***, int **);

        
        int nsym, nnp;
        int ntran, natmin;
        int maxsym;

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
    };
}
#endif