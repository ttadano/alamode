/*
 symmetry.h

 Copyright (c) 2014 Terumasa Tadano

 This file is distributed under the terms of the MIT license.
 Please see the file 'LICENCE.txt' in the root directory 
 or http://opensource.org/licenses/mit-license.php for information.
*/

#pragma once

#include "pointers.h"
#include <string>
#include <fstream>
#include <vector>

#ifdef _USE_EIGEN
#include <Eigen/Core>
#endif

namespace ALM_NS {

    class SymmetryOperation {
    public:
        int rot[3][3];
        double tran[3];

        SymmetryOperation();

        // Declaration construction

        SymmetryOperation(const int rot_in[3][3], const double tran_in[3])
        {
            for (int i = 0; i < 3; ++i){
                for (int j = 0; j < 3; ++j){
                    rot[i][j] = rot_in[i][j];
                }
            }
            for (int i = 0; i < 3; ++i){
                tran[i] = tran_in[i];
            }
        }
    };

    class RotationMatrix {
    public:
        int mat[3][3];

        RotationMatrix();
        RotationMatrix(const int rot[3][3]) {
            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                    mat[i][j] = rot[i][j];
                }
            }
        }
    };

    class Symmetry: protected Pointers {
    public:
        Symmetry(class ALM *);
        ~Symmetry();

        void init();

        unsigned int nsym, ntran, natmin;
        int is_printsymmetry;
        int multiply_data;
        int *symnum_tran;

        double tolerance;
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
        bool *sym_available;

    private:

        void setup_symmetry_operation(int, unsigned int&, double[3][3], double[3][3], 
				      double **, int *);
        void genmaps(int, double **, int **, int **, class Symmetry::Maps *);
        void findsym(int, double [3][3], double **, std::vector<SymmetryOperation> &);
        bool is_translation(int **);
        void symop_in_cart(double [3][3], double[3][3]);
        void pure_translations();
        void print_symmetrized_coordinate(double **);
        void symop_availability_check(double ***, bool *, const int, int &);
        void find_lattice_symmetry(double [3][3], std::vector<RotationMatrix> &);
        void find_crystal_symmetry(int, int, std::vector<unsigned int> *, double **x,
            std::vector<RotationMatrix>, std::vector<SymmetryOperation> &);
	//        void find_nnp_for_translation(unsigned int &, std::vector<SymmetryOperationTransFloat>);

        std::string file_sym;
        int ***symrel_int;
        std::vector<SymmetryOperation> SymmList;
    };
}

