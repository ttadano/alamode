/*
 symmetry.h

 Copyright (c) 2014--2017 Terumasa Tadano

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

namespace ALM_NS
{
    class SymmetryOperation
    {
    public:
        int rotation[3][3]; // in lattice basis
        double tran[3]; // in Cartesian basis
        double rotation_cart[3][3]; // in Cartesian basis
        bool compatible_with_lattice;
        bool compatible_with_cartesian;
        bool is_translation;

        SymmetryOperation();

        SymmetryOperation(const int rot_in[3][3],
                          const double tran_in[3],
                          const double rot_cart_in[3][3],
                          const bool compatibility_lat,
                          const bool compatibility_cart,
                          const bool is_trans_in)
        {
            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                    rotation[i][j] = rot_in[i][j];
                    rotation_cart[i][j] = rot_cart_in[i][j];
                }
            }
            for (int i = 0; i < 3; ++i) {
                tran[i] = tran_in[i];
            }
            compatible_with_lattice = compatibility_lat;
            compatible_with_cartesian = compatibility_cart;
            is_translation = is_trans_in;
        }

        // Operator definition to sort
        bool operator<(const SymmetryOperation &a) const
        {
            std::vector<double> v1, v2;
            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                    v1.push_back(static_cast<double>(rotation[i][j]));
                    v2.push_back(static_cast<double>(a.rotation[i][j]));
                }
            }
            for (int i = 0; i < 3; ++i) {
                if (tran[i] < 0.0) {
                    v1.push_back(1.0 + tran[i]);
                } else {
                    v1.push_back(tran[i]);
                }
                if (a.tran[i] < 0.0) {
                    v2.push_back(1.0 + a.tran[i]);
                } else {
                    v2.push_back(a.tran[i]);
                }
            }
            return std::lexicographical_compare(v1.begin(), v1.end(),
                                                v2.begin(), v2.end());
        }
    };

    class RotationMatrix
    {
    public:
        int mat[3][3];

        RotationMatrix();

        RotationMatrix(const int rot[3][3])
        {
            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                    mat[i][j] = rot[i][j];
                }
            }
        }
    };

    class Symmetry: protected Pointers
    {
    public:
        Symmetry(class ALM *);
        ~Symmetry();

        void init();

        unsigned int nsym, ntran, nat_prim;
        int is_printsymmetry;
        int multiply_data;
        int *symnum_tran;

        double tolerance;

        int **map_sym;
        int **map_p2s;

        class Maps
        {
        public:
            int atom_num;
            int tran_num;
        };

        Maps *map_s2p;

        int trev_sym_mag;

        std::vector<SymmetryOperation> SymmData;

    private:

        void setup_symmetry_operation(int, unsigned int &,
                                      double [3][3], double [3][3],
                                      double **, int *);
        void genmaps(int, double **,
                     int **, int **,
                     class Symmetry::Maps *);

        void findsym(int, double [3][3], double **,
                     std::vector<SymmetryOperation> &);

        bool is_translation(const int [3][3]);
        bool is_proper(const double [3][3]);

        void symop_in_cart(double [3][3], const int [3][3],
                           const double [3][3], const double [3][3]);
        void pure_translations();
        void print_symmetrized_coordinate(double **);

        template <typename T>
        bool is_compatible(const T [3][3], const double tolerance_zero = 1.0e-5);

        void find_lattice_symmetry(double [3][3], std::vector<RotationMatrix> &);

        void find_crystal_symmetry(int, int,
                                   std::vector<unsigned int> *, double **x,
                                   std::vector<RotationMatrix>,
                                   std::vector<SymmetryOperation> &);

        std::string file_sym = "SYMM_INFO";
    };
}
